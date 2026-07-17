/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include <algorithm>
#include <cfloat>

#include "BLI_listbase_iterator.hh"
#include "BLI_set.hh"
#include "BLI_string_utf8.hh"

#include "DNA_node_types.h"

#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_runtime.hh"
#include "BKE_node_tree_update.hh"

#include "NOD_closure_zone.hh"
#include "NOD_geo_bundle.hh"
#include "NOD_mask_stack.hh"
#include "NOD_node_placement.hh"
#include "NOD_socket.hh"
#include "NOD_socket_items.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"

namespace blender::nodes::texture_stack {

bool is_stack(const bNode &node)
{
  return node.is_type("ShaderNodeTextureLayerStack"_ustr);
}

NodeShaderLayerStackItem &StackLayer::item() const
{
  NodeShaderLayerStack &data = layer_stack::storage(*this->stack);
  BLI_assert(this->index >= 0 && this->index < data.items_num);
  return data.items[this->index];
}

bNodeSocket *StackLayer::socket(const LayerSocketRole role) const
{
  const NodeShaderLayerStack &data = layer_stack::storage(*this->stack);
  if (this->index < 0 || this->index >= data.items_num) {
    return nullptr;
  }
  const NodeShaderLayerStackItem &item = data.items[this->index];
  /* Both stack node kinds share the Opacity socket naming; the item socket
   * itself is named per kind, and a mask layer has no mask of its own. */
  const bool is_mask_layer = mask_stack::is_stack(*this->stack);
  std::string id;
  switch (role) {
    case LayerSocketRole::Bundle:
      id = is_mask_layer ? mask_stack::ItemsAccessor::socket_identifier_for_item(item) :
                           ItemsAccessor::socket_identifier_for_item(item);
      break;
    case LayerSocketRole::Opacity:
      id = ItemsAccessor::opacity_socket_identifier_for_item(item);
      break;
    case LayerSocketRole::Mask:
      if (is_mask_layer) {
        return nullptr;
      }
      id = ItemsAccessor::mask_socket_identifier_for_item(item);
      break;
  }
  return bke::node_find_socket(*this->stack, SOCK_IN, UString(id));
}

bNode *StackLayer::source() const
{
  bNodeSocket *sock = this->bundle_socket();
  if (sock == nullptr) {
    return nullptr;
  }
  return bke::node_find_source_node(*this->ntree, *sock);
}

bNode *StackLayer::mask_stack_node() const
{
  bNodeSocket *sock = this->mask_socket();
  if (sock == nullptr) {
    return nullptr;
  }
  bNode *source = bke::node_find_source_node(*this->ntree, *sock);
  if (source && mask_stack::is_stack(*source)) {
    return source;
  }
  return nullptr;
}

bNode *StackLayer::nested_stack() const
{
  bNode *src = this->source();
  if (src == nullptr) {
    return nullptr;
  }
  /* A generator group: the layer's Bundle input is fed directly by a nested
   * stack's Result. */
  if (is_stack(*src)) {
    return src;
  }
  /* An adjustment group: the layer is a closure whose zone body is a nested
   * stack (as opposed to a single adjustment, whose body is a plain node
   * group). */
  if (closure_zone::is_output_node(*src)) {
    if (bNode *body = closure_zone::body_node(*this->ntree, *src)) {
      if (is_stack(*body)) {
        return body;
      }
    }
  }
  return nullptr;
}

bNode *StackLayer::display_source() const
{
  bNode *src = this->source();
  if (src && closure_zone::is_output_node(*src)) {
    if (bNode *body = closure_zone::body_node(*this->ntree, *src)) {
      return body;
    }
  }
  return src;
}

void roots_for_bsdf(const bNodeTree &ntree, bNode &bsdf, Vector<bNode *> &r_roots)
{
  Set<bNode *> seen;
  for (bNodeSocket &sock : bsdf.inputs) {
    bNode *source = bke::node_find_source_node(ntree, sock);
    if (source == nullptr || !texture_channel::is_separate_bundle(*source)) {
      continue;
    }
    bNodeSocket *bundle_in = first_bundle_input(*source);
    if (bundle_in == nullptr) {
      continue;
    }
    bNode *bundle_source = bke::node_find_source_node(ntree, *bundle_in);
    if (bundle_source == nullptr || !is_stack(*bundle_source)) {
      continue;
    }
    if (seen.add(bundle_source)) {
      r_roots.append(bundle_source);
    }
  }
}

bNode *active(bNodeTree &ntree)
{
  bNode *active_node = bke::node_get_active(ntree);
  if (active_node == nullptr) {
    return nullptr;
  }
  /* The active node is usually a layer's source rather than the stack itself
   * (see #layer_for_node), so resolve it to the stack holding that layer. */
  if (const std::optional<StackLayer> layer = layer_for_node(ntree, *active_node)) {
    active_node = layer->stack;
  }
  /* A Mask Stack (a selected mask row) resolves to the Texture Layer Stack its
   * result ultimately feeds, possibly through wedged mask stacks. */
  Set<const bNode *> visited;
  while (active_node && mask_stack::is_stack(*active_node) && visited.add(active_node)) {
    bNode *consumer = nullptr;
    for (const bNodeLink &link : ntree.links) {
      if (link.fromnode == active_node && link.tonode) {
        consumer = link.tonode;
        break;
      }
    }
    active_node = consumer;
  }
  if (active_node == nullptr) {
    return nullptr;
  }
  if (is_stack(*active_node)) {
    return active_node;
  }
  if (!texture_channel::is_bsdf(*active_node)) {
    return nullptr;
  }
  Vector<bNode *> roots;
  roots_for_bsdf(ntree, *active_node, roots);
  return roots.is_empty() ? nullptr : roots[0];
}

static bNode *bsdf_for_node_recursive(const bNodeTree &ntree,
                                      bNode &node,
                                      Set<const bNode *> &visited)
{
  if (!visited.add(&node)) {
    return nullptr;
  }
  if (texture_channel::is_bsdf(node)) {
    return &node;
  }
  for (const bNodeLink &link : ntree.links) {
    if (link.fromnode != &node || link.tonode == nullptr) {
      continue;
    }
    if (bNode *bsdf = bsdf_for_node_recursive(ntree, *link.tonode, visited)) {
      return bsdf;
    }
  }
  return nullptr;
}

bNode *bsdf_for(const bNodeTree &ntree, bNode &stack_node)
{
  Set<const bNode *> visited;
  return bsdf_for_node_recursive(ntree, stack_node, visited);
}

std::optional<StackLayer> parent_group_layer(bNodeTree &ntree, bNode &stack_node)
{
  for (bNode &node : ntree.nodes) {
    if (&node == &stack_node || !is_stack(node)) {
      continue;
    }
    const NodeShaderLayerStack &storage = layer_stack::storage(node);
    for (const int i : IndexRange(storage.items_num)) {
      const StackLayer layer{&ntree, &node, i};
      if (layer.nested_stack() == &stack_node) {
        return layer;
      }
    }
  }
  return std::nullopt;
}

/* The layer of #stack (a Texture Layer Stack or Mask Stack node) that #socket
 * belongs to, or -1. Any of a texture layer's per-item sockets identify it: a
 * node feeding its Opacity or Mask is part of the layer as much as its
 * generator is. */
static int layer_index_for_socket(bNode &stack, const bNodeSocket &socket)
{
  const NodeShaderLayerStack &storage = layer_stack::storage(stack);
  const StringRef identifier = socket.identifier;
  for (const int i : IndexRange(storage.items_num)) {
    const NodeShaderLayerStackItem &item = storage.items[i];
    if (is_stack(stack)) {
      if (ELEM(identifier,
               StringRef(ItemsAccessor::socket_identifier_for_item(item)),
               StringRef(ItemsAccessor::opacity_socket_identifier_for_item(item)),
               StringRef(ItemsAccessor::mask_socket_identifier_for_item(item))))
      {
        return i;
      }
    }
    else if (identifier == StringRef(mask_stack::ItemsAccessor::socket_identifier_for_item(item)))
    {
      return i;
    }
  }
  return -1;
}

/* Follow links downstream from #node, appending every layer whose item socket
 * is reached to #r_layers. The walk never passes through a stack node: what
 * lies beyond it belongs to the layers of that stack. */
static void collect_layers_for_node(bNodeTree &ntree,
                                    bNode &node,
                                    Set<const bNode *> &visited,
                                    Vector<StackLayer> &r_layers)
{
  if (!visited.add(&node)) {
    return;
  }
  for (const bNodeLink &link : ntree.links) {
    if (link.fromnode != &node || link.tonode == nullptr) {
      continue;
    }
    bNode &consumer = *link.tonode;
    if (is_stack(consumer) || mask_stack::is_stack(consumer)) {
      const int index = layer_index_for_socket(consumer, *link.tosock);
      if (index != -1) {
        r_layers.append({&ntree, &consumer, index});
      }
      continue;
    }
    collect_layers_for_node(ntree, consumer, visited, r_layers);
  }
}

std::optional<StackLayer> layer_for_node(bNodeTree &ntree, bNode &node)
{
  /* A stack node stands for its own active layer: a group's nested stack means
   * its children, not the group layer it is the source of. */
  if (is_stack(node) || mask_stack::is_stack(node)) {
    const int index = layer_stack::active_index_in_range(node);
    return (index == -1) ? std::nullopt : std::optional<StackLayer>({&ntree, &node, index});
  }
  Set<const bNode *> visited;
  Vector<StackLayer> layers;
  collect_layers_for_node(ntree, node, visited, layers);
  if (layers.is_empty()) {
    return std::nullopt;
  }
  /* A node shared by several layers (one generator feeding two of them) resolves
   * to the layer its stack has active, i.e. the last explicit pick. */
  for (const StackLayer &layer : layers) {
    if (layer_stack::storage(*layer.stack).active_index == layer.index) {
      return layer;
    }
  }
  return layers[0];
}

bool is_adjustment_group(bNodeTree &ntree, bNode &stack_node)
{
  const NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
  if (storage.items_num == 0) {
    return false;
  }
  bNode *base_src = StackLayer{&ntree, &stack_node, storage.items_num - 1}.source();
  return base_src != nullptr && base_src->is_type("NodeClosureInput"_ustr);
}

static bool contains_impl(bNodeTree &ntree,
                          bNode &outer,
                          const bNode &inner,
                          Set<const bNode *> &visited)
{
  if (&outer == &inner) {
    return true;
  }
  /* Blender keeps cycle-forming stack links in the tree (marked invalid), so guard against
   * unbounded recursion. */
  if (!visited.add(&outer)) {
    return false;
  }
  const NodeShaderLayerStack &storage = layer_stack::storage(outer);
  for (const int i : IndexRange(storage.items_num)) {
    bNode *src = StackLayer{&ntree, &outer, i}.source();
    if (src && is_stack(*src) && contains_impl(ntree, *src, inner, visited)) {
      return true;
    }
  }
  return false;
}

bool contains(bNodeTree &ntree, bNode &outer, const bNode &inner)
{
  Set<const bNode *> visited;
  return contains_impl(ntree, outer, inner, visited);
}

/* Output that accompanies another output rather than carrying the node's own
 * value: an Image Texture's Alpha next to its Color. Driving a channel from it
 * is rarely what the user wants, so it is only used as a last resort. */
static bool is_companion_output(const bNodeSocket &socket)
{
  return StringRef(socket.name) == "Alpha";
}

/* First output socket of the given type on #node, or null. */
static bNodeSocket *first_output_of_type(bNode &node,
                                         const eNodeSocketDatatype type,
                                         const bool skip_companion)
{
  for (bNodeSocket &sock : node.outputs) {
    if (sock.type == type && !(skip_companion && is_companion_output(sock))) {
      return &sock;
    }
  }
  return nullptr;
}

bNodeSocket *preferred_source_output(bNode &node, const eNodeSocketDatatype channel_type)
{
  for (const bool skip_companion : {true, false}) {
    if (bNodeSocket *sock = first_output_of_type(node, channel_type, skip_companion)) {
      return sock;
    }
    for (const eNodeSocketDatatype type : {SOCK_RGBA, SOCK_FLOAT, SOCK_VECTOR}) {
      if (type == channel_type) {
        continue;
      }
      if (bNodeSocket *sock = first_output_of_type(node, type, skip_companion)) {
        return sock;
      }
    }
  }
  return nullptr;
}

/* Nodes that consume any output of #node. */
static Set<bNode *> node_output_consumers(bNodeTree &ntree, const bNode &node)
{
  Set<bNode *> users;
  for (const bNodeLink &link : ntree.links) {
    if (link.fromnode == &node && link.tonode) {
      users.add(link.tonode);
    }
  }
  return users;
}

/* Walk the nodes making up a layer, starting at its source node and appending
 * each to #r_nodes, which doubles as the visited set: a nested stack's
 * per-layer sources, a closure zone's body and paired Closure Input, the
 * generators feeding a Combine Bundle, and so on.
 *
 * With #owned_only the walk stops at nodes something outside it still consumes,
 * so it collects only what the layer owns (what deleting the layer may remove).
 * Without it the whole subtree feeding the layer is collected, including nodes
 * shared with other layers. */
static void collect_source_nodes(bNodeTree &ntree,
                                 bNode *source,
                                 const bool owned_only,
                                 VectorSet<bNode *> &r_nodes)
{
  if (source == nullptr || r_nodes.contains(source)) {
    return;
  }
  if (is_stack(*source)) {
    r_nodes.add(source);
    const NodeShaderLayerStack &storage = layer_stack::storage(*source);
    for (const int i : IndexRange(storage.items_num)) {
      const StackLayer layer{&ntree, source, i};
      bNodeSocket *bundle_sock = layer.bundle_socket();
      collect_source_nodes(ntree,
                           bundle_sock ? bke::node_find_source_node(ntree, *bundle_sock) : nullptr,
                           owned_only,
                           r_nodes);
      bNodeSocket *mask_sock = layer.mask_socket();
      collect_source_nodes(ntree,
                           mask_sock ? bke::node_find_source_node(ntree, *mask_sock) : nullptr,
                           owned_only,
                           r_nodes);
    }
    return;
  }
  if (mask_stack::is_stack(*source)) {
    r_nodes.add(source);
    const NodeShaderLayerStack &storage = layer_stack::storage(*source);
    for (const int i : IndexRange(storage.items_num)) {
      collect_source_nodes(
          ntree, mask_stack::layer_source(ntree, *source, i), owned_only, r_nodes);
    }
    return;
  }
  if (closure_zone::is_output_node(*source)) {
    /* Keep the zone while something outside the walk still uses its closure. A link to a node
     * already collected (e.g. the nested stack being removed around this layer) does not count,
     * otherwise the zone would leak. */
    if (owned_only) {
      for (const bNodeLink &link : ntree.links) {
        if (link.fromnode == source && link.tonode && !r_nodes.contains(link.tonode)) {
          return;
        }
      }
    }
    /* A closure layer owns its zone: the body nodes feeding the zone output,
     * the paired Closure Input, and the output node itself. */
    r_nodes.add(source);
    for (bNodeSocket &sock : source->inputs) {
      collect_source_nodes(ntree, bke::node_find_source_node(ntree, sock), owned_only, r_nodes);
    }
    if (bNode *zone_input = closure_zone::input_node(ntree, *source)) {
      r_nodes.add(zone_input);
    }
    return;
  }
  /* Generic node: keep it if a consumer outside the walk still uses it. Counting all consumers
   * (> 1) is unreliable: at the top level the removed layer's own link is already gone, and in the
   * nested case the sole consumer is the stack being deleted. */
  if (owned_only) {
    for (bNode *consumer : node_output_consumers(ntree, *source)) {
      if (!r_nodes.contains(consumer)) {
        return;
      }
    }
  }
  r_nodes.add(source);
  /* The nodes feeding it belong to the layer too: a generator behind the layer's Combine Bundle,
   * or the asset group behind a mask's Separate Bundle. */
  for (bNodeSocket &sock : source->inputs) {
    collect_source_nodes(ntree, bke::node_find_source_node(ntree, sock), owned_only, r_nodes);
  }
}

void nodes_for_layer(const StackLayer &layer, VectorSet<bNode *> &r_nodes)
{
  bNodeTree &ntree = *layer.ntree;
  collect_source_nodes(ntree, layer.source(), /*owned_only=*/false, r_nodes);
  if (bNodeSocket *mask_sock = layer.mask_socket()) {
    collect_source_nodes(
        ntree, bke::node_find_source_node(ntree, *mask_sock), /*owned_only=*/false, r_nodes);
  }
}

/* Vertical span the nodes of #layer occupy, or false when it has none yet. */
static bool layer_node_span(const StackLayer &layer, float &r_top, float &r_bottom)
{
  VectorSet<bNode *> nodes;
  nodes_for_layer(layer, nodes);
  if (nodes.is_empty()) {
    return false;
  }
  r_top = -FLT_MAX;
  r_bottom = FLT_MAX;
  for (const bNode *node : nodes) {
    r_top = std::max(r_top, node->location[1]);
    r_bottom = std::min(r_bottom, node->location[1] - estimated_node_height(*node));
  }
  return true;
}

static void move_layer_nodes(const StackLayer &layer, const float dy)
{
  VectorSet<bNode *> nodes;
  nodes_for_layer(layer, nodes);
  for (bNode *node : nodes) {
    node->location[1] += dy;
  }
}

void place_layer_source(bNode &node, const StackLayer &layer, const float dx)
{
  const float margin = 40.0f;
  bNodeTree &ntree = *layer.ntree;
  bNode &stack = *layer.stack;
  const int layers_num = layer_stack::storage(stack).items_num;
  const float width = node.typeinfo->default_width;
  const float x = stack.location[0] - dx;

  /* Below the nodes of the layers above this one, so the column reads in the
   * same order as the stack's item sockets and the layers tree view. */
  float top = stack.location[1];
  for (const int i : IndexRange(std::min(layer.index, layers_num))) {
    float other_top, other_bottom;
    if (layer_node_span({&ntree, &stack, i}, other_top, other_bottom)) {
      top = std::min(top, other_bottom - margin);
    }
  }

  /* Nodes sharing the column push the slot down as well, so nothing lands on
   * top of them — including this layer's own nodes, so a second node for it
   * (a mask stack next to the generator) sits below the first. Only the layers
   * below are exempt: they make room further down instead. */
  for (bNode &other : ntree.nodes) {
    if (&other == &node || &other == &stack) {
      continue;
    }
    const float other_width = (other.width > 0.0f) ? other.width : other.typeinfo->default_width;
    if (other.location[0] + other_width < x - margin || other.location[0] > x + width + margin) {
      continue;
    }
    if (other.location[1] > stack.location[1] + margin) {
      continue;
    }
    const std::optional<StackLayer> other_layer = layer_for_node(ntree, other);
    if (other_layer && other_layer->stack == &stack && other_layer->index > layer.index) {
      continue;
    }
    top = std::min(top, other.location[1] - estimated_node_height(other) - margin);
  }

  node.location[0] = x;
  node.location[1] = top;

  /* Push the layers below down by however much the new node overlaps them,
   * moving each layer's nodes as a whole so their arrangement is kept. */
  const float bottom = top - estimated_node_height(node) - margin;
  const IndexRange layers_below = IndexRange(layers_num).drop_front(layer.index + 1);
  float shift = 0.0f;
  for (const int i : layers_below) {
    float other_top, other_bottom;
    if (layer_node_span({&ntree, &stack, i}, other_top, other_bottom)) {
      shift = std::max(shift, other_top - bottom);
    }
  }
  if (shift > 0.0f) {
    for (const int i : layers_below) {
      move_layer_nodes({&ntree, &stack, i}, -shift);
    }
  }
}

void delete_layer_source_recursive(Main &bmain,
                                   bNodeTree &ntree,
                                   bNode *source,
                                   VectorSet<bNode *> &deleted)
{
  const int already_deleted = deleted.size();
  collect_source_nodes(ntree, source, /*owned_only=*/true, deleted);
  for (bNode *node : deleted.as_span().drop_front(already_deleted)) {
    bke::node_remove_node(&bmain, ntree, *node, true);
  }
}

static constexpr LayerSocketRole all_layer_socket_roles[] = {
    LayerSocketRole::Bundle, LayerSocketRole::Opacity, LayerSocketRole::Mask};

LayerState capture_layer_state(const StackLayer &layer)
{
  const NodeShaderLayerStackItem &item = layer.item();

  LayerState state;
  state.name = item.name ? item.name : "";
  state.blend_type = item.blend_type;
  state.flag = item.flag;
  state.item_type = item.item_type;
  for (const int i : IndexRange(item.disabled_channels_num)) {
    if (item.disabled_channels[i].name) {
      state.disabled_channels.append(item.disabled_channels[i].name);
    }
  }

  for (const LayerSocketRole role : all_layer_socket_roles) {
    bNodeSocket *sock = layer.socket(role);
    if (sock == nullptr) {
      continue;
    }
    LayerState::SocketState &socket_state = state.sockets[int(role)];
    if (const bNodeLink *link = bke::node_find_incoming_link(*layer.ntree, *sock)) {
      socket_state.linked = true;
      socket_state.from_node_id = link->fromnode->identifier;
      socket_state.from_socket_id = link->fromsock->identifier;
    }
    else if (role != LayerSocketRole::Bundle) {
      /* Opacity/Mask are plain floats: keep the user's value. */
      socket_state.has_default = true;
      socket_state.default_value = sock->default_value_typed<bNodeSocketValueFloat>()->value;
    }
  }
  return state;
}

void restore_layer_state(const StackLayer &layer, const LayerState &state)
{
  bNodeTree &ntree = *layer.ntree;
  NodeShaderLayerStackItem &item = layer.item();
  item.blend_type = state.blend_type;
  item.flag = state.flag;
  /* Fresh items are created as bundle layers. Restore the captured type and rebuild the node's
   * sockets to match before reconnecting, so a closure (adjustment) layer's source reconnects to a
   * closure input rather than a type-mismatched bundle input. */
  if (item.item_type != state.item_type) {
    item.item_type = state.item_type;
    update_node_declaration_and_sockets(ntree, *layer.stack);
  }
  for (const std::string &channel : state.disabled_channels) {
    layer_stack::set_channel_disabled(item, channel, true);
  }

  for (const LayerSocketRole role : all_layer_socket_roles) {
    bNodeSocket *sock = layer.socket(role);
    if (sock == nullptr) {
      continue;
    }
    const LayerState::SocketState &socket_state = state.sockets[int(role)];
    if (socket_state.linked) {
      bNode *from_node = bke::node_find_node_by_identifier(ntree, socket_state.from_node_id);
      if (from_node == nullptr) {
        continue;
      }
      bNodeSocket *from_sock = bke::node_find_socket(
          *from_node, SOCK_OUT, UString(socket_state.from_socket_id));
      if (from_sock) {
        bke::node_add_link(ntree, *from_node, *from_sock, *layer.stack, *sock);
      }
    }
    else if (socket_state.has_default) {
      sock->default_value_typed<bNodeSocketValueFloat>()->value = socket_state.default_value;
    }
  }
}

void connect_bundle(const StackLayer &layer, bNode &source, const StringRefNull preferred_output)
{
  bNodeSocket *layer_sock = layer.bundle_socket();
  if (layer_sock == nullptr) {
    return;
  }
  bNodeSocket *out = bke::node_find_socket(source, SOCK_OUT, UString(preferred_output));
  if (out == nullptr || out->type != SOCK_BUNDLE) {
    out = first_bundle_output(source);
  }
  if (out) {
    bke::node_add_link(*layer.ntree, source, *out, *layer.stack, *layer_sock);
  }
}

int active_layer_index(bNodeTree &ntree, bNode &stack_node)
{
  if (bNode *active = bke::node_get_active(ntree)) {
    if (const std::optional<StackLayer> layer = layer_for_node(ntree, *active)) {
      if (layer->stack == &stack_node) {
        return layer->index;
      }
    }
  }
  return layer_stack::active_index_in_range(stack_node);
}

int add_layer(Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name)
{
  /* Insert above the active layer, which the active node identifies. */
  layer_stack::storage(stack_node).active_index = std::max(active_layer_index(ntree, stack_node),
                                                           0);
  return layer_stack::add_layer_above_active(bmain, ntree, stack_node, name);
}

bNode *create_root(bNodeTree &ntree, const bNode *bsdf)
{
  bNode *stack = bke::node_add_node(nullptr, ntree, "ShaderNodeTextureLayerStack"_ustr);
  if (stack == nullptr) {
    return nullptr;
  }
  if (bsdf) {
    place_node_left_of(*stack, *bsdf, 300.0f, 200.0f);
  }
  bke::node_set_active(ntree, *stack);
  return stack;
}

void wire_root_into_bsdf(Main &bmain,
                         bNodeTree &ntree,
                         bNode &stack,
                         const int layer_index,
                         bNode &bsdf,
                         bNodeSocket *bsdf_socket)
{
  if (bsdf_socket == nullptr) {
    bsdf_socket = texture_channel::first_fillable_input(bsdf);
  }
  if (bsdf_socket == nullptr) {
    return;
  }

  bNode *combine = bke::node_add_node(nullptr, ntree, "NodeCombineBundle"_ustr);
  bNode *separate = bke::node_add_node(nullptr, ntree, "NodeSeparateBundle"_ustr);
  if (combine == nullptr || separate == nullptr) {
    return;
  }
  /* The chain reads left to right: Combine Bundle, stack, Separate Bundle, BSDF.
   * Each node is placed left of the one it feeds, by more than its own width so
   * they do not overlap. */
  STRNCPY_UTF8(combine->label, "Fill");
  place_node_left_of(*separate, bsdf, 220.0f, 200.0f);
  const float stack_width = (stack.width > 0.0f) ? stack.width : stack.typeinfo->default_width;
  place_node_left_of(stack, *separate, stack_width + 60.0f, 0.0f);
  place_layer_source(*combine, StackLayer{&ntree, &stack, layer_index}, 280.0f);

  /* The Separate Bundle's items are the channels routed into the BSDF: the
   * authoritative bundle signature of the stack. */
  auto &separate_storage = *static_cast<NodeSeparateBundle *>(separate->storage);
  separate_storage.flag = NodeSeparateBundleFlag(separate_storage.flag |
                                                 NODE_SEPARATE_BUNDLE_FLAG_DEFINE_SIGNATURE);
  BKE_main_ensure_invariants(bmain, ntree.id);

  connect_bundle(StackLayer{&ntree, &stack, layer_index}, *combine, "Bundle");

  bNodeSocket *stack_result = bke::node_find_socket(stack, SOCK_OUT, "Result"_ustr);
  bNodeSocket *separate_in = first_bundle_input(*separate);
  if (stack_result && separate_in) {
    bke::node_add_link(ntree, stack, *stack_result, *separate, *separate_in);
  }
  texture_channel::link(ntree, *separate, bsdf, *bsdf_socket);
  BKE_main_ensure_invariants(bmain, ntree.id);
}

}  // namespace blender::nodes::texture_stack
