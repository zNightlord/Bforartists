/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include "BLI_listbase_iterator.hh"

#include "DNA_node_types.h"

#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_tree_update.hh"

#include "NOD_geo_bundle.hh"
#include "NOD_mask_stack.hh"
#include "NOD_node_placement.hh"
#include "NOD_texture_stack.hh"

namespace blender::nodes::mask_stack {

bool is_stack(const bNode &node)
{
  return node.is_type("ShaderNodeMaskStack"_ustr);
}

bNodeSocket *layer_socket(bNode &mask_stack_node, const int layer_index)
{
  const NodeShaderLayerStack &data = layer_stack::storage(mask_stack_node);
  if (layer_index < 0 || layer_index >= data.items_num) {
    return nullptr;
  }
  const std::string id = ItemsAccessor::socket_identifier_for_item(data.items[layer_index]);
  return bke::node_find_socket(mask_stack_node, SOCK_IN, UString(id));
}

bNode *layer_source(const bNodeTree &ntree, bNode &mask_stack_node, const int layer_index)
{
  bNodeSocket *sock = layer_socket(mask_stack_node, layer_index);
  if (sock == nullptr) {
    return nullptr;
  }
  return bke::node_find_source_node(ntree, *sock);
}

/* Create the Mask Stack node for layer #layer_index of #stack_node, placed in
 * that layer's slot. */
static bNode *add_mask_stack_node(bNodeTree &ntree, bNode &stack_node, const int layer_index)
{
  bNode *mask_stack_node = bke::node_add_node(nullptr, ntree, "ShaderNodeMaskStack"_ustr);
  if (mask_stack_node) {
    texture_stack::place_layer_source(
        *mask_stack_node, {&ntree, &stack_node, layer_index}, 360.0f);
  }
  return mask_stack_node;
}

/* Wedge a new Mask Stack between #mask_sock and its current source, so the
 * existing source becomes the base (bottom) mask layer. Returns the new stack. */
static bNode *wedge_mask_stack_before_socket(Main &bmain,
                                             bNodeTree &ntree,
                                             bNode &stack_node,
                                             const int layer_index,
                                             const StringRef mask_socket_id,
                                             const bNodeLink &existing_link)
{
  bNode *existing_source = existing_link.fromnode;
  bNodeSocket *existing_socket = existing_link.fromsock;
  bNode *mask_stack_node = add_mask_stack_node(ntree, stack_node, layer_index);
  if (mask_stack_node == nullptr) {
    return nullptr;
  }
  NodeShaderLayerStackItem *base_item = socket_items::add_item_with_name<ItemsAccessor>(
      *mask_stack_node, "Mask");
  if (base_item == nullptr) {
    return nullptr;
  }
  BKE_ntree_update_tag_node_property(&ntree, mask_stack_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
  const std::string base_socket_id = ItemsAccessor::socket_identifier_for_item(*base_item);
  if (bNodeSocket *base_in = bke::node_find_socket(
          *mask_stack_node, SOCK_IN, UString(base_socket_id)))
  {
    bke::node_add_link(ntree, *existing_source, *existing_socket, *mask_stack_node, *base_in);
  }
  if (bNodeSocket *result_out = bke::node_find_socket(*mask_stack_node, SOCK_OUT, "Result"_ustr)) {
    /* Re-fetch the (possibly reallocated) Mask socket before relinking. */
    bNodeSocket *mask_sock = bke::node_find_socket(stack_node, SOCK_IN, UString(mask_socket_id));
    if (mask_sock) {
      if (const bNodeLink *link = bke::node_find_incoming_link(ntree, *mask_sock)) {
        bke::node_remove_link(&ntree, const_cast<bNodeLink &>(*link));
      }
      bke::node_add_link(ntree, *mask_stack_node, *result_out, stack_node, *mask_sock);
    }
  }
  return mask_stack_node;
}

bNode *ensure_for_layer(Main &bmain,
                        bNodeTree &ntree,
                        bNode &stack_node,
                        const NodeShaderLayerStackItem &item)
{
  const std::string mask_id = texture_stack::ItemsAccessor::mask_socket_identifier_for_item(item);
  bNodeSocket *mask_sock = bke::node_find_socket(stack_node, SOCK_IN, UString(mask_id));
  if (mask_sock == nullptr) {
    return nullptr;
  }

  const NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
  const int layer_index = &item - storage.items;

  if (const bNodeLink *existing_link = bke::node_find_incoming_link(ntree, *mask_sock)) {
    if (existing_link->fromnode && is_stack(*existing_link->fromnode)) {
      /* Already driven by a Mask Stack: reuse it. */
      return existing_link->fromnode;
    }
    if (existing_link->fromnode == nullptr || existing_link->fromsock == nullptr) {
      return nullptr;
    }
    return wedge_mask_stack_before_socket(
        bmain, ntree, stack_node, layer_index, mask_id, *existing_link);
  }

  /* Nothing connected yet: create an empty Mask Stack feeding the Mask socket. */
  bNode *mask_stack_node = add_mask_stack_node(ntree, stack_node, layer_index);
  if (mask_stack_node == nullptr) {
    return nullptr;
  }
  if (bNodeSocket *result_out = bke::node_find_socket(*mask_stack_node, SOCK_OUT, "Result"_ustr)) {
    bke::node_add_link(ntree, *mask_stack_node, *result_out, stack_node, *mask_sock);
  }
  return mask_stack_node;
}

int add_layer(Main &bmain,
              bNodeTree &ntree,
              bNode &mask_stack_node,
              const char *name,
              const std::optional<float> default_value)
{
  /* Insert above the active mask, which the active node identifies. */
  layer_stack::storage(mask_stack_node).active_index = std::max(
      texture_stack::active_layer_index(ntree, mask_stack_node), 0);
  const int new_index = layer_stack::add_layer_above_active(bmain, ntree, mask_stack_node, name);
  if (default_value.has_value()) {
    if (bNodeSocket *sock = layer_socket(mask_stack_node, new_index)) {
      sock->default_value_typed<bNodeSocketValueFloat>()->value = *default_value;
    }
  }
  return new_index;
}

int add_layer_from_source(Main &bmain,
                          bNodeTree &ntree,
                          bNode &mask_stack_node,
                          const char *name,
                          bNode &source_node,
                          bNodeSocket &source_socket)
{
  const int index = add_layer(bmain, ntree, mask_stack_node, name, std::nullopt);
  if (bNodeSocket *mask_in = layer_socket(mask_stack_node, index)) {
    /* Use the caller's node reference rather than source_socket.owner_node(): add_layer just
     * mutated the tree, so the topology cache owner_node() relies on may be dirty. */
    bke::node_add_link(ntree, source_node, source_socket, mask_stack_node, *mask_in);
  }
  return index;
}

void set_layer_blend(bNode &mask_stack_node, const int layer_index, const int blend_type)
{
  NodeShaderLayerStack &data = layer_stack::storage(mask_stack_node);
  if (blend_type == MA_RAMP_BLEND || layer_index < 0 || layer_index >= data.items_num - 1) {
    return;
  }
  data.items[layer_index].blend_type = blend_type;
}

/* The output of #node that best stands in for a mask when it has no float
 * output of its own: a color or vector output, implicitly converted (an Image
 * Texture's Color), and only then a companion Alpha. */
static bNodeSocket *direct_source_output(bNode &node, bNode **r_source_node)
{
  bNodeSocket *sock = texture_stack::preferred_source_output(node, SOCK_FLOAT);
  if (sock == nullptr) {
    return nullptr;
  }
  *r_source_node = &node;
  return sock;
}

bNodeSocket *find_source_output(Main &bmain,
                                bNodeTree &ntree,
                                bNode &group_node,
                                bNode **r_source_node)
{
  *r_source_node = nullptr;
  /* A float output is the mask value itself. An Alpha output is not: it
   * accompanies the node's Color rather than being its value, so it comes after
   * the channel and color outputs below (see #texture_stack::preferred_source_output). */
  for (bNodeSocket &sock : group_node.outputs) {
    if (sock.type == SOCK_FLOAT && !texture_stack::is_companion_output(sock)) {
      *r_source_node = &group_node;
      return &sock;
    }
  }
  /* No float output — extract a channel via Separate Bundle. */
  bNodeSocket *bundle_out = first_bundle_output(group_node);
  if (!bundle_out) {
    return direct_source_output(group_node, r_source_node);
  }
  bNode *separate = bke::node_add_node(nullptr, ntree, "NodeSeparateBundle"_ustr);
  if (!separate) {
    return direct_source_output(group_node, r_source_node);
  }
  separate->location[0] = group_node.location[0] + 200.0f;
  separate->location[1] = group_node.location[1];
  bNodeSocket *separate_in = first_bundle_input(*separate);
  if (!separate_in) {
    bke::node_remove_node(&bmain, ntree, *separate, true);
    return direct_source_output(group_node, r_source_node);
  }
  bke::node_add_link(ntree, group_node, *bundle_out, *separate, *separate_in);
  /* Triggering an update populates the Separate Bundle's output sockets from
   * the source bundle's items. */
  BKE_ntree_update_tag_node_property(&ntree, separate);
  BKE_main_ensure_invariants(bmain, ntree.id);
  /* Prefer "Mask", "Alpha", or any FLOAT output. */
  bNodeSocket *fallback = nullptr;
  for (bNodeSocket &sock : separate->outputs) {
    if (sock.type != SOCK_FLOAT) {
      continue;
    }
    if (sock.identifier == StringRef("Mask") || sock.identifier == StringRef("Alpha")) {
      *r_source_node = separate;
      return &sock;
    }
    if (!fallback) {
      fallback = &sock;
    }
  }
  if (!fallback) {
    /* No float channel to extract a mask from: don't leave the Separate Bundle orphaned. */
    bke::node_remove_node(&bmain, ntree, *separate, true);
    return direct_source_output(group_node, r_source_node);
  }
  *r_source_node = separate;
  return fallback;
}

}  // namespace blender::nodes::mask_stack
