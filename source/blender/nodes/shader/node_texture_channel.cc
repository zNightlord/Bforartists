/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include <algorithm>

#include "BLI_listbase_iterator.hh"

#include "DNA_array_utils.hh"
#include "DNA_node_types.h"

#include "BKE_node.hh"
#include "BKE_node_tree_update.hh"

#include "NOD_geo_bundle.hh"
#include "NOD_layer_stack.hh"
#include "NOD_socket.hh"
#include "NOD_socket_items.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"

namespace blender::nodes::texture_channel {

bool is_bsdf(const bNode &node)
{
  /* Raw socket lists, so this is also usable while the tree is being edited. */
  bool has_shader_out = false;
  for (const bNodeSocket &sock : node.outputs) {
    if (sock.type == SOCK_SHADER) {
      has_shader_out = true;
      break;
    }
  }
  if (!has_shader_out) {
    return false;
  }
  for (const bNodeSocket &sock : node.inputs) {
    if (sock.type == SOCK_SHADER) {
      return false;
    }
  }
  return true;
}

void collect_bsdfs(bNodeTree &ntree, Vector<bNode *> &r_bsdfs)
{
  for (bNode &node : ntree.nodes) {
    if (is_bsdf(node)) {
      r_bsdfs.append(&node);
    }
  }
}

bool is_fillable_input(const bNodeSocket &socket)
{
  return ELEM(socket.type, SOCK_FLOAT, SOCK_VECTOR, SOCK_RGBA) && !(socket.flag & SOCK_UNAVAIL);
}

bNodeSocket *first_fillable_input(bNode &bsdf)
{
  for (bNodeSocket &sock : bsdf.inputs) {
    if (is_fillable_input(sock)) {
      return &sock;
    }
  }
  return nullptr;
}

bool is_separate_bundle(const bNode &node)
{
  return node.is_type("NodeSeparateBundle"_ustr);
}

/* The Texture Layer Stack driving #separate's bundle input, or null. */
static bNode *stack_for_separate_bundle(const bNodeTree &ntree, bNode &separate)
{
  bNodeSocket *bundle_in = first_bundle_input(separate);
  if (bundle_in == nullptr) {
    return nullptr;
  }
  bNode *source = bke::node_find_source_node(ntree, *bundle_in);
  if (source && texture_stack::is_stack(*source)) {
    return source;
  }
  return nullptr;
}

State state(const bNodeTree &ntree, const bNodeSocket &bsdf_input)
{
  const bNodeLink *link = bke::node_find_incoming_link(ntree, bsdf_input);
  if (link == nullptr) {
    return State::Free;
  }
  if (link->fromnode && is_separate_bundle(*link->fromnode) &&
      stack_for_separate_bundle(ntree, *link->fromnode))
  {
    return State::Textured;
  }
  return State::Linked;
}

bNode *separate_bundle_for_bsdf(const bNodeTree &ntree, const bNode &bsdf)
{
  /* A Separate Bundle directly feeding one of the BSDF's inputs wins. */
  for (const bNodeSocket &sock : bsdf.inputs) {
    const bNodeLink *link = bke::node_find_incoming_link(ntree, sock);
    bNode *source = link ? link->fromnode : nullptr;
    if (source && is_separate_bundle(*source) && stack_for_separate_bundle(ntree, *source)) {
      return source;
    }
  }
  /* No channel linked right now: reuse a stack-driven Separate Bundle that
   * feeds no other BSDF, so toggling the last channel off and on again does
   * not create a second stack. */
  for (const bNode &node : ntree.nodes) {
    bNode &candidate = const_cast<bNode &>(node);
    if (!is_separate_bundle(candidate) || !stack_for_separate_bundle(ntree, candidate)) {
      continue;
    }
    bool feeds_other_bsdf = false;
    for (const bNodeLink &link : ntree.links) {
      if (link.fromnode == &candidate && link.tonode && link.tonode != &bsdf &&
          is_bsdf(*link.tonode))
      {
        feeds_other_bsdf = true;
        break;
      }
    }
    if (!feeds_other_bsdf) {
      return &candidate;
    }
  }
  return nullptr;
}

eNodeSocketDatatype socket_type(const bNodeTree &ntree, bNode &stack_node, const StringRef channel)
{
  bNode *bsdf = texture_stack::bsdf_for(ntree, stack_node);
  if (bsdf == nullptr) {
    return SOCK_RGBA;
  }
  for (bNodeSocket &sock : bsdf->inputs) {
    if (channel == sock.name && ELEM(sock.type, SOCK_RGBA, SOCK_FLOAT, SOCK_VECTOR)) {
      return eNodeSocketDatatype(sock.type);
    }
  }
  return SOCK_RGBA;
}

bNode *base_layer_combine_bundle(const bNodeTree &ntree, bNode &stack)
{
  const NodeShaderLayerStack &storage = layer_stack::storage(stack);
  if (storage.items_num == 0) {
    return nullptr;
  }
  const std::string base_id = texture_stack::ItemsAccessor::socket_identifier_for_item(
      storage.items[storage.items_num - 1]);
  bNodeSocket *base_socket = bke::node_find_socket(stack, SOCK_IN, UString(base_id));
  if (base_socket == nullptr) {
    return nullptr;
  }
  bNode *source = bke::node_find_source_node(ntree, *base_socket);
  if (source && source->is_type("NodeCombineBundle"_ustr)) {
    return source;
  }
  return nullptr;
}

bNodeSocket *base_layer_channel_socket(const bNodeTree &ntree,
                                       bNode &stack,
                                       const StringRef channel)
{
  bNode *combine = base_layer_combine_bundle(ntree, stack);
  if (combine == nullptr) {
    return nullptr;
  }
  const auto &storage = *static_cast<const NodeCombineBundle *>(combine->storage);
  for (const int i : IndexRange(storage.items_num)) {
    const NodeCombineBundleItem &item = storage.items[i];
    if (item.name && channel == item.name) {
      const std::string id = CombineBundleItemsAccessor::socket_identifier_for_item(item);
      return bke::node_find_socket(*combine, SOCK_IN, UString(id));
    }
  }
  return nullptr;
}

/* Copy a plain value between two sockets of the same float/vector/color type
 * and tag #to for updates. Only the value is copied, not the UI metadata
 * (subtype, soft range), which must stay the target socket's own. Returns
 * false on type mismatch. */
static bool copy_socket_value(bNodeTree &ntree, const bNodeSocket &from, bNodeSocket &to)
{
  if (from.type != to.type) {
    return false;
  }
  switch (from.type) {
    case SOCK_FLOAT:
      to.default_value_typed<bNodeSocketValueFloat>()->value =
          from.default_value_typed<bNodeSocketValueFloat>()->value;
      break;
    case SOCK_VECTOR:
      std::copy_n(from.default_value_typed<bNodeSocketValueVector>()->value,
                  4,
                  to.default_value_typed<bNodeSocketValueVector>()->value);
      break;
    case SOCK_RGBA:
      std::copy_n(from.default_value_typed<bNodeSocketValueRGBA>()->value,
                  4,
                  to.default_value_typed<bNodeSocketValueRGBA>()->value);
      break;
    default:
      return false;
  }
  BKE_ntree_update_tag_socket_property(&ntree, &to);
  return true;
}

/* True when #combine already has an item for #channel. */
static bool combine_has_channel(const bNode &combine, const StringRef channel)
{
  const auto &storage = *static_cast<const NodeCombineBundle *>(combine.storage);
  for (const int i : IndexRange(storage.items_num)) {
    if (storage.items[i].name && channel == storage.items[i].name) {
      return true;
    }
  }
  return false;
}

/* Keep every non-base layer's Combine Bundle carrying the full channel set. When #channel joins
 * the stack, add it as a disabled channel on each non-base layer (so the layer does not start
 * affecting it, leaving the render unchanged) and, for Combine Bundle sources, add the matching
 * socket so the bundle stays the exact channel set the join expects. The base layer is seeded
 * separately. */
static void spread_channel_to_non_base_layers(bNodeTree &ntree,
                                              bNode &stack,
                                              const StringRef channel,
                                              const bNodeSocket &ui_source)
{
  const eNodeSocketDatatype type = eNodeSocketDatatype(ui_source.type);
  NodeShaderLayerStack &storage = layer_stack::storage(stack);
  const int base_index = storage.items_num - 1;
  for (const int i : IndexRange(storage.items_num)) {
    if (i == base_index) {
      continue;
    }
    layer_stack::set_channel_disabled(storage.items[i], channel, true);
    bNode *source = texture_stack::StackLayer{&ntree, &stack, i}.source();
    if (source == nullptr || !source->is_type("NodeCombineBundle"_ustr) ||
        combine_has_channel(*source, channel))
    {
      continue;
    }
    const std::string name(channel);
    if (NodeCombineBundleItem *new_item =
            socket_items::add_item_with_socket_type_and_name<CombineBundleItemsAccessor>(
                ntree, *source, type, name.c_str()))
    {
      combine_bundle_item_copy_socket_data(*new_item, ui_source);
      update_node_declaration_and_sockets(ntree, *source);
      BKE_ntree_update_tag_node_property(&ntree, source);
    }
  }
}

/* Symmetric to #spread_channel_to_non_base_layers: drop #channel from every layer when it leaves
 * the stack, so no layer's Combine Bundle keeps a stale channel the join no longer carries. */
static void remove_channel_from_all_layers(bNodeTree &ntree, bNode &stack, const StringRef channel)
{
  NodeShaderLayerStack &storage = layer_stack::storage(stack);
  for (const int i : IndexRange(storage.items_num)) {
    layer_stack::set_channel_disabled(storage.items[i], channel, false);
    bNode *source = texture_stack::StackLayer{&ntree, &stack, i}.source();
    if (source == nullptr || !source->is_type("NodeCombineBundle"_ustr)) {
      continue;
    }
    auto &combine_storage = *static_cast<NodeCombineBundle *>(source->storage);
    for (const int j : IndexRange(combine_storage.items_num)) {
      if (combine_storage.items[j].name && channel == combine_storage.items[j].name) {
        socket_items::SocketItemsRef ref = CombineBundleItemsAccessor::get_items_from_node(
            *source);
        dna::array::remove_index(ref.items,
                                 ref.items_num,
                                 ref.active_index,
                                 j,
                                 CombineBundleItemsAccessor::destruct_item);
        update_node_declaration_and_sockets(ntree, *source);
        BKE_ntree_update_tag_node_property(&ntree, source);
        break;
      }
    }
  }
}

void seed_base_value(bNodeTree &ntree,
                     bNode &stack,
                     const StringRef channel,
                     const bNodeSocket &bsdf_socket)
{
  bNode *combine = base_layer_combine_bundle(ntree, stack);
  if (combine == nullptr) {
    return;
  }
  bNodeSocket *channel_socket = base_layer_channel_socket(ntree, stack, channel);
  if (channel_socket == nullptr) {
    const std::string name(channel);
    NodeCombineBundleItem *new_item =
        socket_items::add_item_with_socket_type_and_name<CombineBundleItemsAccessor>(
            ntree, *combine, eNodeSocketDatatype(bsdf_socket.type), name.c_str());
    if (new_item == nullptr) {
      return;
    }
    /* Mirror the BSDF input's slider on the base channel's value. */
    combine_bundle_item_copy_socket_data(*new_item, bsdf_socket);
    update_node_declaration_and_sockets(ntree, *combine);
    BKE_ntree_update_tag_node_property(&ntree, combine);
    channel_socket = base_layer_channel_socket(ntree, stack, channel);
  }
  if (channel_socket && bke::node_find_incoming_link(ntree, *channel_socket) == nullptr) {
    copy_socket_value(ntree, bsdf_socket, *channel_socket);
  }

  /* Keep the other layers' Combine Bundles carrying the full channel set too, so none is left an
   * out-of-date subset once this channel joins the stack. */
  spread_channel_to_non_base_layers(ntree, stack, channel, bsdf_socket);
}

bool link(bNodeTree &ntree, bNode &separate, bNode &bsdf, bNodeSocket &bsdf_socket)
{
  bNode *stack = stack_for_separate_bundle(ntree, separate);
  if (stack == nullptr) {
    return false;
  }
  const StringRef channel = bsdf_socket.name;

  auto &storage = *static_cast<NodeSeparateBundle *>(separate.storage);
  const NodeSeparateBundleItem *item = nullptr;
  for (const int i : IndexRange(storage.items_num)) {
    if (storage.items[i].name && channel == storage.items[i].name) {
      item = &storage.items[i];
      break;
    }
  }
  if (item == nullptr) {
    const std::string name(channel);
    item = socket_items::add_item_with_socket_type_and_name<SeparateBundleItemsAccessor>(
        ntree, separate, eNodeSocketDatatype(bsdf_socket.type), name.c_str());
    if (item == nullptr) {
      return false;
    }
    update_node_declaration_and_sockets(ntree, separate);
    BKE_ntree_update_tag_node_property(&ntree, &separate);
  }
  const std::string out_id = SeparateBundleItemsAccessor::socket_identifier_for_item(*item);
  bNodeSocket *out = bke::node_find_socket(separate, SOCK_OUT, UString(out_id));
  if (out == nullptr) {
    return false;
  }
  if (const bNodeLink *existing = bke::node_find_incoming_link(ntree, bsdf_socket)) {
    bke::node_remove_link(&ntree, const_cast<bNodeLink &>(*existing));
  }
  bke::node_add_link(ntree, separate, *out, bsdf, bsdf_socket);
  seed_base_value(ntree, *stack, channel, bsdf_socket);
  return true;
}

bool unlink(bNodeTree &ntree, bNodeSocket &bsdf_socket)
{
  const bNodeLink *link = bke::node_find_incoming_link(ntree, bsdf_socket);
  if (link == nullptr || link->fromnode == nullptr || link->fromsock == nullptr) {
    return false;
  }
  bNode &separate = *link->fromnode;
  if (!is_separate_bundle(separate)) {
    return false;
  }
  bNode *stack = stack_for_separate_bundle(ntree, separate);
  if (stack == nullptr) {
    return false;
  }

  /* The channel is the Separate Bundle item feeding the socket. */
  const auto &storage = *static_cast<const NodeSeparateBundle *>(separate.storage);
  const NodeSeparateBundleItem *item =
      socket_items::find_item_by_identifier<SeparateBundleItemsAccessor>(
          separate, link->fromsock->identifier);
  const int item_index = item ? int(item - storage.items) : -1;
  const std::string channel = (item && item->name) ? item->name : bsdf_socket.name;

  /* Keep the render unchanged: the base layer's plain value becomes the
   * socket's value again. */
  if (bNodeSocket *base_socket = base_layer_channel_socket(ntree, *stack, channel)) {
    if (bke::node_find_incoming_link(ntree, *base_socket) == nullptr) {
      copy_socket_value(ntree, *base_socket, bsdf_socket);
    }
  }

  bNodeSocket *item_output = link->fromsock;
  bke::node_remove_link(&ntree, const_cast<bNodeLink &>(*link));

  /* Drop the item unless another socket still uses its output. */
  bool still_used = false;
  for (const bNodeLink &other : ntree.links) {
    if (other.fromsock == item_output) {
      still_used = true;
      break;
    }
  }
  if (!still_used && item_index != -1) {
    socket_items::SocketItemsRef ref = SeparateBundleItemsAccessor::get_items_from_node(separate);
    dna::array::remove_index(ref.items,
                             ref.items_num,
                             ref.active_index,
                             item_index,
                             SeparateBundleItemsAccessor::destruct_item);
    BKE_ntree_update_tag_node_property(&ntree, &separate);
    /* The channel has left the stack: drop it from every layer's Combine Bundle so none keeps a
     * stale channel no longer in the join's signature. */
    remove_channel_from_all_layers(ntree, *stack, channel);
  }
  return true;
}

void handle_manual_link(bNodeTree &ntree, bNode &separate, const bNodeLink &link)
{
  if (link.tonode == nullptr || link.tosock == nullptr || link.fromsock == nullptr) {
    return;
  }
  if (!is_bsdf(*link.tonode) || !is_fillable_input(*link.tosock)) {
    return;
  }
  bNode *stack = stack_for_separate_bundle(ntree, separate);
  if (stack == nullptr) {
    return;
  }
  const NodeSeparateBundleItem *item =
      socket_items::find_item_by_identifier<SeparateBundleItemsAccessor>(
          separate, link.fromsock->identifier);
  /* Only name-matched links enable a channel; anything else is intentional
   * cross-wiring that must not touch other nodes. */
  if (item && item->name && StringRef(link.tosock->name) == item->name) {
    seed_base_value(ntree, *stack, item->name, *link.tosock);
  }
}

}  // namespace blender::nodes::texture_channel
