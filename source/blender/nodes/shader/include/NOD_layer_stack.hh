/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * Storage shared by the Texture Layer Stack and Mask Stack shader nodes. Both
 * nodes store a #NodeShaderLayerStack and differ only in the socket types they
 * declare and in how their layers are composited during shader node inlining.
 * This header provides typed storage access, the common part of their
 * socket-item accessors, the node-type callbacks they share, and layer-array
 * operations that work on either node kind.
 *
 * The per-node-kind APIs build on this: #NOD_texture_stack.hh and
 * #NOD_mask_stack.hh.
 */

#include <algorithm>

#include "MEM_guardedalloc.h"

#include "BLI_index_range.hh"
#include "BLI_string.hh"
#include "BLI_string_ref.hh"

#include "DNA_material_types.h"
#include "DNA_node_types.h"

#include "NOD_socket_items.hh"
#include "NOD_socket_items_blend.hh"

struct Main;

namespace blender::nodes::layer_stack {

/** Typed access to the storage shared by the layer stack shader nodes. */
inline NodeShaderLayerStack &storage(bNode &node)
{
  return *static_cast<NodeShaderLayerStack *>(node.storage);
}
inline const NodeShaderLayerStack &storage(const bNode &node)
{
  return *static_cast<const NodeShaderLayerStack *>(node.storage);
}

/** True when #channel is on #item's disabled channels list, i.e. the
 * layer does not contribute to it during stack composition. */
inline bool channel_disabled(const NodeShaderLayerStackItem &item, const StringRef channel)
{
  for (const int i : IndexRange(item.disabled_channels_num)) {
    const char *name = item.disabled_channels[i].name;
    if (name && channel == name) {
      return true;
    }
  }
  return false;
}

/** Add or remove #channel from #item's disabled channels list. The caller
 * is responsible for tagging the node tree for updates. */
inline void set_channel_disabled(NodeShaderLayerStackItem &item,
                                 const StringRef channel,
                                 const bool disabled)
{
  if (disabled == channel_disabled(item, channel)) {
    return;
  }
  if (disabled) {
    NodeShaderLayerStackChannel *channels = MEM_new_array_zeroed<NodeShaderLayerStackChannel>(
        item.disabled_channels_num + 1, __func__);
    if (item.disabled_channels) {
      std::copy_n(item.disabled_channels, item.disabled_channels_num, channels);
      MEM_delete(item.disabled_channels);
    }
    channels[item.disabled_channels_num].name = BLI_strdupn(channel.data(), channel.size());
    item.disabled_channels = channels;
    item.disabled_channels_num++;
    return;
  }
  int kept = 0;
  for (const int i : IndexRange(item.disabled_channels_num)) {
    NodeShaderLayerStackChannel &entry = item.disabled_channels[i];
    if (entry.name && channel == entry.name) {
      MEM_delete(entry.name);
      continue;
    }
    item.disabled_channels[kept++] = entry;
  }
  item.disabled_channels_num = kept;
  if (kept == 0) {
    MEM_SAFE_DELETE(item.disabled_channels);
  }
}

/**
 * Common part of the item accessors for both layer stack nodes. The derived
 * accessor provides what differs per node: #node_idname, #operator_idnames,
 * #rna_names, #item_srna, #socket_prefix (per-item socket identifier prefix)
 * and #default_item_name (name for items created via the extend socket).
 */
template<typename Derived>
struct ItemsAccessorBase : public socket_items::SocketItemsAccessorDefaults {
  using ItemT = NodeShaderLayerStackItem;
  static constexpr bool has_type = false;
  static constexpr bool has_name = true;

  static socket_items::SocketItemsRef<NodeShaderLayerStackItem> get_items_from_node(bNode &node)
  {
    NodeShaderLayerStack &data = storage(node);
    return {&data.items, &data.items_num, &data.active_index};
  }

  static void copy_item(const NodeShaderLayerStackItem &src, NodeShaderLayerStackItem &dst)
  {
    dst = src;
    dst.name = BLI_strdup_null(dst.name);
    if (dst.disabled_channels) {
      NodeShaderLayerStackChannel *channels = MEM_new_array_zeroed<NodeShaderLayerStackChannel>(
          dst.disabled_channels_num, __func__);
      for (const int i : IndexRange(dst.disabled_channels_num)) {
        channels[i].name = BLI_strdup_null(dst.disabled_channels[i].name);
      }
      dst.disabled_channels = channels;
    }
  }

  static void destruct_item(NodeShaderLayerStackItem *item)
  {
    MEM_SAFE_DELETE(item->name);
    for (const int i : IndexRange(item->disabled_channels_num)) {
      MEM_SAFE_DELETE(item->disabled_channels[i].name);
    }
    MEM_SAFE_DELETE(item->disabled_channels);
  }

  static void blend_write_item(BlendWriter *writer, const NodeShaderLayerStackItem &item)
  {
    writer->write_string(item.name);
    writer->write_struct_array_by_id(dna::sdna_struct_id_get<NodeShaderLayerStackChannel>(),
                                     item.disabled_channels_num,
                                     item.disabled_channels);
    for (const int i : IndexRange(item.disabled_channels_num)) {
      writer->write_string(item.disabled_channels[i].name);
    }
  }

  static void blend_read_data_item(BlendDataReader *reader, NodeShaderLayerStackItem &item)
  {
    BLO_read_string(reader, &item.name);
    BLO_read_array_and_validate_size(reader, &item.disabled_channels, &item.disabled_channels_num);
    for (const int i : IndexRange(item.disabled_channels_num)) {
      BLO_read_string(reader, &item.disabled_channels[i].name);
    }
  }

  static char **get_name(NodeShaderLayerStackItem &item)
  {
    return &item.name;
  }

  static void init_with_name(bNode &node, NodeShaderLayerStackItem &item, const char *name)
  {
    NodeShaderLayerStack &data = storage(node);
    item.identifier = data.next_identifier++;
    item.blend_type = MA_RAMP_BLEND;
    item.item_type = SHADER_LAYER_STACK_ITEM_BUNDLE;
    socket_items::set_item_name_and_make_unique<Derived>(node, item, name);
  }

  static std::string socket_identifier_for_item(const NodeShaderLayerStackItem &item)
  {
    return Derived::socket_prefix + std::to_string(item.identifier);
  }

  static std::string opacity_socket_identifier_for_item(const NodeShaderLayerStackItem &item)
  {
    return "Opacity_" + std::to_string(item.identifier);
  }
};

/* -------------------------------------------------------------------- */
/* Node-type callbacks shared by both layer stack nodes. */

inline void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->storage = MEM_new<NodeShaderLayerStack>(__func__);
}

template<typename Accessor> void node_free_storage(bNode *node)
{
  socket_items::destruct_array<Accessor>(*node);
  MEM_delete(static_cast<NodeShaderLayerStack *>(node->storage));
}

template<typename Accessor>
void node_copy_storage(bNodeTree * /*dst_tree*/, bNode *dst_node, const bNode *src_node)
{
  const NodeShaderLayerStack &src_storage = storage(*src_node);
  dst_node->storage = MEM_new<NodeShaderLayerStack>(__func__, dna::shallow_copy(src_storage));
  socket_items::copy_array<Accessor>(*src_node, *dst_node);
}

template<typename Accessor>
void node_blend_write(const bNodeTree & /*tree*/, const bNode &node, BlendWriter &writer)
{
  socket_items::blend_write<Accessor>(&writer, node);
}

template<typename Accessor>
void node_blend_read(bNodeTree & /*tree*/, bNode &node, BlendDataReader &reader)
{
  socket_items::blend_read_data<Accessor>(&reader, node);
}

/**
 * Dragging a link onto the trailing extend socket appends a new item named
 * after #Accessor::default_item_name ("Layer", "Layer.001", ...) and makes it
 * active. The new item lives at the end, so it becomes the new base and the
 * previously-last item gains its blend-weight sockets.
 *
 * For accessors with typed items, a link dropped onto an untyped (UNSET)
 * item's hollow socket types that item from the linked socket instead of
 * creating a new one.
 */
template<typename Accessor> bool node_insert_link(bke::NodeInsertLinkParams &params)
{
  NodeShaderLayerStackItem *new_item = nullptr;
  const bool keep_link = socket_items::try_add_item_via_any_extend_socket<Accessor>(
      params.ntree, params.node, params.node, params.link, StringRef("__extend__"), &new_item);
  if (new_item) {
    socket_items::set_item_name_and_make_unique<Accessor>(
        params.node, *new_item, Accessor::default_item_name);
    storage(params.node).active_index = storage(params.node).items_num - 1;
    return true;
  }
  if constexpr (Accessor::has_type) {
    const bNodeSocket *target = (params.link.tonode == &params.node) ? params.link.tosock :
                                                                       nullptr;
    if (target && STREQ(target->idname, "NodeSocketVirtual") &&
        !STREQ(target->identifier, "__extend__"))
    {
      NodeShaderLayerStackItem *item = socket_items::find_item_by_identifier<Accessor>(
          params.node, target->identifier);
      if (item && item->item_type == SHADER_LAYER_STACK_ITEM_UNSET) {
        const eNodeSocketDatatype src_type = eNodeSocketDatatype(params.link.fromsock->type);
        if (!Accessor::supports_socket_type(src_type, params.ntree.type)) {
          return false;
        }
        item->item_type = (src_type == SOCK_CLOSURE) ? SHADER_LAYER_STACK_ITEM_CLOSURE :
                                                       SHADER_LAYER_STACK_ITEM_BUNDLE;
        update_node_declaration_and_sockets(params.ntree, params.node);
        const std::string identifier = Accessor::socket_identifier_for_item(*item);
        params.link.tosock = bke::node_find_socket(params.node, SOCK_IN, UString(identifier));
        BKE_ntree_update_tag_node_property(&params.ntree, &params.node);
        return params.link.tosock != nullptr;
      }
    }
  }
  return keep_link;
}

/* -------------------------------------------------------------------- */
/* Layer-array operations that work on either layer stack node kind. */

/** The active layer index when it is a valid item, or -1. */
int active_index_in_range(const bNode &stack_node);

/** Append a new layer named #name and move it to #target (clamped).
 * Returns the layer's final index. Does not change the active index.
 * #stack_node may be a Texture Layer Stack or a Mask Stack node. */
int add_layer_at(Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name, int target);

/** Add a new layer above the active one on either stack node kind and make it
 * active. Returns the new layer's index. */
int add_layer_above_active(Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name);

/** Remove layer #index from #stack_node, destructing it and updating the
 * tree. #stack_node may be either layer stack node kind. */
void remove_layer(Main &bmain, bNodeTree &ntree, bNode &stack_node, int index);

}  // namespace blender::nodes::layer_stack
