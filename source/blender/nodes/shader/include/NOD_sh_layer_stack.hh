/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * Shared storage access, socket-item accessors and node callbacks for the
 * Texture Layer Stack and Mask Stack shader nodes. Both nodes store a
 * #NodeShaderLayerStack and differ only in the socket types they declare
 * and in how their layers are composited during shader node inlining.
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

namespace blender::nodes {

/** Typed access to the storage shared by the layer stack shader nodes. */
inline NodeShaderLayerStack &layer_stack_storage(bNode &node)
{
  return *static_cast<NodeShaderLayerStack *>(node.storage);
}
inline const NodeShaderLayerStack &layer_stack_storage(const bNode &node)
{
  return *static_cast<const NodeShaderLayerStack *>(node.storage);
}

/** True when #channel is on #item's disabled channels list, i.e. the
 * layer does not contribute to it during stack composition. */
inline bool layer_stack_channel_disabled(const NodeShaderLayerStackItem &item,
                                         const StringRef channel)
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
inline void layer_stack_channel_set_disabled(NodeShaderLayerStackItem &item,
                                             const StringRef channel,
                                             const bool disabled)
{
  if (disabled == layer_stack_channel_disabled(item, channel)) {
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
struct LayerStackItemsAccessorBase : public socket_items::SocketItemsAccessorDefaults {
  using ItemT = NodeShaderLayerStackItem;
  static constexpr bool has_type = false;
  static constexpr bool has_name = true;

  static socket_items::SocketItemsRef<NodeShaderLayerStackItem> get_items_from_node(bNode &node)
  {
    NodeShaderLayerStack &storage = layer_stack_storage(node);
    return {&storage.items, &storage.items_num, &storage.active_index};
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
    NodeShaderLayerStack &storage = layer_stack_storage(node);
    item.identifier = storage.next_identifier++;
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

struct TextureLayerStackItemsAccessor
    : public LayerStackItemsAccessorBase<TextureLayerStackItemsAccessor> {
  static StructRNA **item_srna;
  static constexpr StringRefNull node_idname = "ShaderNodeTextureLayerStack";
  static constexpr const char *socket_prefix = "Layer_";
  static constexpr const char *default_item_name = "Layer";
  static constexpr bool has_default_item_name = true;
  /** Layers are Bundle-typed (generators) or Closure-typed (adjustments); the
   * type is set from the linked socket when extending the stack. */
  static constexpr bool has_type = true;
  struct operator_idnames {
    static constexpr StringRefNull add_item = "NODE_OT_texture_layer_stack_item_add";
    static constexpr StringRefNull remove_item = "NODE_OT_texture_layer_stack_item_remove";
    static constexpr StringRefNull move_item = "NODE_OT_texture_layer_stack_item_move";
  };
  struct ui_idnames {
    static constexpr StringRefNull list = "NODE_UL_texture_layer_stack_items";
  };
  struct rna_names {
    static constexpr StringRefNull items = "layer_items";
    static constexpr StringRefNull active_index = "active_index";
  };

  static std::string mask_socket_identifier_for_item(const NodeShaderLayerStackItem &item)
  {
    return "Mask_" + std::to_string(item.identifier);
  }

  static eNodeSocketDatatype get_socket_type(const NodeShaderLayerStackItem &item)
  {
    return item.item_type == SHADER_LAYER_STACK_ITEM_CLOSURE ? SOCK_CLOSURE : SOCK_BUNDLE;
  }

  static bool supports_socket_type(const eNodeSocketDatatype socket_type, const int /*ntree_type*/)
  {
    return ELEM(socket_type, SOCK_BUNDLE, SOCK_CLOSURE);
  }

  static void init_with_socket_type_and_name(bNode &node,
                                             NodeShaderLayerStackItem &item,
                                             const eNodeSocketDatatype socket_type,
                                             const char *name)
  {
    init_with_name(node, item, name);
    item.item_type = (socket_type == SOCK_CLOSURE) ? SHADER_LAYER_STACK_ITEM_CLOSURE :
                                                     SHADER_LAYER_STACK_ITEM_BUNDLE;
  }
};

struct MaskStackItemsAccessor : public LayerStackItemsAccessorBase<MaskStackItemsAccessor> {
  static StructRNA **item_srna;
  static constexpr StringRefNull node_idname = "ShaderNodeMaskStack";
  static constexpr const char *socket_prefix = "Mask_";
  static constexpr const char *default_item_name = "Mask";
  static constexpr bool has_default_item_name = true;
  struct operator_idnames {
    static constexpr StringRefNull add_item = "NODE_OT_mask_stack_item_add";
    static constexpr StringRefNull remove_item = "NODE_OT_mask_stack_item_remove";
    static constexpr StringRefNull move_item = "NODE_OT_mask_stack_item_move";
  };
  struct ui_idnames {
    static constexpr StringRefNull list = "NODE_UL_mask_stack_items";
  };
  struct rna_names {
    static constexpr StringRefNull items = "mask_items";
    static constexpr StringRefNull active_index = "active_index";
  };
};

/** Node-type callbacks shared by both layer stack nodes. */
namespace layer_stack {

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
  const NodeShaderLayerStack &src_storage = layer_stack_storage(*src_node);
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
    layer_stack_storage(params.node).active_index = layer_stack_storage(params.node).items_num - 1;
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

}  // namespace layer_stack
}  // namespace blender::nodes
