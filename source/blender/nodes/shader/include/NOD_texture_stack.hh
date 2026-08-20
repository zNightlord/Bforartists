/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * The Texture Layer Stack node and its layers: the single source of truth for
 * how the shader node graph is interpreted as a stack of texture layers, used
 * by the Shader Layers tree view, the layer operators and the asset add-menus,
 * plus the tree edits that maintain that interpretation.
 *
 * The canonical topology the feature builds and expects is (left to right,
 * as displayed in the node editor):
 * \code{.unparsed}
 *   generator / group ── Layer_<id>   ┐
 *                        Opacity_<id> ├ Texture Layer Stack ── Separate Bundle ── in ── BSDF
 *   Mask Stack ───────── Mask_<id>    ┘
 * \endcode
 * A "BSDF" is a node that outputs a shader and takes no shader as input (see
 * #texture_channel::is_bsdf). Each Texture Layer Stack item exposes a Bundle
 * input plus, for non-base items, Opacity and Mask inputs; the per-item socket
 * identifiers are derived from the item's stable id (see #ItemsAccessor).
 * When a layer's Bundle input is fed by another Texture Layer Stack node, that
 * layer is treated as a group and the nested stack's items become its children.
 *
 * All queries walk the raw node/link lists (following through reroutes, see
 * #bke::node_find_source_node), so they are safe to use while the tree is
 * being edited and do not require the topology cache.
 *
 * The related APIs: #NOD_mask_stack.hh (per-layer mask stacks),
 * #NOD_texture_channel.hh (routing channels into the BSDF) and
 * #NOD_closure_zone.hh (adjustment layers).
 */

#include <optional>
#include <string>

#include "BLI_set.hh"
#include "BLI_string_ref.hh"
#include "BLI_vector.hh"
#include "BLI_vector_set.hh"

#include "NOD_layer_stack.hh"

struct Main;

namespace blender::nodes::texture_stack {

struct ItemsAccessor : public layer_stack::ItemsAccessorBase<ItemsAccessor> {
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

/** Texture Layer Stack node. */
bool is_stack(const bNode &node);

/** The per-layer input sockets a layer has on its stack node. */
enum class LayerSocketRole {
  Bundle = 0,
  Opacity = 1,
  Mask = 2,
};

/**
 * Reference to one layer on a Texture Layer Stack node, or to a mask layer when
 * #stack is a Mask Stack node (the two node kinds share their storage and only
 * name their per-item socket differently). A thin handle (tree, stack node,
 * layer index) whose methods answer how the layer is wired; it holds no derived
 * state, so it stays valid across tree edits as long as the node and index do.
 */
struct StackLayer {
  bNodeTree *ntree;
  bNode *stack;
  int index;

  /** The layer's item. Asserts that the index is in range. */
  NodeShaderLayerStackItem &item() const;

  /** The input socket for #role, or null when the index is out of range or
   * the socket is missing (the base layer has no Opacity). */
  bNodeSocket *socket(LayerSocketRole role) const;
  bNodeSocket *bundle_socket() const
  {
    return this->socket(LayerSocketRole::Bundle);
  }
  bNodeSocket *opacity_socket() const
  {
    return this->socket(LayerSocketRole::Opacity);
  }
  bNodeSocket *mask_socket() const
  {
    return this->socket(LayerSocketRole::Mask);
  }

  /** The node feeding the layer's Bundle input, or null when unlinked. */
  bNode *source() const;

  /** The Mask Stack node feeding the layer's Mask input, or null. */
  bNode *mask_stack_node() const;

  /** When the layer is a group, its nested Texture Layer Stack; else null. A
   * generator group's Bundle input is fed by the nested stack directly; an
   * adjustment group's layer is a closure whose zone body is the nested
   * stack. */
  bNode *nested_stack() const;

  /** The node whose properties represent the layer in the UI: the body node
   * inside the closure zone for closure-typed layers, otherwise the direct
   * bundle source. */
  bNode *display_source() const;
};

/** Append the unique Texture Layer Stack roots feeding #bsdf via Separate
 * Bundle nodes to #r_roots. */
void roots_for_bsdf(const bNodeTree &ntree, bNode &bsdf, Vector<bNode *> &r_roots);

/** The stack node that layer edits should target: the active node itself when
 * it is a stack, otherwise the first stack feeding the active BSDF. Null when
 * there is no such stack. */
bNode *active(bNodeTree &ntree);

/** The BSDF ultimately consuming #stack_node's result (following links
 * downstream through Separate Bundle nodes, nested stacks, ...), or null when
 * the stack does not reach one. */
bNode *bsdf_for(const bNodeTree &ntree, bNode &stack_node);

/** The group layer that #stack_node is the nested stack of, or none when the
 * stack is not inside a group. */
std::optional<StackLayer> parent_group_layer(bNodeTree &ntree, bNode &stack_node);

/**
 * The layer #node belongs to, found by following links downstream until a
 * stack's per-item socket is reached; none when the node feeds no layer. This
 * is the inverse of #nodes_for_layer, and lets the active node be a layer's
 * source rather than the stack itself.
 *
 * A stack node resolves to its own active layer, so a group's nested stack
 * stands for its children rather than for the group layer above it. A node
 * feeding several layers resolves to the one its stack has active.
 */
std::optional<StackLayer> layer_for_node(bNodeTree &ntree, bNode &node);

/**
 * Place #node, the node about to feed #layer, in the column #dx to the left of
 * its stack node: below the nodes of the layers above it and above those of the
 * layers below, so the column reads in the same order as the stack's item
 * sockets and the layers tree view (and its links do not cross). The layers
 * below move down to make room; nodes that are not part of the stack stay put.
 */
void place_layer_source(bNode &node, const StackLayer &layer, float dx);

/**
 * The nodes making up a layer, appended to #r_nodes: its source and everything
 * that source pulls in (a nested stack's per-layer sources, a closure zone's
 * body, the generators feeding a Combine Bundle), plus its mask stack and mask
 * sources. Empty when the layer has no source. The stack node itself is not
 * included; the layer is one item of it.
 */
void nodes_for_layer(const StackLayer &layer, VectorSet<bNode *> &r_nodes);

/** Whether #stack_node is the nested stack of an adjustment group: its base
 * (last) layer is fed by the group's Closure Input (the incoming "below"
 * bundle) rather than by a real layer source. */
bool is_adjustment_group(bNodeTree &ntree, bNode &stack_node);

/** True if #inner lives anywhere in #outer's nested groups (cycle guard). */
bool contains(bNodeTree &ntree, bNode &outer, const bNode &inner);

/** Output socket of #node best matching a channel of type #channel_type:
 * an exact type match first, then any color/float/vector output. Outputs that
 * accompany another one (see #is_companion_output) are only used when nothing
 * else is available, so e.g. a float channel takes the image's Color. */
bNodeSocket *preferred_source_output(bNode &node, eNodeSocketDatatype channel_type);

/** True for an output that accompanies another output rather than carrying the
 * node's own value: an Image Texture's Alpha next to its Color. Driving a
 * channel or a mask from one is rarely what the user wants. */
bool is_companion_output(const bNodeSocket &socket);

/** Delete #source and any nodes it transitively owns: #nodes_for_layer's walk
 * restricted to nodes nothing outside the layer uses (a shared Combine Bundle,
 * Group, ... is kept). #deleted tracks already-removed nodes so repeated calls
 * do not recurse over them again. */
void delete_layer_source_recursive(Main &bmain,
                                   bNodeTree &ntree,
                                   bNode *source,
                                   VectorSet<bNode *> &deleted);

/** Snapshot of a layer's data + incoming links, so it can be recreated on
 * another stack (used by convert-to-group, ungroup and cross-stack reparent). */
struct LayerState {
  std::string name;
  int blend_type = 0;
  int flag = 0;
  /** #eShaderLayerStackItemType, so a closure (adjustment) layer keeps its type when recreated on
   * another stack, rather than defaulting to a bundle layer. */
  int item_type = 0;
  Vector<std::string> disabled_channels;
  struct SocketState {
    bool linked = false;
    /** #bNode::identifier of the linked source node. */
    int32_t from_node_id = 0;
    std::string from_socket_id;
    bool has_default = false;
    float default_value = 0.0f;
  };
  SocketState sockets[3];
};

LayerState capture_layer_state(const StackLayer &layer);
/** Apply a snapshot to a freshly created layer. */
void restore_layer_state(const StackLayer &layer, const LayerState &state);

/** Link #source's Bundle output (preferring #preferred_output by
 * identifier) into #layer's Bundle input. */
void connect_bundle(const StackLayer &layer, bNode &source, StringRefNull preferred_output);

/** The active layer of #stack_node: the layer the tree's active node belongs to
 * (see #layer_for_node), or the stack's stored active layer when the active
 * node is none of its nodes. -1 when the stack has no valid active layer. */
int active_layer_index(bNodeTree &ntree, bNode &stack_node);

/** Add a new layer to a Texture Layer Stack, inserted above the active layer,
 * and make it active. Returns the new layer's index. */
int add_layer(Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name);

/** Add a Texture Layer Stack node for a material that has none yet, placed
 * near #bsdf (when given), and made active so it resolves as the target of the
 * layer operators. */
bNode *create_root(bNodeTree &ntree, const bNode *bsdf);

/**
 * Give a freshly bootstrapped root stack's first layer real content and wire
 * it into #bsdf, so the stack is reachable from a BSDF and shows up in the
 * Shader Layers tree view (which only lists stacks connected to one). The
 * layer gets a single channel matching #bsdf_socket (by default the BSDF's
 * first fillable input, e.g. Base Color on a Principled BSDF), via a Combine
 * Bundle feeding the layer and a Separate Bundle feeding the BSDF:
 *
 *   Combine Bundle -> layer -> Stack -> Separate Bundle -> BSDF
 */
void wire_root_into_bsdf(Main &bmain,
                         bNodeTree &ntree,
                         bNode &stack,
                         int layer_index,
                         bNode &bsdf,
                         bNodeSocket *bsdf_socket = nullptr);

}  // namespace blender::nodes::texture_stack
