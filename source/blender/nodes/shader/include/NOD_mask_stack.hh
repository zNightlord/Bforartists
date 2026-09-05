/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * The Mask Stack node: a stack of float mask layers composited into a single
 * mask, used as the source of a texture layer's Mask socket (see
 * #NOD_texture_stack.hh). Queries walk the raw link lists, so they are safe
 * while the tree is being edited.
 */

#include <optional>
#include <string>

#include "BLI_string_ref.hh"

#include "NOD_layer_stack.hh"

struct Main;

namespace blender::nodes::mask_stack {

struct ItemsAccessor : public layer_stack::ItemsAccessorBase<ItemsAccessor> {
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

/** Mask Stack node. */
bool is_stack(const bNode &node);

/** The Mask input socket for layer #layer_index on #mask_stack_node, or
 * null when the index is out of range. */
bNodeSocket *layer_socket(bNode &mask_stack_node, int layer_index);

/** The node feeding mask layer #layer_index on #mask_stack_node, or null. */
bNode *layer_source(const bNodeTree &ntree, bNode &mask_stack_node, int layer_index);

/** Return the Mask Stack feeding #item's Mask socket on #stack_node (a
 * Texture Layer Stack node), creating one (and folding any existing source
 * into it as the base item) when the socket is not already driven by a Mask
 * Stack. Null on failure. */
bNode *ensure_for_layer(Main &bmain,
                        bNodeTree &ntree,
                        bNode &stack_node,
                        const NodeShaderLayerStackItem &item);

/** Add a new layer to a Mask Stack above the active one, make it active, and
 * optionally set its Mask default value. Returns the new layer's index. */
int add_layer(Main &bmain,
              bNodeTree &ntree,
              bNode &mask_stack_node,
              const char *name,
              std::optional<float> default_value);

/** Add a mask layer above the active one on #mask_stack_node, fed by
 * #source_socket (an output socket) of #source_node. Returns the new
 * layer's index. */
int add_layer_from_source(Main &bmain,
                          bNodeTree &ntree,
                          bNode &mask_stack_node,
                          const char *name,
                          bNode &source_node,
                          bNodeSocket &source_socket);

/** Set a freshly added mask layer's blend mode to #blend_type (skipping the
 * base layer, whose blend mode is unused). */
void set_layer_blend(bNode &mask_stack_node, int layer_index, int blend_type);

/** Pick #group_node's output that should feed a mask: prefer a SOCK_FLOAT
 * output (single-channel mask), else fall back to a SOCK_BUNDLE output that
 * is split via a Separate Bundle node, picking the most likely channel.
 * #r_source_node is set to the node owning the returned socket (the group
 * node or the created Separate Bundle), or null when none is found. */
bNodeSocket *find_source_output(Main &bmain,
                                bNodeTree &ntree,
                                bNode &group_node,
                                bNode **r_source_node);

}  // namespace blender::nodes::mask_stack
