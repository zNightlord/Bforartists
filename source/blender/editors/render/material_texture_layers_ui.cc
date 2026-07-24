/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edrend
 *
 * Shader Layers UI template: the tree view and inline layer details shown in
 * the Material properties Texture Layers panel.
 */

#include <fmt/format.h>

#include "BLI_listbase_iterator.hh"
#include "BLI_map.hh"
#include "BLI_set.hh"
#include "BLI_string_ref.hh"
#include "BLI_string_utf8.hh"
#include "BLI_utildefines.hh"
#include "BLI_vector.hh"
#include "BLI_vector_set.hh"

#include "DNA_image_types.h"
#include "DNA_material_types.h"
#include "DNA_node_types.h"

#include "BKE_asset.hh"
#include "BKE_context.hh"
#include "BKE_idprop.hh"
#include "BKE_lib_id.hh"
#include "BKE_main.hh"
#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_runtime.hh"

#include "BLT_translation.hh"

#include "DEG_depsgraph.hh"

#include "GPU_material.hh"

#include "AS_asset_representation.hh"

#include "NOD_layer_stack.hh"
#include "NOD_mask_stack.hh"
#include "NOD_node_declaration.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"

#include "ED_asset_menu_utils.hh"
#include "ED_image.hh"
#include "ED_undo.hh"

#include "RNA_access.hh"
#include "RNA_prototypes.hh"

#include "UI_interface.hh"
#include "UI_interface_c.hh"
#include "UI_interface_layout.hh"
#include "UI_resources.hh"
#include "UI_tree_view.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "render_intern.hh"

namespace blender::ed::render {

eShaderNodeTreeUsage asset_texture_layer_usage(const asset_system::AssetRepresentation &asset)
{
  const AssetMetaData &meta = asset.get_metadata();
  const IDProperty *tree_type = BKE_asset_metadata_idprop_find(&meta, "type");
  if (tree_type == nullptr || tree_type->type != IDP_INT || IDP_int_get(tree_type) != NTREE_SHADER)
  {
    return eShaderNodeTreeUsage(0);
  }
  const IDProperty *usage = BKE_asset_metadata_idprop_find(&meta, "shader_usage");
  if (usage == nullptr || usage->type != IDP_INT) {
    return eShaderNodeTreeUsage(0);
  }
  return eShaderNodeTreeUsage(IDP_int_get(usage));
}

}  // namespace blender::ed::render

namespace blender::ui {

/* The tree view, its items and the panel drawing below it are private to this
 * file: an anonymous namespace so their generic names (#ActiveLayer, ...) do
 * not leak into #blender::ui. */
namespace {

namespace layer_stack = nodes::layer_stack;
namespace texture_stack = nodes::texture_stack;
namespace mask_stack = nodes::mask_stack;
namespace texture_channel = nodes::texture_channel;

/* Row for a socket that already has an incoming link: shows its label in the
 * property-split label column and a "Linked" button that removes the link. */
void draw_linked_socket_row(Layout &layout,
                            bNodeTree &ntree,
                            bNodeSocket &socket,
                            const StringRef label)
{
  Layout &row = layout.row(true);
  uiItemL_respect_property_split(&row, label, ICON_NONE);
  Button *but = uiDefIconTextBut(row.block(),
                                 ButtonType::But,
                                 ICON_LINKED,
                                 IFACE_("Linked"),
                                 0,
                                 0,
                                 UI_UNIT_X * 6,
                                 UI_UNIT_Y,
                                 nullptr,
                                 TIP_("Click to remove the link to this socket"));
  bNodeTree *tree_ptr = &ntree;
  bNodeSocket *sock_ptr = &socket;
  button_func_set(but, [tree_ptr, sock_ptr](bContext &C) {
    if (!ed::render::active_material_editable(C)) {
      return;
    }
    bke::node_remove_socket_links(*tree_ptr, *sock_ptr);
    BKE_main_ensure_invariants(*CTX_data_main(&C), tree_ptr->id);
    WM_event_add_notifier(&C, NC_NODE | NA_EDITED, &tree_ptr->id);
  });
}

/* Append the material-level channel toggle: routes the BSDF input through
 * the texture layer stack or restores it to a plain value. */
void draw_channel_toggle(Layout &layout,
                         const bNode &bsdf,
                         const bNodeSocket &socket,
                         const bool textured)
{
  /* Small button on the right: (+) to start texturing a free input, (-) to
   * stop texturing and restore the plain value. Drawn without emboss so it
   * reads as an inline sub-control rather than a full button. */
  Layout &sub = layout.row(true);
  sub.emboss_set(EmbossType::None);
  PointerRNA op_ptr = sub.op("MATERIAL_OT_texture_channel_toggle",
                             "",
                             textured ? ICON_REMOVE : ICON_ADD,
                             wm::OpCallContext::InvokeDefault,
                             UI_ITEM_NONE);
  RNA_string_set(&op_ptr, "node", bsdf.name);
  RNA_string_set(&op_ptr, "socket", socket.identifier);
}

/* Row for a fillable input on the BSDF. A free input shows its value
 * plus a toggle that routes the channel through the texture layer stack; a
 * textured input shows a wide "Textured" button that restores the plain
 * value on click, mirroring the "Linked" row. Returns false for inputs
 * linked to anything other than the stack's Separate Bundle, which use the
 * generic "Linked" row. */
bool draw_bsdf_channel_row(
    Layout &layout, bNodeTree &ntree, bNode &bsdf, bNodeSocket &socket, const StringRef label)
{
  const texture_channel::State state = texture_channel::state(ntree, socket);
  if (state == texture_channel::State::Linked) {
    return false;
  }

  Layout &row = layout.row(true);
  if (state == texture_channel::State::Textured) {
    /* A grayed, centered "Textured" label replaces the value field, indicating
     * the channel is driven by the layer stack. A (-) button on the right stops
     * texturing. */
    uiItemL_respect_property_split(&row, label, ICON_NONE);
    Layout &value = row.row(true);
    value.alignment_set(LayoutAlign::Center);
    value.active_set(false);
    value.label(IFACE_("Textured"), ICON_NONE);
    draw_channel_toggle(row, bsdf, socket, true);
    return true;
  }

  if (!(socket.flag & SOCK_HIDE_VALUE)) {
    PointerRNA socket_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_NodeSocket, &socket);
    row.prop(&socket_ptr, "default_value", UI_ITEM_NONE, label, ICON_NONE);
  }
  else {
    /* No editable value (e.g. Normal, only meaningful linked): just the
     * label and the toggle. */
    uiItemL_respect_property_split(&row, label, ICON_NONE);
    row.label("", ICON_NONE);
  }
  draw_channel_toggle(row, bsdf, socket, false);
  return true;
}

/* Draw one of a node's input sockets using property-split-friendly layouts.
 * Unlinked sockets use #Layout::prop on their "default_value" so the layout's
 * property-split cooperates; linked sockets show a "Linked" row with an unlink
 * button so the user can see and remove them.
 * With #bsdf_channels set (the node is the consuming BSDF), fillable inputs
 * get the channel-row treatment with the texture-channel toggle. */
void draw_node_input(Layout &layout,
                     bNode &node,
                     PointerRNA *node_ptr,
                     bNodeSocket &socket,
                     const bool bsdf_channels = false)
{
  if (!socket.is_available()) {
    return;
  }
  if (socket.typeinfo == nullptr || socket.typeinfo->draw == nullptr) {
    return;
  }
  if (ELEM(socket.type, SOCK_GEOMETRY, SOCK_MATRIX, SOCK_SHADER, SOCK_BUNDLE, SOCK_CLOSURE)) {
    return;
  }
  if (node.is_reroute()) {
    return;
  }
  if (socket.idname == StringRef("NodeSocketVirtual")) {
    return;
  }

  bNodeTree &ntree = *reinterpret_cast<bNodeTree *>(node_ptr->owner_id);
  const StringRef label = CTX_IFACE_(bke::node_socket_translation_context(socket),
                                     bke::node_socket_label(socket));

  if (bsdf_channels && texture_channel::is_fillable_input(socket)) {
    if (draw_bsdf_channel_row(layout, ntree, node, socket, label)) {
      return;
    }
  }

  if (socket.is_directly_linked()) {
    draw_linked_socket_row(layout, ntree, socket, label);
    return;
  }
  if (socket.flag & SOCK_HIDE_VALUE) {
    return;
  }

  PointerRNA socket_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_NodeSocket, &socket);
  layout.prop(&socket_ptr, "default_value", UI_ITEM_NONE, label, ICON_NONE);
}

using nodes::ItemDeclaration;
using nodes::LayoutDeclaration;
using nodes::NodeDeclaration;
using nodes::PanelDeclaration;
using nodes::SocketDeclaration;

/* Render a node's draw_buttons + input sockets inline, with property split. */
void draw_node_inputs(Layout &layout,
                      bContext *C,
                      PointerRNA *node_ptr,
                      const bool bsdf_channels = false)
{
  bNodeTree &tree = *reinterpret_cast<bNodeTree *>(node_ptr->owner_id);
  bNode &node = *static_cast<bNode *>(node_ptr->data);
  tree.ensure_topology_cache();

  /* Prefer the compact draw_buttons (an Image Texture's interpolation and
   * extension, etc.) over draw_buttons_ex, whose extended templates (e.g. the
   * full image template) add their own sections and separators. Group assets
   * hide their datablock selector: the layer already identifies the asset. */
  if (!node.is_group()) {
    if (node.typeinfo->draw_buttons) {
      node.typeinfo->draw_buttons(layout, C, node_ptr);
    }
    else if (node.typeinfo->draw_buttons_ex) {
      node.typeinfo->draw_buttons_ex(layout, C, node_ptr);
    }
  }

  if (node.declaration()) {
    const NodeDeclaration &node_decl = *node.declaration();
    for (const ItemDeclaration *item_decl : node_decl.root_items) {
      if (const auto *panel_decl = dynamic_cast<const PanelDeclaration *>(item_decl)) {
        /* Reuse the generic node-inputs traversal, drawing each socket as a texture-layer
         * channel row. */
        ui::draw_node_inputs_recursive(
            C, layout, node, node_ptr, *panel_decl, [&](Layout &row, bNodeSocket &socket) {
              draw_node_input(row, node, node_ptr, socket, bsdf_channels);
            });
      }
      else if (const auto *socket_decl = dynamic_cast<const SocketDeclaration *>(item_decl)) {
        if (socket_decl->in_out == SOCK_IN) {
          draw_node_input(
              layout, node, node_ptr, node.socket_by_decl(*socket_decl), bsdf_channels);
        }
      }
      else if (const auto *layout_decl = dynamic_cast<const LayoutDeclaration *>(item_decl)) {
        if (!layout_decl->is_default) {
          layout_decl->draw(layout, C, node_ptr);
        }
      }
    }
  }
  else {
    for (bNodeSocket *input : node.runtime->inputs) {
      draw_node_input(layout, node, node_ptr, *input, bsdf_channels);
    }
  }
}

/* The generator nodes feeding #combine's items, in item order and without
 * duplicates (one generator may drive several channels). */
VectorSet<bNode *> combine_bundle_generators(bNodeTree &ntree, const bNode &combine)
{
  ntree.ensure_topology_cache();
  VectorSet<bNode *> generators;
  for (bNodeSocket *input : combine.runtime->inputs) {
    if (bNode *generator = bke::node_find_source_node(ntree, *input)) {
      generators.add(generator);
    }
  }
  return generators;
}

/*
 * Layer details for a layer whose bundle comes from a Combine Bundle: the
 * unlinked item values are the layer's channel constants (e.g. a Fill layer's
 * colors), and the nodes feeding the linked items are the generators whose
 * inputs are the layer's main options (e.g. a Noise Texture). The bundle
 * node's own item-list UI stays in the node editor; which channels the bundle
 * carries is managed through the Channels panel and the material's channel
 * toggles.
 */
void draw_combine_bundle_source(Layout &layout,
                                bContext *C,
                                bNodeTree &ntree,
                                bNode &combine,
                                const NodeShaderLayerStackItem *item)
{
  ntree.ensure_topology_cache();

  /* A channel disabled in the Influence panel does not contribute, so its input
   * draws inactive (grayed out) but stays editable. The combine input's name is
   * the channel key. */
  auto channel_disabled = [&](const bNodeSocket &input) {
    return item && layer_stack::channel_disabled(*item, input.name);
  };

  PointerRNA combine_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_Node, &combine);
  bool any_values = false;
  for (bNodeSocket *input : combine.runtime->inputs) {
    if (input->is_directly_linked()) {
      continue;
    }
    /* Skip the trailing hollow extend socket: it draws nothing, so counting it
     * as a value would insert a spurious separator above the generator below. */
    if (!input->is_available() || input->idname == StringRef("NodeSocketVirtual")) {
      continue;
    }
    Layout &col = layout.column(true);
    col.active_set(!channel_disabled(*input));
    draw_node_input(col, combine, &combine_ptr, *input);
    any_values = true;
  }

  /* Generators feeding a channel: gray the whole node section when every
   * channel it drives is disabled. */
  const VectorSet<bNode *> generators = combine_bundle_generators(ntree, combine);
  for (const int i : generators.index_range()) {
    bNode *generator = generators[i];
    bool all_disabled = true;
    for (bNodeSocket *input : combine.runtime->inputs) {
      if (bke::node_find_source_node(ntree, *input) == generator && !channel_disabled(*input)) {
        all_disabled = false;
        break;
      }
    }
    if (any_values || i > 0) {
      layout.separator(1.0f, LayoutSeparatorType::Line);
    }
    Layout &col = layout.column(false);
    col.active_set(!all_disabled);
    if (generators.size() > 1) {
      col.label(bke::node_label(ntree, *generator), ICON_NONE);
    }
    PointerRNA generator_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_Node, generator);
    draw_node_inputs(col, C, &generator_ptr);
  }
}

/* Draw an item's header: blend mode, mute and opacity. Shared by texture
 * layer and mask items via #item_srna. The base (last) item has no blend
 * mode or mute, but still shows Opacity when it has an Opacity socket (a base
 * mask). Masks are not mentioned here: the tree view shows and manages each
 * layer's mask stack. */
void draw_layer_header(Layout &layout,
                       bNodeTree &ntree,
                       bNode &stack_node,
                       const int item_index,
                       StructRNA *item_srna)
{
  NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
  if (item_index < 0 || item_index >= storage.items_num) {
    return;
  }
  const bool is_base = (item_index == storage.items_num - 1);
  NodeShaderLayerStackItem &item = storage.items[item_index];

  /* Blend mode, mute and opacity in a single aligned column so the fields
   * connect visually (like col.column(align=True) in the Python API). */
  Layout &col = layout.column(true);
  PointerRNA item_ptr = RNA_pointer_create_discrete(&ntree.id, item_srna, &item);
  if (!is_base) {
    col.prop(&item_ptr, "blend_type", UI_ITEM_NONE, IFACE_("Blend Mode"), ICON_NONE);
    col.prop(&item_ptr, "mute", UI_ITEM_NONE, IFACE_("Mute"), ICON_NONE);
  }

  bNodeSocket *opacity_sock =
      texture_stack::StackLayer{&ntree, &stack_node, item_index}.opacity_socket();
  if (opacity_sock) {
    if (opacity_sock->is_directly_linked()) {
      draw_linked_socket_row(col, ntree, *opacity_sock, IFACE_("Opacity"));
    }
    else {
      PointerRNA sock_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_NodeSocket, opacity_sock);
      col.prop(&sock_ptr, "default_value", UI_ITEM_NONE, IFACE_("Opacity"), ICON_NONE);
    }
  }
}

/* Toggle a channel on the item identified by #item_identifier on
 * #stack_node, with tree update and notifier (the button's automatic
 * BUT_UNDO handles the undo push). Resolved by identifier at click time: the
 * items array may have been reallocated between drawing and the button callback
 * firing. */
void toggle_layer_channel(bContext &C,
                          bNodeTree &ntree,
                          bNode &stack_node,
                          const int item_identifier,
                          const StringRef channel,
                          const bool disable)
{
  if (!ed::render::active_material_editable(C)) {
    return;
  }
  NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
  for (NodeShaderLayerStackItem &item : MutableSpan(storage.items, storage.items_num)) {
    if (item.identifier != item_identifier) {
      continue;
    }
    layer_stack::set_channel_disabled(item, channel, disable);
    BKE_ntree_update_tag_node_property(&ntree, &stack_node);
    BKE_main_ensure_invariants(*CTX_data_main(&C), ntree.id);
    WM_event_add_notifier(&C, NC_NODE | NA_EDITED, &ntree.id);
    /* No ED_undo_push here: this runs from a ButtonType::But callback, which the button
     * handler already flags with BUT_UNDO and pushes after the callback returns (like the
     * unlink button above). Pushing here too would create a redundant second undo step. */
    return;
  }
}

/*
 * True when a group layer enclosing #stack_node does not contribute #channel.
 * Its nested layers then have no way to reach the material through it either,
 * so their toggle for the channel is shown inactive.
 *
 * TODO: an adjustment layer inside the group may still affect the channel
 * through crosstalk (deriving one channel from another), so this does not
 * strictly hold for those.
 */
bool channel_disabled_by_group(bNodeTree &ntree, bNode &stack_node, const StringRef channel)
{
  Set<const bNode *> visited;
  std::optional<texture_stack::StackLayer> parent = texture_stack::parent_group_layer(ntree,
                                                                                      stack_node);
  /* Cyclic stack links are kept in the tree (marked invalid): guard the walk. */
  while (parent && visited.add(parent->stack)) {
    if (layer_stack::channel_disabled(parent->item(), channel)) {
      return true;
    }
    parent = texture_stack::parent_group_layer(ntree, *parent->stack);
  }
  return false;
}

/*
 * Influence sub-panel for the active layer or mask row: its blend mode, mute
 * and opacity, and for a texture layer one row per channel the material
 * textures (BSDF inputs routed through the stack), with a checkbox for whether
 * the layer contributes to it. Channels the layer's bundle does not carry draw
 * grayed out, as the layer cannot affect them regardless. The layer's values
 * are edited in the inputs list above. A mask has no channels: it is one float
 * feeding the layer it belongs to.
 */
void draw_influence_panel(Layout &layout,
                          bContext *C,
                          bNodeTree &ntree,
                          bNode &stack_node,
                          const int item_index,
                          bNode &bsdf,
                          const bool is_mask)
{
  NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
  if (item_index < 0 || item_index >= storage.items_num) {
    return;
  }
  NodeShaderLayerStackItem &item = storage.items[item_index];

  PanelLayout panel = layout.panel(
      C, is_mask ? "mask_layer_influence" : "texture_layer_influence", false);
  panel.header->label(IFACE_("Influence"), ICON_NONE);
  if (!panel.body) {
    return;
  }
  panel.body->use_property_split_set(true);

  /* Blend mode, mute and opacity for the layer, above the per-channel toggles. */
  draw_layer_header(*panel.body,
                    ntree,
                    stack_node,
                    item_index,
                    is_mask ? RNA_ShaderMaskStackItem : RNA_ShaderTextureLayerStackItem);
  if (is_mask) {
    return;
  }

  /* Aligned column so the channel on/off buttons connect visually
   * (like col.column(align=True) in the Python API). Every layer accounts for all of the stack's
   * channels; whether it contributes to each is the per-layer toggle below. */
  Layout &channels_col = panel.body->column(true);
  for (bNodeSocket &sock : bsdf.inputs) {
    if (!texture_channel::is_fillable_input(sock)) {
      continue;
    }
    /* Only the channels the material textures; the others are toggled from
     * the BSDF row. */
    if (texture_channel::state(ntree, sock) != texture_channel::State::Textured) {
      continue;
    }
    const StringRef channel = sock.name;
    const bool enabled = !layer_stack::channel_disabled(item, channel);

    Layout &row = channels_col.row(true);
    row.active_set(!channel_disabled_by_group(ntree, stack_node, channel));
    uiItemL_respect_property_split(&row, IFACE_(sock.name), ICON_NONE);
    Layout &check = row.row(true);
    check.emboss_set(EmbossType::None);
    Button *but = uiDefIconTextBut(check.block(),
                                   ButtonType::But,
                                   enabled ? ICON_CHECKBOX_HLT : ICON_CHECKBOX_DEHLT,
                                   "",
                                   0,
                                   0,
                                   UI_UNIT_X,
                                   UI_UNIT_Y,
                                   nullptr,
                                   TIP_("Whether this layer contributes to the channel"));
    bNodeTree *tree_ptr = &ntree;
    bNode *stack_ptr = &stack_node;
    const int item_identifier = item.identifier;
    button_func_set(
        but,
        [tree_ptr, stack_ptr, item_identifier, name = std::string(channel), enabled](bContext &C) {
          toggle_layer_channel(C, *tree_ptr, *stack_ptr, item_identifier, name, enabled);
        });
  }
}

/* The texture nodes in a layer's source that carry the built-in texture
 * coordinate mapping (Image Texture, Noise Texture, ...). For a generator
 * wrapped in a Combine Bundle these are the nodes feeding its items; otherwise
 * it is the source node itself. */
void collect_texture_mapping_nodes(bNodeTree &ntree, bNode &source, Vector<bNode *> &r_nodes)
{
  VectorSet<bNode *> nodes;
  if (source.is_type("NodeCombineBundle"_ustr)) {
    nodes = combine_bundle_generators(ntree, source);
  }
  else {
    nodes.add(&source);
  }
  for (bNode *node : nodes) {
    PointerRNA node_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_Node, node);
    if (RNA_struct_find_property(&node_ptr, "texture_mapping")) {
      r_nodes.append(node);
    }
  }
}

/*
 * Texture Mapping panel: the built-in texture coordinate mapping (vector type,
 * projection, location/rotation/scale) that shader texture nodes carry, shown
 * for the layer's generator texture nodes. Mirrors the node editor sidebar's
 * NODE_PT_texture_mapping. Collapsed by default. Nothing is drawn when the
 * layer has no mapping-carrying texture nodes.
 */
void draw_texture_mapping_panel(Layout &layout,
                                bContext *C,
                                bNodeTree &ntree,
                                const Span<bNode *> nodes)
{
  if (nodes.is_empty()) {
    return;
  }

  PanelLayout panel = layout.panel(C, "texture_layer_mapping", true);
  panel.header->label(IFACE_("Texture Mapping"), ICON_NONE);
  if (!panel.body) {
    return;
  }
  Layout &body = *panel.body;
  body.use_property_split_set(true);
  body.use_property_decorate_set(false);

  for (const int i : nodes.index_range()) {
    bNode *node = nodes[i];
    if (i > 0) {
      body.separator(1.0f, LayoutSeparatorType::Line);
      body.label(bke::node_label(ntree, *node), ICON_NONE);
    }
    PointerRNA node_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_Node, node);
    PointerRNA mapping = RNA_pointer_get(&node_ptr, "texture_mapping");

    body.prop(&mapping, "vector_type", UI_ITEM_NONE, std::nullopt, ICON_NONE);
    body.separator();
    Layout &col = body.column(true);
    col.prop(&mapping, "mapping_x", UI_ITEM_NONE, IFACE_("Projection X"), ICON_NONE);
    col.prop(&mapping, "mapping_y", UI_ITEM_NONE, IFACE_("Y"), ICON_NONE);
    col.prop(&mapping, "mapping_z", UI_ITEM_NONE, IFACE_("Z"), ICON_NONE);
    body.separator();
    body.prop(&mapping, "translation", UI_ITEM_NONE, std::nullopt, ICON_NONE);
    body.prop(&mapping, "rotation", UI_ITEM_NONE, std::nullopt, ICON_NONE);
    body.prop(&mapping, "scale", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  }
}

/* The nodes in a layer's source whose own properties the layer layout leaves
 * out: the Combine Bundle wrapping a Fill or generator layer (its item list
 * defines which channels the layer carries) and node group nodes (their
 * datablock selector). */
void collect_advanced_nodes(bNodeTree &ntree, bNode &source, Vector<bNode *> &r_nodes)
{
  if (source.is_type("NodeCombineBundle"_ustr)) {
    r_nodes.append(&source);
    for (bNode *generator : combine_bundle_generators(ntree, source)) {
      if (generator->is_group()) {
        r_nodes.append(generator);
      }
    }
    return;
  }
  if (source.is_group()) {
    r_nodes.append(&source);
  }
}

/*
 * Advanced panel: the node properties #collect_advanced_nodes gathers, which
 * describe the nodes backing the layer rather than the layer itself. Collapsed
 * by default, at the bottom. Nothing is drawn when the layer has no such nodes.
 */
void draw_advanced_panel(Layout &layout, bContext *C, bNodeTree &ntree, const Span<bNode *> nodes)
{
  if (nodes.is_empty()) {
    return;
  }

  PanelLayout panel = layout.panel(C, "texture_layer_advanced", true);
  panel.header->label(IFACE_("Advanced"), ICON_NONE);
  if (!panel.body) {
    return;
  }
  Layout &body = *panel.body;
  body.use_property_split_set(true);
  body.use_property_decorate_set(false);

  for (const int i : nodes.index_range()) {
    bNode *node = nodes[i];
    if (i > 0) {
      body.separator(1.0f, LayoutSeparatorType::Line);
      body.label(bke::node_label(ntree, *node), ICON_NONE);
    }
    PointerRNA node_ptr = RNA_pointer_create_discrete(&ntree.id, RNA_Node, node);
    if (node->is_group()) {
      if (node->typeinfo->draw_buttons) {
        node->typeinfo->draw_buttons(body, C, &node_ptr);
      }
    }
    else if (node->typeinfo->draw_buttons_ex) {
      node->typeinfo->draw_buttons_ex(body, C, &node_ptr);
    }
  }
}

/* Tree item kinds for the Shader Layers view. */
enum class EntryKind {
  /* BSDF node at the root of the tree. */
  Bsdf,
  /* A texture layer item from a #ShaderNodeTextureLayerStack. */
  Layer,
  /* A mask item from a #ShaderNodeMaskStack, child of its layer's row. */
  Mask,
};

class ShaderLayerItem : public AbstractTreeViewItem {
 private:
  bNodeTree &ntree_;
  /* Entry kind. */
  EntryKind kind_;
  /* For #Bsdf: the BSDF node. For #Layer: the stack node containing the item. */
  bNode *node_;
  /* For #Layer: index of the layer item in the stack. */
  int item_index_;
  /* Cached: when this layer is a group, the nested stack node. Null otherwise. */
  bNode *nested_stack_;

 public:
  /* Constructor for a BSDF entry. */
  ShaderLayerItem(bNodeTree &ntree, bNode &bsdf)
      : ntree_(ntree),
        kind_(EntryKind::Bsdf),
        node_(&bsdf),
        item_index_(-1),
        nested_stack_(nullptr)
  {
    this->label_ = bke::node_label(ntree_, *node_);
  }

  /* Constructor for a layer or mask entry; #stack_node is the Texture
   * Layer Stack respectively Mask Stack node holding the item. */
  ShaderLayerItem(bNodeTree &ntree, const EntryKind kind, bNode &stack_node, const int item_index)
      : ntree_(ntree),
        kind_(kind),
        node_(&stack_node),
        item_index_(item_index),
        nested_stack_(
            kind == EntryKind::Layer ?
                texture_stack::StackLayer{&ntree, &stack_node, item_index}.nested_stack() :
                nullptr)
  {
    const NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
    BLI_assert(item_index >= 0 && item_index < storage.items_num);
    const NodeShaderLayerStackItem &item = storage.items[item_index];
    this->label_ = item.name ? StringRef(item.name) : StringRef(IFACE_("Layer"));
  }

  /* RNA type of the item this row represents. */
  StructRNA *item_srna() const
  {
    return kind_ == EntryKind::Mask ? RNA_ShaderMaskStackItem : RNA_ShaderTextureLayerStackItem;
  }

  bool is_group() const
  {
    return nested_stack_ != nullptr;
  }

  /* Row icon showing what kind of layer this is; the hook for per-layer
   * previews and further layer kinds (decals, ...) later on. */
  int layer_icon(const NodeShaderLayerStackItem &item) const
  {
    if (kind_ == EntryKind::Mask) {
      return ICON_MOD_MASK;
    }
    if (is_group()) {
      return ICON_FILE_FOLDER;
    }
    if (item.item_type == SHADER_LAYER_STACK_ITEM_CLOSURE) {
      return ICON_SHADERFX;
    }
    return this->generator_icon();
  }

  /* Icon for a generator/Fill layer, reflecting its source node: a Fill
   * (Combine Bundle fed only by constant values) uses the add menu's Fill
   * icon, an Image Texture uses the image icon, and other generators fall back
   * to the generic Generator icon the add menu uses. */
  int generator_icon() const
  {
    bNode *source = texture_stack::StackLayer{&ntree_, node_, item_index_}.source();
    if (source && source->is_type("ShaderNodeTexImage"_ustr)) {
      /* A Paint layer: an Image Texture bundle feeding the layer directly. */
      return ICON_IMAGE_DATA;
    }
    if (source && source->is_type("NodeCombineBundle"_ustr)) {
      const VectorSet<bNode *> generators = combine_bundle_generators(ntree_, *source);
      for (const bNode *generator : generators) {
        if (generator->is_type("ShaderNodeTexImage"_ustr)) {
          return ICON_IMAGE_DATA;
        }
      }
      if (generators.is_empty()) {
        /* No generator feeding it: a Fill layer's constant channel values. */
        return ICON_SHADING_SOLID;
      }
    }
    return ICON_NODE_TEXTURE;
  }

  bNode *stack_node() const
  {
    return kind_ == EntryKind::Layer ? node_ : nullptr;
  }

  bNode *nested_stack() const
  {
    return nested_stack_;
  }

  int item_index() const
  {
    return item_index_;
  }

  void build_row(Layout &row) override
  {
    switch (kind_) {
      case EntryKind::Bsdf:
        row.label(this->label_, ICON_NODE_MATERIAL);
        break;
      case EntryKind::Layer:
      case EntryKind::Mask: {
        NodeShaderLayerStack &storage = layer_stack::storage(*node_);
        NodeShaderLayerStackItem &item = storage.items[item_index_];

        row.label(this->label_, this->layer_icon(item));

        Layout &sub = row.row(true);
        sub.alignment_set(LayoutAlign::Right);
        sub.use_property_decorate_set(false);

        PointerRNA item_ptr = RNA_pointer_create_discrete(&ntree_.id, this->item_srna(), &item);

        /* Show the opacity as a compact, editable percentage (no emboss) when
         * the item has an Opacity socket that isn't driven by a link. The base
         * layer has none, but a base mask does. Right aligned so it sits
         * against the mute button. */
        bNodeSocket *opacity_sock =
            texture_stack::StackLayer{&ntree_, node_, item_index_}.opacity_socket();
        if (opacity_sock && !opacity_sock->is_directly_linked()) {
          Layout &opacity_row = sub.row(true);
          opacity_row.alignment_set(LayoutAlign::Right);
          opacity_row.emboss_set(EmbossType::None);
          opacity_row.prop(&item_ptr, "opacity", UI_ITEM_NONE, "", ICON_NONE);
        }

        /* Mute toggle. */
        Layout &mute_row = sub.row(true);
        mute_row.emboss_set(EmbossType::None);
        const int icon = (item.flag & SHADER_LAYER_STACK_ITEM_MUTED) ? ICON_HIDE_ON :
                                                                       ICON_HIDE_OFF;
        mute_row.prop(&item_ptr, "mute", ITEM_R_ICON_ONLY, std::nullopt, icon);
        break;
      }
    }
  }

  std::optional<bool> should_be_active() const override
  {
    bNode *active = bke::node_get_active(ntree_);
    if (active == nullptr) {
      return false;
    }
    if (kind_ == EntryKind::Bsdf) {
      return active == node_;
    }
    /* The active node is one of a layer's nodes (or the stack node for a layer
     * that has none): either way it resolves to the layer it belongs to, so
     * selecting a node in the node editor highlights its row here. */
    const std::optional<texture_stack::StackLayer> layer = texture_stack::layer_for_node(ntree_,
                                                                                         *active);
    return layer && layer->stack == node_ && layer->index == item_index_;
  }

  /* The collapse state is stored on the layer item so it survives undo and
   * save/load (transient view state is keyed on node pointers, which undo
   * reallocates). #should_be_collapsed runs after the pointer-based state
   * transfer, so the stored value wins. */
  std::optional<bool> should_be_collapsed() const override
  {
    if (kind_ == EntryKind::Bsdf) {
      return std::nullopt;
    }
    const NodeShaderLayerStack &storage = layer_stack::storage(*node_);
    return (storage.items[item_index_].flag & SHADER_LAYER_STACK_ITEM_COLLAPSED) != 0;
  }

  bool set_collapsed(const bool collapsed) override
  {
    if (!AbstractTreeViewItem::set_collapsed(collapsed)) {
      return false;
    }
    if (kind_ != EntryKind::Bsdf) {
      NodeShaderLayerStackItem &item = layer_stack::storage(*node_).items[item_index_];
      SET_FLAG_FROM_TEST(item.flag, collapsed, SHADER_LAYER_STACK_ITEM_COLLAPSED);
    }
    return true;
  }

  void on_collapse_change(bContext &C, const bool is_collapsed) override
  {
    if (kind_ == EntryKind::Bsdf) {
      return;
    }
    /* Route the change through RNA so the notifier fires and the state is
     * tracked consistently (matching e.g. the grease pencil layer tree). */
    NodeShaderLayerStackItem &item = layer_stack::storage(*node_).items[item_index_];
    PointerRNA item_ptr = RNA_pointer_create_discrete(&ntree_.id, this->item_srna(), &item);
    PropertyRNA *prop = RNA_struct_find_property(&item_ptr, "show_expanded");
    RNA_property_boolean_set(&item_ptr, prop, !is_collapsed);
    RNA_property_update(&C, &item_ptr, prop);
  }

  /* The Image Texture node holding the row's paintable image: the layer's
   * direct bundle source (a Paint layer), an image feeding the layer's
   * Combine Bundle (single-channel image generators), or a paint mask's
   * source. Null when the row has no image to paint. */
  bNode *paint_image_node() const
  {
    bNode *source = (kind_ == EntryKind::Mask) ?
                        mask_stack::layer_source(ntree_, *node_, item_index_) :
                        texture_stack::StackLayer{&ntree_, node_, item_index_}.source();
    if (source == nullptr) {
      return nullptr;
    }
    if (source->is_type("ShaderNodeTexImage"_ustr)) {
      return source;
    }
    if (source->is_type("NodeCombineBundle"_ustr)) {
      for (bNode *generator : combine_bundle_generators(ntree_, *source)) {
        if (generator->is_type("ShaderNodeTexImage"_ustr)) {
          return generator;
        }
      }
    }
    return nullptr;
  }

  /* When the row is backed by an Image Texture node, make it the active paint
   * canvas and sync the texture paint slot and image editors, matching the node
   * editor's active-texture handling (see #ED_node_set_active). Runs before the
   * row's own activation, which assigns the active node (and with it the paint
   * canvas flag this checks). */
  void activate_paint_image(bContext &C) const
  {
    bNode *image_node = this->paint_image_node();
    if (image_node == nullptr || (image_node->flag & NODE_ACTIVE_PAINT_CANVAS)) {
      return;
    }
    bke::node_set_active_texture(ntree_, *image_node);

    Main &bmain = *CTX_data_main(&C);
    Image *image = (image_node->id && GS(image_node->id->name) == ID_IM) ?
                       reinterpret_cast<Image *>(image_node->id) :
                       nullptr;
    for (Material &ma : bmain.materials) {
      if (ma.nodetree && bke::node_tree_contains_tree(*ma.nodetree, ntree_)) {
        /* The viewport's texture display follows the active texture. */
        GPU_material_free(&ma.gpumaterial);
        if (image && ma.texpaintslot) {
          for (int i = 0; i < ma.tot_slots; i++) {
            if (ma.texpaintslot[i].ima == image) {
              ma.paint_active_slot = i;
              DEG_id_tag_update(&ma.id, ID_RECALC_SYNC_TO_EVAL);
            }
          }
        }
      }
    }
    if (image) {
      ED_space_image_sync(&bmain, image, true);
    }
    WM_event_add_notifier(&C, NC_IMAGE, nullptr);
  }

  void on_activate(bContext &C) override
  {
    if (!ed::render::active_material_editable(C)) {
      return;
    }
    if (kind_ == EntryKind::Bsdf) {
      bke::node_set_active(ntree_, *node_);
      for (bNode &node : ntree_.nodes) {
        bke::node_set_selected(node, &node == node_);
      }
    }
    else {
      this->activate_paint_image(C);
      /* Activating a layer points the node editor at the nodes making it up,
       * rather than at the stack node, so its active node and selection show
       * what this row's panels edit. */
      ed::render::set_active_layer(ntree_, *node_, item_index_);
    }
    WM_event_add_notifier(&C, NC_NODE | NA_SELECTED, &ntree_.id);
    ED_undo_push(&C, "Activate Texture Layer");
  }

  void build_context_menu(bContext &C, Layout &column) const override
  {
    if (kind_ == EntryKind::Bsdf) {
      return;
    }
    const char *menu_id = (kind_ == EntryKind::Mask) ? "MATERIAL_MT_texture_layer_mask_context" :
                                                       "MATERIAL_MT_texture_layer_context";
    MenuType *mt = WM_menutype_find(menu_id, true);
    if (!mt) {
      return;
    }
    menutype_draw(&C, mt, &column);
  }

  void delete_item(bContext *C) override
  {
    if (kind_ == EntryKind::Bsdf) {
      return;
    }
    /* Make this row the active one so the remove operator targets it; the
     * operator handles the undo push. */
    ed::render::set_active_layer(ntree_, *node_, item_index_);
    const char *op_id = (kind_ == EntryKind::Mask) ? "MATERIAL_OT_texture_layer_mask_remove" :
                                                     "MATERIAL_OT_texture_layer_remove";
    WM_operator_name_call(C, op_id, wm::OpCallContext::InvokeDefault, nullptr, nullptr);
  }

  std::unique_ptr<AbstractViewItemDragController> create_drag_controller() const override;
  std::unique_ptr<TreeViewItemDropTarget> create_drop_target() override;

  /* Layer and mask rows are renameable via double-click; the BSDF row is not. */
  bool supports_renaming() const override
  {
    return kind_ != EntryKind::Bsdf;
  }

  StringRef get_rename_string() const override
  {
    if (kind_ == EntryKind::Bsdf) {
      return "";
    }
    const NodeShaderLayerStack &storage = layer_stack::storage(*node_);
    const NodeShaderLayerStackItem &item = storage.items[item_index_];
    return item.name ? StringRef(item.name) : StringRef("");
  }

  bool rename(const bContext &C, StringRefNull new_name) override
  {
    if (kind_ == EntryKind::Bsdf || !ed::render::active_material_editable(C)) {
      return false;
    }
    NodeShaderLayerStack &storage = layer_stack::storage(*node_);
    NodeShaderLayerStackItem &item = storage.items[item_index_];
    PointerRNA item_ptr = RNA_pointer_create_discrete(&ntree_.id, this->item_srna(), &item);
    PropertyRNA *prop = RNA_struct_find_property(&item_ptr, "name");
    RNA_property_string_set(&item_ptr, prop, new_name.c_str());
    RNA_property_update(&const_cast<bContext &>(C), &item_ptr, prop);
    ED_undo_push(&const_cast<bContext &>(C), "Rename Texture Layer");
    return true;
  }

 protected:
  bool matches_single(const AbstractTreeViewItem &other) const override
  {
    const ShaderLayerItem *o = dynamic_cast<const ShaderLayerItem *>(&other);
    return o && o->kind_ == kind_ && o->node_ == node_ && o->item_index_ == item_index_;
  }
};

/* Drag controller for a Texture Layer item. Captures the source stack node
 * and item index so the drop target can move the item between stacks. */
class ShaderLayerDragController : public AbstractViewItemDragController {
  bNodeTree &ntree_;
  bNode &stack_node_;
  int item_index_;

 public:
  ShaderLayerDragController(AbstractTreeView &tree_view,
                            bNodeTree &ntree,
                            bNode &stack_node,
                            const int item_index)
      : AbstractViewItemDragController(tree_view),
        ntree_(ntree),
        stack_node_(stack_node),
        item_index_(item_index)
  {
  }

  std::optional<eWM_DragDataType> get_drag_type() const override
  {
    return WM_DRAG_TEXTURE_LAYER;
  }

  void *create_drag_data() const override
  {
    wmDragTextureLayer *data = MEM_new_zeroed<wmDragTextureLayer>(__func__);
    data->ntree = &ntree_;
    data->stack_node = &stack_node_;
    data->item_index = item_index_;
    return data;
  }
};

/* Drop target for a Texture Layer item. Reorders around the row in its parent
 * stack (Before/After) and, when the layer is a group, drops INTO the nested
 * stack. The actual move is delegated to #material.texture_layer_reparent. */
class ShaderLayerDropTarget : public TreeViewItemDropTarget {
  bNodeTree &ntree_;
  /* Parent stack and index, used for Before/After reorder. */
  bNode *parent_stack_;
  int parent_index_;
  /* Nested stack when this row is a group, for INTO. Null otherwise. */
  bNode *nested_stack_;

 public:
  ShaderLayerDropTarget(AbstractTreeViewItem &item,
                        DropBehavior behavior,
                        bNodeTree &ntree,
                        bNode *parent_stack,
                        const int parent_index,
                        bNode *nested_stack)
      : TreeViewItemDropTarget(item, behavior),
        ntree_(ntree),
        parent_stack_(parent_stack),
        parent_index_(parent_index),
        nested_stack_(nested_stack)
  {
  }

  /* Asset drag payload when #drag is a texture-layer-capable shader node
   * group asset, else null. */
  static const wmDragAsset *asset_drag_data(const wmDrag &drag)
  {
    if (drag.type != WM_DRAG_ASSET) {
      return nullptr;
    }
    const wmDragAsset *asset_drag = WM_drag_get_asset_data(&drag, ID_NT);
    if (asset_drag == nullptr || asset_drag->asset == nullptr) {
      return nullptr;
    }
    if (ed::render::asset_texture_layer_usage(*asset_drag->asset) == 0) {
      return nullptr;
    }
    return asset_drag;
  }

  bool can_drop(const wmDrag &drag, const char ** /*r_disabled_hint*/) const override
  {
    if (parent_stack_ == nullptr) {
      return false;
    }
    if (drag.type == WM_DRAG_ASSET) {
      /* Generator/adjustment/mask assets drop onto layer rows (not mask rows;
       * their own stack management is by layer). */
      return !mask_stack::is_stack(*parent_stack_) && asset_drag_data(drag) != nullptr;
    }
    if (drag.type != WM_DRAG_TEXTURE_LAYER) {
      return false;
    }
    const wmDragTextureLayer *data = static_cast<const wmDragTextureLayer *>(drag.poin);
    if (data->ntree != &ntree_ || data->stack_node == nullptr) {
      return false;
    }
    /* Masks reorder between mask rows, layers between layer rows — a drag
     * never crosses the two kinds. */
    if (mask_stack::is_stack(*data->stack_node) != mask_stack::is_stack(*parent_stack_)) {
      return false;
    }
    /* Don't drop onto self. */
    if (data->stack_node == parent_stack_ && data->item_index == parent_index_) {
      return false;
    }
    return true;
  }

  std::string drop_tooltip(const DragInfo &drag_info) const override
  {
    if (const wmDragAsset *asset_drag = asset_drag_data(drag_info.drag_data)) {
      const StringRefNull name = asset_drag->asset->get_name();
      if (this->asset_drops_as_mask(*asset_drag, drag_info.drop_location)) {
        return fmt::format(fmt::runtime(TIP_("Add {} as mask")), name);
      }
      if (drag_info.drop_location == DropLocation::After) {
        return fmt::format(fmt::runtime(TIP_("Insert {} layer below")), name);
      }
      return fmt::format(fmt::runtime(TIP_("Insert {} layer above")), name);
    }
    const wmDragTextureLayer *data = static_cast<const wmDragTextureLayer *>(
        drag_info.drag_data.poin);
    StringRef drag_name = IFACE_("Layer");
    if (data && data->stack_node) {
      const NodeShaderLayerStack &storage = layer_stack::storage(*data->stack_node);
      if (data->item_index >= 0 && data->item_index < storage.items_num) {
        const NodeShaderLayerStackItem &item = storage.items[data->item_index];
        if (item.name) {
          drag_name = StringRef(item.name);
        }
      }
    }
    switch (drag_info.drop_location) {
      case DropLocation::Into:
        if (nested_stack_) {
          return fmt::format(fmt::runtime(TIP_("Move {} into group")), drag_name);
        }
        ATTR_FALLTHROUGH;
      case DropLocation::Before:
        return fmt::format(fmt::runtime(TIP_("Move {} above")), drag_name);
      case DropLocation::After:
        return fmt::format(fmt::runtime(TIP_("Move {} below")), drag_name);
    }
    return "";
  }

  /* Dropping ONTO a layer row adds mask-capable assets to that layer's mask
   * stack; the base (last) layer has no mask socket, so it only takes layer
   * inserts. */
  bool asset_drops_as_mask(const wmDragAsset &asset_drag, const DropLocation location) const
  {
    if (location != DropLocation::Into) {
      return false;
    }
    const NodeShaderLayerStack &storage = layer_stack::storage(*parent_stack_);
    if (parent_index_ >= storage.items_num - 1) {
      return false;
    }
    return (ed::render::asset_texture_layer_usage(*asset_drag.asset) &
            SHADER_NODE_TREE_USAGE_MASK_GENERATOR) != 0;
  }

  bool drop_asset(bContext *C, const wmDragAsset &asset_drag, const DragInfo &drag_info) const
  {
    wmOperatorType *ot = WM_operatortype_find("MATERIAL_OT_texture_layer_add_node_group", false);
    if (ot == nullptr) {
      return false;
    }
    PointerRNA props = WM_operator_properties_create_ptr(ot);
    ed::asset::operator_asset_reference_props_set(*asset_drag.asset, props);
    if (this->asset_drops_as_mask(asset_drag, drag_info.drop_location)) {
      /* The mask goes to the active layer: make this row's layer active. */
      ed::render::set_active_layer(ntree_, *parent_stack_, parent_index_);
      RNA_boolean_set(&props, "as_mask", true);
    }
    else {
      /* Activate the target stack so the operator poll (which resolves the active stack)
       * succeeds even when an unrelated node was active at drop time. */
      bke::node_set_active(ntree_, *parent_stack_);
      const int to_index = (drag_info.drop_location == DropLocation::After) ? parent_index_ + 1 :
                                                                              parent_index_;
      RNA_string_set(&props, "to_stack", parent_stack_->name);
      RNA_int_set(&props, "to_index", to_index);
    }
    WM_operator_name_call_ptr(C, ot, wm::OpCallContext::ExecDefault, &props, nullptr);
    WM_operator_properties_free(&props);
    return true;
  }

  bool on_drop(bContext *C, const DragInfo &drag_info) const override
  {
    if (const wmDragAsset *asset_drag = asset_drag_data(drag_info.drag_data)) {
      return this->drop_asset(C, *asset_drag, drag_info);
    }
    const wmDragTextureLayer *data = static_cast<const wmDragTextureLayer *>(
        drag_info.drag_data.poin);
    if (data == nullptr || data->stack_node == nullptr) {
      return false;
    }

    bNode *to_stack = nullptr;
    int to_index = 0;
    if (drag_info.drop_location == DropLocation::Into && nested_stack_) {
      /* Drop INTO this group: place at the top of the nested stack. */
      to_stack = nested_stack_;
      to_index = 0;
    }
    else {
      /* Reorder around this row in the parent stack. */
      to_stack = parent_stack_;
      to_index = parent_index_;
      if (drag_info.drop_location == DropLocation::After) {
        to_index = parent_index_ + 1;
      }
    }
    /* Within the same stack: removing the source first shifts later indices
     * down by one — compensate so the row lands at the visually intended spot. */
    if (to_stack == data->stack_node && data->item_index < to_index) {
      to_index -= 1;
    }

    wmOperatorType *ot = WM_operatortype_find("MATERIAL_OT_texture_layer_reparent", false);
    if (ot == nullptr) {
      return false;
    }
    PointerRNA props = WM_operator_properties_create_ptr(ot);
    RNA_string_set(&props, "from_stack", data->stack_node->name);
    RNA_int_set(&props, "from_index", data->item_index);
    RNA_string_set(&props, "to_stack", to_stack->name);
    RNA_int_set(&props, "to_index", to_index);
    WM_operator_name_call_ptr(C, ot, wm::OpCallContext::ExecDefault, &props, nullptr);
    WM_operator_properties_free(&props);
    return true;
  }
};

std::unique_ptr<AbstractViewItemDragController> ShaderLayerItem::create_drag_controller() const
{
  if (kind_ == EntryKind::Bsdf) {
    return nullptr;
  }
  /* Layer and mask rows share the drag payload; the stack node's type tells
   * the drop target which kind is being dragged. */
  return std::make_unique<ShaderLayerDragController>(get_tree_view(), ntree_, *node_, item_index_);
}

std::unique_ptr<TreeViewItemDropTarget> ShaderLayerItem::create_drop_target()
{
  if (kind_ == EntryKind::Bsdf) {
    return nullptr;
  }
  if (kind_ == EntryKind::Layer) {
    /* Reorder before/after the row, AND drop INTO it: into the nested stack
     * for groups (layer drags), or as a mask on the layer (asset drags). */
    return std::make_unique<ShaderLayerDropTarget>(
        *this, DropBehavior::ReorderAndInsert, ntree_, node_, item_index_, nested_stack_);
  }
  return std::make_unique<ShaderLayerDropTarget>(
      *this, DropBehavior::Reorder, ntree_, node_, item_index_, nullptr);
}

class ShaderLayersTreeView : public AbstractTreeView {
 private:
  bNodeTree &ntree_;
  Vector<bNode *> bsdfs_;
  Map<bNode *, Vector<bNode *>> bsdf_roots_;

 public:
  ShaderLayersTreeView(bNodeTree &ntree,
                       Vector<bNode *> bsdfs,
                       Map<bNode *, Vector<bNode *>> bsdf_roots)
      : ntree_(ntree), bsdfs_(std::move(bsdfs)), bsdf_roots_(std::move(bsdf_roots))
  {
  }

  /* Add a row for each layer in #stack_node with its mask rows as children,
   * recursing into group layers whose Bundle input is fed by another
   * #ShaderNodeTextureLayerStack. */
  void add_stack_items(TreeViewOrItem &parent, bNode &stack_node)
  {
    const NodeShaderLayerStack &storage = layer_stack::storage(stack_node);
    /* An adjustment group's nested stack has a hidden base "input" item that
     * carries the bundle from below the group; it is not a user layer. */
    const bool adjustment_group = texture_stack::is_adjustment_group(ntree_, stack_node);
    for (const int i : IndexRange(storage.items_num)) {
      if (adjustment_group && i == storage.items_num - 1) {
        continue;
      }
      ShaderLayerItem &child = parent.add_tree_item<ShaderLayerItem>(
          ntree_, EntryKind::Layer, stack_node, i);
      if (bNode *mask_stack_node =
              texture_stack::StackLayer{&ntree_, &stack_node, i}.mask_stack_node())
      {
        child.uncollapse_by_default();
        const NodeShaderLayerStack &mask_storage = layer_stack::storage(*mask_stack_node);
        for (const int mask_i : IndexRange(mask_storage.items_num)) {
          child.add_tree_item<ShaderLayerItem>(ntree_, EntryKind::Mask, *mask_stack_node, mask_i);
        }
      }
      if (child.is_group()) {
        child.uncollapse_by_default();
        add_stack_items(child, *child.nested_stack());
      }
    }
  }

  void build_tree() override
  {
    for (bNode *bsdf : bsdfs_) {
      ShaderLayerItem &bsdf_item = this->add_tree_item<ShaderLayerItem>(ntree_, *bsdf);
      bsdf_item.uncollapse_by_default();

      const Vector<bNode *> &roots = bsdf_roots_.lookup(bsdf);
      for (bNode *stack_node : roots) {
        add_stack_items(bsdf_item, *stack_node);
      }
    }
  }

  bool listen(const wmNotifier &notifier) const override
  {
    if (notifier.category == NC_NODE) {
      return true;
    }
    return notifier.category == NC_MATERIAL && notifier.data == ND_SHADING;
  }
};

/* Find the (BSDF, stack_node, item_index) for the active layer or mask
 * entry, if any. */
struct ActiveLayer {
  bNode *bsdf = nullptr;
  /* Texture Layer Stack node, or Mask Stack node when #is_mask is set. */
  bNode *stack_node = nullptr;
  int item_index = -1;
  bool is_mask = false;
};

ActiveLayer resolve_active_layer(bNodeTree &ntree, const Span<bNode *> bsdfs, bNode *active)
{
  ActiveLayer result;
  if (active == nullptr) {
    return result;
  }
  if (bsdfs.contains(active)) {
    result.bsdf = active;
    return result;
  }
  /* The layer the active node belongs to: one of its nodes, or the stack node
   * itself. The stack may be any reachable one — top-level, nested inside a
   * group, or a layer's mask stack. */
  if (const std::optional<texture_stack::StackLayer> layer = texture_stack::layer_for_node(
          ntree, *active))
  {
    result.stack_node = layer->stack;
    result.is_mask = mask_stack::is_stack(*layer->stack);
    result.item_index = layer->index;
    result.bsdf = texture_stack::bsdf_for(ntree, *layer->stack);
    if (result.bsdf == nullptr && !bsdfs.is_empty()) {
      /* Unreachable stack (still being wired up): fall back to the first
       * BSDF so the panel has something to associate with the layer. */
      result.bsdf = bsdfs[0];
    }
  }
  return result;
}

}  // namespace

void template_shader_layers(Layout *layout, bContext *C, PointerRNA *ptr)
{
  if (ptr == nullptr || ptr->data == nullptr) {
    return;
  }
  bNodeTree *ntree = static_cast<bNodeTree *>(ptr->data);
  ntree->ensure_topology_cache();

  Vector<bNode *> bsdfs;
  texture_channel::collect_bsdfs(*ntree, bsdfs);
  if (bsdfs.is_empty()) {
    layout->label(IFACE_("No BSDF nodes"), ICON_NONE);
    return;
  }

  Map<bNode *, Vector<bNode *>> bsdf_roots;
  for (bNode *bsdf : bsdfs) {
    Vector<bNode *> roots;
    texture_stack::roots_for_bsdf(*ntree, *bsdf, roots);
    bsdf_roots.add(bsdf, std::move(roots));
  }

  /* Tree view spanning the full width, with the add/remove/menu buttons in its
   * bottom bar. The layer details below also span the full layout width. */
  Layout &tree_col = layout->column(false);
  Block *block = tree_col.block();

  /* With no stack yet, the (+) button creates one immediately with a default
   * Fill layer rather than opening the add menu. */
  bool has_stack = false;
  for (const Vector<bNode *> &roots : bsdf_roots.values()) {
    if (!roots.is_empty()) {
      has_stack = true;
      break;
    }
  }

  AbstractTreeView *tree_view = block_add_view(
      *block,
      "Shader Layers Tree View",
      std::make_unique<ShaderLayersTreeView>(*ntree, bsdfs, bsdf_roots));
  tree_view->set_default_rows(5);
  tree_view->set_flat_box(true);
  tree_view->set_footer_fn([has_stack](Layout &layout) {
    if (has_stack) {
      layout.menu("MATERIAL_MT_texture_layer_add", "", ICON_ADD);
    }
    else {
      layout.op("MATERIAL_OT_texture_layer_add_default", "", ICON_ADD);
    }
    layout.op("MATERIAL_OT_texture_layer_remove", "", ICON_REMOVE);
    layout.separator();
    layout.menu("MATERIAL_MT_texture_layer_context", "", ICON_DOWNARROW_HLT);
  });
  TreeViewBuilder::build_tree_view(*C, *tree_view, tree_col);

  /* Inline details for the selected tree entry, rendered below the tree. */
  bNode *active = bke::node_get_active(*ntree);
  const ActiveLayer active_layer = resolve_active_layer(*ntree, bsdfs, active);
  if (active_layer.bsdf == nullptr) {
    return;
  }

  layout->separator();
  layout->use_property_split_set(true);
  layout->use_property_decorate_set(false);

  if (active_layer.stack_node == nullptr) {
    /* The BSDF is selected: show its full input list, with the channel
     * toggles that route inputs through the texture layer stack. */
    PointerRNA bsdf_ptr = RNA_pointer_create_discrete(&ntree->id, RNA_Node, active_layer.bsdf);
    draw_node_inputs(*layout, C, &bsdf_ptr, /*bsdf_channels=*/true);
    return;
  }

  /* A layer or mask item is selected. Both draw the same sections: the source
   * node's properties, its texture mapping, an Influence panel with blend mode,
   * mute and opacity, and the node properties left out of the sections above. */
  bNode *stack_node = active_layer.stack_node;
  const bool is_mask = active_layer.is_mask;
  const int item_index = active_layer.item_index;
  const NodeShaderLayerStack &storage = layer_stack::storage(*stack_node);

  /* The node whose properties stand for the row: a mask's source (e.g. an Image
   * Texture for a paint mask), or the layer's source. For closure layers the
   * latter resolves through the zone to the adjustment group, so the user edits
   * the group's inputs directly. A group's source is its nested stack; its
   * contents are redundant with the Influence settings on the group's children
   * in the tree, so it has none. */
  bNode *source;
  if (is_mask) {
    source = mask_stack::layer_source(*ntree, *stack_node, item_index);
  }
  else {
    const texture_stack::StackLayer layer{ntree, stack_node, item_index};
    source = (layer.nested_stack() == nullptr) ? layer.display_source() : nullptr;
  }

  if (source) {
    if (!is_mask && source->is_type("NodeCombineBundle"_ustr)) {
      /* Generator wrapped in a Combine Bundle: channel constants + the
       * generator node's inputs. */
      const NodeShaderLayerStackItem *item = (item_index >= 0 && item_index < storage.items_num) ?
                                                 &storage.items[item_index] :
                                                 nullptr;
      draw_combine_bundle_source(*layout, C, *ntree, *source, item);
    }
    else {
      PointerRNA src_ptr = RNA_pointer_create_discrete(&ntree->id, RNA_Node, source);
      draw_node_inputs(*layout, C, &src_ptr);
    }

    /* Texture Mapping: the built-in texture coordinate mapping of the texture
     * nodes behind the row, above the Influence panel. */
    Vector<bNode *> mapping_nodes;
    collect_texture_mapping_nodes(*ntree, *source, mapping_nodes);
    draw_texture_mapping_panel(*layout, C, *ntree, mapping_nodes);
  }

  draw_influence_panel(*layout, C, *ntree, *stack_node, item_index, *active_layer.bsdf, is_mask);

  /* Advanced: the node properties left out above, at the very bottom. */
  if (source) {
    Vector<bNode *> advanced_nodes;
    collect_advanced_nodes(*ntree, *source, advanced_nodes);
    draw_advanced_panel(*layout, C, *ntree, advanced_nodes);
  }
}

}  // namespace blender::ui
