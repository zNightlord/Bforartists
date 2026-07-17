/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edrend
 *
 * Add menus + operator for inserting Texture Layer Stack layers driven by
 * shader node group assets, filtered by their #eShaderNodeTreeUsage flag.
 */

#include "AS_asset_library.hh"
#include "AS_asset_representation.hh"

#include "BLI_listbase_iterator.hh"
#include "BLI_string_utf8.hh"
#include "BLI_vector.hh"

#include <algorithm>

#include "DNA_material_types.h"
#include "DNA_node_types.h"
#include "DNA_screen_types.h"

#include "BKE_asset.hh"
#include "BKE_context.hh"
#include "BKE_global.hh"
#include "BKE_idprop.hh"
#include "BKE_lib_id.hh"
#include "BKE_main.hh"
#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_runtime.hh"
#include "BKE_node_tree_update.hh"
#include "BKE_report.hh"
#include "BKE_screen.hh"

#include "BLT_translation.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"
#include "RNA_prototypes.hh"

#include "ED_asset.hh"
#include "ED_asset_menu_utils.hh"
#include "ED_screen.hh"

#include "NOD_closure_zone.hh"
#include "NOD_geo_bundle.hh"
#include "NOD_layer_stack.hh"
#include "NOD_mask_stack.hh"
#include "NOD_node_placement.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"

#include "UI_interface.hh"
#include "UI_interface_layout.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "render_intern.hh"

namespace blender::ed::render {

using asset_system::AssetRepresentation;
using nodes::CombineBundleItemsAccessor;

namespace layer_stack = nodes::layer_stack;
namespace texture_stack = nodes::texture_stack;
namespace mask_stack = nodes::mask_stack;
namespace texture_channel = nodes::texture_channel;
namespace closure_zone = nodes::closure_zone;

/* -------------------------------------------------------------------- */
/** \name Asset Iteration / Menu Drawing
 * \{ */

/* #asset_texture_layer_usage is defined with the tree-view drop code in
 * material_texture_layers_ui.cc, which also consumes it. */

/* The blend mode (MA_RAMP_*) a mask driven by this asset gets by default,
 * from #bNodeTree::default_mask_blend baked into the asset metadata. */
static int asset_mask_blend_type(const AssetRepresentation &asset)
{
  const IDProperty *blend = BKE_asset_metadata_idprop_find(&asset.get_metadata(),
                                                           "mask_blend_type");
  if (blend == nullptr || blend->type != IDP_INT) {
    return MA_RAMP_BLEND;
  }
  return IDP_int_get(blend);
}

static void draw_assets_for_usage(const bContext *C,
                                  ui::Layout *layout,
                                  const int usage_flag,
                                  const bool as_mask)
{
  const AssetLibraryReference library_ref = asset_system::all_library_reference();
  asset::list::storage_fetch(&library_ref, C);
  if (!asset::list::library_get_once_available(library_ref)) {
    layout->label(IFACE_("Loading Asset Libraries"), ICON_INFO);
    return;
  }

  wmOperatorType *ot = WM_operatortype_find("MATERIAL_OT_texture_layer_add_node_group", true);
  if (!ot) {
    return;
  }

  /* Collect and name-sort the matching assets, so the menu order is stable across library
   * reloads (matching the built-in node entries above). */
  Vector<AssetRepresentation *> assets;
  asset::list::iterate(library_ref, [&](AssetRepresentation &asset) {
    if (asset.get_id_type() != ID_NT) {
      return true;
    }
    if ((asset_texture_layer_usage(asset) & usage_flag) == 0) {
      return true;
    }
    assets.append(&asset);
    return true;
  });
  std::sort(assets.begin(),
            assets.end(),
            [](const AssetRepresentation *a, const AssetRepresentation *b) {
              return StringRef(a->get_name()) < StringRef(b->get_name());
            });
  for (AssetRepresentation *asset : assets) {
    PointerRNA props = layout->op(
        ot, IFACE_(asset->get_name()), ICON_NONE, wm::OpCallContext::InvokeDefault, UI_ITEM_NONE);
    asset::operator_asset_reference_props_set(*asset, props);
    /* Mask menu invocations: tell the operator to route the asset's output
     * into the active layer's mask stack instead of inserting a new layer. */
    RNA_boolean_set(&props, "as_mask", as_mask);
    /* Adjustment menu: a dual-usage asset should be inserted as an adjustment, not a generator. */
    RNA_boolean_set(
        &props, "as_adjustment", usage_flag == SHADER_NODE_TREE_USAGE_TEXTURE_ADJUSTMENT);
  }
  if (assets.is_empty()) {
    layout->label(IFACE_("No assets available"), ICON_INFO);
  }
}

/* Entries for built-in shader node types whose #bNodeType::texture_layer_usage
 * carries #usage_flag, drawn above the node group assets. */
static void draw_builtin_nodes_for_usage(ui::Layout *layout,
                                         const int usage_flag,
                                         const bool as_mask)
{
  wmOperatorType *ot = WM_operatortype_find("MATERIAL_OT_texture_layer_add_node", true);
  if (!ot) {
    return;
  }
  Vector<const bke::bNodeType *> types;
  for (const bke::bNodeType *ntype : bke::node_types_get()) {
    if (ntype->texture_layer_usage & usage_flag) {
      types.append(ntype);
    }
  }
  std::sort(types.begin(), types.end(), [](const bke::bNodeType *a, const bke::bNodeType *b) {
    return a->ui_name < b->ui_name;
  });
  for (const bke::bNodeType *ntype : types) {
    PointerRNA props = layout->op(ot,
                                  IFACE_(ntype->ui_name.c_str()),
                                  ntype->ui_icon,
                                  wm::OpCallContext::InvokeDefault,
                                  UI_ITEM_NONE);
    RNA_string_set(&props, "type", ntype->idname.c_str());
    RNA_boolean_set(&props, "as_mask", as_mask);
  }
  if (!types.is_empty()) {
    layout->separator();
  }
}

static void generator_assets_draw(const bContext *C, Menu *menu)
{
  draw_builtin_nodes_for_usage(menu->layout, SHADER_NODE_TREE_USAGE_TEXTURE_GENERATOR, false);
  draw_assets_for_usage(C, menu->layout, SHADER_NODE_TREE_USAGE_TEXTURE_GENERATOR, false);
}

static void adjustment_assets_draw(const bContext *C, Menu *menu)
{
  draw_builtin_nodes_for_usage(menu->layout, SHADER_NODE_TREE_USAGE_TEXTURE_ADJUSTMENT, false);
  draw_assets_for_usage(C, menu->layout, SHADER_NODE_TREE_USAGE_TEXTURE_ADJUSTMENT, false);
}

static void mask_assets_draw(const bContext *C, Menu *menu)
{
  draw_builtin_nodes_for_usage(menu->layout, SHADER_NODE_TREE_USAGE_MASK_GENERATOR, true);
  draw_assets_for_usage(C, menu->layout, SHADER_NODE_TREE_USAGE_MASK_GENERATOR, true);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add Layer From Asset Operator
 * \{ */

static bool add_node_group_poll(bContext *C)
{
  return resolve_active_stack(*C).has_value();
}

static wmOperatorStatus add_mask_node_group_exec(bContext *C,
                                                 bNodeTree &ntree,
                                                 bNodeTree &group,
                                                 const StringRefNull name,
                                                 const int default_blend_type)
{
  Main &bmain = *CTX_data_main(C);

  /* The selected mask's own stack, or the active layer's (created on demand). */
  bNode *mask_stack_node = ensure_active_layer_mask_stack(*C, ntree);
  if (!mask_stack_node) {
    return OPERATOR_CANCELLED;
  }

  /* Add the asset group node. */
  bNode *group_node = bke::node_add_node(C, ntree, "ShaderNodeGroup"_ustr);
  if (!group_node) {
    return OPERATOR_CANCELLED;
  }
  group_node->id = &group.id;
  id_us_plus(&group.id);
  STRNCPY_UTF8(group_node->label, name.c_str());
  BKE_ntree_update_tag_node_property(&ntree, group_node);
  BKE_main_ensure_invariants(bmain, ntree.id);

  /* Add a new mask item above the active one, fed by the asset's mask output,
   * with the blend mode the asset declares (e.g. Multiply for occlusion). Name it after the asset
   * (the group datablock name may carry a ".001" suffix from import). */
  int new_index = -1;
  bNode *mask_source = nullptr;
  if (bNodeSocket *mask_out = mask_stack::find_source_output(
          bmain, ntree, *group_node, &mask_source))
  {
    new_index = mask_stack::add_layer_from_source(
        bmain, ntree, *mask_stack_node, name.c_str(), *mask_source, *mask_out);
  }
  else {
    new_index = mask_stack::add_layer(bmain, ntree, *mask_stack_node, name.c_str(), std::nullopt);
  }
  mask_stack::set_layer_blend(*mask_stack_node, new_index, default_blend_type);
  BKE_main_ensure_invariants(bmain, ntree.id);
  /* Placed once the mask layer exists, so it lands in that layer's slot. */
  texture_stack::place_layer_source(*group_node, {&ntree, mask_stack_node, new_index}, 320.0f);
  if (mask_source && mask_source != group_node) {
    /* The Separate Bundle extracting the mask channel sits between them. */
    nodes::place_node_left_of(*mask_source, *group_node, -200.0f, 0.0f);
  }
  set_active_layer(ntree, *mask_stack_node, new_index);
  return OPERATOR_FINISHED;
}

static wmOperatorStatus add_node_group_exec(bContext *C, wmOperator *op)
{
  Material *mat = active_material(*C);
  bNodeTree *ntree_ptr = mat ? mat->nodetree : nullptr;
  if (mat == nullptr || ntree_ptr == nullptr) {
    BKE_report(op->reports, RPT_ERROR, "No active material node tree");
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ntree_ptr;

  if (G.background) {
    /* In background/test mode, asset library loading is asynchronous by default.
     * Force a blocking fetch so the operator can resolve the asset. */
    asset::list::storage_fetch_blocking(asset_system::all_library_reference(), *C);
  }

  const AssetRepresentation *asset =
      asset::operator_asset_reference_props_get_asset_from_all_library(*C, *op->ptr, op->reports);
  if (!asset) {
    return OPERATOR_CANCELLED;
  }

  const bool as_mask = RNA_boolean_get(op->ptr, "as_mask");

  /* Resolve the target stack for a layer insert *before* importing the asset. A tree-view drop
   * passes an explicit target stack/index; a menu invocation uses the active stack. Doing this
   * first means a stale target (node renamed/deleted mid-drop) fails without importing a 0-user
   * group, and no longer requires an active stack when the target is explicit. */
  bNode *stack_node = nullptr;
  if (!as_mask) {
    char to_stack_name[MAX_NAME];
    RNA_string_get(op->ptr, "to_stack", to_stack_name);
    if (to_stack_name[0]) {
      stack_node = bke::node_find_node_by_name(ntree, to_stack_name);
      if (stack_node == nullptr || !texture_stack::is_stack(*stack_node)) {
        BKE_report(op->reports, RPT_ERROR, "Target texture layer stack not found");
        return OPERATOR_CANCELLED;
      }
    }
    else {
      const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
      if (!ctx) {
        BKE_report(op->reports, RPT_ERROR, "No active Texture Layer Stack");
        return OPERATOR_CANCELLED;
      }
      stack_node = ctx->stack;
    }
    /* A mask-only asset (no generator/adjustment bundle output) has nothing to feed a new layer's
     * Bundle input: reject the insert rather than create a dead empty layer plus an orphan node.
     */
    if ((asset_texture_layer_usage(*asset) &
         (SHADER_NODE_TREE_USAGE_TEXTURE_GENERATOR | SHADER_NODE_TREE_USAGE_TEXTURE_ADJUSTMENT)) ==
        eShaderNodeTreeUsage(0))
    {
      BKE_report(op->reports, RPT_ERROR, "Asset cannot be inserted as a texture layer");
      return OPERATOR_CANCELLED;
    }
  }

  Main &bmain = *CTX_data_main(C);
  bNodeTree *group = reinterpret_cast<bNodeTree *>(
      asset::asset_local_id_ensure_imported(bmain, *asset));
  if (!group) {
    return OPERATOR_CANCELLED;
  }
  if (group->type != NTREE_SHADER) {
    BKE_report(op->reports, RPT_ERROR, "Asset is not a shader node group");
    return OPERATOR_CANCELLED;
  }

  /* Mask-routing path: don't insert a new layer; route the asset's output
   * through the active layer's mask stack. */
  if (as_mask) {
    const wmOperatorStatus status = add_mask_node_group_exec(
        C, ntree, *group, asset->get_name(), asset_mask_blend_type(*asset));
    if (status == OPERATOR_FINISHED) {
      tree_changed(*C, ntree, mat);
    }
    return status;
  }

  const int to_index = RNA_int_get(op->ptr, "to_index");

  /* Add the new layer; its dynamic Bundle input socket is rebuilt by the tree
   * update inside the helper before we link to it. */
  int new_item_index;
  if (to_index >= 0) {
    new_item_index = layer_stack::add_layer_at(
        bmain, ntree, *stack_node, asset->get_name().c_str(), to_index);
    layer_stack::storage(*stack_node).active_index = new_item_index;
    bke::node_set_active(ntree, *stack_node);
  }
  else {
    new_item_index = texture_stack::add_layer(
        bmain, ntree, *stack_node, asset->get_name().c_str());
  }

  /* Add the asset's node group to the material tree and link its bundle output
   * to the new layer's Bundle input. */
  bNode *group_node = bke::node_add_node(C, ntree, "ShaderNodeGroup"_ustr);
  if (!group_node) {
    return OPERATOR_CANCELLED;
  }
  group_node->id = &group->id;
  id_us_plus(&group->id);
  STRNCPY_UTF8(group_node->label, asset->get_name().c_str());
  BKE_ntree_update_tag_node_property(&ntree, group_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
  /* Placed once its sockets exist, so its height is known. */
  texture_stack::place_layer_source(*group_node, {&ntree, stack_node, new_item_index}, 280.0f);

  const eShaderNodeTreeUsage usage_flags = asset_texture_layer_usage(*asset);
  /* Wrap in a closure zone (adjustment semantics) when the asset acts as an adjustment: either the
   * user invoked it from the Adjustment menu, or it is adjustment-only (no generator role). A
   * dual-usage asset added from the Generator menu is treated as a generator. */
  const bool as_adjustment = RNA_boolean_get(op->ptr, "as_adjustment");
  const bool use_adjustment = (usage_flags & SHADER_NODE_TREE_USAGE_TEXTURE_ADJUSTMENT) &&
                              (as_adjustment ||
                               !(usage_flags & SHADER_NODE_TREE_USAGE_TEXTURE_GENERATOR));
  if (use_adjustment) {
    /* Adjustments transform the stack accumulated below their layer: wrap the
     * group in a closure zone whose Closure output feeds the layer. */
    closure_zone::wrap_group(bmain, ntree, *group_node, *stack_node, new_item_index);
  }
  else {
    bNodeSocket *layer_in =
        texture_stack::StackLayer{&ntree, stack_node, new_item_index}.bundle_socket();
    bNodeSocket *bundle_out = nodes::first_bundle_output(*group_node);
    if (layer_in && bundle_out) {
      bke::node_add_link(ntree, *group_node, *bundle_out, *stack_node, *layer_in);
    }
  }
  set_active_layer(ntree, *stack_node, new_item_index);
  tree_changed(*C, ntree, mat);
  return OPERATOR_FINISHED;
}

static std::string add_node_group_get_description(bContext *C,
                                                  wmOperatorType * /*ot*/,
                                                  PointerRNA *ptr)
{
  const AssetRepresentation *asset =
      asset::operator_asset_reference_props_get_asset_from_all_library(*C, *ptr, nullptr);
  if (!asset || !asset->get_metadata().description) {
    return "";
  }
  return TIP_(asset->get_metadata().description);
}

static void MATERIAL_OT_texture_layer_add_node_group(wmOperatorType *ot)
{
  ot->name = "Add Texture Layer from Asset";
  ot->description = "Add a new texture layer driven by a shader node group asset";
  ot->idname = "MATERIAL_OT_texture_layer_add_node_group";

  ot->exec = add_node_group_exec;
  ot->poll = add_node_group_poll;
  ot->get_description = add_node_group_get_description;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  asset::operator_asset_reference_props_register(*ot->srna);

  PropertyRNA *prop = RNA_def_boolean(
      ot->srna,
      "as_mask",
      false,
      "As Mask",
      "Route the asset's output through the active layer's mask stack instead "
      "of inserting a new layer");
  RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);

  /* Set by the Adjustment menu so a dual-usage (generator + adjustment) asset is wrapped in a
   * closure zone rather than inserted as a plain generator. */
  prop = RNA_def_boolean(ot->srna,
                         "as_adjustment",
                         false,
                         "As Adjustment",
                         "Insert the asset as an adjustment layer (closure zone) rather than a "
                         "generator");
  RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);

  /* Placement for drag & drop onto the layers tree view. */
  prop = RNA_def_string(
      ot->srna, "to_stack", nullptr, MAX_NAME, "Target Stack", "Stack node to insert into");
  RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);
  prop = RNA_def_int(ot->srna,
                     "to_index",
                     -1,
                     -1,
                     INT32_MAX,
                     "Target Index",
                     "Layer index to insert at",
                     -1,
                     INT32_MAX);
  RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add Layer From Built-in Node Operator
 * \{ */

static wmOperatorStatus add_node_exec(bContext *C, wmOperator *op)
{
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    BKE_report(op->reports, RPT_ERROR, "No active Texture Layer Stack");
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode *stack_node = ctx->stack;
  Main &bmain = *CTX_data_main(C);

  char type_idname[128];
  RNA_string_get(op->ptr, "type", type_idname);

  /* Mask-routing path: no new layer; feed the node's float output through the
   * active layer's mask stack. */
  if (RNA_boolean_get(op->ptr, "as_mask")) {
    bNode *mask_stack_node = ensure_active_layer_mask_stack(*C, ntree);
    if (mask_stack_node == nullptr) {
      return OPERATOR_CANCELLED;
    }
    bNode *node = bke::node_add_node(C, ntree, UString(type_idname));
    if (node == nullptr) {
      return OPERATOR_CANCELLED;
    }
    BKE_ntree_update_tag_node_property(&ntree, node);
    BKE_main_ensure_invariants(bmain, ntree.id);
    bNode *mask_source = nullptr;
    if (bNodeSocket *mask_out = mask_stack::find_source_output(bmain, ntree, *node, &mask_source))
    {
      const int new_index = mask_stack::add_layer_from_source(bmain,
                                                              ntree,
                                                              *mask_stack_node,
                                                              node->typeinfo->ui_name.c_str(),
                                                              *mask_source,
                                                              *mask_out);
      /* The node type's default mask blend, e.g. Multiply for AO. */
      mask_stack::set_layer_blend(
          *mask_stack_node, new_index, node->typeinfo->texture_layer_mask_blend);
      /* Placed once the mask layer exists, so it lands in that layer's slot. */
      texture_stack::place_layer_source(*node, {&ntree, mask_stack_node, new_index}, 320.0f);
      if (mask_source && mask_source != node) {
        /* The Separate Bundle extracting the mask channel sits between them. */
        nodes::place_node_left_of(*mask_source, *node, -200.0f, 0.0f);
      }
      set_active_layer(ntree, *mask_stack_node, new_index);
    }
    tree_changed(*C, ntree, ctx->material);
    return OPERATOR_FINISHED;
  }

  /* Texture layer path: the node's output feeds a Combine Bundle keyed to the
   * chosen channel, which becomes the new layer's bundle. */
  char channel[MAX_NAME];
  RNA_string_get(op->ptr, "channel", channel);
  const eNodeSocketDatatype channel_type = texture_channel::socket_type(
      ntree, *stack_node, channel);

  bNode *node = bke::node_add_node(C, ntree, UString(type_idname));
  if (node == nullptr) {
    return OPERATOR_CANCELLED;
  }
  const int new_item_index = texture_stack::add_layer(
      bmain, ntree, *stack_node, node->typeinfo->ui_name.c_str());

  bNode *combine = bke::node_add_node(C, ntree, "NodeCombineBundle"_ustr);
  if (combine == nullptr) {
    return OPERATOR_CANCELLED;
  }
  NodeCombineBundleItem *item =
      nodes::socket_items::add_item_with_socket_type_and_name<CombineBundleItemsAccessor>(
          ntree, *combine, channel_type, channel);
  BKE_ntree_update_tag_node_property(&ntree, node);
  BKE_ntree_update_tag_node_property(&ntree, combine);
  BKE_main_ensure_invariants(bmain, ntree.id);
  /* Placed once their sockets exist, so their heights are known. */
  texture_stack::place_layer_source(*combine, {&ntree, stack_node, new_item_index}, 280.0f);
  nodes::place_node_left_of(*node, *combine, 250.0f, 0.0f);

  texture_stack::connect_bundle({&ntree, stack_node, new_item_index}, *combine, "Bundle");
  if (item) {
    const std::string in_id = CombineBundleItemsAccessor::socket_identifier_for_item(*item);
    bNodeSocket *combine_in = bke::node_find_socket(*combine, SOCK_IN, UString(in_id));
    bNodeSocket *source_out = texture_stack::preferred_source_output(*node, channel_type);
    if (combine_in && source_out) {
      bke::node_add_link(ntree, *node, *source_out, *combine, *combine_in);
    }
  }

  /* Route the channel out of the stack to the consuming BSDF when it is
   * not already: add the Separate Bundle output for the channel and link it to
   * the BSDF input. Without this the new layer feeds the stack, but the channel
   * never reaches the BSDF (e.g. adding a Base Color generator to a stack that
   * only routed Metallic). */
  if (bNode *bsdf = texture_stack::bsdf_for(ntree, *stack_node)) {
    for (bNodeSocket &bsdf_socket : bsdf->inputs) {
      if (!texture_channel::is_fillable_input(bsdf_socket) ||
          StringRef(bsdf_socket.name) != channel)
      {
        continue;
      }
      /* Inherit the BSDF channel's slider/range onto the combine item. */
      if (item) {
        nodes::combine_bundle_item_copy_socket_data(*item, bsdf_socket);
        BKE_ntree_update_tag_node_property(&ntree, combine);
      }
      if (texture_channel::state(ntree, bsdf_socket) == texture_channel::State::Free) {
        if (bNode *separate = texture_channel::separate_bundle_for_bsdf(ntree, *bsdf)) {
          texture_channel::link(ntree, *separate, *bsdf, bsdf_socket);
        }
      }
      break;
    }
  }

  set_active_layer(ntree, *stack_node, new_item_index);
  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

/* Generator adds ask which channel the node output should drive: a popup menu
 * listing the consuming BSDF's channels (its color/float/vector inputs), each
 * entry re-invoking the operator with the channel set. Mask adds and calls
 * with an explicit channel skip the menu. */
static wmOperatorStatus add_node_invoke(bContext *C, wmOperator *op, const wmEvent * /*event*/)
{
  if (RNA_boolean_get(op->ptr, "as_mask") || RNA_struct_property_is_set(op->ptr, "channel")) {
    return add_node_exec(C, op);
  }
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNode *bsdf = texture_stack::bsdf_for(*ctx->ntree, *ctx->stack);
  if (bsdf == nullptr) {
    /* No BSDF consuming the stack yet: fall back to the default channel. */
    return add_node_exec(C, op);
  }

  char type_idname[128];
  RNA_string_get(op->ptr, "type", type_idname);

  /* Only offer the channels currently routed through the stack (the Textured
   * subset), not every fillable BSDF input. Adding a generator layers onto an
   * existing channel, so the wider universe would only clutter the menu. */
  Vector<bNodeSocket *> channels;
  for (bNodeSocket &sock : bsdf->inputs) {
    if (texture_channel::is_fillable_input(sock) &&
        texture_channel::state(*ctx->ntree, sock) == texture_channel::State::Textured)
    {
      channels.append(&sock);
    }
  }
  if (channels.is_empty()) {
    /* No channel routed yet: fall back to the default channel. */
    return add_node_exec(C, op);
  }

  ui::PopupMenu *pup = ui::popup_menu_begin(C, IFACE_("Channel"), ICON_NONE);
  ui::Layout &layout = *ui::popup_menu_layout(pup);
  for (bNodeSocket *sock : channels) {
    PointerRNA props = layout.op(
        op->type, IFACE_(sock->name), ICON_NONE, wm::OpCallContext::ExecDefault, UI_ITEM_NONE);
    RNA_string_set(&props, "type", type_idname);
    RNA_string_set(&props, "channel", sock->name);
  }
  ui::popup_menu_end(C, pup);
  return OPERATOR_INTERFACE;
}

static std::string add_node_get_description(bContext * /*C*/,
                                            wmOperatorType * /*ot*/,
                                            PointerRNA *ptr)
{
  char type_idname[128];
  RNA_string_get(ptr, "type", type_idname);
  const bke::bNodeType *ntype = bke::node_type_find(UString(StringRef(type_idname)));
  if (ntype == nullptr || ntype->ui_description.empty()) {
    return "";
  }
  return TIP_(ntype->ui_description.c_str());
}

static void MATERIAL_OT_texture_layer_add_node(wmOperatorType *ot)
{
  ot->name = "Add Texture Layer from Node";
  ot->description = "Add a new texture layer driven by a built-in shader node";
  ot->idname = "MATERIAL_OT_texture_layer_add_node";

  ot->invoke = add_node_invoke;
  ot->exec = add_node_exec;
  ot->poll = add_node_group_poll;
  ot->get_description = add_node_get_description;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  PropertyRNA *prop = RNA_def_string(
      ot->srna, "type", nullptr, 128, "Node Type", "Shader node type to add as the layer source");
  RNA_def_property_flag(prop, PROP_HIDDEN);
  prop = RNA_def_string(ot->srna,
                        "channel",
                        "Base Color",
                        MAX_NAME,
                        "Channel",
                        "Channel of the BSDF that the node output drives");
  /* Transient per-invocation: without PROP_SKIP_SAVE the OPTYPE_REGISTER last-used value keeps
   * "channel" set, so add_node_invoke's RNA_struct_property_is_set() check skips the channel popup
   * on every later generator add. */
  RNA_def_property_flag(prop, PROP_SKIP_SAVE);
  prop = RNA_def_boolean(ot->srna,
                         "as_mask",
                         false,
                         "As Mask",
                         "Route the node's output through the active layer's mask stack instead "
                         "of inserting a new layer");
  RNA_def_property_flag(prop, PROP_HIDDEN | PROP_SKIP_SAVE);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Registration
 * \{ */

static MenuType make_menu(const char *idname, void (*draw)(const bContext *, Menu *))
{
  MenuType type{};
  STRNCPY_UTF8(type.idname, idname);
  type.draw = draw;
  type.listener = asset::list::asset_reading_region_listen_fn;
  return type;
}

void material_texture_layer_assets_register()
{
  WM_menutype_add(MEM_new<MenuType>(
      __func__,
      make_menu("MATERIAL_MT_texture_layer_add_generator_assets", generator_assets_draw)));
  WM_menutype_add(MEM_new<MenuType>(
      __func__,
      make_menu("MATERIAL_MT_texture_layer_add_adjustment_assets", adjustment_assets_draw)));
  WM_menutype_add(MEM_new<MenuType>(
      __func__, make_menu("MATERIAL_MT_texture_layer_add_mask_assets", mask_assets_draw)));
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_node_group);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_node);
}

/** \} */

}  // namespace blender::ed::render
