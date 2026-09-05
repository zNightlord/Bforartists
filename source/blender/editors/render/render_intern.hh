/* SPDX-FileCopyrightText: 2008 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edrend
 */

#pragma once

#include <optional>

#include "DNA_node_types.h"

namespace blender {

struct ReportList;
struct Scene;
struct ScrArea;
struct bContext;
struct wmOperatorType;

/* `render_shading.cc` */

void OBJECT_OT_material_slot_add(wmOperatorType *ot);
void OBJECT_OT_material_slot_remove(wmOperatorType *ot);
void OBJECT_OT_material_slot_assign(wmOperatorType *ot);
void OBJECT_OT_material_slot_select(wmOperatorType *ot);
void OBJECT_OT_material_slot_deselect(wmOperatorType *ot);
void OBJECT_OT_material_slot_copy(wmOperatorType *ot);
void OBJECT_OT_material_slot_move(wmOperatorType *ot);
void OBJECT_OT_material_slot_remove_unused(wmOperatorType *ot);
void OBJECT_OT_material_slot_remove_all(wmOperatorType *ot);

void MATERIAL_OT_new(wmOperatorType *ot);
void TEXTURE_OT_new(wmOperatorType *ot);
void WORLD_OT_new(wmOperatorType *ot);

void MATERIAL_OT_copy(wmOperatorType *ot);
void MATERIAL_OT_paste(wmOperatorType *ot);

void SCENE_OT_view_layer_add(wmOperatorType *ot);
void SCENE_OT_view_layer_remove(wmOperatorType *ot);
void SCENE_OT_view_layer_add_aov(wmOperatorType *ot);
void SCENE_OT_view_layer_remove_aov(wmOperatorType *ot);
void SCENE_OT_view_layer_add_lightgroup(wmOperatorType *ot);
void SCENE_OT_view_layer_remove_lightgroup(wmOperatorType *ot);
void SCENE_OT_view_layer_add_used_lightgroups(wmOperatorType *ot);
void SCENE_OT_view_layer_remove_unused_lightgroups(wmOperatorType *ot);

void OBJECT_OT_lightprobe_cache_bake(wmOperatorType *ot);
void OBJECT_OT_lightprobe_cache_free(wmOperatorType *ot);

void SCENE_OT_render_view_add(wmOperatorType *ot);
void SCENE_OT_render_view_remove(wmOperatorType *ot);

#ifdef WITH_FREESTYLE
void SCENE_OT_freestyle_module_add(wmOperatorType *ot);
void SCENE_OT_freestyle_module_remove(wmOperatorType *ot);
void SCENE_OT_freestyle_module_move(wmOperatorType *ot);
void SCENE_OT_freestyle_lineset_add(wmOperatorType *ot);
void SCENE_OT_freestyle_lineset_copy(wmOperatorType *ot);
void SCENE_OT_freestyle_lineset_paste(wmOperatorType *ot);
void SCENE_OT_freestyle_lineset_remove(wmOperatorType *ot);
void SCENE_OT_freestyle_lineset_move(wmOperatorType *ot);
void SCENE_OT_freestyle_linestyle_new(wmOperatorType *ot);
void SCENE_OT_freestyle_color_modifier_add(wmOperatorType *ot);
void SCENE_OT_freestyle_alpha_modifier_add(wmOperatorType *ot);
void SCENE_OT_freestyle_thickness_modifier_add(wmOperatorType *ot);
void SCENE_OT_freestyle_geometry_modifier_add(wmOperatorType *ot);
void SCENE_OT_freestyle_modifier_remove(wmOperatorType *ot);
void SCENE_OT_freestyle_modifier_move(wmOperatorType *ot);
void SCENE_OT_freestyle_modifier_copy(wmOperatorType *ot);
void SCENE_OT_freestyle_stroke_material_create(wmOperatorType *ot);
#endif

void TEXTURE_OT_slot_copy(wmOperatorType *ot);
void TEXTURE_OT_slot_paste(wmOperatorType *ot);
void TEXTURE_OT_slot_move(wmOperatorType *ot);

/* `render_texture_cache.cc` */

void RENDER_OT_generate_texture_cache(wmOperatorType *ot);
void RENDER_OT_clear_texture_cache(wmOperatorType *ot);

/* `material_texture_layer_assets.cc` */

namespace asset_system {
class AssetRepresentation;
}

namespace ed::render {
/** The #eShaderNodeTreeUsage flags of a shader node group asset, from its
 * metadata (0 when it is not a shader node group or carries no usage). */
eShaderNodeTreeUsage asset_texture_layer_usage(const asset_system::AssetRepresentation &asset);
void material_texture_layer_assets_register();
}  // namespace ed::render

/* `material_texture_layers.cc` */

struct Material;

namespace ed::render {
/** The active material, its node tree and the active Texture Layer Stack in
 * it, resolved from the context. All pointers are non-null. */
struct ActiveStackContext {
  Material *material;
  bNodeTree *ntree;
  bNode *stack;
  /** The active layer of #stack, or -1. Resolved from the active node (which
   * points at the layer, see #set_active_layer) rather than read from the
   * stack, so it stays right when the active node changed without going
   * through the layer operators. */
  int layer_index;
};
std::optional<ActiveStackContext> resolve_active_stack(const bContext &C);

/** The active material from the context, or null. */
Material *active_material(const bContext &C);
/** True when the active material's embedded node tree may be structurally edited (not linked or a
 * library override). Gates the layer edit operators and tree-view callbacks. */
bool active_material_editable(const bContext &C);
/** Tag the tree as changed: ensure invariants and send node + material notifiers. */
void tree_changed(bContext &C, bNodeTree &ntree, Material *mat);
/**
 * Make layer #index of #stack (a Texture Layer Stack or Mask Stack node) the
 * active one: stores it on the stack, points the node editor's active node at
 * the layer (see #layer_active_node) and selects the layer's nodes, so the
 * node editor shows what the layer panel edits.
 */
void set_active_layer(bNodeTree &ntree, bNode &stack, int index);
/** The node a layer is represented by in the node editor: the node whose
 * properties the layer panel shows, or the stack node itself for a layer with
 * no source and for a group layer (whose nested stack stands for its
 * children). */
bNode *layer_active_node(bNodeTree &ntree, bNode &stack, int index);
/** Select the nodes making up layer #index of #stack, and deselect the rest. */
void select_layer_nodes(bNodeTree &ntree, bNode &stack, int index);
/** Ensure the active layer has a Mask Stack and return it, or null when there
 * is no active layer to attach one to. */
bNode *ensure_active_layer_mask_stack(bContext &C, bNodeTree &ntree);
void material_texture_layers_register();
}  // namespace ed::render

/* `render_internal.cc` */

/* Base class for all WM_JOB_TYPE_RENDER jobs. */
struct RenderJobBase {
  Scene *scene = nullptr;
  Scene *current_scene = nullptr;
};

/**
 * Contextual render, using current scene, view3d?
 */
void RENDER_OT_render(wmOperatorType *ot);
void RENDER_OT_shutter_curve_preset(wmOperatorType *ot);

/* `render_view.cc` */

/**
 * New window uses x,y to set position.
 */
ScrArea *render_view_open(bContext *C, int mx, int my, ReportList *reports);

void RENDER_OT_view_show(wmOperatorType *ot);
void RENDER_OT_view_cancel(wmOperatorType *ot);

/* `render_opengl.cc` */

void RENDER_OT_opengl(wmOperatorType *ot);

}  // namespace blender
