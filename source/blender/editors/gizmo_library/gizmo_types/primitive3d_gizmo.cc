/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edgizmolib
 *
 * \name Primitive Gizmo
 *
 * 3D Gizmo
 *
 * \brief Gizmo with primitive drawing type (plane, cube, etc.).
 * Currently only plane primitive supported without its own handling, use with operator only.
 */

#include "MEM_guardedalloc.h"

#include "BLI_math_matrix.h"
#include "BLI_math_vector_types.hh"

#include "BKE_context.hh"

#include "DNA_screen_types.h"
#include "DNA_space_types.h"
#include "DNA_view3d_types.h"

#include "GPU_immediate.hh"
#include "GPU_immediate_util.hh"
#include "GPU_matrix.hh"
#include "GPU_select.hh"
#include "GPU_state.hh"

#include "RNA_define.hh"
#include "RNA_access.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "ED_gizmo_library.hh"
#include "ED_screen.hh"
#include "ED_transform_snap_object_context.hh"
#include "ED_view3d.hh"

/* own includes */
#include "../gizmo_library_intern.hh"

#define MVAL_MAX_PX_DIST 12.0f

static float verts_plane[4][3] = {
    {-1, -1, 0},
    {1, -1, 0},
    {1, 1, 0},
    {-1, 1, 0},
};

struct PrimitiveGizmo3D {
  wmGizmo gizmo;

  int draw_style;
  float arc_inner_factor;
  bool draw_inner;
  /* Added to 'matrix_basis' when calculating the matrix. */
  float prop_co[3];
};

struct PrimitiveInteraction {
  struct {
    float mval[2];
    /* Only for when using properties. */
    float prop_co[3];
    float matrix_final[4][4];
  } init;
  struct {
    eWM_GizmoFlagTweak tweak_flag;
  } prev;

  /* We could have other snap contexts, for now only support 3D view. */
  SnapObjectContext *snap_context_v3d;
};

static void gizmo_primitive_matrix_basis_get(const wmGizmo *gz, float r_matrix[4][4])
{
  PrimitiveGizmo3D *gz_prim = (PrimitiveGizmo3D *)gz;

  copy_m4_m4(r_matrix, gz_prim->gizmo.matrix_basis);
  add_v3_v3(r_matrix[3], gz_prim->prop_co);
}

/* -------------------------------------------------------------------- */
/** \name RNA callbacks */

static PrimitiveGizmo3D *gizmo_primitive_rna_find_operator(PointerRNA *ptr)
{
  return (PrimitiveGizmo3D *)gizmo_find_from_properties(
      static_cast<const IDProperty *>(ptr->data), SPACE_TYPE_ANY, RGN_TYPE_ANY);
}

static int gizmo_primitive_rna__draw_style_get_fn(PointerRNA *ptr, PropertyRNA * /*prop*/)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  return gz_prim->draw_style;
}

static void gizmo_primitive_rna__draw_style_set_fn(PointerRNA *ptr,
                                                   PropertyRNA * /*prop*/,
                                                   int value)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  gz_prim->draw_style = value;
}

static float gizmo_primitive_rna__arc_inner_factor_get_fn(PointerRNA *ptr, PropertyRNA * /*prop*/)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  return gz_prim->arc_inner_factor;
}

static void gizmo_primitive_rna__arc_inner_factor_set_fn(PointerRNA *ptr,
                                                         PropertyRNA * /*prop*/,
                                                         float value)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  gz_prim->arc_inner_factor = value;
}

static bool gizmo_primitive_rna__draw_inner_get_fn(PointerRNA *ptr, PropertyRNA * /*prop*/)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  return gz_prim->draw_inner;
}

static void gizmo_primitive_rna__draw_inner_set_fn(PointerRNA *ptr,
                                                   PropertyRNA * /*prop*/,
                                                   bool value)
{
  PrimitiveGizmo3D *gz_prim = gizmo_primitive_rna_find_operator(ptr);
  gz_prim->draw_inner = value;
}

/* -------------------------------------------------------------------- */

static void gizmo_primitive_draw_geom(PrimitiveGizmo3D *gz_prim,
                                      const float col_inner[4],
                                      const float col_outer[4],
                                      const int nsegments,
                                      const bool draw_inner)
{
  uint pos = GPU_vertformat_attr_add(immVertexFormat(), "pos", GPU_COMP_F32, 3, GPU_FETCH_FLOAT);
  const bool use_polyline_shader = gz_prim->gizmo.line_width > 1.0f;

  if (draw_inner || !use_polyline_shader) {
    immBindBuiltinProgram(GPU_SHADER_3D_UNIFORM_COLOR);
  }

  if (draw_inner) {
    if (gz_prim->draw_style == ED_GIZMO_PRIMITIVE_STYLE_PLANE) {
      wm_gizmo_vec_draw(col_inner, verts_plane, ARRAY_SIZE(verts_plane), pos, GPU_PRIM_TRI_FAN);
    }
    else {
      immUniformColor4fv(col_inner);
      if (gz_prim->draw_style == ED_GIZMO_PRIMITIVE_STYLE_CIRCLE) {
        imm_draw_circle_fill_3d(pos, 0.0f, 0.0f, 1.0f, nsegments);
      }
      else {
        BLI_assert(gz_prim->draw_style == ED_GIZMO_PRIMITIVE_STYLE_ANNULUS);
        imm_draw_disk_partial_fill_3d(
            pos, 0.0f, 0.0f, 0.0f, gz_prim->arc_inner_factor, 1.0f, nsegments, 0.0f, 360.0f);
      }
    }
  }

  /* Draw outline. */

  if (use_polyline_shader) {
    if (draw_inner) {
      immUnbindProgram();
    }
    immBindBuiltinProgram(GPU_SHADER_3D_POLYLINE_UNIFORM_COLOR);

    float viewport[4];
    GPU_viewport_size_get_f(viewport);
    immUniform2fv("viewportSize", &viewport[2]);
    immUniform1f("lineWidth", gz_prim->gizmo.line_width * U.pixelsize);
  }

  if (gz_prim->draw_style == ED_GIZMO_PRIMITIVE_STYLE_PLANE) {
    wm_gizmo_vec_draw(col_outer, verts_plane, ARRAY_SIZE(verts_plane), pos, GPU_PRIM_LINE_LOOP);
  }
  else {
    immUniformColor4fv(col_outer);
    if (gz_prim->draw_style == ED_GIZMO_PRIMITIVE_STYLE_CIRCLE) {
      imm_draw_circle_wire_3d(pos, 0.0f, 0.0f, 1.0f, nsegments);
    }
    else {
      imm_draw_circle_wire_3d(pos, 0.0f, 0.0f, gz_prim->arc_inner_factor, nsegments);
      imm_draw_circle_wire_3d(pos, 0.0f, 0.0f, 1.0f, nsegments);
    }
  }
  immUnbindProgram();
}

static void gizmo_primitive_draw_intern(wmGizmo *gz, const bool select, const bool highlight)
{
  PrimitiveGizmo3D *gz_prim = (PrimitiveGizmo3D *)gz;

  float color_inner[4], color_outer[4];
  float matrix_final[4][4];

  gizmo_color_get(gz, highlight, color_outer);
  copy_v4_v4(color_inner, color_outer);
  color_inner[3] *= 0.5f;

  WM_gizmo_calc_matrix_final(gz, matrix_final);

  GPU_blend(GPU_BLEND_ALPHA);
  GPU_matrix_push();
  GPU_matrix_mul(matrix_final);

  gizmo_primitive_draw_geom(gz_prim,
                            color_inner,
                            color_outer,
                            select ? 24 : DIAL_RESOLUTION,
                            gz_prim->draw_inner || select);

  GPU_matrix_pop();

  if (gz->interaction_data) {
    GizmoInteraction *inter = static_cast<GizmoInteraction *>(gz->interaction_data);

    copy_v4_fl(color_inner, 0.5f);
    copy_v3_fl(color_outer, 0.5f);
    color_outer[3] = 0.8f;

    GPU_matrix_push();
    GPU_matrix_mul(inter->init_matrix_final);

    gizmo_primitive_draw_geom(
        gz_prim, color_inner, color_outer, DIAL_RESOLUTION, gz_prim->draw_inner);

    GPU_matrix_pop();
  }
  GPU_blend(GPU_BLEND_NONE);
}

static void gizmo_primitive_draw_select(const bContext * /*C*/, wmGizmo *gz, int select_id)
{
  GPU_select_load_id(select_id);
  gizmo_primitive_draw_intern(gz, true, false);
}

static void gizmo_primitive_draw(const bContext * /*C*/, wmGizmo *gz)
{
  gizmo_primitive_draw_intern(gz, false, (gz->state & WM_GIZMO_STATE_HIGHLIGHT));
}

static void primitive_get_translate(const wmGizmo *gz,
                                 const wmEvent *event,
                                 const ARegion *region,
                                 float co_delta[3])
{
  PrimitiveInteraction *inter = static_cast<PrimitiveInteraction *>(gz->interaction_data);
  const float xy_delta[2] = {
      event->mval[0] - inter->init.mval[0],
      event->mval[1] - inter->init.mval[1],
  };

  RegionView3D *rv3d = static_cast<RegionView3D *>(region->regiondata);
  float co_ref[3];
  mul_v3_mat3_m4v3(co_ref, gz->matrix_space, inter->init.prop_co);
  const float zfac = ED_view3d_calc_zfac(rv3d, co_ref);

  ED_view3d_win_to_delta(region, xy_delta, zfac, co_delta);

  float matrix_space_inv[3][3];
  copy_m3_m4(matrix_space_inv, gz->matrix_space);
  invert_m3(matrix_space_inv);
  mul_m3_v3(matrix_space_inv, co_delta);
}

static int gizmo_primitive_modal(bContext *C,
                            wmGizmo *gz,
                            const wmEvent *event,
                            eWM_GizmoFlagTweak tweak_flag)
{
  PrimitiveInteraction *inter = static_cast<PrimitiveInteraction *>(gz->interaction_data);
  if ((event->type != MOUSEMOVE) && (inter->prev.tweak_flag == tweak_flag)) {
    return OPERATOR_RUNNING_MODAL;
  }
  PrimitiveGizmo3D *prim = (PrimitiveGizmo3D *)gz;
  ARegion *region = CTX_wm_region(C);

  float prop_delta[3];
  if (CTX_wm_area(C)->spacetype == SPACE_VIEW3D) {
    primitive_get_translate(gz, event, region, prop_delta);
  }
  else {
    float mval_proj_init[2], mval_proj_curr[2];
    if ((gizmo_window_project_2d(C, gz, inter->init.mval, 2, false, mval_proj_init) == false) ||
        (gizmo_window_project_2d(
             C, gz, blender::float2(blender::int2(event->mval)), 2, false, mval_proj_curr) ==
         false))
    {
      return OPERATOR_RUNNING_MODAL;
    }
    sub_v2_v2v2(prop_delta, mval_proj_curr, mval_proj_init);
    if ((gz->flag & WM_GIZMO_DRAW_NO_SCALE) == 0) {
      mul_v2_fl(prop_delta, gz->scale_final);
    }
    prop_delta[2] = 0.0f;
  }

  if (tweak_flag & WM_GIZMO_TWEAK_PRECISE) {
    mul_v3_fl(prop_delta, 0.1f);
  }

  add_v3_v3v3(prim->prop_co, inter->init.prop_co, prop_delta);

  if (tweak_flag & WM_GIZMO_TWEAK_SNAP) {
    if (inter->snap_context_v3d) {
      float dist_px = MVAL_MAX_PX_DIST * U.pixelsize;
      const float mval_fl[2] = {float(event->mval[0]), float(event->mval[1])};
      float co[3];
      SnapObjectParams params{};
      params.snap_target_select = SCE_SNAP_TARGET_ALL;
      params.edit_mode_type = SNAP_GEOM_EDIT;
      params.occlusion_test = SNAP_OCCLUSION_AS_SEEM;
      if (ED_transform_snap_object_project_view3d(
              inter->snap_context_v3d,
              CTX_data_ensure_evaluated_depsgraph(C),
              region,
              CTX_wm_view3d(C),
              (SCE_SNAP_TO_VERTEX | SCE_SNAP_TO_EDGE | SCE_SNAP_TO_FACE),
              &params,
              nullptr,
              mval_fl,
              nullptr,
              &dist_px,
              co,
              nullptr))
      {
        float matrix_space_inv[4][4];
        invert_m4_m4(matrix_space_inv, gz->matrix_space);
        mul_v3_m4v3(prim->prop_co, matrix_space_inv, co);
      }
    }
  }

  /* set the property for the operator and call its modal function */
  wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "offset");
  if (WM_gizmo_target_property_is_valid(gz_prop)) {
    WM_gizmo_target_property_float_set_array(C, gz, gz_prop, prim->prop_co);
  }
  else {
    zero_v3(prim->prop_co);
  }

  ED_region_tag_redraw_editor_overlays(region);

  inter->prev.tweak_flag = tweak_flag;

  return OPERATOR_RUNNING_MODAL;
}

static void gizmo_primitive_setup(wmGizmo *gz)
{
  gz->flag |= WM_GIZMO_DRAW_MODAL;

  /* Default Values. */
  PrimitiveGizmo3D *gz_prim = (PrimitiveGizmo3D *)gz;
  gz_prim->draw_style = ED_GIZMO_PRIMITIVE_STYLE_PLANE;
  gz_prim->arc_inner_factor = 1.0f;
  gz_prim->draw_inner = true;
}

static int gizmo_primitive_invoke(bContext * /*C*/, wmGizmo *gz, const wmEvent *event)
{
   PrimitiveInteraction *inter = static_cast<PrimitiveInteraction *>(
      MEM_callocN(sizeof(PrimitiveInteraction), __func__));
    inter->init.mval[0] = event->mval[0];
    inter->init.mval[1] = event->mval[1];

 wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "offset");
  if (WM_gizmo_target_property_is_valid(gz_prop)) {
    WM_gizmo_target_property_float_get_array(gz, gz_prop, inter->init.prop_co);
  }

  WM_gizmo_calc_matrix_final(gz, inter->init.matrix_final);

  gz->interaction_data = inter;

  return OPERATOR_RUNNING_MODAL;
}

static void gizmo_primitive_exit(bContext *C, wmGizmo *gz, const bool cancel)
{
  PrimitiveInteraction *inter = static_cast<PrimitiveInteraction *>(gz->interaction_data);
  bool use_reset_value = false;
  const float *reset_value = nullptr;
  if (cancel) {
    /* Set the property for the operator and call its modal function. */
    wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "offset");
    if (WM_gizmo_target_property_is_valid(gz_prop)) {
      use_reset_value = true;
      reset_value = inter->init.prop_co;
    }
  }

  if (use_reset_value) {
    wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "offset");
    if (WM_gizmo_target_property_is_valid(gz_prop)) {
      WM_gizmo_target_property_float_set_array(C, gz, gz_prop, reset_value);
    }
  }

  if (inter->snap_context_v3d) {
    ED_transform_snap_object_context_destroy(inter->snap_context_v3d);
    inter->snap_context_v3d = nullptr;
  }

  if (!cancel) {
    wmGizmoProperty *gz_prop = WM_gizmo_target_property_find(gz, "offset");
    if (WM_gizmo_target_property_is_valid(gz_prop)) {
      WM_gizmo_target_property_anim_autokey(C, gz, gz_prop);
    }
  }
}


static int gizmo_primitive_test_select(bContext *C, wmGizmo *gz, const int mval[2])
{
  float point_local[2];

  if (gizmo_window_project_2d(C, gz, blender::float2(blender::int2(mval)), 2, true, point_local) ==
      false)
  {
    return -1;
  }

  /* The 'gz->scale_final' is already applied to the projection
   * when #WM_GIZMO_DRAW_NO_SCALE isn't set. */
  const float radius = (gz->flag & WM_GIZMO_DRAW_NO_SCALE) ? gz->scale_final : 1.0f;
  if (len_squared_v2(point_local) < radius) {
    return 0;
  }

  return -1;
}

static void gizmo_move_property_update(wmGizmo *gz, wmGizmoProperty *gz_prop)
{
  PrimitiveGizmo3D *prim = (PrimitiveGizmo3D *)gz;
  if (WM_gizmo_target_property_is_valid(gz_prop)) {
    WM_gizmo_target_property_float_get_array(gz, gz_prop, prim->prop_co);
  }
  else {
    zero_v3(prim->prop_co);
  }
}

static int gizmo_primitive_cursor_get(wmGizmo * /*gz*/)
{
  return WM_CURSOR_NSEW_SCROLL;
}


/* -------------------------------------------------------------------- */
/** \name Primitive Gizmo API
 * \{ */

static void GIZMO_GT_primitive_3d(wmGizmoType *gzt)
{
  /* identifiers */
  gzt->idname = "GIZMO_GT_primitive_3d";

  /* api callbacks */
  gzt->draw = gizmo_primitive_draw;
  gzt->draw_select = gizmo_primitive_draw_select;
  gzt->test_select = gizmo_primitive_test_select;
  gzt->matrix_basis_get = gizmo_primitive_matrix_basis_get;
  gzt->setup = gizmo_primitive_setup;
  gzt->invoke = gizmo_primitive_invoke;
  gzt->property_update = gizmo_move_property_update;
  gzt->modal = gizmo_primitive_modal;
  gzt->exit = gizmo_primitive_exit;
  gzt->cursor_get = gizmo_primitive_cursor_get;

  gzt->struct_size = sizeof(PrimitiveGizmo3D);

  static const EnumPropertyItem rna_enum_draw_style[] = {
      {ED_GIZMO_PRIMITIVE_STYLE_PLANE, "PLANE", 0, "Plane", ""},
      {ED_GIZMO_PRIMITIVE_STYLE_CIRCLE, "CIRCLE", 0, "Circle", ""},
      {ED_GIZMO_PRIMITIVE_STYLE_ANNULUS, "ANNULUS", 0, "Annulus", ""},
      {0, nullptr, 0, nullptr, nullptr},
  };

  PropertyRNA *prop;
  prop = RNA_def_enum(gzt->srna,
                      "draw_style",
                      rna_enum_draw_style,
                      ED_GIZMO_PRIMITIVE_STYLE_PLANE,
                      "Draw Style",
                      "");
  RNA_def_property_enum_funcs_runtime(prop,
                                      gizmo_primitive_rna__draw_style_get_fn,
                                      gizmo_primitive_rna__draw_style_set_fn,
                                      nullptr);

  prop = RNA_def_float_factor(
      gzt->srna, "arc_inner_factor", 0.0f, 0.0f, FLT_MAX, "Arc Inner Factor", "", 0.0f, 1.0f);
  RNA_def_property_float_funcs_runtime(prop,
                                       gizmo_primitive_rna__arc_inner_factor_get_fn,
                                       gizmo_primitive_rna__arc_inner_factor_set_fn,
                                       nullptr);

  prop = RNA_def_boolean(gzt->srna, "draw_inner", true, "Draw Inner", "");
  RNA_def_property_boolean_funcs_runtime(
      prop, gizmo_primitive_rna__draw_inner_get_fn, gizmo_primitive_rna__draw_inner_set_fn);
    
  WM_gizmotype_target_property_def(gzt, "offset", PROP_FLOAT, 3);
}

void ED_gizmotypes_primitive_3d()
{
  WM_gizmotype_append(GIZMO_GT_primitive_3d);
}

/** \} */
