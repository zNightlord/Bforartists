/* SPDX-FileCopyrightText: 2019-2022 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "infos/overlay_paint_infos.hh"

VERTEX_SHADER_CREATE_INFO(overlay_paint_point)

#include "draw_model_lib.glsl"
#include "draw_view_clipping_lib.glsl"
#include "draw_view_lib.glsl"

void main()
{
  float3 world_pos = drw_point_object_to_world(pos);
  gl_Position = drw_point_world_to_homogenous(world_pos);

  bool is_select = (paint_overlay_flag & 1) != 0;
  bool is_hidden = (paint_overlay_flag & 2) != 0;

  if (is_hidden) {
    gl_Position = float4(-2.0f, -2.0f, -2.0f, 1.0f);
  }

  if (use_colored_vertex) {
    float4 vg_col = float4(vgroup_color_blended, 1.0f);
    if (is_select) {
      /* Slight warm tint toward white — not fully white. */
      final_color = mix(vg_col, float4(1.0f, 1.0f, 1.0f, 1.0f), 0.35f);
    }
    else {
      final_color = vg_col;
    }
  }
  else {
    /* Normal vertex display — white when selected. */
    final_color = is_select ? float4(1.0f) : theme.colors.wire;
  }

  /* Opacity controls all vertex visibility including colored mode. */
  final_color.a *= colored_opacity;

  gl_PointSize = theme.sizes.vert * 2.0f;
  view_clipping_distances(world_pos);
}