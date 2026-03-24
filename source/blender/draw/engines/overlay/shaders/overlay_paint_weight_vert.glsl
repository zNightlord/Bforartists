/* SPDX-FileCopyrightText: 2018-2022 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "infos/overlay_paint_infos.hh"

VERTEX_SHADER_CREATE_INFO(overlay_paint_weight)

#include "draw_model_lib.glsl"
#include "draw_view_clipping_lib.glsl"
#include "draw_view_lib.glsl"

/* Integer hash → pseudo-random RGB */
float3 vgroup_index_to_color(int idx)
{
  uint h = uint(idx + vgroup_color_random_id);
  h ^= h >> 16u;
  h *= 0x45d9f3bu;
  h ^= h >> 16u;
  return float3(float( h         & 0xFFu),
                float((h >>  8u) & 0xFFu),
                float((h >> 16u) & 0xFFu)) / 255.0f;
}

void main()
{
  float3 world_pos = drw_point_object_to_world(pos);
  gl_Position = drw_point_world_to_homogenous(world_pos);

  if (vgroup_color_mode == 1) {
    /* SINGLE: whole mesh gets one color based on active group index.
     * vertex_group_index per vertex is ignored — use the active index
     * passed via vgroup_color_random_id offset trick or a separate constant.
     * Simplest: color is uniform, computed from 0 + random_id. */
    vgroup_color = vgroup_index_to_color(0);
  }
  else if (vgroup_color_mode == 2) {
    vgroup_color = vgroup_index_to_color(vertex_group_index);
    weight_interp = max(float2(vertex_group_dominant_weight, -vertex_group_dominant_weight), 0.0f);
  }
  else {
    vgroup_color = float3(0.0f);
  }

  /* Separate actual weight and alerts for independent interpolation */
  weight_interp = max(float2(weight, -weight), 0.0f);

  /* Saturate the weight to give a hint of the geometry behind the weights. */
#ifdef FAKE_SHADING
  float3 view_normal = normalize(drw_normal_object_to_view(nor));
  color_fac = abs(dot(view_normal, light_dir));
  color_fac = color_fac * 0.9f + 0.1f;

#else
  color_fac = 1.0f;
#endif

  view_clipping_distances(world_pos);
}
