/* SPDX-FileCopyrightText: 2019 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma BLENDER_REQUIRE(juniper_material_lib.glsl)

void node_bevel(float radius, vec3 N, out vec3 result)
{
  vec2 screenuv = coordinate_screen(g_data.P).xy;
  ivec2 texSize = textureSize(in_normal_tx, 0);

  vec2 pix_scale = 1.0 / texSize;
  result = get_prepass_normal(screenuv + (N.xy * pix_scale));
}
