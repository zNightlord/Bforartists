/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/**
 * Screen-space ray-tracing routine.
 *
 * Based on "Efficient GPU Screen-Space Ray Tracing"
 * by Morgan McGuire & Michael Mara
 * https://jcgt.org/published/0003/04/04/paper.pdf
 *
 * Many modifications were made for our own usage.
 */

#include "draw_view_lib.glsl"
#include "eevee_bxdf_diffuse_lib.glsl"
#include "eevee_bxdf_microfacet_lib.glsl"
#include "eevee_ray_types_lib.glsl"
#include "eevee_reverse_z_lib.glsl"
#include "eevee_thickness_lib.glsl"
#include "gpu_shader_codegen_lib.glsl"
#include "gpu_shader_math_fast_lib.glsl"

/* Inputs expected to be in view-space. */
void raytrace_clip_ray_to_near_plane(inout Ray ray)
{
  float near_dist = drw_view_near();
  if ((ray.origin.z + ray.direction.z * ray.max_time) > near_dist) {
    ray.max_time = abs((near_dist - ray.origin.z) / ray.direction.z);
  }
}

struct ScreenTraceHitData {
  /* Screen space hit position [0..1]. Last component is the ray depth, not the occluder depth. */
  float3 ss_hit_P;
  /* View space hit position. */
  float3 v_hit_P;
  /* Tracing time in world space. */
  float time;
  /* True if there was a valid intersection. False if went out of screen without intersection. */
  bool valid;
};

/**
 * Ray-trace against the given HIZ-buffer height-field.
 *
 * \param stride_rand: Random number in [0..1] range. Offset along the ray to avoid banding
 *                     artifact when steps are too large.
 * \param roughness: Determine how lower depth mipmaps are used to make the tracing faster. Lower
 *                   roughness will use lower mipmaps.
 * \param discard_backface: If true, ray-trace will return false if we hit a surface from behind.
 * \param allow_self_intersection: If false, ray-trace will return false if the ray is not covering
 *                                 at least one pixel.
 * \param ray: View-space ray. Direction pre-multiplied by maximum length.
 *
 * \return True if there is a valid intersection.
 */
ScreenTraceHitData raytrace_screen(RayTraceData rt_data,
                                   HiZData hiz_data,
                                   sampler2D hiz_tx,
                                   float stride_rand,
                                   float roughness,
                                   const bool discard_backface,
                                   const bool allow_self_intersection,
                                   Ray ray)
{
  /* Clip to near plane for perspective view where there is a singularity at the camera origin. */
  if (drw_view().winmat[3][3] == 0.0f) {
    raytrace_clip_ray_to_near_plane(ray);
  }

  /* NOTE: The 2.0 factor here is because we are applying it in NDC space. */
  ScreenSpaceRay ssray = raytrace_screenspace_ray_create(
      ray, 2.0f * rt_data.full_resolution_inv, rt_data.thickness);

  /* Avoid no iteration. */
  if (!allow_self_intersection && ssray.max_time < 1.1f) {
    /* Still output the clipped ray. */
    float3 hit_ssP = ssray.origin.xyz + ssray.direction.xyz * ssray.max_time;
    float3 hit_P = drw_point_screen_to_world(float3(hit_ssP.xy, saturate(hit_ssP.z)));
    ray.direction = hit_P - ray.origin;

    ScreenTraceHitData no_hit;
    no_hit.time = 0.0f;
    no_hit.valid = false;
    return no_hit;
  }

  ssray.max_time = max(1.1f, ssray.max_time);

  float prev_delta = 0.0f, prev_time = 0.0f;
  float depth_sample = drw_depth_view_to_screen(ray.origin.z);
  float delta = depth_sample - ssray.origin.z;

  float lod_fac = saturate(sqrt_fast(roughness) * 2.0f - 0.4f);

  /* Cross at least one pixel. */
  float t = 1.001f, time = 1.001f;
  bool hit = false;
  constexpr int max_steps = 255;
  for (int iter = 1; !hit && (time < ssray.max_time) && (iter < max_steps); iter++) {
    float stride = 1.0f + float(iter) * rt_data.quality;
    float lod = log2(stride) * lod_fac;

    prev_time = time;
    prev_delta = delta;

    time = min(t + stride * stride_rand, ssray.max_time);
    t += stride;

    float4 ss_p = ssray.origin + ssray.direction * time;
    depth_sample = textureLod(hiz_tx, ss_p.xy * hiz_data.uv_scale, floor(lod)).r;

    delta = depth_sample - ss_p.z;
    /* Check if the ray is below the surface ... */
    hit = (delta < 0.0f);
    /* ... and above it with the added thickness. */
    hit = hit && (delta > ss_p.z - ss_p.w || abs(delta) < abs(ssray.direction.z * stride * 2.0f));
  }
  /* Discard back-face hits. */
  hit = hit && !(discard_backface && prev_delta < 0.0f);
  /* Reject hit if background. */
  hit = hit && (depth_sample != 1.0f);
  /* Refine hit using intersection between the sampled height-field and the ray.
   * This simplifies nicely to this single line. */
  time = mix(prev_time, time, saturate(prev_delta / (prev_delta - delta)));

  ScreenTraceHitData result;
  result.ss_hit_P = ssray.origin.xyz + ssray.direction.xyz * time;
  result.v_hit_P = drw_point_screen_to_view(result.ss_hit_P);
  /* Convert to world space ray time. */
  result.time = length(result.v_hit_P - ray.origin) / length(ray.direction);
  /* Update the validity as ss_hit_P can point to a background sample. */
  result.valid = hit &&
                 (textureLod(hiz_tx, result.ss_hit_P.xy * hiz_data.uv_scale, 0.0f).r != 1.0f);

  return result;
}

#ifdef PLANAR_PROBES

ScreenTraceHitData raytrace_planar(RayTraceData rt_data,
                                   sampler2DArrayDepth planar_depth_tx,
                                   PlanarProbeData planar,
                                   float stride_rand,
                                   Ray ray)
{
  /* Clip to near plane for perspective view where there is a singularity at the camera origin. */
  if (drw_view().winmat[3][3] == 0.0f) {
    raytrace_clip_ray_to_near_plane(ray);
  }

  float2 inv_texture_size = 1.0f / float2(textureSize(planar_depth_tx, 0).xy);
  /* NOTE: The 2.0 factor here is because we are applying it in NDC space. */
  ScreenSpaceRay ssray = raytrace_screenspace_ray_create(
      ray, planar.winmat, 2.0f * inv_texture_size);

  float prev_delta = 0.0f, prev_time = 0.0f;
  float depth_sample = reverse_z::read(
      texture(planar_depth_tx, float3(ssray.origin.xy, planar.layer_id)).r);
  float delta = depth_sample - ssray.origin.z;

  float t = 0.0f, time = 0.0f;
  bool hit = false;
  constexpr int max_steps = 32;
  for (int iter = 1; !hit && (time < ssray.max_time) && (iter < max_steps); iter++) {
    float stride = 1.0f + float(iter) * rt_data.quality;

    prev_time = time;
    prev_delta = delta;

    time = min(t + stride * stride_rand, ssray.max_time);
    t += stride;

    float4 ss_ray = ssray.origin + ssray.direction * time;

    depth_sample = reverse_z::read(texture(planar_depth_tx, float3(ss_ray.xy, planar.layer_id)).r);

    delta = depth_sample - ss_ray.z;
    /* Check if the ray is below the surface. */
    hit = (delta < 0.0f);
  }
  /* Reject hit if background. */
  hit = hit && (depth_sample != 1.0f);
  /* Refine hit using intersection between the sampled height-field and the ray.
   * This simplifies nicely to this single line. */
  time = mix(prev_time, time, saturate(prev_delta / (prev_delta - delta)));

  ScreenTraceHitData result;
  result.ss_hit_P = ssray.origin.xyz + ssray.direction.xyz * time;
  /* Update the validity as ss_hit_P can point to a not loaded sample. */
  result.valid =
      hit &&
      textureLod(planar_depth_tx, float3(result.ss_hit_P.xy, planar.layer_id), 0.0f).r != 0.0;

  /* NOTE: v_hit_P is in planar reflected view space. */
  result.v_hit_P = project_point(planar.wininv, drw_screen_to_ndc(result.ss_hit_P));
  /* Convert to world space ray time. */
  result.time = length(result.v_hit_P - ray.origin) / length(ray.direction);
  return result;
}

#endif

/* Modify the ray origin before tracing it. We must do this because ray origin is implicitly
 * reconstructed from gbuffer depth which we cannot modify. */
Ray raytrace_thickness_ray_amend(Ray ray, ClosureUndetermined cl, float3 V, float thickness)
{
  switch (cl.type) {
    case CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID:
      return bxdf_ggx_ray_amend_transmission(cl, V, ray, thickness);
    case CLOSURE_BSDF_TRANSLUCENT_ID:
      return bxdf_translucent_ray_amend(cl, V, ray, thickness);
    case CLOSURE_NONE_ID:
    case CLOSURE_BSDF_DIFFUSE_ID:
    case CLOSURE_BSDF_MICROFACET_GGX_REFLECTION_ID:
    case CLOSURE_BSSRDF_BURLEY_ID:
      break;
  }
  return ray;
}

/**
 * Adapted from G3D Innovation Engine http://casual-effects.com/g3d
 * Copyright 2000-2019, Morgan McGuire (BSD License)
 * (Updated version of Efficient GPU Screen-Space Ray Tracing -
 * https://jcgt.org/published/0003/04/04/paper.pdf)
 *
 * \param vs_origin: Camera-space ray origin, which must be within the view volume and must have
 *                   z < -0.01 and project within the valid screen rectangle.
 * \param vs_direction: Unit length camera-space ray direction.
 * \param max_raytrace_distance: Maximum camera-space distance to trace before returning a miss.
 * \param hiz_tx: The depth Z buffer.
 * \param thickness: Camera space thickness to ascribe to each pixel in the depth buffer.
 * \param stride: Step in horizontal or vertical pixels between samples. This is a float because
 *                integer math is slow on GPUs, but should be set to an integer >= 1.
 * \param jitter_fraction: Number between 0 and 1 for how far to bump the ray in stride units to
 *                         conceal banding artifacts, plus the stride ray offset. It is recommended
 *                         to set this to at least 1.0 to avoid self-intersection artifacts. Using
 *                         1 + float((int(gl_FragCoord.x) + int(gl_FragCoord.y)) & 1) * 0.5 gives
 *                         a nice dither pattern when stride is > 1.0.
 * \param max_steps: Maximum number of iterations. Higher gives better images but may be slow.
 * \param vs_hit_point: Camera space location of the ray hit.
 */
bool raytrace_screen_2(float3 vs_origin,
                       float3 vs_direction,
                       float max_raytrace_distance,
                       sampler2D hiz_tx,
                       float thickness,
                       float stride,
                       float jitter_fraction,
                       float max_steps,
                       out float3 vs_hit_point,
                       out float3 debug_color)
{
  if (drw_view_is_perspective()) {
    /* Negative number. Doesn't have to be THE actual near plane, just a reasonable value for
     * clipping rays headed towards the camera. */
    const float near_z_plane = drw_view_near();
    /* Clip ray to a near plane in 3D. */
    if ((vs_origin.z + vs_direction.z * max_raytrace_distance) > near_z_plane) {
      max_raytrace_distance = abs((near_z_plane - vs_origin.z) / vs_direction.z);
    }
  }

  float3 vs_end_point = vs_origin + vs_direction * max_raytrace_distance;

  float2 half_extent = float2(uniform_buf.film.render_extent) / 2.0f;

  /* Project into screen (pixel) space. */
  float4 h0 = drw_point_view_to_homogenous(vs_origin);
  h0.xy += half_extent;
  float4 h1 = drw_point_view_to_homogenous(vs_end_point);
  h1.xy += half_extent;

  /* There are a lot of divisions by w that can be turned into multiplications at some minor
   * precision loss... and we need to interpolate these 1/w values anyway.
   * Because the caller was required to clip to the near plane, this homogeneous division
   * (projecting from 4D to 2D) is guaranteed to succeed. */
  float k0 = 1.0 / h0.w;
  float k1 = 1.0 / h1.w;

  /* Switch the original points to values that interpolate linearly in 2D. */
  float3 q0 = vs_origin * k0;
  float3 q1 = vs_end_point * k1;

  /* Screen-space endpoints */
  float2 p0 = h0.xy * k0;
  float2 p1 = h1.xy * k1;

  /* Optional clipping to frustum sides. */
  const bool do_clip = false;
  if (do_clip) {
    float x_min = 0.5;
    float x_max = float(uniform_buf.film.render_extent.x) - 0.5;
    float y_min = 0.5;
    float y_max = float(uniform_buf.film.render_extent.y) - 0.5;
    float alpha = 0.0;

    /* Assume p0 is in the viewport (p1 - p0 is never zero when clipping) */
    if ((p1.y > y_max) || (p1.y < y_min)) {
      alpha = (p1.y - ((p1.y > y_max) ? y_max : y_min)) / (p1.y - p0.y);
    }

    if ((p1.x > x_max) || (p1.x < x_min)) {
      alpha = max(alpha, (p1.x - ((p1.x > x_max) ? x_max : x_min)) / (p1.x - p0.x));
    }

    p1 = mix(p1, p0, alpha);
    k1 = mix(k1, k0, alpha);
    q1 = mix(q1, q0, alpha);
  }

  /* Initialize to off screen. */
  float2 hit_pixel = float2(-1.0, -1.0);

  /* If the line is degenerate,
   * make it cover at least one pixel to avoid handling zero-length segments later. */
  p1 += float2((distance_squared(p0, p1) < 0.0001) ? 0.01 : 0.0);

  float2 delta = p1 - p0;

  /* Primary iteration is in x to reduce large branches later */
  bool permute = (abs(delta.x) < abs(delta.y));
  if (permute) {
    /* More-vertical line.
     * Create a permutation that swaps x and y in the output by directly swizzling the inputs. */
    delta = delta.yx;
    p1 = p1.yx;
    p0 = p0.yx;
  }

  /* From now on, "x" is the primary iteration direction and "y" is the secondary one. */
  float step_direction = sign(delta.x);
  float inv_dx = step_direction / delta.x;
  float2 d_p = float2(step_direction, inv_dx * delta.y);

  /* Track the derivatives of Q and k. */
  float3 d_q = (q1 - q0) * inv_dx;
  float d_k = (k1 - k0) * inv_dx;

  /* Because we test 1/2 a texel forward along the ray, on the very last iteration the
   * interpolation can go past the end of the ray. Use these bounds to clamp it. */
  float z_min = min(vs_end_point.z, vs_origin.z);
  float z_max = max(vs_end_point.z, vs_origin.z);

  /* Scale derivatives by the desired pixel stride. */
  d_p *= stride;
  d_q *= stride;
  d_k *= stride;

  /* Offset the starting values by the jitter fraction. */
  p0 += d_p * jitter_fraction;
  q0 += d_q * jitter_fraction;
  k0 += d_k * jitter_fraction;

  /* Slide p from p0 to p1, (now-homogeneous) q from q0 to q1, and k from k0 to k1. */
  float3 q = q0;
  float k = k0;

  /* We track the ray depth at +/- 1/2 pixel to treat pixels as clip-space solid voxels.
   * Because the depth at -1/2 for a given pixel will be the same as at +1/2 for the previous
   * iteration, we actually only have to compute one value per iteration. */
  float prev_z_max_estimate = vs_origin.z;
  float step_count = 0.0;
  float ray_z_max = prev_z_max_estimate;
  float ray_z_min = prev_z_max_estimate;
  float scene_z_max = ray_z_max + 1e4;

  /* p1.x is never modified after this point,
   * so pre-scale it by the step direction for a signed comparison. */
  float end = p1.x * step_direction;

  /* We only advance the z field of q in the inner loop, since q.xy is never used until after the
   * loop terminates. */
  float2 p;
  for (p = p0; ((p.x * step_direction) <= end) && (step_count < max_steps) &&
               ((ray_z_max < scene_z_max - thickness) || (ray_z_min > scene_z_max)) &&
               (scene_z_max != 0.0);
       p += d_p, q.z += d_q.z, k += d_k, step_count += 1.0)
  {

    /* The depth range that the ray covers within this loop iteration.
     * Assume that the ray is moving in increasing z and swap if backwards.
     * Because one end of the interval is shared between adjacent iterations,
     * we track the previous value and then swap as needed to ensure correct ordering. */
    ray_z_min = prev_z_max_estimate;

    /* Compute the value at 1/2 step into the future */
    ray_z_max = (d_q.z * 0.5 + q.z) / (d_k * 0.5 + k);
    ray_z_max = clamp(ray_z_max, z_min, z_max);
    prev_z_max_estimate = ray_z_max;

    /* Since we don't know if the ray is stepping forward or backward in depth, maybe swap.
     * Note that we preserve our original z "max" estimate first. */
    if (ray_z_min > ray_z_max) {
      /* Swap. */
      float tmp = ray_z_min;
      ray_z_min = ray_z_max;
      ray_z_max = tmp;
    }

    /* Camera-space z of the background. */
    hit_pixel = permute ? p.yx : p;
    scene_z_max = texelFetch(hiz_tx, int2(hit_pixel), 0).r;

    scene_z_max = drw_depth_screen_to_view(scene_z_max);
  }

  /* Undo the last increment, which ran after the test variables were set up. */
  p -= d_p;
  q.z -= d_q.z;
  k -= d_k;
  step_count -= 1.0;

  bool hit = (ray_z_max >= scene_z_max - thickness) && (ray_z_min <= scene_z_max);

  /* If using non-unit stride and we hit a depth surface... */
  if ((stride > 1.0) && hit) {
    /* Refine the hit point within the last large-stride step. */

    /* Retreat one whole stride step from the previous loop so that we can re-run that iteration at
     * finer scale. */
    p -= d_p;
    q.z -= d_q.z;
    k -= d_k;
    step_count -= 1.0;

    /* Take the derivatives back to single-pixel stride. */
    float inv_stride = 1.0 / stride;
    d_p *= inv_stride;
    d_q.z *= inv_stride;
    d_k *= inv_stride;

    /* For this test, we don't bother checking thickness or passing the end, since we KNOW there
     * will be a hit point.
     * As soon as the ray passes behind an object, call it a hit.
     * Advance (stride + 1) steps to fully check this interval (we could skip the very first
     * iteration, but then we'd need identical code to prime the loop). */
    float refinement_step_count = 0.0;

    /* This is the current sample point's z-value, taken back to camera space. */
    prev_z_max_estimate = q.z / k;
    ray_z_min = prev_z_max_estimate;

    /* Ensure that the FOR-loop test passes on the first iteration since we
     * won't have a valid value of scene_z_max to test. */
    scene_z_max = ray_z_min - 1e7;

    for (; (refinement_step_count <= stride * 1.4) && (ray_z_min > scene_z_max) &&
           (scene_z_max != 0.0);
         p += d_p, q.z += d_q.z, k += d_k, refinement_step_count += 1.0)
    {

      ray_z_min = prev_z_max_estimate;

      /* Compute the ray camera-space Z value at 1/2 fine step (pixel) into the future. */
      ray_z_max = (d_q.z * 0.5 + q.z) / (d_k * 0.5 + k);
      ray_z_max = clamp(ray_z_max, z_min, z_max);

      prev_z_max_estimate = ray_z_max;
      ray_z_min = min(ray_z_max, ray_z_min);

      hit_pixel = permute ? p.yx : p;
      scene_z_max = texelFetch(hiz_tx, int2(hit_pixel), 0).r;

      scene_z_max = drw_depth_screen_to_view(scene_z_max);
    }

    /* Undo the last increment, which happened after the test variables were set up. */
    q.z -= d_q.z;
    refinement_step_count -= 1.0;

    /* Count the refinement steps as fractions of the original stride. Save a register
     * by not retaining inv_stride until here. */
    step_count += refinement_step_count / stride;
    // debug_color = float3(refinement_step_count / stride);
  }

  /* note: q was q0.. but we use q variable - adapt to q ensure q.xy is updated. */
  p.xy += d_q.xy * step_count;
  q.xy += d_q.xy * step_count;
  vs_hit_point = q * (1.0 / k);

  /* Support debugging. This will compile away if debug_color is unused */
  if ((p.x * step_direction) > end) {
    /* Hit the max ray distance -> blue. */
    debug_color = float3(0.0, 0.0, 1.0);
  }
  else if (step_count >= max_steps) {
    /* Ran out of steps -> red. */
    debug_color = float3(1.0, 0.0, 0.0);
  }
  else if (scene_z_max == 0.0) {
    /* Went off screen -> yellow. */
    debug_color = float3(1.0, 1.0, 0.0);
  }
  else {
    /* Encountered a valid hit -> green.
     * ((ray_z_max >= scene_z_max - thickness) && (ray_z_min <= scene_z_max))*/
    debug_color = float3(0.0, 1.0, 0.0);
  }

  /* Does the last point discovered represent a valid hit? */
  return hit;
}
