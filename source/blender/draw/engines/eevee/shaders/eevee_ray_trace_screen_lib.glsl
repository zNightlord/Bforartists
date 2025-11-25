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

float4 hs_interpolation(float4 a, float4 b, float t)
{
  return mix(a * a.w, b * b.w, t);
}

bool raytrace_screen_2(float3 vs_origin,
                       float3 vs_direction,
                       float max_distance,
                       sampler2D hiz_tx,
                       float thickness,
                       int max_steps,
                       float jitter,
                       out float3 vs_hit_point)
{
  // TODO: Clip and discard if fully out of frustum.

  float4 start = drw_point_view_to_homogenous(vs_origin);
  float4 end = drw_point_view_to_homogenous(vs_origin + vs_direction * max_distance);

  float2 extent = float2(uniform_buf.film.render_extent);
  float2 half_extent = extent / 2.0f;

  /* Convert xy coordinates to pixel space (cartesian/non-homogeneous). */
  start.xy = (start.xy / start.w) * half_extent + half_extent;
  end.xy = (end.xy / end.w) * half_extent + half_extent;

  float4 delta = end - start;

  /* Number of steps required to trace a fully contiguous line. */
  int steps = int(max(abs(delta.x), abs(delta.y))) + 1;
  /* Limit to max steps. */
  steps = min(steps, max_steps);
  /* And scale steps accordingly. */
  float4 delta_step = delta / float(steps);

  float previous_step_depth = (start.z / start.w) * 0.5f + 0.5f;

  /* Skip the first step to avoid self-occlusion. But iterate at least once. */
  for (int i = 1; i < steps || i == 1; i++) {
    /* Ensure we don't go past ray end. */
    bool is_last_step = i + 1 >= steps;
    float step_jitter = is_last_step ? 0.0f : jitter;

    /* Compute the position of the current step. */
    float4 step = start + (delta_step * (float(i) + step_jitter));
    /* Convert to fragment coordinates. */
    step.z = (step.z / step.w) * 0.5f + 0.5f;

    // TODO: Clip the rays at the function start.
    if (any(lessThan(step.xyz, float3(0.0f))) ||
        any(greaterThan(step.xyz, float3(extent.x, extent.y, 1.0f))))
    {
      previous_step_depth = clamp(step.z, 0.0f, 1.0f);
      continue;
    }

    float screen_depth = texelFetch(hiz_tx, int2(step.xy), 0).r;

    /* Take thickness into account, but ensure it's not lower than the step depth. */
    float min_depth = min(previous_step_depth,
                          drw_depth_view_to_screen(drw_depth_screen_to_view(step.z) + thickness));

    if (step.z > screen_depth && screen_depth > min_depth) {
      /* We have a hit. Compute the hit point. */
      float3 ndc_hit_point;
      ndc_hit_point.xy = (step.xy - half_extent) / half_extent;
      ndc_hit_point.z = screen_depth * 2.0 - 1.0;
      vs_hit_point = drw_point_ndc_to_view(ndc_hit_point);
      /* Project the hit point into the ray. */
      float hit_distance = dot(vs_hit_point - vs_origin, vs_direction);
      vs_hit_point = vs_origin + vs_direction * hit_distance;
      return true;
    }

    previous_step_depth = step.z;
  }

  return false;
}
