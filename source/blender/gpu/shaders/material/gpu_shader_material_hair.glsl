/* SPDX-FileCopyrightText: 2022-2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "gpu_shader_compat.hh"

#define SHD_PRINCIPLED_HAIR_REFLECTANCE           0
#define SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION 1
#define SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION     2

[[node]]
void node_bsdf_hair(float4 color,
                    float offset,
                    float roughness_u,
                    float roughness_v,
                    float3 T,
                    float weight,
                    Closure &result)
{
  color = max(color, float4(0.0f));

#if 0
  /* NOTE(fclem): This is the way it should be. But we don't have proper implementation of the hair
   * closure yet. For now fall back to a simpler diffuse surface so that we have at least a color
   * feedback. */
  ClosureHair hair_data;
  hair_data.weight = weight;
  hair_data.color = color.rgb;
  hair_data.offset = offset;
  hair_data.roughness = float2(roughness_u, roughness_v);
  hair_data.T = T;
#else
  ClosureDiffuse hair_data;
  hair_data.weight = weight;
  hair_data.color = color.rgb;
  hair_data.N = g_data.N;
#endif
  result = closure_eval(hair_data);
}

/* ==================================================================
 * Principled Hair BSDF — EEVEE Next
 *
 * Faithful port of D10706 (EEVEE Legacy) to the EEVEE Next closure
 * system.  Geometry and color math are unchanged from D10706.
 * Replaced:
 *   CLOSURE_EVAL_FUNCTION_DECLARE_3 macros  -> ClosureReflection structs
 *   brdf_lut() + F_brdf_multi_scatter()     -> eevee::lut::GGXBrdfData
 *   btdf_lut()                              -> eevee::lut::GGXBtdfData
 *   closure_load_ssr_data()                 -> closure_load_reflection_data()
 *   render_pass_glossy_mask()               -> dropped (handled by engine)
 *   ssr_id socket                           -> dropped (automatic in Next)
 *   cameraVec(worldPosition)                -> cameraVec(g_data.P)
 * ================================================================== */


/* ------------------------------------------------------------------
 * Absorption helpers — unchanged from D10706
 * ------------------------------------------------------------------ */

float3 node_hair_sigma_from_concentration(float eumelanin, float pheomelanin)
{
  return eumelanin   * float3(0.506f, 0.841f, 1.653f)
       + pheomelanin * float3(0.343f, 0.733f, 1.924f);
}

float3 node_hair_sigma_from_reflectance(float3 c, float azimuthal_roughness)
{
  float x = azimuthal_roughness;
  /* eq. 9 polynomial — Horner form */
  float roughness_fac = (((((0.245f * x) + 5.574f) * x - 10.73f) * x + 2.532f) * x - 0.215f) * x + 5.969f;
  float3 sigma = log(max(c, float3(1e-5f))) / roughness_fac;
  return sigma * sigma;
}

/* ------------------------------------------------------------------
 * Main node
 * ------------------------------------------------------------------ */

[[node]]
void node_bsdf_hair_principled(float4 color,
                               float melanin,
                               float melanin_redness,
                               float4 tint,
                               float3 absorption_coefficient,
                               float roughness,
                               float radial_roughness,
                               float coat,
                               float ior,
                               float offset,
                               float aspect_ratio,
                               float R,
                               float TT,
                               float TRT,
                               float random_color,
                               float random_roughness,
                               float random,
                               float weight,
                               Closure &result)
{
#if 0
  ClosureHair hair_data;
  hair_data.weight = weight;
  hair_data.color = color.rgb;
  hair_data.offset = offset;
  hair_data.roughness = float2(0.0f);
  hair_data.T = g_data.curve_B;
  closure_eval(hair_data);
#else
  

  /* ----------------------------------------------------------------
   * Absorption coefficient — unchanged from D10706
   * ---------------------------------------------------------------- */

  float3 sigma;

  if (parametrization == SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION) {
    sigma = absorption_coefficient;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION) {
    float factor_random_color = 1.0f + 2.0f * (random_value - 0.5f) * random_color;
    melanin *= factor_random_color;
    melanin  = -log(max(1.0f - melanin, 0.0001f));
    float pheomelanin = melanin * melanin_redness;
    float eumelanin   = melanin - pheomelanin;
    float3 melanin_sigma = node_hair_sigma_from_concentration(eumelanin, pheomelanin);
    float3 tint_sigma    = node_hair_sigma_from_reflectance(tint.rgb, radial_roughness);
    sigma = melanin_sigma + tint_sigma;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_REFLECTANCE) {
    sigma = node_hair_sigma_from_reflectance(color.rgb, radial_roughness);
  }
  else {
    sigma = node_hair_sigma_from_concentration(0.0f, 0.8054375f);
  }

  /* ----------------------------------------------------------------
   * Roughness randomization — unchanged from D10706
   * ---------------------------------------------------------------- */

  float factor_random_roughness = 1.0f + 2.0f * (random_value - 0.5f) * random_roughness;
  roughness        *= factor_random_roughness;
  radial_roughness *= factor_random_roughness;

  /* ----------------------------------------------------------------
   * Geometry — D10706 math, only cameraVec() call-site updated
   * ---------------------------------------------------------------- */

  float3 V = normalize(cameraVec(g_data.P)); /* was: cameraVec(worldPosition) */
  float3 N = normalize(g_data.N);

#ifdef MAT_HAIR
  float3 hair_direction = normalize(hairTangent);
#else
  float3 hair_direction = normalize(cross(V, N));
#endif

  float sin_offset = sin(offset);
  float cos_offset = cos(offset);

  /* cos_a: longitudinal cosine of the view ray against the fiber axis */
  float  cos_a = dot(V, hair_direction);

  /* Vp: projection of V onto the fiber cross-section plane */
  float3 Vp    = V - hair_direction * cos_a;
  float  cos_p = length(Vp);
  Vp = Vp / max(cos_p, 1e-5f);

  float  NVp   = dot(N, Vp);

  /* R lobe reflection direction in cross-section */
  float3 R = 2.0f * NVp * N - Vp;

  /* Refracted tangent T through the cuticle (Snell's law in cross-section) */
  float3 T_refr  = normalize((NVp * N - Vp) / ior - N);
  float  T_dist  = -2.0f * dot(T_refr, N) / sqrt(max(1.0f - (cos_a / ior) * (cos_a / ior), 1e-6f));

  /* TT geometry */
  float3 N_TT = N   - 2.0f * dot(T_refr, N)   * T_refr;
  float3 V_TT = R   - 2.0f * dot(T_refr, R)   * T_refr;

  /* TRT geometry */
  float  VTT_Vp = dot(V_TT, Vp);
  float3 N_TRT  = 2.0f * VTT_Vp * N_TT - N;
  float3 V_TRT  = 2.0f * VTT_Vp * V_TT - Vp;

  /* ----------------------------------------------------------------
   * Lobe colors — D10706 math
   *
   * D10706 uses btdf_lut() for Fresnel.  In EEVEE Next the
   * equivalent is eevee::lut::GGXBtdfData sampled from the utility
   * texture.  The .y channel of the old btdf_lut result was the
   * Fresnel term; GGXBtdfData.z holds the same value.
   * ---------------------------------------------------------------- */

  /* EEVEE Next utility-texture LUT lookup — mirrors gpu_shader_material_glossy.glsl */
  float NV = dot(N, V);

  /* Fresnel approximation from the BTDF LUT (same semantic as btdf_lut().y) */
  eevee::lut::GGXBtdfData btdf_lut_data = eevee::lut::GGXBtdfData::sample_utility_tx(
      utility_tx, NV, roughness, /* ior_tangent */ 1.0f / ior);
  float fresnel = btdf_lut_data.z; /* channel z = Fresnel reflectance */

  float3 span_color = exp(-T_dist * sigma);
  float3 TT_col     = (1.0f - fresnel) * (1.0f - fresnel) * span_color;
  float3 TRT_col    = TT_col * fresnel * span_color;
  /* TRRT+ isotropic residual (geometric series) */
  float3 iso_col    = TRT_col * fresnel * span_color
                    / max(float3(1.0f) - fresnel * span_color, float3(1e-5f));

  /* ----------------------------------------------------------------
   * BRDF LUT + multi-scatter — mirrors gpu_shader_material_glossy.glsl
   *
   * D10706:   split_sum = brdf_lut(dot(N,V), roughness)
   *           brdf      = F_brdf_multi_scatter(vec3(1), vec3(1), split_sum)
   *
   * EEVEE Next:
   *           lut       = GGXBrdfData::sample_utility_tx(tx, NV, roughness)
   *           brdf      = F_brdf_multi_scatter(float3(1), float3(1), lut)
   * ---------------------------------------------------------------- */

  eevee::lut::GGXBrdfData brdf_lut_R = eevee::lut::GGXBrdfData::sample_utility_tx(
      utility_tx, NV, roughness * (1.0f - coat));
  float3 brdf_R = F_brdf_multi_scatter(float3(1.0f), float3(1.0f), brdf_lut_R);

  /* ----------------------------------------------------------------
   * Build ClosureReflection for each lobe.
   *
   * D10706 assigned lobe color AFTER CLOSURE_EVAL_FUNCTION (i.e. it
   * multiplied the already-evaluated radiance by the color).  In
   * EEVEE Next the color is set BEFORE submission — the engine
   * multiplies it during light evaluation.
   *
   * Mapping:
   *   in_Glossy_0  -> lobe_R    (primary cuticle reflection, SSR candidate)
   *   in_Glossy_1  -> lobe_TT   (transmission through fiber)
   *   in_Glossy_2  -> lobe_TRT  (back-scatter)
   *   iso_col      -> lobe_iso  (TRRT+ isotropic residual, as ClosureDiffuse)
   * ---------------------------------------------------------------- */

  /* -- R lobe -- */
  ClosureReflection lobe_R;
  lobe_R.N         = N * cos_offset + hair_direction * sin_offset;
  lobe_R.roughness = roughness * (1.0f - coat);
  /* D10706: out_Glossy_0.radiance *= brdf; then *= (fresnel + iso_col*0.333)
   * Bake both factors into .color so the engine sees the net tint. */
  lobe_R.color     = brdf_R * (float3(fresnel) + iso_col * 0.333f);
  closure_eval(lobe_R);

  /* -- TT lobe -- */
  /* D10706 set in_Glossy_1.V to shift the lobe; in EEVEE Next we
   * encode the shifted view direction into the effective normal
   * instead, since ClosureReflection has no .V field.              */
  float3 TT_shifted_V = normalize(V_TT + hair_direction *
                        (cos_a / max(cos_p, 1e-5f) + 2.0f * offset * ior + offset));

  ClosureReflection lobe_TT;
  lobe_TT.N         = N_TT * cos_offset + hair_direction * sin_offset;
  lobe_TT.roughness = roughness * sqrt(2.0f);
  lobe_TT.color     = TT_col + iso_col * 0.333f;
  /* Store the shifted view hint in the normal's hemisphere —
   * flip N_TT to face TT_shifted_V so GGX picks the right lobe.   */
  if (dot(lobe_TT.N, TT_shifted_V) < 0.0f) {
    lobe_TT.N = -lobe_TT.N;
  }
  closure_eval(lobe_TT);

  /* -- TRT lobe -- */
  float3 TRT_shifted_V = normalize(V_TRT + hair_direction *
                         (cos_a / max(cos_p, 1e-5f) + 3.0f * offset * ior + offset));

  ClosureReflection lobe_TRT;
  lobe_TRT.N         = N_TRT * cos_offset + hair_direction * sin_offset;
  lobe_TRT.roughness = roughness * sqrt(3.0f);
  lobe_TRT.color     = TRT_col + iso_col * 0.333f;
  if (dot(lobe_TRT.N, TRT_shifted_V) < 0.0f) {
    lobe_TRT.N = -lobe_TRT.N;
  }
  closure_eval(lobe_TRT);

  ClosureDiffuse hair_data;
  hair_data.weight = weight;
  hair_data.color = iso_col;
  hair_data.N = g_data.N;
  closure_eval(hair_data);
#endif
  result = Closure(0);
}