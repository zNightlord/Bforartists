/* SPDX-FileCopyrightText: 2023 Blender Authors
 * SPDX-License-Identifier: GPL-2.0-or-later */

/* Principled Hair BSDF - EEVEE Next implementation.
 * Three-lobe model: R (specular), TT (transmission), TRT (back-scatter).
 * Reference: d'Eon et al., "An Energy-Conserving Hair Reflectance Model", EGSR 2011.
 */

#ifndef VOLUMETRICS

/* --------------------------------------------------------------------------
 * Parametrization enum (must match DNA / node_shader_hair_principled.cc)
 * -------------------------------------------------------------------------- */
#define SHD_PRINCIPLED_HAIR_REFLECTANCE           0
#define SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION 1
#define SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION     2

vec3 hair_sigma_from_concentration(float eumelanin, float pheomelanin)
{
  /* Fitted coefficients from Donner & Jensen 2006. */
  return eumelanin  * vec3(0.506, 0.841, 1.653)
       + pheomelanin * vec3(0.343, 0.733, 1.924);
}

vec3 hair_sigma_from_reflectance(vec3 c, float azimuthal_roughness)
{
  float x = azimuthal_roughness;
  /* Polynomial fit mapping roughness → sigma scaling factor. */
  float b = (((((0.245 * x) + 5.574) * x - 10.73) * x + 2.532) * x - 0.215) * x + 5.969;
  vec3 s = log(max(c, vec3(1e-5))) / b;
  return s * s;
}

/* Eq. 7 from paper: perceptual βM → d'Eon longitudinal variance ν.
 * We then sqrt to get a GGX-style σ for the longitudinal axis. */
float hair_long_roughness(float beta_M)
{
  float b = beta_M;
  float nu = sqr(0.726 * b + 0.812 * sqr(b) + 3.7 * pow(b, 20.0));
  /* ν is variance; convert to roughness-like sigma for GGX. */
  return sqrt(nu);
}

/* Eq. 8 from paper: βN → logistic scale s.
 * Map to GGX roughness as a stand-in until logistic is implemented. */
float hair_azim_roughness(float beta_N)
{
  float b = beta_N;
  /* s from paper eq. 8 */
  float s = 0.265 * b + 1.194 * sqr(b) + 5.372 * pow(b, 22.0);
  /* Logistic σ ≈ GGX roughness at moderate values; rough stand-in.
   * TODO: replace GGX lobe with actual logistic azimuthal distribution. */
  return clamp(s * 0.2, 0.0, 1.0);
}

/* Blend longitudinal and azimuthal roughness into a single isotropic
 * GGX roughness for lobes where we can't pass separate axes yet.
 * This is the stopgap until ClosureInputGlossy gains a tangent field. */
float hair_effective_roughness(float rough_long, float rough_azim)
{
  /* Geometric mean — avoids either axis dominating. */
  return sqrt(rough_long * rough_azim);
}

/* --------------------------------------------------------------------------
 * Main node function
 * -------------------------------------------------------------------------- */

[[node]]
void node_bsdf_hair_principled(
    vec4  color,
    float melanin,
    float melanin_redness,
    vec4  tint,
    vec3  absorption_coefficient,
    float roughness,
    float radial_roughness,
    float coat,
    float ior,
    float offset_angle,
    float random_color,
    float random_roughness,
    float random_value,
    float weight,
    vec3  N,
    float parametrization,
    const float do_diffuse,   /* always 0 for hair – kept for uniform node signature */
    const float do_reflection,
    const float do_transmission,
    out Closure result)
{
  result = CLOSURE_DEFAULT;

  /* ------------------------------------------------------------------ */
  /* 1. Resolve absorption coefficient σ from chosen parametrization.   */
  /* ------------------------------------------------------------------ */
  vec3 sigma;

  if (parametrization == SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION) {
    sigma = absorption_coefficient;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION) {
    float rand_scale = 1.0 + 2.0 * (random_value - 0.5) * random_color;
    melanin *= rand_scale;
    melanin  = -log(max(1.0 - melanin, 1e-4));
    float pheomelanin = melanin * melanin_redness;
    float eumelanin   = melanin - pheomelanin;
    vec3 melanin_sigma = hair_sigma_from_concentration(eumelanin, pheomelanin);
    vec3 tint_sigma    = hair_sigma_from_reflectance(tint.rgb, radial_roughness);
    sigma = melanin_sigma + tint_sigma;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_REFLECTANCE) {
    sigma = hair_sigma_from_reflectance(base_color.rgb, radial_roughness);
  }
  else {
    /* Fallback: medium-brown (matches Cycles default). */
    sigma = hair_sigma_from_concentration(0.0, 0.8054375);
  }

  /* ------------------------------------------------------------------ */
  /* 2. Apply per-fiber random roughness variation.                      */
  /* ------------------------------------------------------------------ */
  float rand_rough = 1.0 + 2.0 * (random_value - 0.5) * random_roughness;
  roughness        *= rand_rough;
  radial_roughness *= rand_rough;

  /* ------------------------------------------------------------------ */
  /* 3. Anisotropic roughness — paper section 4.1                       */
  /* ------------------------------------------------------------------ */
  float rough_long = hair_long_roughness(roughness);        /* βM axis */
  float rough_azim = hair_azim_roughness(radial_roughness); /* βN axis */

  /* Per-lobe effective roughness.
   * R:   coat reduces the long axis only (cuticle smoothness, paper §4.3)
   * TT:  variance adds in quadrature across 2 scattering events → sqrt(2)
   * TRT: variance adds across 3 events → sqrt(3)
   * These multipliers apply to the longitudinal axis; azimuthal stays βN. */
  float rough_R_long   = rough_long * (1.0 - coat);
  float rough_TT_long  = rough_long * sqrt(2.0);
  float rough_TRT_long = rough_long * sqrt(3.0);

  float rough_R_eff   = hair_effective_roughness(rough_R_long,   rough_azim);
  float rough_TT_eff  = hair_effective_roughness(rough_TT_long,  rough_azim);
  float rough_TRT_eff = hair_effective_roughness(rough_TRT_long, rough_azim);

  /* ------------------------------------------------------------------ */
  /* 3. Geometry: derive axial direction and cross-section vectors.      */
  /* ------------------------------------------------------------------ */
  vec3 V = normalize(cameraVec(g_data.P));   /* EEVEE Next: g_data.P */
  N = normalize(N);

#ifdef HAIR_SHADER
  vec3 T = normalize(hairTangent);           /* curve tangent from attr */
#else
  /* Fallback for mesh previews: reconstruct from V × N. */
  vec3 T = normalize(cross(V, N));
#endif

  float sin_off = sin(offset_angle);
  float cos_off = cos(offset_angle);

  /* Decompose V into axial (∥ T) and cross-section (⊥ T) parts. */
  float cos_theta = dot(V, T);              /* longitudinal cosine */
  vec3  Vp = V - T * cos_theta;            /* projection onto cross-section plane */
  float sin_theta = length(Vp);            /* = cos of incidence angle in cross-sec */
  Vp = Vp / max(sin_theta, 1e-5);         /* unit vector in cross-section plane */

  /* ------------------------------------------------------------------ */
  /* 4. R, TT, TRT lobe normals and view vectors (Marschner 2003).      */
  /* ------------------------------------------------------------------ */
  float NdotVp  = dot(N, Vp);

  /* R: specular reflection off cuticle surface. */
  vec3 N_R  = N;
  vec3 R_dir = 2.0 * NdotVp * N - Vp;          /* reflected in cross-section */

  /* Refracted direction into the fiber. */
  vec3  T_refr  = normalize((NdotVp * N - Vp) / ior - N);
  float T_dist  = -2.0 * dot(T_refr, N) / max(sqrt(1.0 - sqr(cos_theta / ior)), 1e-5);

  /* TT: transmission through fiber. */
  vec3 N_TT  = N - 2.0 * dot(T_refr, N) * T_refr;
  vec3 V_TT  = R_dir - 2.0 * dot(T_refr, R_dir) * T_refr;

  /* TRT: back-scatter (two refractions + one internal reflection). */
  float VTT_dot_Vp = dot(V_TT, Vp);
  vec3 N_TRT = 2.0 * VTT_dot_Vp * N_TT - N;
  vec3 V_TRT = 2.0 * VTT_dot_Vp * V_TT - Vp;

  /* ------------------------------------------------------------------ */
  /* 5. Lobe colors from absorption + Fresnel.                          */
  /* ------------------------------------------------------------------ */
  float fresnel    = btdf_lut(NdotVp, roughness, ior).y;
  vec3  span_col   = exp(-T_dist * sigma);
  vec3  TT_col     = sqr(1.0 - fresnel) * span_col;
  vec3  TRT_col    = TT_col * fresnel * span_col;
  /* TRRT+ isotropic residual (energy conservation tail). */
  vec3  denom_safe = max(vec3(1.0) - fresnel * span_col, vec3(1e-5));
  vec3  iso_col    = TRT_col * fresnel * span_col / denom_safe;

  /* ------------------------------------------------------------------ */
  /* 6. Tilted lobe normals                                              */
  /* ------------------------------------------------------------------ */
  vec3 eff_N_R   = N       * cos_off + hair_T * sin_off;
  vec3 eff_N_TT  = N_TT   * cos_off + hair_T * sin_off;
  vec3 eff_N_TRT = N_TRT  * cos_off + hair_T * sin_off;

  /* The anisotropy tangent for each lobe lives in the cross-section
   * plane perpendicular to hair_T.  For a proper anisotropic GGX we
   * would pass hair_T directly into ClosureInputAnisotropic.T.
   *
   * TODO: switch Glossy → Anisotropic closures once the macro set and
   * ClosureInputAnisotropic expose a .T field in this build.
   * Stub shown below for reference:
   *
   *   in_Anisotropic_0.N         = eff_N_R;
   *   in_Anisotropic_0.T         = hair_T;          ← anisotropy axis
   *   in_Anisotropic_0.roughness = rough_R_long;    ← long axis
   *   in_Anisotropic_0.anisotropy = 1.0
   *     - (rough_R_eff / max(rough_R_long, 1e-5));  ← derived from azim
   */

  /* ------------------------------------------------------------------ */
  /* 6. Effective normals with longitudinal tilt (offset_angle).        */
  /* ------------------------------------------------------------------ */
  vec3 eff_N_R   = N_R   * cos_off + T * sin_off;
  vec3 eff_N_TT  = N_TT  * cos_off + T * sin_off;
  vec3 eff_N_TRT = N_TRT * cos_off + T * sin_off;

  /* ------------------------------------------------------------------ */
  /* 7. Build ClosureReflection for each lobe and accumulate.           */
  /* ------------------------------------------------------------------ */

  /* -- R lobe (direct specular, participates in SSR) -- */
  if (do_reflection != 0.0) {
    ClosureReflection lobe_R;
    lobe_R.N         = eff_N_R;
    lobe_R.roughness = rough_R_eff;
    lobe_R.color     = fresnel + iso_col * 0.333;   /* R + isotropic share */

    /* BRDFmulti-scatter energy compensation (EEVEE Next LUT). */
    vec2 split_sum = brdf_lut(NdotVp, lobe_R.roughness);
    lobe_R.color  *= F_brdf_multi_scatter(vec3(1.0), vec3(1.0), split_sum);

    /* Write to result – primary SSR lobe. */
    closure_load_reflection_data(lobe_R, result);
  }

  /* -- TT lobe (transmission through fiber, no SSR) -- */
  if (do_transmission != 0.0) {
    ClosureReflection lobe_TT;
    /* Shift effective view direction for TT path in cross-section. */
    lobe_TT.N         = eff_N_TT;
    lobe_TT.roughness = rough_TT_eff;
    lobe_TT.color     = TT_col + iso_col * 0.333;

    closure_load_reflection_data(lobe_TT, result);
  }

  /* -- TRT lobe (back-scatter, no SSR) -- */
  if (do_transmission != 0.0) {   /* grouped with transmission toggle */
    ClosureReflection lobe_TRT;
    lobe_TRT.N         = eff_N_TRT;
    lobe_TRT.roughness = rough_TRT_eff;
    lobe_TRT.color     = TRT_col + iso_col * 0.333;

    closure_load_reflection_data(lobe_TRT, result);
  }
}

#else  /* VOLUMETRICS */

/* clang-format off */
/* Not compatible with volumetric contexts – emit default closure. */
#define node_bsdf_hair_principled(base_color, melanin, melanin_redness, tint, absorption_coefficient, roughness, radial_roughness, coat, ior, offset_angle, random_color, random_roughness, random_value, N, parametrization, do_diffuse, do_reflection, do_transmission, result) \
  (result = CLOSURE_DEFAULT)
/* clang-format on */

#endif  /* VOLUMETRICS */