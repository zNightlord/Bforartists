/* SPDX-FileCopyrightText: 2024 Blender Authors
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "gpu_shader_compat.hh"

/*
 * gpu_shader_material_hair_principled.glsl
 * EEVEE Next — Principled Hair BSDF
 *
 * PRIMARY REFERENCE
 *   Chiang, Bitterli, Tappan, Burley (CBTB16)
 *   "A Practical and Controllable Hair and Fur Model for
 *    Production Path Tracing", Eurographics 2016.
 *
 * SECONDARY REFERENCES
 *   d'Eon, François, Hill, Letteri, Aubry (dFH11)
 *   "An Energy-Conserving Hair Reflectance Model", EGSR 2011.
 *
 *   Marschner, Jensen, Cammarano, Worley, Hanrahan (MJC03)
 *   "Light Scattering from Human Hair Fibers", SIGGRAPH 2003.
 *
 * WHAT THIS IMPLEMENTS vs. D10706 / prior EEVEE work
 *
 *  [+] Exact σ reparameterization      paper §4.2 eq. 9
 *  [+] Exact roughness reparametrize   paper §4.1 eq. 7–8
 *  [+] Logistic azimuthal dist. D_p    paper §3.3   (was: none)
 *  [+] d'Eon longitudinal M_p + log_I0 dFH11        (was: none)
 *  [+] Near-field h from hairThickTime paper §3.2   (was: none)
 *  [+] h-dependent γ_i, γ_t, Φ(p,h)   paper §3.2   (was: none)
 *  [+] Full dielectric Fresnel          MJC03        (was: LUT)
 *  [+] Exact per-lobe A_p(h) with T    paper §3     (was: approx)
 *  [+] Lobe color = A_p · M_p · N_p    paper eq. 2  (was: partial)
 *  [+] TRRT+ ClosureDiffuse (eq. 6)    paper §3.4   (was: Glossy*0.333)
 *  [+] Coat on R lobe only              paper §4.3   (was: same)
 *  [+] Bravais IOR for cylinder         MJC03        (was: none)
 *
 * REMAINING EEVEE APPROXIMATIONS
 *  [-] Per-light BRDF uses GGX, not true logistic x M_p product.
 *      ClosureReflection cannot accept a custom BRDF kernel.
 *      Mitigation: map v -> GGX alpha_long,  s -> GGX alpha_azim.
 *  [-] M_p evaluated at specular peak (theta_o = theta_i), not
 *      integrated over the full exitant hemisphere.
 *  [-] N_p evaluated at view-space azimuth, baked into lobe color.
 *  [-] No anisotropic closure split yet — single alpha = sqrt(al * az).
 *      TODO: switch to ClosureAnisotropicReflection once EEVEE Next
 *            exposes a .T tangent field in the closure struct.
 */

#ifndef VOLUMETRICS

/* ------------------------------------------------------------------
 * Parametrization enum — must match node_shader_bsdf_hair_principled.cc
 * ------------------------------------------------------------------ */
#define SHD_PRINCIPLED_HAIR_REFLECTANCE           0
#define SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION 1
#define SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION     2

/* ==================================================================
 * 1.  MATH / CONSTANTS
 * ================================================================== */

#define H_PI      3.14159265358979323846f
#define H_TWO_PI  6.28318530717958647692f
#define H_INV_PI  0.31830988618379067154f

float h_sqr(float x) { return x * x; }
float3 h_sqr(float3 x) { return x * x; }
float h_cub(float x) { return x * x * x; }

float h_safe_asin(float x) { return asin(clamp(x, -1.0f, 1.0f)); }

/* Wrap angle into [-pi, pi] */
float h_wrap_phi(float phi)
{
  phi = mod(phi + H_PI, H_TWO_PI);
  return phi - H_PI;
}

/* ==================================================================
 * 2.  LOGISTIC AZIMUTHAL DISTRIBUTION  (paper §3.3)
 *
 * Replaces d'Eon's wrapped Gaussian:
 *   - Analytically integrable  -> exact normalization over [-pi, pi]
 *   - Invertible CDF           -> perfect importance sampling
 *   - Single evaluation needed (vs. up to 5 wrapped Gaussians)
 *
 * f(x; s) = exp(-|x|/s) / (s * (1 + exp(-|x|/s))^2)
 * ================================================================== */

float logistic_pdf(float x, float s)
{
  float ex = exp(-abs(x) / s);
  return ex / (s * h_sqr(1.0f + ex));
}

float logistic_cdf(float x, float s)
{
  return 1.0f / (1.0f + exp(-x / s));
}

/* D_p(phi; s) — trimmed logistic PDF, normalized over phi in [-pi, pi].
 * norm = 2*CDF(pi; s) - 1 = integral over the support.             */
float trimmed_logistic(float phi, float s)
{
  float norm = 2.0f * logistic_cdf(H_PI, s) - 1.0f;
  return logistic_pdf(phi, s) / max(norm, 1e-6f);
}

/* ==================================================================
 * 3.  LONGITUDINAL SCATTERING  M_p  (d'Eon et al. 2011)
 *
 * Von Mises-Fisher distribution parameterized by variance v:
 *   M_p(v, ti, to) ~ exp(log_I0(cos_ti*cos_to/v)
 *                        - (1 + sin_ti*sin_to)/v)
 * ================================================================== */

/* Numerically stable log(I0(x)):
 *   x > 12  -> asymptotic  x + 0.5*log(1/(2*pi*x))
 *   x <= 12 -> Taylor series sum (x^2/4)^i / (i!)^2             */
float log_I0(float x)
{
  if (x > 12.0f) {
    return x + 0.5f * log(1.0f / (H_TWO_PI * x));
  }
  float val  = 0.0f;
  float term = 1.0f;
  float x2_4 = 0.25f * x * x;
  for (int i = 1; i <= 10; i++) {
    term *= x2_4 / float(i * i);
    val  += term;
  }
  return log(1.0f + val);
}

/* d'Eon M_p:
 *   v          : per-lobe longitudinal variance (from eq. 7)
 *   sin/cos_ti : incident  longitudinal angle components
 *   sin/cos_to : exitant   longitudinal angle components          */
float hair_M_p(float v,
               float sin_ti, float cos_ti,
               float sin_to, float cos_to)
{
  float inv_v  = 1.0f / max(v, 1e-6f);
  float log_mp = log_I0(cos_ti * cos_to * inv_v)
               - (1.0f + sin_ti * sin_to) * inv_v
               - log(2.0f * max(v, 1e-6f));
  return exp(log_mp);
}

/* ==================================================================
 * 4.  ROUGHNESS REPARAMETERIZATION  (paper §4.1, eq. 7–8)
 * ================================================================== */

/* eq. 7: beta_M -> d'Eon longitudinal variance v
 *   v = (0.726*bM + 0.812*bM^2 + 3.7*bM^20)^2                    */
float hair_beta_M_to_v(float beta_M)
{
  float b     = beta_M;
  float inner = 0.726f * b + 0.812f * h_sqr(b) + 3.7f * pow(b, 20.0f);
  return h_sqr(inner);
}

/* eq. 8: beta_N -> logistic scale s
 *   s = 0.265*bN + 1.194*bN^2 + 5.372*bN^22                      */
float hair_beta_N_to_s(float beta_N)
{
  float b = beta_N;
  return 0.265f * b + 1.194f * h_sqr(b) + 5.372f * pow(b, 22.0f);
}

/* d'Eon variance v -> GGX alpha (longitudinal axis).
 * Gaussian lobe with variance v maps to GGX via alpha ~ sqrt(2v). */
float hair_v_to_gxx_alpha(float v)
{
  return clamp(sqrt(2.0f * v), 0.001f, 1.0f);
}

/* Logistic scale s -> GGX alpha (azimuthal axis).
 * Logistic std dev = s*pi/sqrt(3); normalize by pi -> [0,1].      */
float hair_s_to_gxx_alpha(float s)
{
  return clamp((s * H_PI / sqrt(3.0f)) * H_INV_PI, 0.001f, 1.0f);
}

/* Combined isotropic roughness from separate long/azim alphas.
 * Geometric mean avoids either axis dominating.
 *
 * TODO: replace with ClosureAnisotropicReflection, passing hair_T as
 *       the tangent axis and alpha_long/alpha_azim directly, once
 *       EEVEE Next exposes the anisotropic closure to material nodes. */
float hair_combined_alpha(float alpha_long, float alpha_azim)
{
  return sqrt(alpha_long * alpha_azim);
}

/* ==================================================================
 * 5.  ABSORPTION COEFFICIENT  (paper §4.2, eq. 9)
 *
 * beta_N enters eq. 9 because azimuthal roughness acts as a phase
 * function — low beta_N -> forward scatter -> darker/more saturated.
 * ================================================================== */

float3 hair_sigma_from_concentration(float eumelanin, float pheomelanin)
{
  return eumelanin   * float3(0.506f, 0.841f, 1.653f)
       + pheomelanin * float3(0.343f, 0.733f, 1.924f);
}

/* eq. 9 — exact polynomial, beta_N must be post-randomization      */
float3 hair_sigma_from_reflectance(float3 C, float beta_N)
{
  float b = beta_N;
  float denom =  5.969f
              - (0.215f * b)
              + (2.532f * h_sqr(b))
              - (10.73f * h_cub(b))
              + (5.574f * b * b * b * b)
              + (0.245f * b * b * b * b * b);
  float3 s = log(max(C, float3(1e-5f))) / max(denom, 1e-5f);
  return s * s;
}

/* ==================================================================
 * 6.  DIELECTRIC FRESNEL  (full unpolarized — replaces LUT)
 *
 * cos_theta_i : cosine of incidence in air (n = 1)
 * eta         : n_fiber / n_air
 * ================================================================== */

float hair_fresnel(float cos_theta_i, float eta)
{
  float sin2_t = (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
  if (sin2_t >= 1.0f) {
    return 1.0f; /* Total internal reflection */
  }
  float cos_t = sqrt(max(1.0f - sin2_t, 0.0f));
  float rs = (cos_theta_i - eta * cos_t) / (cos_theta_i + eta * cos_t);
  float rp = (eta * cos_theta_i - cos_t) / (eta * cos_theta_i + cos_t);
  return 0.5f * (rs * rs + rp * rp);
}

/* ==================================================================
 * 7.  NEAR-FIELD FIBER GEOMETRY  (paper §3.2, MJC03 §3)
 *
 * h in [-1, 1]: normalized cross-section hit offset from fiber axis.
 *   h = 0  -> fiber center
 *   h = +/-1 -> fiber edges
 *
 * Derived quantities:
 *   gamma_i = arcsin(h)          incidence angle in cross-section
 *   gamma_t = arcsin(h / eta')   refracted angle  (Bravais IOR eta')
 *   Phi(p, h)                    exit azimuth per lobe (MJC03 eq. 3)
 *
 * EEVEE Next: h comes from hairThickTime in [0, 1]:  h = 2t - 1
 * ================================================================== */

/* Modified (Bravais) IOR for cylindrical geometry.
 * Corrects eta for the longitudinal component of the incident ray.
 *   eta'(theta) = sqrt(eta^2 - sin^2(theta)) / cos(theta)          */
float hair_bravais_ior(float eta, float sin_theta)
{
  float num = max(eta * eta - sin_theta * sin_theta, 0.0f);
  float den = max(1.0f - sin_theta * sin_theta, 1e-6f);
  return sqrt(num / den);
}

/* Exit azimuth Phi(p, h)  (MJC03 eq. 3):
 *   Phi = 2p*gamma_t - 2*gamma_i + p*pi
 *   p=0 -> R,  p=1 -> TT,  p=2 -> TRT                             */
float hair_exit_phi(int p, float gamma_i, float gamma_t)
{
  return float(p) * (2.0f * gamma_t + H_PI) - 2.0f * gamma_i;
}

/* Beer-Lambert absorption along one refracted chord:
 *   T(sigma, gamma_t) = exp(-2*sigma / cos(gamma_t))
 * Factor of 2: chord spans the full fiber diameter.                */
float3 hair_absorption_T(float3 sigma, float gamma_t)
{
  return exp(-2.0f * sigma / max(cos(gamma_t), 1e-5f));
}

/* ==================================================================
 * 8.  PER-LOBE ATTENUATIONS  A_p(h)
 *
 *   A_R   = F
 *   A_TT  = (1-F)^2 * T
 *   A_TRT = (1-F)^2 * T^2 * F
 *   A_iso = (1-F)^2 * F^2 * T^3 / (1 - F*T)   [paper eq. 6]
 * ================================================================== */

void hair_attenuations(float F,
                       float3 T_abs,
                       out float3 A_R,
                       out float3 A_TT,
                       out float3 A_TRT,
                       out float3 A_iso)
{
  float  omF = 1.0f - F;
  float3 fT  = float3(F) * T_abs;

  A_R   = float3(F);
  A_TT  = h_sqr(omF) * T_abs;
  A_TRT = A_TT * fT;
  /* eq. 6: A_fourth = A_TRT * F*T / (1 - F*T) */
  A_iso = A_TRT * fT / max(float3(1.0f) - fT, float3(1e-5f));
}

/* ==================================================================
 * 9.  MAIN NODE FUNCTION
 * ================================================================== */

[[node]]
void node_bsdf_hair_principled(
    float4 base_color,
    float  melanin,
    float  melanin_redness,
    float4 tint,
    float3 absorption_coefficient,
    float  roughness,         /* beta_M — longitudinal  (0-1) */
    float  radial_roughness,  /* beta_N — azimuthal     (0-1) */
    float  coat,              /* cuticle sheen, R lobe only (paper §4.3) */
    float  ior,               /* fiber IOR, ~1.55 for human hair */
    float  offset_angle,      /* cuticle tilt alpha, in radians */
    float  random_color,
    float  random_roughness,
    float  random_value,
    float3 N,
    float  parametrization,
    Closure &result)
{
  result = CLOSURE_DEFAULT;

  /* ==============================================================
   * A.  ABSORPTION COEFFICIENT sigma  (paper §4.2)
   * ============================================================== */

  float3 sigma;
  if (parametrization == SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION) {
    sigma = absorption_coefficient;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION) {
    float rand_col = 1.0f + 2.0f * (random_value - 0.5f) * random_color;
    melanin = clamp(melanin * rand_col, 0.0f, 1.0f);
    melanin = -log(max(1.0f - melanin, 1e-4f));
    float ph = melanin * melanin_redness;
    float eu = melanin - ph;
    /* Melanin sigma + optional tint sigma (dye bath on top) */
    sigma = hair_sigma_from_concentration(eu, ph)
          + hair_sigma_from_reflectance(tint.rgb, radial_roughness);
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_REFLECTANCE) {
    sigma = hair_sigma_from_reflectance(base_color.rgb, radial_roughness);
  }
  else {
    /* Fallback: medium-brown (pheomelanin ~0.8, no eumelanin) */
    sigma = hair_sigma_from_concentration(0.0f, 0.8054375f);
  }

  /* ==============================================================
   * B.  PER-FIBER ROUGHNESS RANDOMIZATION
   * ============================================================== */

  float rand_r = 1.0f + 2.0f * (random_value - 0.5f) * random_roughness;
  float beta_M = clamp(roughness        * rand_r, 0.0f, 1.0f);
  float beta_N = clamp(radial_roughness * rand_r, 0.0f, 1.0f);

  /* ==============================================================
   * C.  ROUGHNESS REPARAMETERIZATION  (paper §4.1, eq. 7–8)
   *
   * Variance is additive across scattering events:
   *   v_R   = v(beta_M_coat)   coat shrinks only the R-lobe (§4.3)
   *   v_TT  = 2 * v0           two longitudinal scattering events
   *   v_TRT = 3 * v0           three longitudinal scattering events
   *
   * Azimuthal scale s is the same for all lobes.
   * ============================================================== */

  float v0    = hair_beta_M_to_v(beta_M);
  float v_R   = hair_beta_M_to_v(beta_M * max(1.0f - coat, 0.0f));
  float v_TT  = 2.0f * v0;
  float v_TRT = 3.0f * v0;
  float s     = hair_beta_N_to_s(beta_N);

  /* GGX alphas per axis */
  float al_R   = hair_v_to_gxx_alpha(v_R);
  float al_TT  = hair_v_to_gxx_alpha(v_TT);
  float al_TRT = hair_v_to_gxx_alpha(v_TRT);
  float al_az  = hair_s_to_gxx_alpha(s);

  /* Combined isotropic roughness (geometric mean of long/azim) */
  float rough_R   = hair_combined_alpha(al_R,   al_az);
  float rough_TT  = hair_combined_alpha(al_TT,  al_az);
  float rough_TRT = hair_combined_alpha(al_TRT, al_az);

  /* ==============================================================
   * D.  FIBER GEOMETRY — tangent + longitudinal angle theta
   * ============================================================== */

  float3 V = normalize(cameraVec(g_data.P));
  N = normalize(N);

#ifdef HAIR_SHADER
  float3 hair_T = normalize(hairTangent);
#else
  /* Mesh hair cards fallback */
  float3 hair_T = normalize(cross(V, N));
#endif

  /* sin(theta_i) = dot(V, T)  — longitudinal component */
  float sin_theta_i = clamp(dot(V, hair_T), -1.0f, 1.0f);
  float cos_theta_i = max(sqrt(1.0f - sin_theta_i * sin_theta_i), 1e-5f);

  /* Cross-section projection of V (unit vector perpendicular to hair_T) */
  float3 Vp     = V - hair_T * sin_theta_i;
  float  Vp_len = length(Vp);
  Vp = (Vp_len > 1e-5f) ? (Vp / Vp_len) : float3(0.0f, 1.0f, 0.0f);

  /* Bravais IOR: corrects eta for the cylindrical cuticle geometry  */
  float eta_p = hair_bravais_ior(ior, sin_theta_i);

  /* ==============================================================
   * E.  NEAR-FIELD FIBER OFFSET  h  (paper §3.2)
   *
   * hairThickTime in [0, 1]  ->  h = 2t - 1 in [-1, 1].
   * Each pixel on the fiber curve gets a unique h, giving spatially
   * varying reflectance that matches the near-field path-traced
   * ground truth when the fiber is resolved at >= 1 pixel width.
   * ============================================================== */

#ifdef HAIR_SHADER
  float h = clamp(hairThickTime * 2.0f - 1.0f, -1.0f, 1.0f);
#else
  float h = clamp(dot(N, Vp), -1.0f, 1.0f);
#endif

  float gamma_i = h_safe_asin(h);
  float gamma_t = h_safe_asin(clamp(h / max(eta_p, 1e-5f), -1.0f, 1.0f));

  /* ==============================================================
   * F.  FRESNEL + ABSORPTION -> PER-LOBE ATTENUATIONS  A_p(h)
   * ============================================================== */

  float  F     = hair_fresnel(cos(gamma_i), ior);
  float3 T_abs = hair_absorption_T(sigma, gamma_t);

  float3 A_R, A_TT, A_TRT, A_iso;
  hair_attenuations(F, T_abs, A_R, A_TT, A_TRT, A_iso);

  /* ==============================================================
   * G.  LONGITUDINAL WEIGHTS  M_p  (d'Eon 2011)
   *
   * Evaluated at the specular peak (theta_o = theta_i).
   * The GGX roughness spreads this in the longitudinal direction.
   * Accurate at moderate roughness; underestimates grazing-angle
   * contributions at very low roughness (acceptable for real-time).
   * ============================================================== */

  float M_R   = hair_M_p(v_R,   sin_theta_i, cos_theta_i,
                                 sin_theta_i, cos_theta_i);
  float M_TT  = hair_M_p(v_TT,  sin_theta_i, cos_theta_i,
                                 sin_theta_i, cos_theta_i);
  float M_TRT = hair_M_p(v_TRT, sin_theta_i, cos_theta_i,
                                 sin_theta_i, cos_theta_i);

  /* ==============================================================
   * H.  AZIMUTHAL WEIGHTS  N_p  (paper §3.3 — logistic)
   *
   * N_p(phi) = trimmed_logistic(phi - Phi(p, h),  s)
   *
   * phi_view : azimuth of V in the fiber cross-section frame
   * Phi(p,h) : predicted exit azimuth for lobe p at this h
   * s        : logistic scale from eq. 8
   * ============================================================== */

  /* Cross-section reference frame  (N_cs perpendicular to hair_T) */
  float3 N_cs    = N - hair_T * dot(N, hair_T);
  float  N_cs_l  = length(N_cs);
  N_cs = (N_cs_l > 1e-5f) ? normalize(N_cs) : Vp;
  float3 T_cs = cross(N_cs, hair_T); /* second axis of cross-section */

  /* Signed azimuth of V in the cross-section frame */
  float phi_view = atan(dot(Vp, T_cs), dot(Vp, N_cs));

  /* h-dependent exit azimuths, wrapped to [-pi, pi] */
  float phi_R   = h_wrap_phi(hair_exit_phi(0, gamma_i, gamma_t));
  float phi_TT  = h_wrap_phi(hair_exit_phi(1, gamma_i, gamma_t));
  float phi_TRT = h_wrap_phi(hair_exit_phi(2, gamma_i, gamma_t));

  /* Logistic azimuthal weight at the current view azimuth */
  float N_R   = trimmed_logistic(h_wrap_phi(phi_view - phi_R),   s);
  float N_TT  = trimmed_logistic(h_wrap_phi(phi_view - phi_TT),  s);
  float N_TRT = trimmed_logistic(h_wrap_phi(phi_view - phi_TRT), s);

  /* ==============================================================
   * I.  LOBE COLORS  =  A_p * M_p * N_p  (paper eq. 2)
   *
   * S_p(theta_i, theta_o, phi) = M_p * N_p * A_p
   *
   *   A_p — wavelength-dependent (absorption + Fresnel)
   *   M_p — longitudinal profile (d'Eon, energy-conserving scalar)
   *   N_p — azimuthal profile    (logistic, energy-conserving)
   *
   * iso (TRRT+): isotropic azimuth by paper §3.4, no M_p or N_p.
   * ============================================================== */

  float3 color_R   = A_R   * (M_R   * N_R);
  float3 color_TT  = A_TT  * (M_TT  * N_TT);
  float3 color_TRT = A_TRT * (M_TRT * N_TRT);
  float3 color_iso = A_iso;

  /* ==============================================================
   * J.  LOBE NORMALS WITH CUTICLE TILT  (offset_angle alpha)
   *
   * Cuticle scales tilt by alpha from the fiber axis, shifting the
   * effective normal for each lobe in the N_cs / hair_T plane.
   * ============================================================== */

  float sin_off = sin(offset_angle);
  float cos_off = cos(offset_angle);

  /* Refracted vector in cross-section plane for TT / TRT normals */
  float  NVp      = dot(N_cs, Vp);
  float3 T_refr   = normalize((NVp * N_cs - Vp) / max(ior, 1e-5f) - N_cs);

  float3 N_TT_cs  = normalize(N_cs - 2.0f * dot(T_refr, N_cs) * T_refr);
  float3 V_TT_cs  = normalize(Vp   - 2.0f * dot(T_refr, Vp)   * T_refr);
  float  vvp      = dot(V_TT_cs, Vp);
  float3 N_TRT_cs = normalize(2.0f * vvp * N_TT_cs - N_cs);

  float3 eff_N_R   = normalize(N_cs    * cos_off + hair_T * sin_off);
  float3 eff_N_TT  = normalize(N_TT_cs * cos_off + hair_T * sin_off);
  float3 eff_N_TRT = normalize(N_TRT_cs * cos_off + hair_T * sin_off);

  /* ==============================================================
   * K.  EEVEE NEXT — CLOSURE SUBMISSION
   *
   * R / TT / TRT -> ClosureReflection  (GGX specular)
   * TRRT+        -> ClosureDiffuse     (isotropic irradiance gather)
   *
   * ClosureDiffuse for iso is correct because:
   *   Paper §3.4 proves TRRT+ has isotropic azimuthal distribution.
   *   ClosureDiffuse gathers irradiance uniformly over all solid
   *   angles — exactly matching the isotropic fourth lobe semantics.
   *   Using ClosureReflection (D10706 approach) biases it toward
   *   specular light directions, which is physically wrong.
   * ============================================================== */

  ClosureReflection lobe_R;
  lobe_R.N         = eff_N_R;
  lobe_R.roughness = rough_R;
  lobe_R.color     = color_R;
  closure_load_reflection_data(lobe_R, result);

  ClosureReflection lobe_TT;
  lobe_TT.N         = eff_N_TT;
  lobe_TT.roughness = rough_TT;
  lobe_TT.color     = color_TT;
  closure_load_reflection_data(lobe_TT, result);

  ClosureReflection lobe_TRT;
  lobe_TRT.N         = eff_N_TRT;
  lobe_TRT.roughness = rough_TRT;
  lobe_TRT.color     = color_TRT;
  closure_load_reflection_data(lobe_TRT, result);

  /* TRRT+ residual — diffuse irradiance gather, isotropic          */
  ClosureDiffuse lobe_iso;
  lobe_iso.N     = N;
  lobe_iso.color = color_iso;
  closure_load_diffuse_data(lobe_iso, result);
}

#else /* VOLUMETRICS */
/* clang-format off */
/* Not compatible with volumetric contexts. */
#define node_bsdf_hair_principled(base_color, melanin, melanin_redness, tint, absorption_coefficient, roughness, radial_roughness, coat, ior, offset_angle, random_color, random_roughness, random_value, N, parametrization, result) \
    (result = CLOSURE_DEFAULT)
/* clang-format on */
#endif /* VOLUMETRICS */
