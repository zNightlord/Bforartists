#ifndef VOLUMETRICS

CLOSURE_EVAL_FUNCTION_DECLARE_3(node_bsdf_hair_principled, Glossy, Glossy, Glossy)

/* principled hair parametrization */
#define SHD_PRINCIPLED_HAIR_REFLECTANCE 0
#define SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION 1
#define SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION 2

vec3 sigma_from_concentration(float eumelanin, float pheomelanin)
{
  return eumelanin * vec3(0.506, 0.841, 1.653) + pheomelanin * vec3(0.343, 0.733, 1.924);
}

vec3 sigma_from_reflectance(vec3 c, float azimuthal_roughness)
{
  float x = azimuthal_roughness;
  float roughness_fac = (((((0.245 * x) + 5.574) * x - 10.73) * x + 2.532) * x - 0.215) * x + 5.969;
  vec3 sigma = log(c) / roughness_fac;
  return sigma * sigma;
}

void node_bsdf_hair_principled(vec4 base_color,
                               float melanin,
                               float melanin_redness,
                               vec4 tint_color,
                               vec3 absorption_coefficient,
                               float roughness,
                               float radial_roughness,
                               float coat,
                               float ior,
                               float offset_angle,
                               float random_color,
                               float random_roughness,
                               float random_value,
                               vec3 N,
                               float parametrization,
                               float ssr_id,
                               out Closure result)
{
  // bool do_ssr = (ssrToggle && int(ssr_id) == outputSsrId);

  CLOSURE_VARS_DECLARE_3(Glossy, Glossy, Glossy);

  /* Convert color into absorbtion coefficient. Match Cycles. */
  vec3 sigma;

  if (parametrization == SHD_PRINCIPLED_HAIR_DIRECT_ABSORPTION) {
    sigma = absorption_coefficient;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION) {
    float factor_random_color = 1 + 2 * (random_value - 0.5) * random_color;
    melanin *= factor_random_color;
    melanin = -log(max(1 - melanin, 0.0001));
    float pheomelanin = melanin * melanin_redness;
    float eumelanin = melanin - pheomelanin;
    vec3 melanin_sigma = sigma_from_concentration(eumelanin, pheomelanin);
    vec3 tint_sigma = sigma_from_reflectance(tint_color.rgb, radial_roughness);
    sigma = melanin_sigma + tint_sigma;
  }
  else if (parametrization == SHD_PRINCIPLED_HAIR_REFLECTANCE) {
    sigma = sigma_from_reflectance(base_color.rgb, radial_roughness);
  }
  else {
    /* Fallback to brownish hair, same as defaults for melanin. */
    sigma = sigma_from_concentration(0.0, 0.8054375);
  }

  float factor_random_roughness = 1 + 2 * (random_value - 0.5) * random_roughness;
  roughness *= factor_random_roughness;
  radial_roughness *= factor_random_roughness;

  vec3 V = normalize(cameraVec(worldPosition));
  N = normalize(N);

#ifdef HAIR_SHADER
  vec3 hair_direction = normalize(hairTangent);
#else
  vec3 hair_direction = normalize(cross(V, N));
#endif

  float sin_offset = sin(offset_angle);
  float cos_offset = cos(offset_angle);

  float cos_a = dot(V, hair_direction);
  // vec3 Va = hair_direction * cos_a;
  /* Projection of the view vector V onto axial plane - hair cross section */
  vec3 Vp = V - hair_direction * cos_a;
  float cos_p = length(Vp);
  //Va = Va / cos_p; /* Scale so it can be added to the TT and TRT rays */
  Vp = Vp / cos_p; /* Normalize */

  float NVp = dot(N, Vp);
  vec3 R = 2 * NVp * N - Vp; /* |R|=1 */
  vec3 T = normalize((NVp*N - Vp) / ior - N);
  float T_dist = -2 * dot(T, N) / sqrt(1 - sqr(cos_a / ior)); /* sin_p = cos_a */
  vec3 N_TT = N - 2 * dot(T, N) * T; /* |N_TT|=1 */
  vec3 V_TT = R - 2 * dot(T, R) * T; /* |V_TT|=1 */
  vec3 N_TRT = 2 * dot(V_TT, Vp) * N_TT - N;  /* |N_TRT|=1 */
  vec3 V_TRT = 2 * dot(V_TT, Vp) * V_TT - Vp; /* |V_TRT|=1 */

  /* At each scattering incident brdf is the same, but total rougness adds up. */
  // float fresnel = F_eta(ior, dot(N, V));
  float fresnel = btdf_lut(dot(N, V), roughness, ior).y;
  vec3 span_color = exp(-T_dist * sigma);
  vec3 TT_col = sqr(1 - fresnel) * span_color;
  vec3 TRT_col = TT_col * fresnel * span_color;
  vec3 iso_col = TRT_col * fresnel * span_color / (1 - fresnel * span_color);

  in_Glossy_0.N = N * cos_offset + hair_direction * sin_offset;
  in_Glossy_0.roughness = roughness * (1 - coat);
  in_Glossy_0.ensure_valid = false;

  in_Glossy_1.N = N_TT * cos_offset + hair_direction * sin_offset;
  /* offset_angle is usually small - use linear approximation */
  in_Glossy_1.V = V_TT + hair_direction * (cos_a/cos_p + 2 * offset_angle * ior + offset_angle);
  /* Roughness is defined as somewhat like gaussian variance. Sqr(total variance)
   * is a sum of squares of variances of all scattering events. */
  in_Glossy_1.roughness = roughness * sqrt(2);
  in_Glossy_1.ensure_valid = false;

  in_Glossy_2.N = N_TRT * cos_offset + hair_direction * sin_offset;
  in_Glossy_2.V = V_TRT + hair_direction * (cos_a/cos_p + 3 * offset_angle * ior + offset_angle);
  in_Glossy_2.roughness = roughness * sqrt(3);
  in_Glossy_2.ensure_valid = false;

  CLOSURE_EVAL_FUNCTION_3(node_bsdf_hair_principled, Glossy, Glossy, Glossy);

  result = CLOSURE_DEFAULT;

  vec2 split_sum = brdf_lut(dot(N, V), in_Glossy_0.roughness);
  vec3 brdf = F_brdf_multi_scatter(vec3(1), vec3(1), split_sum);

  out_Glossy_0.radiance = closure_mask_ssr_radiance(out_Glossy_0.radiance, ssr_id);
  out_Glossy_0.radiance *= brdf;
  out_Glossy_0.radiance = render_pass_glossy_mask(vec3(1), out_Glossy_0.radiance);
  /* TODO: compute TRRT+ contribution from axial lighting. */
  out_Glossy_0.radiance *= fresnel + iso_col * 0.333;
  closure_load_ssr_data(
      out_Glossy_0.radiance, in_Glossy_0.roughness, in_Glossy_0.N, ssr_id, result);

  // out_Glossy_1.radiance = closure_mask_ssr_radiance(out_Glossy_1.radiance, ssr_id);
  // out_Glossy_1.radiance *= brdf;
  out_Glossy_1.radiance = render_pass_glossy_mask(vec3(1), out_Glossy_1.radiance);
  out_Glossy_1.radiance *= TT_col + iso_col * 0.333;
  // closure_load_ssr_data(
  //     out_Glossy_1.radiance, in_Glossy_1.roughness, in_Glossy_1.N, ssr_id, result);

  // out_Glossy_2.radiance = closure_mask_ssr_radiance(out_Glossy_2.radiance, ssr_id);
  // out_Glossy_2.radiance *= brdf;
  out_Glossy_2.radiance = render_pass_glossy_mask(vec3(1), out_Glossy_2.radiance);
  out_Glossy_2.radiance *= TRT_col + iso_col * 0.333;
  // closure_load_ssr_data(
  //     out_Glossy_2.radiance, in_Glossy_2.roughness, in_Glossy_2.N, ssr_id, result);

  result.radiance += out_Glossy_1.radiance + out_Glossy_2.radiance;
}

#else
/* clang-format off */
/* Stub because it is not compatible with volumetrics. */
#  define node_bsdf_hair_principled(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, result) (result = CLOSURE_DEFAULT)
/* clang-format on */
#endif
 