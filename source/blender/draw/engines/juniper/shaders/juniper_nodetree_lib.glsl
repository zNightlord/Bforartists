
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(gpu_shader_codegen_lib.glsl)

/* The Closure type is never used. Use float as dummy type. */
#define Closure float
#define CLOSURE_DEFAULT 0.0
#define closure_eval(data) CLOSURE_DEFAULT

float ambient_occlusion_eval(vec3 normal,
float distance,
const float inverted,
const float sample_count)
{
    /* TODO */
    return 1.0;
}

void attrib_load();
Closure nodetree_surface();
Closure nodetree_volume();
vec3 nodetree_displacement();
float nodetree_thickness();
vec4 closure_to_rgba(Closure cl);

void init_globals()
{
    g_data.P = interp.P;
    g_data.Ni = interp.N;
    g_data.N = safe_normalize(interp.N);
    g_data.Ng = g_data.N;
    g_data.is_strand = false;
    g_data.hair_time = 0.0;
    g_data.hair_thickness = 0.0;
    g_data.hair_strand_id = 0;
    g_data.ray_type = RAY_TYPE_CAMERA;
    g_data.ray_depth = 0.0;
    g_data.ray_length = distance(g_data.P, cameraPos);
    g_data.barycentric_coords = vec2(0.0);
    g_data.barycentric_dists = vec3(0.0);

    #ifdef GPU_FRAGMENT_SHADER
    g_data.N = (FrontFacing) ? g_data.N : -g_data.N;
    g_data.Ni = (FrontFacing) ? g_data.Ni : -g_data.Ni;
    g_data.Ng = safe_normalize(cross(dFdx(g_data.P), dFdy(g_data.P)));
    #endif

    #if defined(USE_BARYCENTRICS) && defined(GPU_FRAGMENT_SHADER) && defined(MAT_GEOM_MESH)
    g_data.barycentric_coords = gpu_BaryCoord.xy;
    g_data.barycentric_dists = barycentric_distances_get();
    #endif
}

void output_aov(vec4 color, float value, uint hash)
{
    #if defined(MAT_AOV_SUPPORT) && defined(GPU_FRAGMENT_SHADER)
    for (int i = 0; i < AOV_MAX && i < aov_buf.color_len; i++) {
        if (aov_buf.hash_color[i] == hash) {
            imageStore(aov_color_img, ivec3(gl_FragCoord.xy, i), color);
            return;
        }
    }
    for (int i = 0; i < AOV_MAX && i < aov_buf.value_len; i++) {
        if (aov_buf.hash_value[i] == hash) {
            imageStore(aov_value_img, ivec3(gl_FragCoord.xy, i), vec4(value));
            return;
        }
    }
        #endif
}

    #ifdef MAT_DISPLACEMENT_BUMP
/* Return new shading normal. */
vec3 displacement_bump()
{
    #  if defined(GPU_FRAGMENT_SHADER) && !defined(MAT_GEOM_CURVES)
    vec2 dHd;
    dF_branch(dot(nodetree_displacement(), g_data.N + dF_impl(g_data.N)), dHd);

    vec3 dPdx = dFdx(g_data.P);
    vec3 dPdy = dFdy(g_data.P);

    /* Get surface tangents from normal. */
    vec3 Rx = cross(dPdy, g_data.N);
    vec3 Ry = cross(g_data.N, dPdx);

    /* Compute surface gradient and determinant. */
    float det = dot(dPdx, Rx);

    vec3 surfgrad = dHd.x * Rx + dHd.y * Ry;

    float facing = FrontFacing ? 1.0 : -1.0;
    return normalize(abs(det) * g_data.N - facing * sign(det) * surfgrad);
    #  else
    return g_data.N;
    #  endif
}
    #endif

void fragment_displacement()
{
    #ifdef MAT_DISPLACEMENT_BUMP
    g_data.N = displacement_bump();
    #endif
}

/* -------------------------------------------------------------------- */
/** \name Coordinate implementations
 *
 * Callbacks for the texture coordinate node.
 *
 * \{ */

vec3 coordinate_camera(vec3 P)
{
    vec3 vP;
    if (false /* probe */) {
        /* Unsupported. It would make the probe camera-dependent. */
        vP = P;
    }
    else {
        #ifdef MAT_WORLD
        vP = transform_direction(ViewMatrix, P);
        #else
        vP = transform_point(ViewMatrix, P);
        #endif
    }
    vP.z = -vP.z;
    return vP;
}

vec3 coordinate_screen(vec3 P)
{
    vec3 window = vec3(0.0);
    if (false /* probe */) {
        /* Unsupported. It would make the probe camera-dependent. */
        window.xy = vec2(0.5);
    }
    else {
        /* TODO(fclem): Actual camera transform. */
        window.xy = project_point(ProjectionMatrix, transform_point(ViewMatrix, P)).xy * 0.5 + 0.5;
        // window.xy = window.xy * camera_buf.uv_scale + camera_buf.uv_bias;
    }
    return window;
}

vec3 coordinate_reflect(vec3 P, vec3 N)
{
    #ifdef MAT_WORLD
    return N;
    #else
    return -reflect(cameraVec(P), N);
    #endif
}

vec3 coordinate_incoming(vec3 P)
{
    #ifdef MAT_WORLD
    return -P;
    #else
    return cameraVec(P);
    #endif
}

/** \} */

float attr_load_temperature_post(float attr)
{
    return attr;
}
vec4 attr_load_color_post(vec4 attr)
{
    return attr;
}

vec4 attr_load_uniform(vec4 attr, const uint attr_hash)
{
    #if defined(OBATTR_LIB)
    uint index = floatBitsToUint(ObjectAttributeStart);
    for (uint i = 0; i < floatBitsToUint(ObjectAttributeLen); i++, index++) {
        if (drw_attrs[index].hash_code == attr_hash) {
            return vec4(drw_attrs[index].data_x,
            drw_attrs[index].data_y,
            drw_attrs[index].data_z,
            drw_attrs[index].data_w);
        }
    }
    return vec4(0.0);
    #else
    return attr;
    #endif
}

