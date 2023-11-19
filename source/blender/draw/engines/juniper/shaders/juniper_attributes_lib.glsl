
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(gpu_shader_codegen_lib.glsl)


#ifdef OBINFO_LIB
vec3 attr_load_orco(vec4 orco)
{
    /* We know when there is no orco layer when orco.w is 1.0 because it uses the generic vertex
     * attribute (which is [0,0,0,1]). */
    if (orco.w == 1.0) {
        /* If the object does not have any deformation, the orco layer calculation is done on the fly
         * using the orco_madd factors. */
        return OrcoTexCoFactors[0].xyz + pos * OrcoTexCoFactors[1].xyz;
    }
    return orco.xyz * 0.5 + 0.5;
}
#endif

vec4 attr_load_tangent(vec4 tangent)
{
    tangent.xyz = safe_normalize(normal_object_to_world(tangent.xyz));
    return tangent;
}
vec4 attr_load_vec4(vec4 attr)
{
    return attr;
}
vec3 attr_load_vec3(vec3 attr)
{
    return attr;
}
vec2 attr_load_vec2(vec2 attr)
{
    return attr;
}
float attr_load_float(float attr)
{
    return attr;
}
vec4 attr_load_color(vec4 attr)
{
    return attr;
}
vec3 attr_load_uv(vec3 attr)
{
    return attr;
}
