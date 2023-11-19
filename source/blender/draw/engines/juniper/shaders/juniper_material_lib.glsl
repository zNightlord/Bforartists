
#ifndef JNPR_MATERIAL_LIB
#define JNPR_MATERIAL_LIB

// Utility functions

bool mask_light_groups(LightData ld) {
    return !((ld.light_groups.x & light_groups.x) == 0 &&
    (ld.light_groups.y & light_groups.y) == 0 &&
    (ld.light_groups.z & light_groups.z) == 0 &&
    (ld.light_groups.w & light_groups.w) == 0);
}

#if defined(GPU_FRAGMENT_SHADER) && !defined(JUNIPER_PREPASS)
vec3 get_prepass_normal(vec2 uvs)
{
    return texture(in_normal_tx, uvs).xyz;
}

float get_prepass_depth(vec2 uvs)
{
    return texture(in_depth_tx, uvs).x;
}
#else
vec3 get_prepass_normal(vec2 uvs)
{
    return vec3(0.0);
}

float get_prepass_depth(vec2 uvs)
{
    return 0.0;
}
#endif // GPU_FRAGMENT_SHADER

#endif // JNPR_MATERIAL_LIB
