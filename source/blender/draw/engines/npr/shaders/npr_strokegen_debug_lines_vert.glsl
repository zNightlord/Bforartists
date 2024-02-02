
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)

/* https://www.shadertoy.com/view/7tlXR4 */
vec3 hash32(vec2 p) 
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}
vec3 hsl2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
vec3 rand_col_rgb(uint seed0, uint seed1)
{
    vec3 rnd = hash32(vec2(float(seed0 * 17), float(seed1 * 33))); 

    float hue = rnd.x;
    float saturation = 0.6 + rnd.z*0.4;
    float luminosity  = 0.6 + rnd.y*0.4;

    return hsl2rgb(vec3(hue, saturation, luminosity)); 
}

void main()
{
    uint vid = gl_VertexID;
    uint line_id = vid / 2u; 
    
    uint ld_addr = vid % 2u == 0u ? line_id*6u : line_id*6u+3u;
    uvec3 vpos_enc; /* world space pos */
    vpos_enc.x = ssbo_dbg_lines_[ld_addr+0u]; 
    vpos_enc.y = ssbo_dbg_lines_[ld_addr+1u];
    vpos_enc.z = ssbo_dbg_lines_[ld_addr+2u];

    vec4 vpos = vec4(uintBitsToFloat(vpos_enc).xyz, 1.0f);

	mat4 world_to_view = ubo_view_matrices_.viewmat;
    mat4 proj = ubo_view_matrices_.winmat;
    vec4 pos_hclip = world_to_view * vpos; 
    pos_hclip = proj * vec4(pos_hclip.xyz, 1.0f); 

    gl_Position = pos_hclip;
    // gl_Position.z -= 2.0e-5 * pos_hclip.w;


    color.rgb = rand_col_rgb(line_id, line_id * 7u);
    color.a = 1.0f; 

    // if ((line_id % 2u) == 0u) 
    //     color = vec4(0, 1, 0, 1);
    // else 
    //     color = vec4(1, 0, 0, 1); 

    // if ((line_id % 3u) == 0u) 
    //     color = vec4(0, 1, 0, 1);
    // else if ((line_id % 3u) == 1u)
    //     color = vec4(1, 0, 0, 1); 
    // else
    //     color = vec4(0, 0, 1, 1); 
}
