
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)

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
    
    uint line_id_local = vid / 2u; 
    uint line_id = line_id_local + get_debug_line_offset(pcs_line_type_); 

    DebugVertData vtx_data = load_debug_vtx_data(line_id, vid); 
    vec4 vpos = vec4(vtx_data.pos.xyz, 1.0f);

	mat4 world_to_view = ubo_view_matrices_.viewmat;
    mat4 proj = ubo_view_matrices_.winmat;
    vec4 pos_hclip = world_to_view * vpos; 
    pos_hclip = proj * vec4(pos_hclip.xyz, 1.0f); 

    float z_bias = 5e-4f; 
    if (vtx_data.dbg_data.x == DBG_LINE_SPEC_DATA_X__BOTTOM_LAYER) 
        z_bias = 4e-4f; 

    gl_Position = pos_hclip;
    gl_Position.z -= z_bias * pos_hclip.w;

    
    color.a = 1.0f; 
    
    color.rgb = vtx_data.col; 

	// if (vtx_data.dbg_data.x != 0u)
	// 	color.r = vtx_data.dbg_data.x; 
	// if (vtx_data.dbg_data.y != 0u)
	// 	color.g = vtx_data.dbg_data.y; 
	// if (vtx_data.dbg_data.z != 0u)
	// 	color.b = vtx_data.dbg_data.z; 
	// if (vtx_data.dbg_data.w != 0u)
	// 	color.a = vtx_data.dbg_data.w; 

    // color = vec4(0, 1, 0, 1);
    // color *= .6f; 

}
