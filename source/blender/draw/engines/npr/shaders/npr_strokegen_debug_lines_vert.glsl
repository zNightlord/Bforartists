
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_debug_view_lib.glsl)



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

    float z_bias = 3e-4f; 
    if (vtx_data.dbg_data.x == DBG_LINE_SPEC_DATA_X__BOTTOM_LAYER) 
        z_bias = 2e-4f; 

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
