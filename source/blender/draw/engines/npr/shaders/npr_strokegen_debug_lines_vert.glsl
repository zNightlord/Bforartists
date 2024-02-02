
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)

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


    color = vec4(0, 1, 0, 1);

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
