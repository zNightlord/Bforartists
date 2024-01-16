
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)

void main()
{
    uint vid = gl_VertexID;
    uint line_id = vid / 2u; 
    
    uint ld_addr = vid % 2u == 0u ? line_id*6u : line_id*6u+3u;
    uvec3 vpos_enc; /* world space pos */
    Load3(ssbo_debug_lines_, ld_addr, vpos_enc); 
    vec4 vpos = vec4(uintBitsToFloat(vpos_enc).xyz, 1.0f);

	mat4 world_to_view = ubo_view_matrices_.viewmat;
    mat4 proj = ubo_view_matrices_.winmat;
    vec4 pos_hclip = world_to_view * vpos; 
    pos_hclip = proj * vec4(pos_hclip.xyz, 1.0f); 

    gl_Position = pos_hclip;

    color = vec4(1, 0, 0, 1); 
}
