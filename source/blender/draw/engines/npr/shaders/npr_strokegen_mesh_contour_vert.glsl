
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
    uint vid = gl_VertexID;
    uint contour_edge_id = vid / 2u; 
    
    uint base_addr = vid * 6; 
    uint num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 

    /* Fetch vert pos from ssbo */
    uint addr_vpos = mesh_pool_addr__wpos(contour_edge_id); 
    if (vid % 2u == 1u) 
        addr_vpos += 3u; 
    vec3 pos_ws_curr = vec3(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+1]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+2])
    );
    vec4 pos_ndc = point_world_to_ndc(pos_ws_curr); 

    /* Fetch edge dir from ssbo */
    uint addr_edgedir = mesh_pool_addr__edgedir(contour_edge_id, num_contour_edges); 
    vec2 edgedir_uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+1])
    );
    color = vec4(edgedir_uv.xy * 0.5f + .5f, 0, 1); 


    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    vec4 pos_gl = pos_ndc;
    gl_Position = pos_gl;
    /* Apply depth bias to curve breaks, 
     * which can happen due to Z-fighting artifacts. */
    gl_Position.z -= 0.00005/* 0.00005 */;
}
