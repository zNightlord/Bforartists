
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
    uint vid = gl_VertexID;
    uint contour_edge_id = vid / 2u; 
    
    uint base_addr = vid * 6; 
    uint num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 

    /* Fetch vert pos from ssbo */
    /* uint addr_vpos = mesh_pool_addr__wpos(contour_edge_id); 
    if (vid % 2u == 1u) 
        addr_vpos += 3u; 
    vec3 pos_ws_curr = vec3(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+1]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_vpos+2])
    );
    vec4 pos_ndc = point_world_to_ndc(pos_ws_curr); */

    uint addr_uv = mesh_pool_addr__edgeuv(contour_edge_id, num_contour_edges); 
    if (vid % 2u == 1u) 
        addr_uv += 2u;
    vec2 uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+1])
    );

    uint addr_zhclip = mesh_pool_addr__zwhclip(contour_edge_id, num_contour_edges); 
    if (vid % 2u == 1u) 
        addr_zhclip += 2u; 
    float zhclip = uintBitsToFloat(buf_strokegen_mesh_pool[addr_zhclip+0]); 
    float whclip = uintBitsToFloat(buf_strokegen_mesh_pool[addr_zhclip+1]); 


    vec3 pos_view = get_view_space_from_depth(uv, zhclip);

    vec4 pos_hclip; 
    pos_hclip.xy = uv * 2.0f - 1.0f;
    pos_hclip.z = zhclip;
    pos_hclip.w = whclip; 
    pos_hclip.xy *= pos_hclip.w; 

    /* Fetch edge dir from ssbo */
    uint addr_edgedir = mesh_pool_addr__edgedir(contour_edge_id, num_contour_edges); 
    vec2 edgedir_uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+1])
    );
    vec2 edgenor_uv = vec2(edgedir_uv.y, -edgedir_uv.x); 
    color = vec4(edgedir_uv.xy * 0.5f + .5f, 0, 1); 


    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    gl_Position = pos_hclip;
    /* Apply depth bias to curve breaks, 
     * which can happen due to Z-fighting artifacts. */
    /* TODO: Further imporve this */
    gl_Position.z -= 1.0e-6 * whclip;
    gl_Position.xy += edgenor_uv * whclip * pcs_screen_size_inv_ * 1.2f; 
}
