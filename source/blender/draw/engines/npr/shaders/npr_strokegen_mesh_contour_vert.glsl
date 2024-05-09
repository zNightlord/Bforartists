
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)

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
    uint contour_edge_id = vid / 2u; 
    
    uint base_addr = vid * 6; 
    uint num_contour_lines = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 


    // Debug stroke topo
    uint contour_edge_rank = ssbo_contour_snake_rank_[contour_edge_id]; 
    uint contour_edge_list_len = ssbo_contour_snake_list_len_[contour_edge_id]; 
    float contour_edge_param = float(contour_edge_rank) / float(contour_edge_list_len); 
    uint contour_edge_list_head = ssbo_contour_snake_list_head_[contour_edge_id]; 
    uint contour_seg_rank = ssbo_contour_snake_seg_rank_[contour_edge_id];
    uint contour_seg_len = ssbo_contour_snake_seg_len_[contour_edge_id];
    uint prev_contour_id = ssbo_contour_to_contour_[2*contour_edge_id];
    uint next_contour_id = ssbo_contour_to_contour_[2*contour_edge_id+1u];
    ContourFlags cf = load_contour_flags(contour_edge_id); 
    
    uint end_node_id = contour_edge_list_head; 
    uvec2 link_end = uvec2(
        ssbo_contour_to_contour_[2*end_node_id], 
        ssbo_contour_to_contour_[2*end_node_id+1u]
    ); 
	uint end_node_rank = ssbo_contour_snake_rank_[end_node_id]; 

    




    uint draw_contour_edge_id = contour_edge_id;
    uint addr_uv = mesh_pool_addr__edgeuv(draw_contour_edge_id); 
    uint wedge_id = buf_strokegen_mesh_pool[addr_uv+4u];

    if (vid % 2u == 1u) 
        addr_uv += 2u;
    vec2 uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+1])
    );


    uint addr_uv_another = mesh_pool_addr__edgeuv(draw_contour_edge_id);
    if (vid % 2u == 0u) 
        addr_uv_another += 2u;
    vec2 uv_another_vtx = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv_another+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv_another+1])
    );
    float edge_len = length((uv - uv_another_vtx) / pcs_screen_size_inv_.xy); 


    uint addr_zhclip = mesh_pool_addr__zwhclip(draw_contour_edge_id); 
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
    uint addr_edgedir = mesh_pool_addr__edgedir(draw_contour_edge_id); 
    vec2 edgedir_uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+1])
    );
    // goes CCW on screen. 
    vec2 edgenor_uv = vec2(edgedir_uv.y, -edgedir_uv.x); 


    color = /* cf.occluded ? vec4(1.0f, 0.0f, 0.0f, 0.0f) :  */
        vec4(rand_col_rgb(contour_seg_len / 16, contour_seg_len / 16), contour_edge_id);
        // (contour_edge_list_len < contour_edge_rank || contour_edge_rank == 0) ? 
        // vec4(prev_contour_id, contour_edge_rank, contour_edge_list_len, contour_edge_list_head) : vec4(.0f);


    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    gl_Position = pos_hclip;
    /* Apply depth bias to curve breaks, 
     * which can happen due to Z-fighting artifacts. */
    /* TODO: Further imporve this */
    gl_Position.z -= 8.0e-5 * whclip;
    gl_Position.xy += edgenor_uv * whclip * pcs_screen_size_inv_ * 1.0f; 
	vec2 edge_dir_ext = (vid % 2u == 1u) ? edgedir_uv : -edgedir_uv; 
	gl_Position.xy += edge_dir_ext * whclip * pcs_screen_size_inv_ * 1.5f; 
}
