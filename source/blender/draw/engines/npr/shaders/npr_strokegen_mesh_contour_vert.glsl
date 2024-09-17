
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_brush_toolbox_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_debug_view_lib.glsl)




#if defined(_KERNEL_MULTICOMPILE_DRAW_CONTOUR__EDGES)
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
    uint contour_seg_list_head = move_elem_along_loop(
        contour_edge_id, -int(contour_seg_rank), contour_edge_list_head, contour_edge_list_len
    ); 
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


    color = // cf.occluded ? vec4(.5f, .5f, .5f, 0.0f) : 
        cf.cusp_func_pstv ? vec4(.0f, 1.0f, 0.0f, 1.0f) : vec4(1.0f, 0.0f, 0.0f, 1.0f); 
        // cf.dbg_flag_0 ? vec4(.0f, 1.0f, 0.0f, 1.0f) : vec4(1.0f, 0.0f, 0.0f, 1.0f); 
        // vec4(rand_col_rgb(contour_seg_list_head, contour_seg_list_head), 1.0f); // contour_edge_id);
        // vec4(rand_col_rgb(contour_seg_len / 12, contour_seg_len / 12), 1.0f); // contour_edge_id);
        // contour_seg_len == contour_edge_list_len ? vec4(1.0f, 0.0f, 0.0f, 1.0f) : vec4(0.0f, 1.0f, 0.0f, 1.0f); 
        // (contour_edge_list_len < contour_edge_rank || contour_edge_rank == 0) ?
        // vec4(prev_contour_id, contour_edge_rank, contour_edge_list_len, contour_edge_list_head) : vec4(.0f);


    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    gl_Position = pos_hclip;
    /* Apply depth bias to curve breaks,
     * which can happen due to Z-fighting artifacts. */
    /* TODO: Further imporve this */
    gl_Position.z -= 1.0e-5 * whclip;
    gl_Position.xy += edgenor_uv * whclip * pcs_screen_size_inv_ * 1.0f;
	vec2 edge_dir_ext = (vid % 2u == 1u) ? edgedir_uv : -edgedir_uv;
	gl_Position.xy += edge_dir_ext * whclip * pcs_screen_size_inv_ * 1.5f;
}
#endif


#if defined(_KERNEL_MULTICOMPILE_DRAW_CONTOUR__2D_SAMPLES)
void main()
{
    uint vid = gl_VertexID;
    uint spine_id = vid / POINTS_PER_WING_QUAD; // == sample_id
    uint num_2d_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 

    SpineTopoInfo sti = load_ssbo_stroke_mesh_pool__skeletal_topo_info(vid, num_2d_samples);
    vec2 coord = load_ssbo_stroke_mesh_pool__skeletal_VB(vid, sti); 
    vec4 col = load_ssbo_stroke_mesh_pool__skeletal_color(vid, num_2d_samples);

    vec4 pos_hclip = vec4(coord * 2.0f - 1.0f, 0.0f, 1.0f); 
    gl_Position.xyzw = pos_hclip;

    color = col; // vec4(1.0f, 0.0f, 0.0f, 1.0f); 
    // normal = vec3(0, 0, 1);
    // tangent.xyz = vec3(0, 0, 1);
}
#endif