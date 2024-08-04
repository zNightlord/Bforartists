
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_debug_view_lib.glsl)


#if defined(_KERNEL_MULTICOMPILE__CALC_SUBD_TREE_NODES__PER_EDGE)
void main()
{
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 

#if defined(_KERNEL_MULTICOMPILE__CALC_SUBD_TREE_NODE_FOR_NEW_EDGES) 
    if (!valid_thread) return; 

    uvec2 iverts_cwedge = mark__cwedge_to_verts(1u);
    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge[0]]; 
    VertFlags vf_1 = load_vert_flags(v1);
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge[1]];
    VertFlags vf_3 = load_vert_flags(v3);

    EdgeFlags ef = load_edge_flags(wedge_id); 

    if ((vf_1.new_by_split != vf_3.new_by_split) && ef.new_by_split && !ef.new_by_split_on_old_edge) 
    { // an new edge that is neither an face-edge nor a split-edge, 
      // usually happens at incomplete subdived faces, around the border of subdivided region
      // invalidate the link to prevent writing into the same position of another face/split-edge
      // when we are building the tree. 
        LoopSubdEdgeTreeUpNode node = init_loop_subd_tree_leaf__face_edge(); // point to nothing 
        ssbo_subd_edge_tree_node_up_[wedge_id] = encode_loop_subd_tree_node(node);
    }
    if (!(vf_1.new_by_split && vf_3.new_by_split)) return; // only do this for new edges
    

    LoopSubdEdgeTreeUpNode node = init_loop_subd_tree_leaf__face_edge(); // point to nothing 
    
    // Step 1 - find par edge id & cache 2 opposite verts 
    uint oppo_vert_old = 0xffffffffu; 
    uint iface_oppo_vert_old = 0u; 
    for (uint iface = 0; iface < 2; ++iface) 
    { 
        uint vtx_oppo = ssbo_edge_to_vert_[wedge_id*4u + mark__center_wedge_to_oppo_vert__at_face(iface)];
        VertFlags vf = load_vert_flags(vtx_oppo); 

        if (vf.new_by_split) node.parent_edge_id = ssbo_subd_edge_vert_to_old_edge_[vtx_oppo];
        else {
            oppo_vert_old = vtx_oppo; 
            iface_oppo_vert_old = iface; 
        }
    } // Note: only one vert should be newly generated. (true == new_by_split)

    // Step 2 - calculate sub-edge code & transfer subtree orientation
    if (node.parent_edge_id != LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID && oppo_vert_old != 0xffffffffu)
    { 
        LoopSubdEdgeTreeUpNode par_node = decode_loop_subd_tree_node(ssbo_subd_edge_tree_node_up_[node.parent_edge_id]); 
        uvec2 par_oppo_verts_at_iface = uvec2(
            ssbo_edge_to_vert_[node.parent_edge_id*4u + mark__center_wedge_to_oppo_vert__at_face(0u)], 
            ssbo_edge_to_vert_[node.parent_edge_id*4u + mark__center_wedge_to_oppo_vert__at_face(1u)]
        ); 
        // try to locate new edge at which iface of the parent edge 
        uint iface_new_edge = (
                (oppo_vert_old == par_oppo_verts_at_iface[0])
                // note: rather than oppo_vert_old, 
                // the opposite vert of the parent edge could change to v1/v3 due to edge split
                || (v1 == par_oppo_verts_at_iface[0])
                || (v3 == par_oppo_verts_at_iface[0])
            ) ? 0u : 1u; 
        uint par_subd_face_id = (iface_new_edge == par_node.iface_subd_f0) ? 0u : 1u; 
        if (par_subd_face_id == 0u) 
        { // at parent edge's subd-face#0
            node.iface_subd_f0 = iface_oppo_vert_old; 
            node.code = get_loop_subd_tree_leaf_code__face_edge(0u); 
        }else 
        { // at parent edge's subd-face#1
            node.iface_subd_f0 = iface_oppo_vert_old == 0u ? 1u : 0u; // flip 
            node.code = get_loop_subd_tree_leaf_code__face_edge(1u); 
        }
    }

    if (valid_thread)
        ssbo_subd_edge_tree_node_up_[wedge_id] = encode_loop_subd_tree_node(node);
#endif

// Note: this kernel is specially called for all edges, including the old edges. 
#if defined(_KERNEL_MULTICOMPILE__BUILD_SUBD_TREE_DOWNWARD__INIT)
    if (gl_GlobalInvocationID.x == 0u)
        temporal_record_counter(pc_obj_id_, pc_frame_id_) = 0u; 

    wedge_id = gl_GlobalInvocationID.x; 
    uint num_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_edges; 
    valid_thread = wedge_id < num_edges; 

    if (!valid_thread) return; 
    // Point every edge to itself as 4 sub-nodes in the subd-tree
    uvec4 init_nodes_enc = uvec4(encode_loop_subd_tree_node_dw(LoopSubdEdgeTreeDwNode(wedge_id)));
    Store4(ssbo_subd_edge_tree_node_dw_, wedge_id, init_nodes_enc); 

    // Point each edge to a null record
    store_ssbo_edge_to_new_temporal_record_(wedge_id, PER_EDGE_TEMPORAL_REC_ID_NULL); 
    store_ssbo_edge_to_old_temporal_record_(wedge_id, PER_EDGE_TEMPORAL_REC_ID_NULL); 
#endif

#if defined(_KERNEL_MULTICOMPILE__BUILD_SUBD_TREE_DOWNWARD__TRANSFER_FROM_UPWARD)
    if (!valid_thread) return; 

    LoopSubdEdgeTreeUpNode node = decode_loop_subd_tree_node(ssbo_subd_edge_tree_node_up_[wedge_id]);
    if (node.parent_edge_id == LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID) return;
    if (node.parent_edge_id == wedge_id) return; // root node

    EdgeFlags ef = load_edge_flags(wedge_id);
    if (ef.dupli || ef.sel_border) return; 

    LoopSubdEdgeTreeDwNode par_node_dw; 
    par_node_dw.wedge_id = wedge_id; 
    ssbo_subd_edge_tree_node_dw_[node.parent_edge_id * 4u + node.code] = encode_loop_subd_tree_node_dw(par_node_dw); 
#endif
}
#endif



#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS)

#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
struct TraceDebugContext
{
    bool record_dbg_lines; 

    uvec2 dbg_line_id_beg; // beg dbg line addr on 2 paths 
    bool dbg_show_both_paths_if_matched; 
    bool dbg_show_full_path; 
    bool dbg_show_both_paths; 

    bool dbg_use_opti_0; 
}; 
TraceDebugContext init_trace_dbg_context()
{
    TraceDebugContext ctx;
    ctx.record_dbg_lines = false;  
    ctx.dbg_line_id_beg = uvec2(0u); 
    ctx.dbg_show_both_paths_if_matched = (pc_dbg_matching_line_mode_ & 1u) != 0u;
    ctx.dbg_show_full_path             = (pc_dbg_matching_line_mode_ & 2u) != 0u;
    ctx.dbg_show_both_paths            = (pc_dbg_matching_line_mode_ & 4u) != 0u;
    ctx.dbg_use_opti_0                 = (pc_dbg_matching_line_mode_ & 8u) != 0u;
    return ctx; 
}

vec3 calc_trace_start_pos(uint beg_wedge_id, vec3 cam_pos_ws)
{
    uint v1 = ssbo_edge_to_vert_[beg_wedge_id*4u + 1u]; 
    uint v3 = ssbo_edge_to_vert_[beg_wedge_id*4u + 3u]; 
    vec3 vpos_v1v3[2] = { ld_vpos(v1), ld_vpos(v3) };  

    vec3 p_beg; 
    if (0u < pcs_extract_interpo_contour_) {
        vec3 vnor_v1v3[2] = { ld_vnor(v1), ld_vnor(v3) };
        p_beg = calc_interp_contour_vert_pos(vnor_v1v3, vpos_v1v3, cam_pos_ws); 
    }
    else p_beg = (vpos_v1v3[0] + vpos_v1v3[1]) * .5f;  

    return p_beg; 
}

struct PathMatchingResult
{
    bool matched; 
    uint num_steps_til_match; 
    uint matched_rec_id; 
    TemporalRecordContourData matched_trcd; 
}; 
PathMatchingResult init_path_matching_result()
{
    PathMatchingResult match; 
    match.matched = false; 
    match.num_steps_til_match   = 0u; 
    match.matched_rec_id        = PER_EDGE_TEMPORAL_REC_ID_NULL; 
    match.matched_trcd = init_temporal_record_contour_data(); 

    return match; 
}
void match_to_history_record_on_edge(
    uint trace_iter, AdjWedgeInfo wlk_awi, uint num_history_recs, 
    inout PathMatchingResult match
){
    uint old_rec_id = load_ssbo_edge_to_old_temporal_record_(wlk_awi.wedge_id); 
    if (old_rec_id != PER_EDGE_TEMPORAL_REC_ID_NULL)
    { // try match to this edge
        TemporalRecordFlags trf = 
            load_ssbo_contour_temporal_records_old__flags(old_rec_id); 
        TemporalRecordContourData trcd = 
            load_ssbo_contour_temporal_records_old__contour_data(old_rec_id, num_history_recs); 
        
        bool matched = match_history_rec(trf, trcd); 
        if (matched) 
        { // for now we only use the first match
            // in the future we want to find the best match
            match.num_steps_til_match = trace_iter + 1; 
            match.matched_rec_id      = old_rec_id;
            match.matched_trcd        = trcd; 
        }
        match.matched = match.matched || matched; 
    }
}

struct TracedTriangle
{
    //    D0       
    //     \       
    //     q0     
    //     | \__  
    //     X    \_
    //     |\     \__
    //    Q1 \       Q2 
    //     |  \  CCW    \__     D2 
    // D1  |   D           \_    |
    //   \ |    \            \__ |    
    //    \q1 -- p ------------ q2
    // we are currently walking at p on edge q1q2(wlk_awi)
    uvec3 q; // vertex ids of q0, q1, q2
    vec3 qpos[3]; 
    vec3 Q1; 
    vec3 Q2; 
    uvec2 ibwedges_Q1Q2; 
    vec3 D[3]; 
}; 
vec3 load_trace_gradient_at_vertex(uint vid)
{
    uvec3 Dk_enc; 
    Load3(ssbo_vgrad_contour_, vid, Dk_enc);   
    return normalize(uintBitsToFloat(Dk_enc)); 
}
void load_curr_walked_triangle_geometry(
    AdjWedgeInfo wlk_awi, out TracedTriangle tri
)
{
    tri.q = uvec3( // winding: CCW
        ssbo_edge_to_vert_[
            wlk_awi.wedge_id*4u + mark__center_wedge_to_oppo_vert__at_face(wlk_awi.iface_adj)], // q0: oppo vert
        ssbo_edge_to_vert_[wlk_awi.wedge_id*4u + mark__cwedge_to_beg_vert(wlk_awi.iface_adj)], // q1: edge vert
        ssbo_edge_to_vert_[wlk_awi.wedge_id*4u + mark__cwedge_to_end_vert(wlk_awi.iface_adj)]  // q2: edge vert
    );       
    
    for (uint k = 0; k < 3; ++k) tri.qpos[k] = ld_vpos(tri.q[k]); 
    
    tri.Q1 = tri.qpos[1] - tri.qpos[0]; 
    tri.Q2 = tri.qpos[2] - tri.qpos[0]; 
    tri.ibwedges_Q1Q2 = uvec2( // iwedge ids of Q1, Q2 
        mark__cwedge_rotate_back(wlk_awi.iface_adj), 
        mark__cwedge_rotate_next(wlk_awi.iface_adj) 
    ); 

    for (uint k = 0; k < 3; ++k) {
        tri.D[k] = load_trace_gradient_at_vertex(tri.q[k]); 
    }
}
void move_triangle_to_adjacent_wedge(
    bool walk_to_Q1, bool walk_to_Q2, bool walk_to_adj_iface,  
    AdjWedgeInfo wlk_awi_old, AdjWedgeInfo wlk_awi_new, TracedTriangle tri_old, 
    out TracedTriangle tri_new
){
    // load new vtx q0: opposite vertex to the new wedge(wlk_awi_new)
    uint ivert_oppo = mark__center_wedge_to_oppo_vert__at_face(wlk_awi_new.iface_adj); 
    tri_new.q[0]    = ssbo_edge_to_vert_[wlk_awi_new.wedge_id*4u + ivert_oppo]; 
    tri_new.qpos[0] = ld_vpos(tri_new.q[0]);
    tri_new.D[0]    = load_trace_gradient_at_vertex(tri_new.q[0]); 
    
    // reuse data on old triangle to init q1, q2
    uvec2 new_q1q2_to_old; 
    if (walk_to_Q1 || walk_to_Q2) {
        new_q1q2_to_old = walk_to_Q1 ? uvec2(1, 0) : uvec2(0, 2); // yes this works
    }else if (walk_to_adj_iface) {
        new_q1q2_to_old = uvec2(2 , 1); // reverse q1q2 order
    }
    for (uint iq = 1; iq <= 2; ++iq)
    { 
        uint iq_old = new_q1q2_to_old[iq - 1u]; 
        tri_new.q[iq]    = tri_old.q[iq_old]; 
        tri_new.qpos[iq] = tri_old.qpos[iq_old]; 
        tri_new.D[iq]    = tri_old.D[iq_old]; 
    }

    tri_new.Q1 = tri_new.qpos[1] - tri_new.qpos[0]; 
    tri_new.Q2 = tri_new.qpos[2] - tri_new.qpos[0]; 
    tri_new.ibwedges_Q1Q2 = uvec2( // iwedge ids of Q1, Q2 
        mark__cwedge_rotate_back(wlk_awi_new.iface_adj), 
        mark__cwedge_rotate_next(wlk_awi_new.iface_adj) 
    ); 
}



struct PathTraceContext
{
    bool auto_init_D; 
    vec3 D_init; 

    vec3 D_prev; 
    vec3 D; 

    AdjWedgeInfo wlk_awi; 
    vec3 p; 
}; 
void calc_tracing_direction(
    inout PathTraceContext trace_ctx, TracedTriangle tri, 
    uint path, uint trace_iter, vec3 cam_pos_ws, TraceDebugContext dbg_ctx
){
    
    vec3 D; {
        if (trace_iter == 0u) 
        {
            if (!trace_ctx.auto_init_D) 
                D = trace_ctx.D_init; 
            else {
                // when we start from a interpolated contour point, 
                // since the point is a local minima of the scalar field ||dot(N,V)||
                // gradients around it - tri.D[1] and tri.D[2], or even tri.D[0] may bot be reliable              
                // So we choose the initial gradient from the most "un-contour" vtx on the triangle.
                vec3 vnor_tri[3]; 
                vec3 ndv_tri; 
                for (uint i = 0; i < 3; ++i) {
                    vnor_tri[i] = ld_vnor(tri.q[i]); 
                    ndv_tri[i] = dot(vnor_tri[i], normalize(cam_pos_ws - tri.qpos[i])); 
                }
                
                uint ivert_max = 0; 
                float abs_ndv_max = .0f; 
                for (uint i = 0; i < 3; ++i) {
                    if (abs_ndv_max < abs(ndv_tri[i])) {
                        abs_ndv_max = abs(ndv_tri[i]); 
                        ivert_max = i; 
                    }
                }
                D = tri.D[ivert_max]; 
            }
        }
        else
        {
            float w_sum = .0f; 
            vec3 dir_sum = vec3(.0f); 
            for (uint k = 0; k < 3; ++k)
                if (.0f < dot(tri.D[k], trace_ctx.D_prev)) {
                    w_sum += 1.0f; 
                    dir_sum += tri.D[k]; 
                }
            if (.0f < w_sum) D = normalize(dir_sum / w_sum); 
            else D = - normalize((tri.D[0] + tri.D[1] + tri.D[2]) / 3.0f); 
        }
    }

    if (trace_ctx.auto_init_D && path == 0u && trace_iter == 0u) trace_ctx.D_init = D; // cache initial direction
    if (trace_ctx.auto_init_D && path == 1u && trace_iter == 0u && dot(D, trace_ctx.D_init) > .0f) D = -D; // reverse direction for the second path

    trace_ctx.D = D;
}


#define TRACE_STEPS ((pc_dbg_history_trace_steps_)) // 32u
void trace_path(
    uint path, uint num_trace_steps, 
    vec3 cam_pos_ws, uint num_history_recs, 
    TraceDebugContext dbg_ctx, 
    inout PathTraceContext trace_ctx, 
    bool record_path_matching_info, inout PathMatchingResult path_matching_info
){
    uint curr_dbg_line_id = dbg_ctx.dbg_line_id_beg[path];  

    TracedTriangle tri; 
    for (uint trace_iter = 0u; trace_iter < num_trace_steps; ++trace_iter) 
    { 
        match_to_history_record_on_edge( 
            trace_iter, trace_ctx.wlk_awi, num_history_recs, 
            /*inout*/ path_matching_info
        ); 
        //    D0       
        //     \       
        //     q0     
        //     | \__  
        //     X    \_
        //     |\     \__
        //    Q1 \       Q2 
        //     |  \  CCW    \__     D2 
        // D1  |   D           \_    |
        //   \ |    \            \__ |    
        //    \q1 -- p ------------ q2
        if (trace_iter == 0u) load_curr_walked_triangle_geometry(trace_ctx.wlk_awi, /*out*/tri); 
        
        calc_tracing_direction(/*inout*/trace_ctx, tri, path, trace_iter, cam_pos_ws, dbg_ctx); 

        // Project p, D onto the plane of the triangle q012
        vec3 D_proj; vec2 uvD;
        proj_vec_to_triangle_plane(trace_ctx.D, tri.Q1, tri.Q2, /*out*/D_proj, /*out*/uvD); 
        trace_ctx.D = D_proj; 
        #define u_D ((uvD.x))
        #define v_D ((uvD.y))


        vec2 uvP; { // Easy to solve since p is on the Q1Q2 edge
            float factor = length(trace_ctx.p - tri.qpos[1]) / length(tri.qpos[2] - tri.qpos[1]); 
            uvP = vec2(1.0f - factor, factor);
        }
        #define u_p ((uvP.x))
        #define v_p ((uvP.y))

        // Intersect ray p + tD with two edges q0 + s1*Q1, q0 + s2*Q2
        float t1 = -v_p / v_D; 
        float s1 = u_p + t1 * u_D; 
        float t2 = -u_p / u_D;
        float s2 = v_p + t2 * v_D; 
        bool walk_to_Q1 = -0.000000f < s1 && s1 < 1.000000f && .0f <= t1;
        bool walk_to_Q2 = -0.000000f < s2 && s2 < 1.000000f && .0f <= t2; 
        if (walk_to_Q1 && walk_to_Q2) { // pick further one
            if (t1 < t2) walk_to_Q2 = false;
            else walk_to_Q1 = false; 
        }
        float min_t = walk_to_Q1 ? t1 : walk_to_Q2 ? t2 : .0f;
        vec3 p_next = trace_ctx.p + min_t * D_proj; 

        // the ray p+tD is leaving face q012, we can try flip to another iface
        //    q0         q0                       
        //   / \        / \                      
        //  /  iface0  /   \                      
        //q1==p==q2  q1==p==q2                    
        //     \       \ iface1                   
        //      D       \ /                       
        //               *                        
        bool walk_to_adj_iface = 
            (t1 < .0f && (1.0f <= s2 || s2 <= .0f)) 
            || (t2 < .0f && (1.0f <= s1 || s1 <= .0f)); 
        #undef u_D
        #undef v_D
        #undef u_p
        #undef v_p

        // Determine the next iface & edge 
        AdjWedgeInfo wlk_next_wedge = trace_ctx.wlk_awi; 
        if (walk_to_Q1 || walk_to_Q2)
        { // walk to new edge - either Q1 or Q2
            uint wlk_next_iwedge = walk_to_Q1 ? tri.ibwedges_Q1Q2.x : tri.ibwedges_Q1Q2.y; 
            wlk_next_wedge = decode_adj_wedge_info(
                ssbo_edge_to_edges_[trace_ctx.wlk_awi.wedge_id*4u + wlk_next_iwedge]
            ); 
        }else if (walk_to_adj_iface)
        { // keep position in current edge, but flip to the other iface
            wlk_next_wedge.iface_adj = 1u - trace_ctx.wlk_awi.iface_adj;
        }

        // Determine the next triangle 
        TracedTriangle tri_next; 
        if (walk_to_Q1 || walk_to_Q2 || walk_to_adj_iface) {
            move_triangle_to_adjacent_wedge(
                walk_to_Q1, walk_to_Q2, walk_to_adj_iface, 
                trace_ctx.wlk_awi, wlk_next_wedge, tri, 
                /*out*/ tri_next
            );
        }else {
            tri_next = tri; 
        }

        // Output debug lines  
        if (dbg_ctx.record_dbg_lines)
        {
            vec4 vpos_ws_10 = vec4(trace_ctx.p, 1.0f);
            vec4 vpos_ws_11 = vec4(p_next, 1.0f);
            bool dead_end = !(walk_to_Q1 || walk_to_Q2 || walk_to_adj_iface); 
            vec3 dvd_col = vec3(float(trace_iter) / float(TRACE_STEPS).xx, 1.0); 
            DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dvd_col, uvec4(0u)); 
            DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dvd_col, uvec4(0u));
            if (path_matching_info.matched && !dbg_ctx.dbg_show_full_path) dvd_11 = dvd_10; 
            store_debug_line_data(curr_dbg_line_id, dvd_10, dvd_11); 
            curr_dbg_line_id++; 
        }

        // Update contexts for next iteration
        if (walk_to_Q1 || walk_to_Q2 || walk_to_adj_iface) { 
            trace_ctx.wlk_awi = wlk_next_wedge;
            trace_ctx.D_prev  = trace_ctx.D;
            tri = tri_next; 
        } 
        if (walk_to_Q1 || walk_to_Q2) {
            trace_ctx.p       = p_next; 
        }
    }
}
#endif


void main()
{
#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__COMPACT)
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 

    EdgeFlags ef = load_edge_flags(wedge_id); 

    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + 1u]; 
    VertFlags vf_1 = load_vert_flags(v1); 
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + 3u]; 
    VertFlags vf_3 = load_vert_flags(v3); 


    // TODO: also support for cotour edges. 
    bool is_contour = false; 
    if (0 < pcs_extract_interpo_contour_)
        is_contour = valid_thread && is_interp_contour_edge__before_tessellation(vf_1, vf_3, ef); 
    else
    {
        mat4 view_to_world = ubo_view_matrices_.viewinv; // transform matrices, see "common_view_lib.glsl"
        vec3 cam_pos_ws = view_to_world[3].xyz; // see "#define cameraPos ViewMatrixInverse[3].xyz" 

        uint v0 = ssbo_edge_to_vert_[wedge_id*4u + 0u];
        uint v2 = ssbo_edge_to_vert_[wedge_id*4u + 2u];
        vec3 vpos[4] = { ld_vpos(v0), ld_vpos(v1), ld_vpos(v2), ld_vpos(v3) }; // yes it is bad but fuck this I dont care

        float face_orient_123, face_orient_013; 
        is_contour = valid_thread
            && is_contour_edge(
                vpos[0], vpos[1], vpos[2], vpos[3], cam_pos_ws, ef, 
                /*out*/face_orient_123, face_orient_013
            );
    }

	const uint groupId = gl_LocalInvocationID.x;
    uint rec_id = compact_temporal_contour_record(is_contour, groupId); 
    
    if (valid_thread && is_contour)
        ssbo_contour_temporal_records_new_[rec_id] = wedge_id; // store temp data, will be replaced by final output
    if (valid_thread && is_contour)
        store_ssbo_edge_to_new_temporal_record_(wedge_id, rec_id); 
#endif






#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
    uint rec_id = gl_GlobalInvocationID.x; 
    uint groupIdx = gl_LocalInvocationID.x; 
    uint num_recs = temporal_record_counter(pc_obj_id_, pc_frame_id_); 
    uint num_history_recs = temporal_record_counter(pc_obj_id_, pc_frame_id_history_); 
    bool valid_thread = rec_id < num_recs; 
    
    uint rec_wedge_id = ssbo_contour_temporal_records_new_[rec_id]; // load temp data 


    // Contract subd tree nodes on the path to root, compress the path to form a code
    // ----------------------------------------------------------------------------------------------
    uint curr_edge_id = rec_wedge_id; 
    uint code_combine = 0u; 
    for (uint tree_depth = 0; tree_depth < pcs_loop_subd_iters_; tree_depth++)
    {
        LoopSubdEdgeTreeUpNode node = decode_loop_subd_tree_node(
            ssbo_subd_edge_tree_node_up_[curr_edge_id]
        ); 
        if (node.parent_edge_id == curr_edge_id) break; // root node at curr_edge_id  

        append_subd_tree_node_to_code(node.code, /*inout*/code_combine); 
        
        curr_edge_id = node.parent_edge_id; 
    }

    TemporalRecordFlags trf = init_temporal_record_flags(code_combine);
    if (valid_thread) {
        store_ssbo_contour_temporal_records_new__flags(rec_id, trf); 
        store_ssbo_contour_temporal_records_new__subd_root_edge_id(rec_id, curr_edge_id, num_recs); 
    }


    // Match to history records 
    // Trace along the gradient field by walking on the triangular mesh
    // ----------------------------------------------------------------------------------------------
    PathMatchingResult path_matching_info[2]; 
    path_matching_info[0] = path_matching_info[1] = init_path_matching_result(); 
    
    mat4 view_to_world = ubo_view_matrices_.viewinv; // transform matrices, see "common_view_lib.glsl"
    vec3 cam_pos_ws = view_to_world[3].xyz; // see "#define cameraPos ViewMatrixInverse[3].xyz" 
    
    uint dbg_line_id = compact_general_dbg_lines(valid_thread, groupIdx, 2u * TRACE_STEPS);
    dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 
    TraceDebugContext dbg_ctx = init_trace_dbg_context(); 

    if (valid_thread)
    {
        vec3 p_beg = calc_trace_start_pos(rec_wedge_id, cam_pos_ws); 
        uint beg_wedge_id = rec_wedge_id; 
        PathTraceContext trace_ctx;

        for (uint path = 0; path < 2; ++path)
        {
            uint beg_iface = path; // TODO: classify 2 probing paths(fwd/bck) from contour curve orientation
 
            {
                trace_ctx.wlk_awi = AdjWedgeInfo(beg_wedge_id, beg_iface); 
                trace_ctx.p = p_beg; 
                trace_ctx.auto_init_D = true; 

                dbg_ctx.dbg_line_id_beg[path] = dbg_line_id + TRACE_STEPS * path; 
                dbg_ctx.record_dbg_lines = true; 
                
                trace_path(
                    path, TRACE_STEPS, 
                    cam_pos_ws, num_history_recs, dbg_ctx,
                    /*inout*/trace_ctx,  
                    /*record_path_matching_info*/true, 
                    /*inout*/path_matching_info[path]
                );
            } 
        }
    } /* if (valid_thread) */

    if (valid_thread) {
        // find the best matched history record from both search paths 
        uint best_matched_path = 0u; 
        if (path_matching_info[0].matched != path_matching_info[1].matched) 
            best_matched_path = path_matching_info[0].matched ? 0 : 1;
        else if (path_matching_info[0].matched && path_matching_info[1].matched) 
            best_matched_path = path_matching_info[0].num_steps_til_match < path_matching_info[1].num_steps_til_match ? 0 : 1;

        for (uint path = 0; path < 2; ++path) 
        {
            bool rmv_curr_path = (!path_matching_info[path].matched) && (!dbg_ctx.dbg_show_both_paths); 
            if (path_matching_info[0].matched
                && path_matching_info[1].matched 
                && !(dbg_ctx.dbg_show_both_paths_if_matched || dbg_ctx.dbg_show_both_paths)
            )
                rmv_curr_path = rmv_curr_path || (
                    path_matching_info[path].num_steps_til_match > 
                        path_matching_info[path == 0u ? 1u : 0u].num_steps_til_match
                    ); 
            
            if (dbg_ctx.record_dbg_lines)
            {
                uint update_dbg_line_id = dbg_ctx.dbg_line_id_beg[path]; 
                for (uint trace_iter = 0u; trace_iter < TRACE_STEPS; ++trace_iter)
                {
                    bool rmv_curr_line = rmv_curr_path; 
                    if (rmv_curr_line) {
                        DebugVertData dvd_rmv = DebugVertData(vec3(.0f), vec3(.0f), uvec4(0u)); 
                        store_debug_line_data(update_dbg_line_id, dvd_rmv, dvd_rmv); 
                    }else{
                        // update the color of the line
                        uint seg_key = path_matching_info[path].matched_trcd.seg_key; 
                        vec3 dvd_col = rand_col_rgb(seg_key / 8, seg_key / 8); 
                        store_debug_line_color(update_dbg_line_id, dvd_col);
                    }
                    
                    update_dbg_line_id++; 
                }
            }
        }
    }
    
    // Initialize temporal contour data 
    // --------------------------------------------------------------------------
    if (valid_thread) {
        TemporalRecordContourData trcd = init_temporal_record_contour_data(); 
        store_ssbo_contour_temporal_records_new__contour_data(rec_id, trcd, num_recs); 
    }
#endif
}

#endif




#if defined(_KERNEL_MULTICOMPILE__RECONSTRUCT_HISTORY_TEMPORAL_RECORDS)
void ld_edge_vpos(uint edge_id, out vec3 vpos_1, out vec3 vpos_3)
{
    // Generate debug lines ---
    uint v1 = ssbo_edge_to_vert_[edge_id * 4u + 1u]; 
    vpos_1 = ld_vpos(v1); 
    uint v3 = ssbo_edge_to_vert_[edge_id * 4u + 3u]; 
    vpos_3 = ld_vpos(v3); 
}
void main()
{
    uint rec_id = gl_GlobalInvocationID.x;  // history record index
    uint groupIdx = gl_LocalInvocationID.x; 
    uint num_recs = temporal_record_counter(pc_obj_id_, pc_frame_id_history_); 
    bool valid_thread = rec_id < num_recs;

    TemporalRecordFlags trf = load_ssbo_contour_temporal_records_old__flags(rec_id); 
    uint base_edge_id = load_ssbo_contour_temporal_records_old__subd_root_edge_id(rec_id, num_recs); 

    uint tree_code_chain = trf.subd_tree_code_chain; 
    uint par_edge_id = base_edge_id; 
uvec2 dbg_par_edge_ids = uvec2(0u); 
    bool hit_leaf_node = false; 
    for (uint subd_level = 0; subd_level < pc_loop_subd_iters_; ++subd_level)
    {
        uint tree_code; bool null_node;
        pop_subd_tree_code_chain(tree_code_chain, /*out*/tree_code, null_node); 
        if (null_node) break; 



        LoopSubdEdgeTreeDwNode sub_node = decode_loop_subd_tree_node_dw(
            ssbo_subd_edge_tree_node_dw_[4u * par_edge_id + tree_code]
        ); 
        hit_leaf_node = (sub_node.wedge_id == par_edge_id) || (subd_level + 1 == pc_loop_subd_iters_); 
        par_edge_id = sub_node.wedge_id;

dbg_par_edge_ids[subd_level] = par_edge_id;

        if (hit_leaf_node) break; /* arrived at a leaf node */
    }

    EdgeFlags ef = load_edge_flags(par_edge_id); 
    hit_leaf_node = hit_leaf_node && valid_thread && !(ef.del_by_split || ef.dupli); 
 

#define DRW_DEBUG_LINES 1u
#if defined(DRW_DEBUG_LINES)
    uint dbg_line_id = compact_general_dbg_lines(hit_leaf_node, groupIdx, 1u);
    dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 
#endif

    if (hit_leaf_node)
    {
        store_ssbo_edge_to_old_temporal_record_(par_edge_id, rec_id); 

#if defined(DRW_DEBUG_LINES)
        uint dbg_edge_id = par_edge_id; 
        uint v1 = ssbo_edge_to_vert_[dbg_edge_id * 4u + 1u]; 
        vec3 vpos_1 = ld_vpos(v1); 
        uint v3 = ssbo_edge_to_vert_[dbg_edge_id * 4u + 3u]; 
        vec3 vpos_3 = ld_vpos(v3);

        TemporalRecordContourData trcd = load_ssbo_contour_temporal_records_old__contour_data(rec_id, num_recs);

        uint dbg_line_id = dbg_line_id; 
        {
            vec3 vpos_ws_0 = vpos_1;
            vec3 vpos_ws_1 = vpos_3;
            vec3 dvd_col = vec3(1.0f, 1.0f, 1.0f); 
            if (trf.valid_contour_data)
                dvd_col = rand_col_rgb(trcd.seg_key / 8, trcd.seg_key / 8); 
			uvec4 dvd_data = uvec4(base_edge_id, dbg_par_edge_ids[0], dbg_par_edge_ids[1], 0); 
            DebugVertData dvd_0 = DebugVertData(vpos_ws_0.xyz, dvd_col, dvd_data); 
            DebugVertData dvd_1 = DebugVertData(vpos_ws_1.xyz, dvd_col, dvd_data);

            store_debug_line_data(dbg_line_id, dvd_0, dvd_1); 
            dbg_line_id++; 
        }
#endif
    }

}
    
#endif








#if defined(_KERNEL_MULTICOMPILE__RECORD_TEMPORAL_CONTOUR_DATA)
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
    uint groupIdx = gl_LocalInvocationID.x; 
	
    const uint contour_id = idx; 
	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	bool valid_thread = contour_id < num_contours; 

    uint rec_id = ssbo_contour_snake_to_temporal_record_[contour_id]; 
    if (rec_id == PER_CONTOUR_TEMPORAL_REC_ID_NULL) valid_thread = false; 
    
    ContourFlags cf = load_contour_flags(contour_id);
    ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf);
    uint curve_rank = ssbo_contour_snake_rank_[contour_id]; 
    uint seg_rank = ssbo_contour_snake_seg_rank_[contour_id]; 
    uint seg_len = ssbo_contour_snake_seg_len_[contour_id]; 
    uint seg_head_id = move_contour_id_along_loop(cct, contour_id, -float(seg_len)); 

    TemporalRecordContourData trcd; 
    trcd.curve_key  = cct.head_id & 0x00ffffffu; // 24 bits  
    trcd.seg_key    = seg_len/* seg_head_id */ & 0x00ffffffu;
    trcd.curve_rank = curve_rank & 0x00ffffffu;
    trcd.seg_rank   = seg_rank & 0x00ffffffu;
    trcd.cf         = cf;
    
    if (valid_thread) {
        // Inject contour data to record
        // TODO: use total count of records from all objects
        uint num_recs = temporal_record_counter(/* pc_obj_id_ */0u, pc_frame_id_); 
        store_ssbo_contour_temporal_records_new__contour_data(rec_id, trcd, num_recs); 
        // and mark it as valid
        TemporalRecordFlags trf = load_ssbo_contour_temporal_records_new__flags(rec_id); 
        trf.valid_contour_data = true; 
        store_ssbo_contour_temporal_records_new__flags(rec_id, trf); 
    }

#if defined(DRW_DEBUG_LINES)
#undef DRW_DEBUG_LINES
#endif
// #define DRW_DEBUG_LINES 1u
#if defined(DRW_DEBUG_LINES)
    uint dbg_line_id = compact_general_dbg_lines(valid_thread, groupIdx, 1u);
    dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 
    
    // TODO: encapsulate this
	uint next_contour_id = move_contour_id_along_loop(cct, contour_id, 1.0f); 
	bool non_looped_curve_tail = (cct.tail_id == contour_id) && (!cf.looped_curve); 
	vec3 vpos_0, vpos_1;
	{ // read vertex pos
		uvec3 vpos_0_enc, vpos_1_enc;
		Load3(ssbo_contour_snake_vpos_, contour_id,    	vpos_0_enc);
		Load3(ssbo_contour_snake_vpos_, next_contour_id, vpos_1_enc);
		vpos_0 = uintBitsToFloat(vpos_0_enc);
		vpos_1 = uintBitsToFloat(vpos_1_enc);
        if (cct.tail_id == contour_id && !cct.looped_curve) vpos_1 = vpos_0;
	}

    if (valid_thread) {
        vec3 vpos_ws_0 = vpos_0;
        vec3 vpos_ws_1 = vpos_1;
        vec3 dvd_col = vec3(1.0f, .0f, .0f); 
        uvec4 dvd_data = uvec4(0u);
        dvd_data.a = rec_id; 
        DebugVertData dvd_0 = DebugVertData(vpos_ws_0.xyz, dvd_col, uvec4(0u)); 
        DebugVertData dvd_1 = DebugVertData(vpos_ws_1.xyz, dvd_col, uvec4(0u));
        store_debug_line_data(dbg_line_id, dvd_0, dvd_1); 
        dbg_line_id++; 
    }
#endif
}
#endif