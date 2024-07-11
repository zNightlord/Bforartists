
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)


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

#if defined(_KERNEL_MULTICOMPILE__BUILD_SUBD_TREE_DOWNWARD__INIT)
    if (gl_GlobalInvocationID.x == 0u)
        temporal_record_counter(pc_obj_id_, pc_frame_id_) = 0u; 

    if (!valid_thread) return; 
    // Point every mesh-edge to itself
    uvec4 init_nodes_enc = uvec4(encode_loop_subd_tree_node_dw(LoopSubdEdgeTreeDwNode(wedge_id)));
    Store4(ssbo_subd_edge_tree_node_dw_, wedge_id, init_nodes_enc); 
#endif

#if defined(_KERNEL_MULTICOMPILE__BUILD_SUBD_TREE_DOWNWARD__TRANSFER_FROM_UPWARD)
    if (!valid_thread) return; 

    LoopSubdEdgeTreeUpNode node = decode_loop_subd_tree_node(ssbo_subd_edge_tree_node_up_[wedge_id]);
    if (node.parent_edge_id == LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID) return;
    if (node.parent_edge_id == wedge_id) return; // root node

    LoopSubdEdgeTreeDwNode par_node_dw;
    par_node_dw.wedge_id = wedge_id; 
    ssbo_subd_edge_tree_node_dw_[node.parent_edge_id * 4u + node.code] = encode_loop_subd_tree_node_dw(par_node_dw); 
#endif
}
#endif



#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS)

/* Inputs: 
 * _KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS
 * INCLUDE_VERTEX_POSITION
 * int pcs_contour_detect_type_
 * int pcs_loop_subd_iters_
 * ubo_view_matrices_
 * + ssbo_bnpr_mesh_pool_counters_.num_temporal_recs
 * ssbo_contour_new_temporal_records_[]
 * ssbo_edge_to_temporal_record_[]
 * ssbo_vbo_full_[]
*/

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
        ssbo_contour_new_temporal_records_[rec_id] = wedge_id; // store temp data, will be replaced by final output
    if (valid_thread)
        store_ssbo_edge_to_new_temporal_record_(wedge_id, is_contour ? rec_id : TEMPORAL_REC_ID_NULL); 
#endif

// note: fill indirect dispatch args from ssbo_bnpr_mesh_pool_counters_.num_temporal_recs

#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
    uint rec_id = gl_GlobalInvocationID.x; 
    uint groupIdx = gl_LocalInvocationID.x; 
    uint num_recs = temporal_record_counter(pc_obj_id_, pc_frame_id_); 
    bool valid_thread = rec_id < num_recs; 
    
    uint rec_wedge_id = ssbo_contour_new_temporal_records_[rec_id]; // load temp data 


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

        // Traverse tree upwards and concatenate node code, 3 bits each level
        code_combine <<= 3u; 
        code_combine |= ((node.code << 1u) | 1u/*valid code*/); 
    }

    TemporalRecordFlags trf = init_temporal_record_flags(code_combine);
    if (valid_thread)
    {
        store_ssbo_contour_temporal_new_records__flags(rec_id, trf); 
        store_ssbo_contour_temporal_new_records__subd_root_edge_id(rec_id, curr_edge_id, num_recs); 
    }


    // Walking on the triangulation
    // ----------------------------------------------------------------------------------------------
    mat4 view_to_world = ubo_view_matrices_.viewinv; // transform matrices, see "common_view_lib.glsl"
    vec3 cam_pos_ws = view_to_world[3].xyz; // see "#define cameraPos ViewMatrixInverse[3].xyz" 
    
#define TRACE_STEPS 8u

    uint dbg_line_id = compact_general_dbg_lines(valid_thread, groupIdx, 2u * TRACE_STEPS);
    dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 

    if (valid_thread)
    {
        vec3 p_beg; { // Init as interpolated contour point
            uint v1 = ssbo_edge_to_vert_[rec_wedge_id*4u + 1u]; 
            uint v3 = ssbo_edge_to_vert_[rec_wedge_id*4u + 3u]; 
            vec3 vpos_v1v3[2] = { ld_vpos(v1), ld_vpos(v3) };  
            
            if (0u < pcs_extract_interpo_contour_) {
                vec3 vnor_v1v3[2] = { ld_vnor(v1), ld_vnor(v3) };
                p_beg = calc_interp_contour_vert_pos(vnor_v1v3, vpos_v1v3, cam_pos_ws); 
            }
            else p_beg = (vpos_v1v3[0] + vpos_v1v3[1]) * .5f; 
        }


        for (uint path = 0; path < 2; ++path)
        {
            uint beg_iface = path; // TODO: classify 2 probing paths(fwd/bck) from contour curve orientation
            uint dbg_line_id_beg = dbg_line_id + TRACE_STEPS * path; 
            
            // inputs: rec_wedge_id, beg_iface, p_beg, dbg_line_id_beg, cam_pos_ws    
            AdjWedgeInfo wlk_awi = AdjWedgeInfo(rec_wedge_id, beg_iface); 
            vec3 p = p_beg;
            uint curr_dbg_line_id = dbg_line_id_beg; 
            //  q0     
            //  | \__  
            //  X    \_
            //  |\     \__
            // Q1 \       Q2 
            //  |  \  CCW    \__  
            //  |   D           \_    
            //  |    \            \__      
            //  q1 -- p ------------ q2
            //
            vec3 D_init; vec3 D_prev; 
            for (uint trace_iter = 0u; trace_iter < TRACE_STEPS; ++trace_iter)
            {
                uvec3 q = uvec3( // winding: CCW
                    ssbo_edge_to_vert_[
                        wlk_awi.wedge_id*4u + mark__center_wedge_to_oppo_vert__at_face(wlk_awi.iface_adj)], // q0: oppo vert
                    ssbo_edge_to_vert_[wlk_awi.wedge_id*4u + mark__cwedge_to_beg_vert(wlk_awi.iface_adj)], // q1: edge vert
                    ssbo_edge_to_vert_[wlk_awi.wedge_id*4u + mark__cwedge_to_end_vert(wlk_awi.iface_adj)]  // q2: edge vert
                );       

                vec3 qpos[3]; { 
                    for (uint k = 0; k < 3; ++k) qpos[k] = ld_vpos(q[k]); 
                }
                vec3 Q1 = qpos[1] - qpos[0];
                vec3 Q2 = qpos[2] - qpos[0];
                uvec2 ibwedges_Q1Q2 = uvec2( // iwedge ids of Q1, Q2 
                    mark__cwedge_rotate_back(wlk_awi.iface_adj), 
                    mark__cwedge_rotate_next(wlk_awi.iface_adj)
                ); 

                vec3 D; {
                    uvec3 D0_enc, D1_enc, D2_enc;  
                    Load3(ssbo_vgrad_contour_, q[0], D0_enc);   
                    
                    vec3 D_[3]; 
                    D_[0] = normalize(uintBitsToFloat(D0_enc)); 
                    D_[1] = D_[2] = vec3(.0f); 
                    if (trace_iter == 0u) D = D_[0]; 
                    // when we start from a interpolated contour point, 
                    // since the point is a local minima of the scalar field ||dot(N,V)||
                    // gradients around it - D_[1] and D_[2], or even D_[0] may bot be reliable  
                    else{
                        Load3(ssbo_vgrad_contour_, q[1], D1_enc); 
                        D_[1] = normalize(uintBitsToFloat(D1_enc)); 
                        Load3(ssbo_vgrad_contour_, q[2], D2_enc); 
                        D_[2] = normalize(uintBitsToFloat(D2_enc)); 
                        
                        float w_sum = .0f; 
                        vec3 dir_sum = vec3(.0f); 
                        for (uint k = 0; k < 3; ++k)
                            if (.0f < dot(D_[k], D_prev)){
                                w_sum += 1.0f; 
                                dir_sum += D_[k]; 
                            }
                        if (.0f < w_sum) D = normalize(dir_sum / w_sum); 
                        else D = - normalize((D_[0] + D_[1] + D_[2]) / 3.0f); 
                    }
                }
                
                if (path == 0u && trace_iter == 0u) D_init = D; // save initial direction at the first path
                if (path == 1u && trace_iter == 0u) D = -D_init; // reverse direction for the second path



                // Project p, D onto the plane of the triangle q012
                vec2 uvD = proj_vec_to_triangle_plane(D, Q1, Q2); 
                #define u_D ((uvD.x))
                #define v_D ((uvD.y))
                
                vec2 uvP; { // Easy to solve since p is on the Q1Q2 edge
                    float factor = length(p - qpos[1]) / length(qpos[2] - qpos[1]); 
                    uvP = vec2(1.0f - factor, factor);
                }
                #define u_p ((uvP.x))
                #define v_p ((uvP.y))

                // Intersect ray p + tD with two edges q0 + s1*Q1, q0 + s2*Q2
                float t1 = -v_p / v_D; 
                float s1 = u_p + t1 * u_D; 
                float t2 = -u_p / u_D;
                float s2 = v_p + t2 * v_D; 
                bool walk_to_Q1 = -0.0000001f < s1 && s1 < 1.0000001f && .0f <= t1;
                bool walk_to_Q2 = -0.0000001f < s2 && s2 < 1.0000001f && .0f <= t2; 
				if (walk_to_Q1 && walk_to_Q2) { // pick further one
					if (t1 < t2) walk_to_Q2 = false;
					else walk_to_Q1 = false; 
				}
                float min_t = walk_to_Q1 ? t1 : walk_to_Q2 ? t2 : .0f;
                vec3 p_next = p + min_t * D; 

                // the ray p+tD is leaving current triangle, we can try flip to another iface
                //    *          *                        
                //   / \        / \                      
                //  /  iface0  /   \                      
                // *==p==*    *==p==*                     
                //     \       \ iface1                   
                //      D       \ /                       
                //               *                        
                bool walk_to_adj_face = (t1 < .0f && 1.0f < s2) || (t2 < .0f && 1.0f < s1); 
                #undef u_D
                #undef v_D
                #undef u_p
                #undef v_p
                
                // Find the next iface & edges
                AdjWedgeInfo wlk_next_wedge = wlk_awi; 
                if (walk_to_Q1 || walk_to_Q2)
                { // walk to new edges Q1 or Q2
                    uint wlk_next_iwedge = walk_to_Q1 ? ibwedges_Q1Q2.x : ibwedges_Q1Q2.y; 
                    wlk_next_wedge = decode_adj_wedge_info(
                        ssbo_edge_to_edges_[wlk_awi.wedge_id*4u + wlk_next_iwedge]
                    ); 
                }else if (walk_to_adj_face)
                { // keep position in current edge, but flip to the other iface
                    wlk_next_wedge.iface_adj = 1u - wlk_awi.iface_adj;
                }

                
                // Output debug lines            
                vec4 vpos_ws_10 = vec4(p, 1.0f);
                vec4 vpos_ws_11 = vec4(p_next, 1.0f);
                vec3 dvd_col = vec3(float(trace_iter) / float(TRACE_STEPS).xx, 1.0); 
                DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dvd_col); 
                DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dvd_col); 
                store_debug_line_data(curr_dbg_line_id, dvd_10, dvd_11); 
                curr_dbg_line_id++; 

                // Update contexts for next iteration
                if (walk_to_Q1 || walk_to_Q2 || walk_to_adj_face) { 
                    wlk_awi = wlk_next_wedge;
                    D_prev  = D;
                } 
                if (walk_to_Q1 || walk_to_Q2) {
                    p       = p_next; 
                }
            }
        }
    }
#endif

}



#endif







