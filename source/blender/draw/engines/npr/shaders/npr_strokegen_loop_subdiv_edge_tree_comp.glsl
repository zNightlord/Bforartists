
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)


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
                // note: the opposite vert of the parent edge could change to v1/v3 due to edge split
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








