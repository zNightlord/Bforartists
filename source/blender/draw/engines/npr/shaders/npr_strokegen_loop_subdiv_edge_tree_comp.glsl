
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)


#if defined(_KERNEL_MULTICOMPILE__CALC_SUBD_TREE_NODE_FOR_NEW_EDGES) 
/*
uint pcs_edge_count_

uint ssbo_dyn_mesh_counters_out_[]
?    ssbo_bnpr_mesh_pool_counters_
uint ssbo_selected_edge_to_edge_[]
uint ssbo_edge_to_vert_[]
uint ssbo_vert_flags_[]
uint ssbo_subd_edge_tree_node_[]

uint ssbo_subd_edge_vert_to_old_edge_[]
*/
void main()
{
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 

    if (!valid_thread) return; 


    uvec2 iverts_cwedge = mark__cwedge_to_verts(1u);
    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge[0]]; 
    VertFlags vf_1 = load_vert_flags(v1);
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge[1]];
    VertFlags vf_3 = load_vert_flags(v3);

    if (!(vf_1.new_by_split && vf_3.new_by_split)) return; // only do this for new edges
    

    LoopSubdEdgeTreeNode node = init_loop_subd_tree_leaf(); // point to nothing 
    for (uint iface = 0; iface < 2; ++iface) 
    { 
        uint vtx_oppo = ssbo_edge_to_vert_[wedge_id*4u + mark__center_wedge_to_oppo_vert__at_face(iface)];
        VertFlags vf = load_vert_flags(vtx_oppo); 

        if (vf.new_by_split){
            node.parent_edge_id = ssbo_subd_edge_vert_to_old_edge_[vtx_oppo];
            node.code = get_loop_subd_tree_leaf_code(iface); 
            break; 
        }
    } // Note: only one vert should be newly generated. (true == new_by_split)

    if (valid_thread)
        ssbo_subd_edge_tree_node_[wedge_id] = encode_loop_subd_tree_ptr(node);
}
#endif





