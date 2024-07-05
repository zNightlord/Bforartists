#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)


#ifndef NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H
#define NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H



#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__SCAN)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    dispatch_args.x = ssbo_tree_scan_infos_.num_thread_groups;
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_scan_dispatch_args_.num_groups_x = args.x;
    ssbo_scan_dispatch_args_.num_groups_y = args.y;
    ssbo_scan_dispatch_args_.num_groups_z = args.z;
    ssbo_scan_dispatch_args_._pad0 = 0; 
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__SEGLOOPCONV1D)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = ssbo_segloopconv1d_info_.num_conv_items; 
    dispatch_args.x = compute_num_groups(num_work_items, pc_segloopconv1d_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_segloopconv1d_dispatch_args_.num_groups_x = args.x;
    ssbo_segloopconv1d_dispatch_args_.num_groups_y = args.y;
    ssbo_segloopconv1d_dispatch_args_.num_groups_z = args.z;
    ssbo_segloopconv1d_dispatch_args_._pad0 = 0; 
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__LIST_RANKING_ANCHORS)

#define DISTPATCH_GRANULARITY__PER_ANCHOR 0u
#define DISTPATCH_GRANULARITY__PER_SPLICED 1u

void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = 0;
    if (pc_listranking_indirect_arg_granularity_ == DISTPATCH_GRANULARITY__PER_ANCHOR)
    {
        if (pc_listranking_counter_buffer_slot_id_ == 0)
        { /* in practice this could be from other GPU buffer */
            num_work_items = (pc_listranking_custom_ == 0) ? 
                NUM_ITEMS_BNPR_LIST_RANK_TEST : ssbo_list_ranking_inputs_.num_nodes; 
        }else{
            num_work_items = ssbo_list_ranking_anchor_counters_[pc_listranking_counter_buffer_slot_id_];
        }
    }
    if (pc_listranking_indirect_arg_granularity_ == DISTPATCH_GRANULARITY__PER_SPLICED)
    {
        num_work_items = ssbo_list_ranking_splice_counters_[pc_listranking_counter_buffer_slot_id_];
    }

    dispatch_args.x = compute_num_groups(num_work_items, pc_listranking_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    if (pc_listranking_indirect_arg_granularity_ == DISTPATCH_GRANULARITY__PER_ANCHOR)
    {
        ssbo_list_ranking_indirect_dispatch_args_per_anchor.num_groups_x = args.x;
        ssbo_list_ranking_indirect_dispatch_args_per_anchor.num_groups_y = args.y;
        ssbo_list_ranking_indirect_dispatch_args_per_anchor.num_groups_z = args.z;
        ssbo_list_ranking_indirect_dispatch_args_per_anchor._pad0 = 0; 
    }
    if (pc_listranking_indirect_arg_granularity_ == DISTPATCH_GRANULARITY__PER_SPLICED)
    {
        ssbo_list_ranking_indirect_dispatch_args_per_spliced.num_groups_x = args.x;
        ssbo_list_ranking_indirect_dispatch_args_per_spliced.num_groups_y = args.y;
        ssbo_list_ranking_indirect_dispatch_args_per_spliced.num_groups_z = args.z;
        ssbo_list_ranking_indirect_dispatch_args_per_spliced._pad0 = 0;
    }
}

#endif

#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING)
/* In: 
 * int pc_meshing_dispatch_group_size_
 * ssbo_bnpr_mesh_pool_counters_
 * ssbo_indirect_dispatch_args_per_filtered_edge_
 * ssbo_indirect_dispatch_args_per_filtered_vert_
*/
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = 0; 
    
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_EDGE)
    num_work_items = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_VERT)
    num_work_items = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts;
#endif
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_meshing_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_EDGE)
    ssbo_indirect_dispatch_args_per_filtered_edge_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_filtered_edge_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_filtered_edge_.num_groups_z = args.z; 
    ssbo_indirect_dispatch_args_per_filtered_edge_._pad0 = 0;
#endif
 
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_VERT)
    ssbo_indirect_dispatch_args_per_filtered_vert_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_filtered_vert_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_filtered_vert_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_filtered_vert_._pad0 = 0;
#endif
}

#endif



#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING)
/* new params: ssbo_bnpr_mesh_pool_counters_.num_filtered_edges, pcs_only_selected_elems_ */
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = 0; 
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_EDGE)
    num_work_items = ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges_pass_1;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_edge_split_dispatch_group_size_);
#endif
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_COLLAPSED_EDGE)
    num_work_items = ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapsed_edges_pass_1;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_edge_collapse_dispatch_group_size_);
#endif
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_FLIP_EDGE)
    num_work_items = ssbo_edge_flip_counters_[pcs_flip_iter_].num_flip_edges_pass_1;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_edge_flip_dispatch_group_size_);
#endif
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_FACE)
    num_work_items = ssbo_face_split_counters_[pcs_split_iter_].num_split_faces;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_face_split_dispatch_group_size_);
#endif
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_EDGE)
    uint num_edges_dynamesh = ssbo_dyn_mesh_counters_.num_edges; 
    uint num_edges_static = (0 < pcs_only_selected_elems_) ? ssbo_bnpr_mesh_pool_counters_.num_filtered_edges : pcs_edge_count_; 
    num_work_items = num_edges_static + num_edges_dynamesh;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_remeshed_edges_dispatch_group_size_);
#endif
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_VERT)
    uint num_verts_dynamesh = ssbo_dyn_mesh_counters_.num_verts; 
    uint num_verts_static = (0 < pcs_only_selected_elems_) ? ssbo_bnpr_mesh_pool_counters_.num_filtered_verts : pcs_vert_count_; 
    num_work_items = num_verts_static + num_verts_dynamesh;
    dispatch_args.x = compute_num_groups(num_work_items, pcs_remeshed_verts_dispatch_group_size_);
#endif
    dispatch_args.y = 1;
    dispatch_args.z = 1; 
}

void FillDispatchArgsBuffer(uvec3 args)
{
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_EDGE)
    ssbo_indirect_dispatch_args_per_split_edge_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_split_edge_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_split_edge_.num_groups_z = args.z; 
    ssbo_indirect_dispatch_args_per_split_edge_._pad0 = 0;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_COLLAPSED_EDGE)
    ssbo_indirect_dispatch_args_per_collapsed_edge_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_collapsed_edge_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_collapsed_edge_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_collapsed_edge_._pad0 = 0;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_FLIP_EDGE)
    ssbo_indirect_dispatch_args_per_flip_edge_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_flip_edge_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_flip_edge_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_flip_edge_._pad0 = 0;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_FACE)
    ssbo_indirect_dispatch_args_per_split_face_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_split_face_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_split_face_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_split_face_._pad0 = 0;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_EDGE)
    ssbo_indirect_dispatch_args_per_remeshed_edges_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_remeshed_edges_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_remeshed_edges_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_remeshed_edges_._pad0 = 0;
#endif

#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_VERT)
    ssbo_indirect_dispatch_args_per_remeshed_verts_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_remeshed_verts_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_remeshed_verts_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_remeshed_verts_._pad0 = 0;
#endif
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FILL_DRAW_ARGS__REMESHING)
void GetDrawArgs(out DrawCommand draw_args)
{
#if defined(_KERNEL_MULTICOMPILE__FILL_DRAW_ARGS__REMESHING__DBG_LINES)
    draw_args.vertex_len = 2u * get_debug_line_counter(pcs_line_type_);
#endif
    draw_args.instance_len 	= 1;  		/*#instances*/
    draw_args.vertex_first 	= 0;  		/*ibo offset*/
    draw_args.base_index 	= 0;  		/*vbo offset*/
    draw_args.instance_first_indexed = 0; 	/*instance offset*/
    draw_args._pad0 = 0;
    draw_args._pad1 = 0;
    draw_args._pad2 = 0; 
}
void FillDrawArgsBuffer(DrawCommand draw_args)
{
    ssbo_bnpr_vert_debug_draw_args_.vertex_len             = draw_args.vertex_len; 
    ssbo_bnpr_vert_debug_draw_args_.instance_len           = draw_args.instance_len; 
    ssbo_bnpr_vert_debug_draw_args_.vertex_first           = draw_args.vertex_first; 
    ssbo_bnpr_vert_debug_draw_args_.base_index             = draw_args.base_index; 
    ssbo_bnpr_vert_debug_draw_args_.instance_first_indexed = draw_args.instance_first_indexed; 
    ssbo_bnpr_vert_debug_draw_args_._pad0                  = draw_args._pad0; 
    ssbo_bnpr_vert_debug_draw_args_._pad1                  = draw_args._pad1; 
    ssbo_bnpr_vert_debug_draw_args_._pad2                  = draw_args._pad2; 
}
#endif



#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_EDGE)
/* In: 
 * pc_per_contour_edge_dispatch_group_size_
 * ssbo_bnpr_mesh_pool_counters_
 * ssbo_indirect_dispatch_args_per_contour_edge_
*/
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = 0; 
    if (pc_dispatch_for_all_edges_ > 0)
        num_work_items = ssbo_bnpr_mesh_pool_counters_.num_contour_edges;
    else
        num_work_items = ssbo_bnpr_mesh_pool_counters_.num_contour_edges - ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges; 
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_per_contour_edge_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_indirect_dispatch_args_per_contour_edge_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_contour_edge_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_contour_edge_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_contour_edge_._pad0 = 0;
}
#endif


#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_FRAG)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items;
    if (0 < pc_dispatch_for_all_frags_)
        num_work_items = ssbo_bnpr_mesh_pool_counters_.num_frags; 
    else
        num_work_items = ssbo_bnpr_mesh_pool_counters_.num_frags - ssbo_bnpr_mesh_pool_counters_prev_.num_frags;
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_per_contour_frag_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_bnpr_mesh_contour_frag_dispatch_args_.num_groups_x = args.x;
    ssbo_bnpr_mesh_contour_frag_dispatch_args_.num_groups_y = args.y;
    ssbo_bnpr_mesh_contour_frag_dispatch_args_.num_groups_z = args.z;
    ssbo_bnpr_mesh_contour_frag_dispatch_args_._pad0 = 0;
}
#endif


#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_VERT)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = ssbo_bnpr_mesh_pool_counters_.num_contour_verts;
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_per_contour_vert_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_indirect_dispatch_args_per_contour_vert_.num_groups_x = args.x;
    ssbo_indirect_dispatch_args_per_contour_vert_.num_groups_y = args.y;
    ssbo_indirect_dispatch_args_per_contour_vert_.num_groups_z = args.z;
    ssbo_indirect_dispatch_args_per_contour_vert_._pad0 = 0;
}
#endif


#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_2D_SAMPLE)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items;
    num_work_items = ssbo_bnpr_mesh_pool_counters_.num_2d_samples;
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_per_contour_2d_sample_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_.num_groups_x = args.x;
    ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_.num_groups_y = args.y;
    ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_.num_groups_z = args.z;
    ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_._pad0 = 0;
}
#endif


#if defined (_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_TEMPORAL_RECORD)
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items;
    num_work_items = temporal_record_counter(pc_obj_id_, pc_frame_id_);
    
    dispatch_args.x = compute_num_groups(num_work_items, pc_temporal_records_dispatch_group_size_);
    dispatch_args.y = 1;
    dispatch_args.z = 1;
}

void FillDispatchArgsBuffer(uvec3 args)
{
    ssbo_bnpr_temporal_record_dispatch_args_.num_groups_x = args.x;
    ssbo_bnpr_temporal_record_dispatch_args_.num_groups_y = args.y;
    ssbo_bnpr_temporal_record_dispatch_args_.num_groups_z = args.z;
    ssbo_bnpr_temporal_record_dispatch_args_._pad0 = 0;
}
#endif


#endif
