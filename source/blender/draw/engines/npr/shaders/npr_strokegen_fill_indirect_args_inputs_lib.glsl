#ifndef NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H
#define NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H



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
/* In: 
 * int pc_edge_split_dispatch_group_size_
 * int pcs_split_iter_
 * ssbo_edge_split_counters_[]
 * ssbo_indirect_dispatch_args_per_split_edge_
*/
void GetDispatchArgs(out uvec3 dispatch_args)
{
    uint num_work_items = 0; 
#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_EDGE)
    num_work_items = ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges_pass_1;
#endif
    dispatch_args.x = compute_num_groups(num_work_items, pc_edge_split_dispatch_group_size_);
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
    {
        num_work_items = ssbo_bnpr_mesh_pool_counters_.num_contour_edges - ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges; 
        ssbo_bnpr_mesh_pool_counters_.num_contour_edges_curr = num_work_items; /* we record newly genderated contour edge count here */
    }
    
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







#endif
