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
            num_work_items = NUM_ITEMS_BNPR_LIST_RANK_TEST;
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



#endif

