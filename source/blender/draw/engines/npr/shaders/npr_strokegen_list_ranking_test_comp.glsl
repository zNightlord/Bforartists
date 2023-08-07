
#pragma BLENDER_REQUIRE(npr_strokegen_scan_no_subgroup_codegen_lib.glsl)

#pragma BLENDER_REQUIRE(npr_strokegen_list_ranking_lib.glsl)



/* -----------------------------------------------
 uint pc_listranking_tagging_iter_; 
 ubo_list_ranking_tagging_; 
* ----------------------------------------------- */

#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_UPDATE_ANCHORS
void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;
    
    uint node_id = idx; 
    uint tagging_iter = pc_listranking_tagging_iter_;
    uint tag_subbuff_size = ubo_list_ranking_tagging_.subbuff_size;
    uint tag_new = 0u; 
    _FUNC_UPDATE_NODE_TAG(node_id, tagging_iter, tag_subbuff_size, tag_new);

    if (idx == 0u && pc_listranking_tagging_iter_ == 0u) /* clear the atomics */ 
    {
        ssbo_list_ranking_atomics_.counter_anchors = 0; 
        ssbo_list_ranking_atomics_.counter_lists = 0; 
    }
}
#endif




#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS
shared uint LDS_scan_block_offset; 

#define _TEST_ATOMIC_COMPACTION 1u 
#ifdef _TEST_ATOMIC_COMPACTION
 shared uint LDS_hist_per_lane_slot[32u]; 
 shared uint LDS_hist_blk; 
 shared uint LDS_offset_per_lane_slot[32u]; 
#endif

void main()
{
    uint tagging_iter = pc_listranking_tagging_iter_;
    uint tag_subbuff_size = ubo_list_ranking_tagging_.subbuff_size;

    const uint groupId = gl_LocalInvocationID.x;

    TreeScanIndices scan_ids = GetTreeScanIndices(groupId, gl_WorkGroupID.x);

    bool valid_item_A = scan_ids.global_x2.x < ubo_list_ranking_tagging_.num_nodes; 
    bool valid_item_B = scan_ids.global_x2.y < ubo_list_ranking_tagging_.num_nodes; 

    ListRankingLink links_A; 
    ListRankingLink links_B; 
    bool is_anchor_A = _FUNC_IS_NODE_ANCHOR(scan_ids.global_x2.x, tagging_iter, tag_subbuff_size, links_A);
    bool is_anchor_B = _FUNC_IS_NODE_ANCHOR(scan_ids.global_x2.y, tagging_iter, tag_subbuff_size, links_B); 
    bool is_head_or_tail_A = FUNC_IS_NODE_HEAD_OR_TAIL(links_A, scan_ids.global_x2.x); 
    bool is_head_or_tail_B = FUNC_IS_NODE_HEAD_OR_TAIL(links_B, scan_ids.global_x2.y); 
    bvec2 is_tail_AB = bvec2(
        FUNC_IS_NODE_TAIL(links_A, scan_ids.global_x2.x), 
        FUNC_IS_NODE_TAIL(links_B, scan_ids.global_x2.y) 
    ); 
    is_anchor_A = is_anchor_A && valid_item_A; 
    is_anchor_B = is_anchor_B && valid_item_B;


#ifdef _TEST_ATOMIC_COMPACTION
    const uint wave_id = groupId >> 5u; 
    const uint lane_id = groupId % 32u;
    if (wave_id == 0u) 
    { /* Clear LDS counters */
        if (lane_id == 0u)
            LDS_hist_blk = 0u; 
        LDS_hist_per_lane_slot[lane_id] = 0u; 
    }
    LDS_hist_per_lane_slot[lane_id] = 0u;
    barrier(); 

    /* Alloc for all 32 set of lanes */ 
    uint lds_compact_input = uint(is_anchor_A) + uint(is_anchor_B); 
    uint lds_compact_output = atomicAdd(LDS_hist_per_lane_slot[lane_id], lds_compact_input); 
    barrier(); 

    /* Prefix sum on lane sums */
    if (wave_id == 0u)
    {
        LDS_offset_per_lane_slot[lane_id] = 
            atomicAdd(LDS_hist_blk, LDS_hist_per_lane_slot[lane_id]);  
    }
    barrier(); 

    /* Compute final offset */
    uint local_offset = LDS_offset_per_lane_slot[lane_id] + lds_compact_output; 

    /* Add block sum to global counter. */
    if (groupId == gl_WorkGroupSize.x - 1u)
    {
        LDS_scan_block_offset = atomicAdd(
            ssbo_list_ranking_atomics_.counter_anchors, 
            LDS_hist_blk
        ); 
    }
    barrier(); 
#endif

#ifndef _TEST_ATOMIC_COMPACTION
    uint scanval_A, scanval_B;
    scanval_A = is_anchor_A ? 1u : 0u; 
    scanval_B = is_anchor_B ? 1u : 0u; 
    /* avoid invalid loads */
    _FUNC_CLEAN_SCAN_DATA(
        scan_ids, ubo_list_ranking_scan_infos_.num_scan_items,
        scanval_A, scanval_B /* <- inout */
    );

    /* execute block-wise exlusive scan */
	uint scanres_A, scanres_B;
    scanres_A = scanres_B = 0u;  
    _FUNC_TREE_SCAN_BLOCK(
		groupId,
		gl_WorkGroupID.x,
		scanval_A,
		scanval_B,
		/* -out- */
		scanres_A,
		scanres_B
	);

    /* Add block sum to global counter. */
    if (groupId == gl_WorkGroupSize.x - 1u)
    {
        LDS_scan_block_offset = atomicAdd(
            ssbo_list_ranking_atomics_.counter_anchors, 
            scanres_B + scanval_B
        ); 
    }
    barrier(); 

    uint blk_offset = LDS_scan_block_offset; 

    scanres_A += blk_offset; 
    scanres_B += blk_offset; 
#endif

#ifdef _TEST_ATOMIC_COMPACTION
    uint scanres_A = local_offset;
    uint scanres_B = local_offset + (is_anchor_A ? 1u : 0u); 
    
    uint blk_offset = LDS_scan_block_offset; 
    scanres_A += blk_offset; 
    scanres_B += blk_offset; 
#endif

    uvec2 anchor_id_AB = uvec2(scanres_A, scanres_B); 
    uvec2 node_id_AB = scan_ids.global_x2.xy; 
    if (is_anchor_A)
    { /* addressing already secured by is_anchor_A/B */
        FUNC_DEVICE_STORE_PER_ANCHOR_NODEID(anchor_id_AB.x, node_id_AB.x); 
    }
    if (valid_item_A)
    {
        FUNC_DEVICE_STORE_PER_NODE_CONTRACTION_INFO(
            node_id_AB.x, links_A.next_node_id, is_anchor_A, anchor_id_AB.x
        ); 
    }

    if (is_anchor_B)
    {
        FUNC_DEVICE_STORE_PER_ANCHOR_NODEID(anchor_id_AB.y, node_id_AB.y); 
    }
    if (valid_item_B)
    {
        FUNC_DEVICE_STORE_PER_NODE_CONTRACTION_INFO(
            node_id_AB.y, links_B.next_node_id, is_anchor_B, anchor_id_AB.y
        ); 
    } 

}
#endif



#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_RANKING
void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    uint anchor_id = idx; 
    uint num_anchors = FUNC_GET_NUM_ANCHORS(); 
    bool b_valid_anchor = FUNC_DEVICE_VALIDATE_THREAD_PER_ANCHOR(anchor_id, num_anchors); 
    if (!b_valid_anchor) return; 

    uint anchor_node_id = FUNC_DEVICE_LOAD_PER_ANCHOR_NODEID(anchor_id); 


    /* Traverse the sublist, starting from the succssor node.
    *  we ensure that #nodes<=BNPR_LIST_RANKING_MAX_SUBLIST_LEN between two anchors (See ComputeTaggingIters) 
    */
    ListRankingLink links = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(anchor_node_id);
    uint curr_node_id = links.next_node_id; 
    uint curr_node_rank_in_sublist = 1u; 
    uint sublist_len = 1u; 
    ContractionInfo ci;
    bool reached_end = false; 
    for (uint i = 0; i < BNPR_LIST_RANKING_MAX_SUBLIST_LEN; i++)
    {
        FUNC_DEVICE_LOAD_PER_NODE_CONTRACTION_INFO(curr_node_id, ci);
        reached_end = reached_end || ci.is_list_end;  
        if (ci.is_anchor || reached_end) /* hit the next sublist's header, exit */ 
            break; 

        /* Link non-anchor node to its sub-list head(i.e this anchor) */
        FUNC_DEVICE_STORE_PER_NODE_TO_ANCHOR_ID(curr_node_id, anchor_id); 
        /* Store sub-list rank for non-anchor nodes */
        FUNC_DEVICE_STORE_PER_NODE_SUBLIST_RANK(curr_node_id, curr_node_rank_in_sublist); 
        

        /* Move to next non-anchor node */
        curr_node_id = ci.id; 
        curr_node_rank_in_sublist++; 
        sublist_len++; 
    }
    /* Store sub-list rank for this anchor as sub-list length */
    FUNC_DEVICE_STORE_PER_NODE_SUBLIST_RANK(anchor_node_id, sublist_len); 
    

    FUNC_DEVICE_STORE_PER_ANCHOR_NEXT_ANCHORID(anchor_id, ci.id); 
}


#endif







#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    uint anchor_id = idx; 
    uint num_anchors = FUNC_GET_NUM_ANCHORS(); 
    bool b_valid_anchor = FUNC_DEVICE_VALIDATE_THREAD_PER_ANCHOR(anchor_id, num_anchors); 
    if (!b_valid_anchor) return; /* invalid thread, do nothing */

    uint subbuff_size = (((num_anchors + 3u) >> 2u) << 2u); /* ping-pong */

    uint curr_iter = pc_listranking_jumping_iter_; 
    uint num_iters = FUNC_GET_NUM_JUMPS();

    JumpingInfo ji; 
    if (0u == curr_iter)
    {
        uint anchor_node_id = FUNC_DEVICE_LOAD_PER_ANCHOR_NODEID(anchor_id); 
        ji = FUNC_DEVICE_INIT_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, anchor_node_id);  

        FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, curr_iter, ji, subbuff_size);
         
        return; /* 0th iter to init jumping info */
    }
    if (num_iters < curr_iter) return; /* more than needed, do nothing */


    /* Pointer-Jumping */
    JumpingInfo ji_next, ji_updated;
    ji      = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id,              curr_iter, subbuff_size); 
    ji_next = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(ji.jump_next_anchor_id, curr_iter, subbuff_size);
    
    bool jumped_to_end = false; 
    ji_updated = FUNC_DEVICE_UPDATE_ANCHOR_LIST_JUMPING_INFO(ji, ji_next, jumped_to_end); 


    FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, curr_iter, ji_updated, subbuff_size);


    /* In the last iter, anchors which is/are list header(s) put 
     * neccessary data into the tail, so that every anchors in the same list
     * can access it using ji.jump_next_anchor_id
     */
    if ((num_iters == pc_listranking_jumping_iter_) && ji_updated.is_list_head)
    {
        JumpingListInfo list_info; 
        list_info.list_len = ji_updated.data + 1; /* +1 for tail node */  
        list_info.list_addr = atomicAdd(ssbo_list_ranking_atomics_.counter_lists, list_info.list_len); 
        FUNC_DEVICE_STORE_TAIL_ANCHOR_LIST_INFO(ji_updated.jump_next_anchor_id, list_info, subbuff_size); 
    }
}

#endif






#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_STITCHING

#endif