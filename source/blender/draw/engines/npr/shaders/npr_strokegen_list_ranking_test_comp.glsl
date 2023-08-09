
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
        ssbo_list_ranking_atomics_.counter_spliced = 0;
        ssbo_list_ranking_atomics_.counter_lists = 0; 
    }
}
#endif




#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS

shared uint LDS_digit_per_lane[32u]; 
shared uint LDS_offset_per_lane_slot[32u][2u]; 

shared uint LDS_hist_blk[2];  
shared uint LDS_scan_block_offset[2]; 

void main()
{
    const uint tagging_iter = pc_listranking_tagging_iter_;
    const uint tag_subbuff_size = ubo_list_ranking_tagging_.subbuff_size;
    const uint num_nodes = ubo_list_ranking_tagging_.num_nodes; 

    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    const uint node_id = idx; 
    bool valid_item = node_id < num_nodes; 

    ListRankingLink links; 
    bool is_head_or_tail = FUNC_IS_NODE_HEAD_OR_TAIL(links, node_id); 
    bool is_tail         = FUNC_IS_NODE_TAIL(links, node_id); 
    bool is_anchor = _FUNC_IS_NODE_ANCHOR(node_id, tagging_iter, tag_subbuff_size, links);
    is_anchor = is_anchor && valid_item; 


    const uint wave_id = groupId >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */
    const uint lane_id = groupId % 32u;
    const uint num_waves = gl_WorkGroupSize.x >> 5u; 
    if (wave_id == 0u) 
    { /* Clear LDS counters */
        if (lane_id == 0u)
            LDS_hist_blk = 0u; 
        LDS_digit_per_lane[lane_id] = 0u; 
    }
    barrier(); 

    /* Mark 1/0 at bit #wave_id */ 
    uint compact_bitval = uint(is_anchor); 
    uint lds_compact_input = compact_bitval << wave_id; 
    atomicOr(LDS_digit_per_lane[lane_id], lds_compact_input); 
    barrier(); 

    /* Prefix sum on lane sums */
    uint lane_digit = LDS_digit_per_lane[lane_id]; 
#define DIGIT_BITS_BEHIND wave_id
    uint wave_mask = (0xffffffffu >> (32u - DIGIT_BITS_BEHIND));
    uint lane_digit_masked = lane_digit & wave_mask; 
    uint num_1_bits_low = bitCount(lane_digit_masked);
    uint num_0_bits_low = DIGIT_BITS_BEHIND - lane_offset_1; 
#undef DIGIT_BITS_BEHIND
    uint lane_offset = is_anchor ? num_1_bits_low : num_0_bits_low; 

    if (wave_id <= 1u)
    {
        uint num_x_bits = bitCount(lane_digit); 
        if (wave_id == 0u) num_x_bits = num_waves/*#valid bits*/ - lane_digit; 
        LDS_offset_per_lane_slot[lane_id][wave_id] = atomicAdd(LDS_hist_blk[wave_id], num_x_bits);  
    }
    barrier(); 

    /* Add block sum to global counter. */
    if (groupId == gl_WorkGroupSize.x - 2u)
    {
        LDS_scan_block_offset[0] = atomicAdd(
            ssbo_list_ranking_atomics_.counter_spliced, 
            LDS_hist_blk[0]
        ); 
    }
    if (groupId == gl_WorkGroupID.x - 1u)
    {
        LDS_scan_block_offset[1] = atomicAdd(
            ssbo_list_ranking_atomics_.counter_anchors, 
            LDS_hist_blk[1]
        ); 
    }
    barrier(); 

    /* Compute final offset */
    uint local_offset = LDS_offset_per_lane_slot[lane_id][compact_bitval] + lane_offset; 
    uint blk_offset   = LDS_scan_block_offset[compact_bitval]; 
    
    uint scanres = local_offset + blk_offset; 


    uint output_id = scanres; /* anchor_id or spliced_node_id */
    if (is_anchor)
    { /* addressing already secured by is_anchor */
        FUNC_DEVICE_STORE_PER_ANCHOR_NODEID(output_id, node_id); 
    }else{
        /* FUNC_DEVICE_STORE_PER_SPLICED_NODEID(output_id, node_id); */
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