

#pragma BLENDER_REQUIRE(npr_strokegen_list_ranking_lib.glsl)

#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_FILL_CPU_DATA

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;
    
    uint node_id = idx; 

    ssbo_list_ranking_links_out_[node_id*2] = ssbo_list_ranking_links_in_[node_id*2];
    ssbo_list_ranking_links_out_[node_id*2+1] = ssbo_list_ranking_links_in_[node_id*2+1]; 
}

#endif

#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING
void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;
    
    const uint splicing_iter = pc_listranking_splice_iter_; 
    const uint tagging_iter = pc_listranking_tagging_iter_; 
    const uint num_nodes_total = NUM_NODES_TOTAL; 
    const uint num_anchors = splicing_iter == 0 ? 
        num_nodes_total : ssbo_list_ranking_anchor_counters_[splicing_iter]; 
    
    /* compute node tag */
    uint anchor_id = idx; 
    uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter); 
    if (anchor_id < num_anchors && node_id < num_nodes_total)
    {
        _FUNC_UPDATE_NODE_TAG(node_id, tagging_iter); 
    }

    /* bootstraping */
    if (splicing_iter == 0u && tagging_iter == 0u)
    { 
        FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK(node_id, 1u); 
        if (idx == 0u)
        {
            ssbo_list_ranking_anchor_counters_[0] = num_nodes_total;
            ssbo_list_ranking_splice_counters_[0] = 0; 
        }
    }

    /* clear the atomics */ 
    if (idx == 0u && tagging_iter == 0u) 
    { 
        ssbo_list_ranking_anchor_counters_[splicing_iter + 1] = 0;
        ssbo_list_ranking_splice_counters_[splicing_iter + 1] = 0; 
    }
}
#endif




#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS

shared uint LDS_digit_per_lane[32u]; 
shared uint LDS_offset_per_lane_slot[32u][2u]; 

shared uint LDS_hist_blk[2];  
shared uint LDS_scan_block_offset[2]; 
/* 6 loads + 2 atomic_adds/tg + 1 store */
void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    const uint tag_subbuff_size = SUB_BUFF_SIZE;
    const uint num_nodes_total = NUM_NODES_TOTAL; 

    const uint splicing_iter = pc_listranking_splice_iter_; 
    const uint num_anchors = ssbo_list_ranking_anchor_counters_[splicing_iter]; 

    const uint anchor_id = idx; 
    const uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter); 
    bool valid_item = (node_id < num_nodes_total) && (anchor_id < num_anchors); 

    const uint num_full_groups = (num_anchors / blockSize) * blockSize; 
    const uint num_valid_threads_in_group = 
        blockIdx > num_full_groups ? num_anchors % blockSize : blockSize;

    ListRankingLink links; 
    bool is_anchor = (false == _FUNC_IS_NODE_INSIDE_INDEPENDENT_SET(node_id, /*out*/links));
    bool is_head_or_tail = FUNC_IS_NODE_HEAD_OR_TAIL(links, node_id); 
    bool is_tail         = FUNC_IS_NODE_TAIL(links.next_node_id, node_id); 
    is_anchor = is_anchor && valid_item; 

    const uint wave_id = groupIdx >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */
    const uint lane_id = groupIdx % 32u;
    const uint num_waves = gl_WorkGroupSize.x >> 5u; 
    uint valid_digit_mask = 0xffffffff; 
    const uint num_valid_bits = 
        _FUNC_CALC_LDS_DIGIT_MASK(lane_id, blockIdx, blockSize, num_anchors, /*out*/valid_digit_mask);
    
    if (wave_id == 0u) 
    { /* Clear LDS counters */
        if (lane_id <= 1u)
            LDS_hist_blk[lane_id] = 0u; 
        LDS_digit_per_lane[lane_id] = 0u; 
    }
    barrier(); 

    /* Mark 1/0 at bit #wave_id */ 
    uint compact_bitval = uint(is_anchor); 
    uint lds_compact_input = compact_bitval << wave_id; 
    atomicOr(LDS_digit_per_lane[lane_id], lds_compact_input); 
    barrier(); 

    /* Prefix sum on lane sums */
    uint lane_digit = (LDS_digit_per_lane[lane_id] & valid_digit_mask); 
    uint wave_mask = (~(0xffffffffu << wave_id));
    uint lane_digit_masked = lane_digit & wave_mask; 
    uint num_1_bits_low = bitCount(lane_digit_masked);
    uint num_0_bits_low = bitCount((~lane_digit_masked) & (valid_digit_mask & wave_mask)); 
    uint lane_offset = is_anchor ? num_1_bits_low : num_0_bits_low; 

    if (wave_id <= 1u)
    {
        uint num_x_bits =               bitCount(lane_digit); /*#1-bits*/ 
        if (wave_id == 0u) num_x_bits = bitCount((~lane_digit) & (valid_digit_mask)); /*#valid 0-bits*/
        LDS_offset_per_lane_slot[lane_id][wave_id] = atomicAdd(LDS_hist_blk[wave_id], num_x_bits);  
    }
    barrier(); 

    /* Add block sum to global counter. */
    if (groupIdx == gl_WorkGroupSize.x - 2u)
    {
        LDS_scan_block_offset[0] = atomicAdd(
            ssbo_list_ranking_splice_counters_[splicing_iter + 1],  
            LDS_hist_blk[0]
        ); 
    }
    if (groupIdx == gl_WorkGroupSize.x - 1u)
    {
        LDS_scan_block_offset[1] = atomicAdd(
            ssbo_list_ranking_anchor_counters_[splicing_iter + 1], 
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
        FUNC_DEVICE_STORE_PER_SPLICED_NODEID(output_id, node_id); 
    }
}
#endif






#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES
/* 4 loads + 3 stores */
void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    const uint splicing_iter = pc_listranking_splice_iter_; 
    const uint num_nodes_total = NUM_NODES_TOTAL; 
    const uint num_spliced = ssbo_list_ranking_splice_counters_[splicing_iter]; 

    const uint splice_id = idx; 
    const uint node_id = FUNC_GET_NODE_ID_FOR_SPLICED(splice_id); 
    bool valid_item = (node_id < num_nodes_total) && (splice_id < num_spliced); 

    ListRankingLink links = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(node_id); 

    uint node_rank = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(node_id); 
    _FUNC_SPLICE_NODE_OUT(splice_id, node_id, node_rank, links.prev_node_id, links.next_node_id); 
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
        uint anchor_node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter); 
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
}

#endif






#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_STITCHING

#endif