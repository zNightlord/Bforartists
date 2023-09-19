

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
        uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id);
        bool is_tail_node = node_id == next_node_id; 
        FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK(node_id, is_tail_node ? 0u : 1u); 
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
        ssbo_list_ranking_addressing_counters_[0] = 0; 
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

    const uint num_nodes_total = NUM_NODES_TOTAL; 

    const uint splicing_iter = pc_listranking_splice_iter_; 
    const uint num_splicing_iters = pc_num_splice_iters_; 
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
        if (splicing_iter == num_splicing_iters - 1u)
        { /* save back pointers for later pointer jumping */
            FUNC_DEVICE_STORE_PER_NODE_ANCHORID(node_id, output_id); 
        }
    }else if (valid_item)
    {
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
    const uint num_spliced = ssbo_list_ranking_splice_counters_[splicing_iter + 1]; 

    const uint splice_id = idx; 
    const uint node_id = FUNC_GET_NODE_ID_FOR_SPLICED(splice_id); 
    bool valid_item = (node_id < num_nodes_total) && (splice_id < num_spliced); 

    if (!valid_item) return; 

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

    const uint curr_jump_iter = pc_listranking_jumping_iter_; 
    const uint splicing_iter = pc_listranking_splice_iter_; 

    const uint num_nodes_total = NUM_NODES_TOTAL; 
    const uint num_anchors = ssbo_list_ranking_anchor_counters_[splicing_iter]; 

    const uint num_jump_iters = FUNC_GET_NUM_JUMPS(num_anchors); 
    const uint iter_bootstrap   = 0; 
    const uint iter_jump_beg    = 1; 
    const uint iter_jump_end    = num_jump_iters; 
    const uint iter_finalize    = num_jump_iters + 1; 
    const uint num_iters = 2 + num_jump_iters; 
    /* +1 iter for initialization */
    /* +1 iter for finalization */

    const uint anchor_id = idx; 
    bool b_valid_anchor = (anchor_id < num_anchors); 
    if (!b_valid_anchor) return; /* invalid thread, do nothing */

    JumpingInfo ji; 

    /* Extra iters - initialization & finalization */
    if (iter_bootstrap == curr_jump_iter)
    { /* The first iteration sets up anchor-to-anchor links and per-anchor ranks
           from per-node buffers. This is for strided mem access. */
        ji = FUNC_DEVICE_INIT_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, splicing_iter);  
        FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, ji);
         
        return; /* !!! EXIT !!! ------------------------- */
    }
    if (iter_finalize == curr_jump_iter)
    { /* The last iteration sets up the topology for each anchor */
        uint curr_node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter); 
        uint curr_broadcast_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(curr_node_id); 
        /* retrieve broadcasted data */
        uint list_len, list_addr; 
        FUNC_DEVICE_RETRIEVE_LIST_TOPOLOGY(curr_broadcast_node_id, /*out*/list_len, list_addr);
        FUNC_DEVICE_STORE_LIST_TOPOLOGY(curr_node_id, list_len, list_addr); 

        return; /* !!! EXIT !!! -------------------------- */
    }
    if (num_iters <= curr_jump_iter) return; /* more than needed, do nothing */



    /* Pointer-Jumping */
    JumpingInfo ji_next, ji_updated;
    ji      = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id); 
    ji_next = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(ji.jump_next_anchor_id);
    
    bool jumped_to_end = false; 
    ji_updated = FUNC_DEVICE_UPDATE_ANCHOR_LIST_JUMPING_INFO(ji, ji_next, /*out*/jumped_to_end); 

    if (curr_jump_iter < iter_jump_end)
    { /* iteratively update jumping pointer & data */
        FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, ji_updated);
    } /* don't need to store at last jump iter (curr_jump_iter == iter_jump_end) */

    /* At last jumping iter, use per-anchor data to update per-node buffers */
    if (curr_jump_iter == iter_jump_end)
    { 
        uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter); 

        /* *) Output anchor rank ------------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
        if (pc_listranking_ranking_pass_with_broken_loops_ > 0)
        { /* cache head flag so that each node can simutaneously know the Rank&Head_Flag of its next node  */
            FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK_WITH_HEAD_TAIL_BITS(
                node_id, ji_updated.data, ji_updated.is_list_head, ji_updated.is_list_tail
            );
        }
#else 
        FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK(node_id, ji_updated.data); 
#endif

        /* *) Store broadcast position ------------------------------------------------ */ 
        /* for non-loop lists, this is the tail; for loops this is the last anchor. */
        /* in either way, we get a uniform position for all nodes in a list */
        uint list_broadcast_node_id = ji_updated.jump_next_anchor_id; 
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
        if (pc_listranking_ranking_pass_with_broken_loops_ > 0)
#endif
            FUNC_DEVICE_STORE_LISTRANKING_NODE_NEXT_NODE_ID(node_id, list_broadcast_node_id); 

        /* Also put list topology at the tail node, so that all list nodes can access it */ 
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
        if (pc_listranking_ranking_pass_with_broken_loops_ > 0)
#endif
        if (ji_updated.is_list_head)
        { /* head node is responsible for this */
            /* for each anchor, rank == dist to list tail */
            /* => for head, rank+1 == list length */
            uint list_len = 1 + ji_updated.data; 
            uint list_addr = FUNC_DEVICE_ALLOC_LIST_ADDR(list_len); 
            FUNC_DEVICE_BROADCAST_LIST_TOPOLOGY(list_broadcast_node_id, list_len, list_addr); 
        }
    }
}

#endif


#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    const uint splicing_iter = pc_listranking_splice_iter_; 

    const uint num_nodes_total = NUM_NODES_TOTAL; 
    const uint num_anchors = ssbo_list_ranking_anchor_counters_[splicing_iter]; 

    const uint anchor_id = idx; 
    bool b_valid_anchor = (anchor_id < num_anchors); 
    if (!b_valid_anchor) return; /* invalid thread, do nothing */


    JumpingInfo ji;
    ji = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id); 
    
    bool is_loop_head, is_loop_tail; 
    /* Non-loop head/tail flags are already there. */
    is_loop_head = is_loop_tail = false; 

    /* Find loop head */
    uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id); 
    const uint anchor_code = node_id; 
    if (ji.data == anchor_code && !ji.is_list_head/*exclude non-loop cases*/) 
        is_loop_head = true;

    /* Find loop tail */
    /* TODO: this can be rough, consider cache next anchor id for performance. */
    uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id); 
    uint next_anchor_id = FUNC_DEVICE_LOAD_PER_NODE_ANCHORID(next_node_id); 
    JumpingInfo ji_next = FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO(next_anchor_id); 
    
    const uint next_anchor_code = next_node_id; 
    if (ji_next.data == next_node_id && !ji.is_list_tail)
        is_loop_tail = true;     

    /* Update linkage to break the loop */
    if (is_loop_head)
    {
        FUNC_DEVICE_STORE_LISTRANKING_NODE_PREV_NODE_ID(node_id, node_id);             
        ji.is_list_head = true; 
    }
    if (is_loop_tail)
    {
        FUNC_DEVICE_STORE_LISTRANKING_NODE_NEXT_NODE_ID(node_id, node_id); 
        ji.is_list_tail = true; 
    }

    /* Also setup the jumping info for next jumping pass (which does ranking) */
    ji.data = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(node_id); 
    ji.next_anchor_id = is_loop_tail ? anchor_id : next_anchor_id; 
    FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO(anchor_id, ji);  
}

#endif



#ifdef _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx      = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    const uint blockSize = gl_WorkGroupSize.x;

    uint num_relink_iters = pc_listranking_num_relink_iters_; 
    uint curr_relink_iter = pc_listranking_relink_iter_; 

    const uint num_nodes_total = NUM_NODES_TOTAL; 
    const uint num_spliced_nodes = ssbo_list_ranking_splice_counters_[num_relink_iters - curr_relink_iter]; 

    const uint spliced_id = idx; 
    if (spliced_id >= num_spliced_nodes) return; /* invalid thread, do nothing */ 



    const uint node_id = FUNC_GET_NODE_ID_FOR_SPLICED(spliced_id); 
    uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id);



    /* Update Rank ------------------------------------------------------------------ */
    uint node_rank = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(node_id); 
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
    bool is_next_node_head, is_next_node_tail; 
    uint next_node_rank = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK_WITH_HEAD_TAIL_BITS(
        next_node_id, /*out*/is_next_node_head, is_next_node_tail
    );
#else
    uint next_node_rank = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(next_node_id); 
#endif

#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
    if (false == is_next_node_head)
        node_rank += next_node_rank; 
#else
    node_rank += next_node_rank; 
#endif

    FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK(node_id, node_rank);   




    /* diffuse broadcast-node-id from top spliced layer to the bottom ----------------- */
    uint broadcast_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(next_node_id); 

    uint list_len, list_addr; /* retrieve broadcasted data */
    FUNC_DEVICE_RETRIEVE_LIST_TOPOLOGY(broadcast_node_id, /*out*/list_len, list_addr); 
    FUNC_DEVICE_STORE_LIST_TOPOLOGY(node_id, list_len, list_addr); 

    /* Since each layer of spliced nodes form a independent set
     * read/write simutaneuosly to this buffer is safe here. */ 
    FUNC_DEVICE_STORE_LISTRANKING_NODE_NEXT_NODE_ID(node_id, broadcast_node_id); 

}

#endif