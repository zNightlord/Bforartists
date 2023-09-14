

#pragma BLENDER_REQUIRE(npr_strokegen_list_ranking_inputs_lib.glsl)

#ifndef NPR_STROKEGEN_LIST_RANKING_LIB_H
#define NPR_STROKEGEN_LIST_RANKING_LIB_H

/* Macro expansion */
#ifndef CAT
	#define CAT_(x, y) x ## y
	#define CAT(x, y) CAT_(x, y)
#endif

/* Inputs ///////////////////////////////////////////////////////////////////// */
/* TAG_LIST_RANKING */
#define tag_list_ranking TAG_LIST_RANKING
/* Inputs ///////////////////////////////////////////////////////////////////// */
/* See npr_strokegen_list_ranking_inputs_lib.glsl */


/* Implementations ////////////////////////////////////////////////////////////////// */

/** Deterministic Coin Tossing.
    https://algo2.iti.kit.edu/download/mem_hierarchy_11_unedited.pdf
    https://algo2.iti.kit.edu/download/mem_hierarchy_12_unedited.pdf
    https://algo2.iti.kit.edu/download/mem_hierarchy_06.pdf */
uint calc_node_tag(uint id_curr, uint id_next)
{
    uint xor = id_curr ^ id_next; 
    /** findLSB returns the bit number of the least significant bit 
    * that is set to 1 in the binary representation of value.  */
    int lsb_diff = findLSB(xor); 
    /* If value is zero, -1 will be returned. */
    if (lsb_diff == -1)
        return 0xffffffff; /* error: id_curr == id_next */ 
    uint lsb_diff_u = uint(max(0, lsb_diff)); 

    uint bitval_diff = ((id_curr >> lsb_diff_u) & 1u);

    uint tag = 2u * lsb_diff_u + bitval_diff; 

    /* tag_curr != tag_next is always true from this contruction */
    return tag; 
}
bool is_node_in_independent_set(uint tag_curr, uint tag_prev, uint tag_next, bool head_or_tail)
{ /* independent set node has the minimal ID among its neighbours 
    * anchors breaks the chains into small segments, 
    * each seg len < 2(log_2(N)), where N is #tagged nodes*/
    return ((tag_curr < tag_prev) && (tag_curr < tag_next)) && (!head_or_tail/*protect two ends*/);
}

bool is_head_node(uint prev_node_id, uint node_id)
{
    return (prev_node_id == node_id); 
}
#define FUNC_IS_NODE_HEAD is_head_node
bool is_tail_node(uint next_node_id, uint node_id)
{
    return (next_node_id == node_id); 
}
#define FUNC_IS_NODE_TAIL is_tail_node
bool is_head_or_tail_node(ListRankingLink link, uint node_id)
{
    return is_head_node(link.prev_node_id, node_id) || is_tail_node(link.next_node_id, node_id); 
}
#define FUNC_IS_NODE_HEAD_OR_TAIL is_head_or_tail_node



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING)

    /* Public Interfaces ///////////////////////////////// */
    #define _FUNC_UPDATE_NODE_TAG CAT(CalcNewTag_, tag_list_ranking)
    #define _FUNC_IS_NODE_INSIDE_INDEPENDENT_SET CAT(IsNodeAnchor_, tag_list_ranking)

	/**
	 * \brief Iteratively calculate & store tag for each node  
	 */
    void _FUNC_UPDATE_NODE_TAG(uint node_id, uint tagging_iter)
    {
        uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id); 

        uint tag, tag_next;
        if (tagging_iter == 0u)
        {
            tag      = node_id; 
            tag_next = next_node_id; 
        }
        else /*tagging_iter > 0u*/
        {
            tag      = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(node_id); 
            tag_next = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(next_node_id); 
        }

        if (FUNC_IS_NODE_TAIL(next_node_id, node_id))
            /* for tail this can be arbitrary unequal value
             * but we give larger tag to increase the possibility
             * to mark this into the independent set */
            tag_next = tag + 1; 

        uint tag_new = calc_node_tag(tag, tag_next);

        FUNC_DEVICE_STORE_LISTRANKING_NODE_TAG(node_id, tag_new); 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)
    /**
	 * \brief Determine if a node is anchor, head/tail node must be anchor for code to work
	*/
    bool _FUNC_IS_NODE_INSIDE_INDEPENDENT_SET(
        uint node_id, out ListRankingLink link
    ){
        link = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(node_id);
        bool b_list_head_or_tail = FUNC_IS_NODE_HEAD_OR_TAIL(link, node_id); 

        uint tag =      FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(node_id);
        uint tag_prev = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(link.prev_node_id);
        uint tag_next = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(link.next_node_id); 

        return is_node_in_independent_set(tag, tag_prev, tag_next, b_list_head_or_tail); 
    }

    /* In parallel split, each lane maps to a tgsm u32 counter.
     * Each node votes 1/0 based on which side they would be split into.
     * We use counter to collect votes from threads in the same lane but different waves.
     * But there can be invalid thread (does not map to a node) that votes 0 and corrupts the counter.
     * hence we need a mask to filter them out. 
     * 
     * This function returns an integer value to indicate the number of valid bits 
     * of the corresponding per-lane u32 counter. 
     */
    uint _FUNC_CALC_LDS_DIGIT_MASK(uint lane_id, uint blk_id, uint blk_size, uint num_nodes, out uint digit_mask)
    {
        const uint num_blks = (num_nodes + blk_size - 1) / blk_size; 
        const uint num_nodes_in_blk = 
            (blk_id == num_blks - 1u) ? num_nodes % blk_size : blk_size;

        int max_wave_id_valid = /* max wave_id that has lane_id mapped to a valid node */
            lane_id >= num_nodes_in_blk ? -1 : int(((num_nodes_in_blk - 1u) - lane_id) >> 5u);
        uint num_valid_waves = uint(max_wave_id_valid + 1); 
        digit_mask = num_valid_waves == 32u ? 0xffffffff : ~(0xffffffff << num_valid_waves); 

        return num_valid_waves; 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES)

/* update linkage + store rank(no load) : 0.95ms */
/* load rank + store rank : 0.53ms */
    void _FUNC_SPLICE_NODE_OUT(
        uint splice_id, uint node_id, uint node_rank, uint prev_node_id, uint next_node_id)
    {
        /* update linkage */
        if (prev_node_id != node_id)
            FUNC_DEVICE_STORE_LISTRANKING_NODE_NEXT_NODE_ID(prev_node_id, next_node_id); 
        if (next_node_id != node_id)
            FUNC_DEVICE_STORE_LISTRANKING_NODE_PREV_NODE_ID(next_node_id, prev_node_id); 
        /* update rank */
        uint prev_node_rank = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(prev_node_id); 
        prev_node_rank += node_rank; 

        FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK(prev_node_id, prev_node_rank); 
    }
    
#endif




#endif
