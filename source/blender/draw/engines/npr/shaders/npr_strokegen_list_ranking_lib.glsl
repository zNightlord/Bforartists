
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
bool is_anchor_node(uint tag_curr, uint tag_prev, uint tag_next, bool head_or_tail)
{ /* "Anchor" node has the minimal ID among its neighbours 
    * anchors breaks the chains into small segments, 
    * each seg len < 2(log_2(N)), where N is #tagged nodes*/
    return head_or_tail || ((tag_curr < tag_prev) && (tag_curr < tag_next));
}

bool is_head_node(ListRankingLink link, uint node_id)
{
    return (link.prev_node_id == node_id); 
}
#define FUNC_IS_NODE_HEAD is_head_node
bool is_tail_node(ListRankingLink link, uint node_id)
{
    return (link.next_node_id == node_id); 
}
#define FUNC_IS_NODE_TAIL is_tail_node
bool is_head_or_tail_node(ListRankingLink link, uint node_id)
{
    return is_head_node(link, node_id) || is_tail_node(link, node_id); 
}
#define FUNC_IS_NODE_HEAD_OR_TAIL is_head_or_tail_node



#if _KERNEL_MULTICOMPILE__TEST_LIST_RANKING_UPDATE_ANCHORS

    /* Public Interfaces ///////////////////////////////// */
    #define _FUNC_UPDATE_NODE_TAG CAT(CalcNewTag_, tag_list_ranking)
    #define _FUNC_IS_NODE_ANCHOR CAT(IsNodeAnchor_, tag_list_ranking)

	/**
	 * \brief Iteratively calculate & store tag for each node  
	 */
    void _FUNC_UPDATE_NODE_TAG(uint node_id, uint iter, uint tag_subbuff_size, out uint tag_new)
    {
        ListRankingLink link = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(node_id);
        uint tag =      FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(node_id,           iter, tag_subbuff_size);
        uint tag_next = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(link.next_node_id, iter, tag_subbuff_size); 
        
        tag_new = calc_node_tag(tag, tag_next);

        FUNC_DEVICE_STORE_LISTRANKING_NODE_TAG(tag_new, node_id, iter, tag_subbuff_size); 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)
    /**
	 * \brief Determine if a node is anchor, head/tail node must be anchor for code to work
	*/
    bool _FUNC_IS_NODE_ANCHOR(
        uint node_id, uint iter, uint tag_subbuff_size, 
        out ListRankingLink link
    ){
        link = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(node_id);
        bool b_list_head_or_tail = FUNC_IS_NODE_HEAD_OR_TAIL(link, node_id); 

        uint tag =      FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(node_id,           iter, tag_subbuff_size);
        uint tag_prev = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(link.prev_node_id, iter, tag_subbuff_size);
        uint tag_next = FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG(link.next_node_id, iter, tag_subbuff_size); 

        return is_anchor_node(tag, tag_prev, tag_next, b_list_head_or_tail); 
    }
#endif




#endif












