#ifndef NPR_STROKEGEN_LIST_RANKING_INPUTS_LIB_H
#define NPR_STROKEGEN_LIST_RANKING_INPUTS_LIB_H


struct ListRankingLink
{
    uint prev_node_id; 
    uint next_node_id; 
};

struct JumpingInfo
{
    uint jump_next_anchor_id; 
    uint data; 
    bool is_list_head; 
    bool is_list_tail; 
};


uvec2 EncodePointerJumpingInfo(JumpingInfo ji)
{
    uint packed_x = (uint(ji.is_list_head) | (ji.jump_next_anchor_id << 1));
    uint packed_y = (uint(ji.is_list_tail) | (ji.data << 1)); 

    return uvec2(packed_x, packed_y); 
}
JumpingInfo DecodePointerJumpingInfo(uvec2 encoded)
{
    JumpingInfo ji; 
    ji.is_list_head = (1u == (encoded.x & 1u)); 
    ji.jump_next_anchor_id = (encoded.x >> 1);
    ji.is_list_tail = (1u == (encoded.y & 1u));  
    ji.data = (encoded.y >> 1); 

    return ji; 
}


/* Global Configs */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SETUP) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)

#define NUM_NODES_TOTAL ((ssbo_list_ranking_inputs_.num_nodes))

#define ANCHOR_COUNTER_PREV_ITER(splice_iter) ((ssbo_list_ranking_anchor_counters_[(splice_iter)])) 
#define ANCHOR_COUNTER(splice_iter) ((ssbo_list_ranking_anchor_counters_[(splice_iter) + 1]))  
#define SPLICE_COUNTER(splice_iter) ((ssbo_list_ranking_splice_counters_[(splice_iter) + 1]))  
#define SPLICE_COUNTER_RELINK(num_relink_iters, relink_iter) \
    ((ssbo_list_ranking_splice_counters_[(num_relink_iters) - (relink_iter)]))

#endif



/* Links */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING) || defined (_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)
uint func_device_load_listranking_test_next_node_id(uint node_id)
{ 
    uint node_offset = node_id * 2; 
    return ssbo_list_ranking_links_[node_offset + 1]; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID func_device_load_listranking_test_next_node_id

uint func_device_load_listranking_test_prev_node_id(uint node_id)
{ 
    uint node_offset = node_id * 2; 
    return ssbo_list_ranking_links_[node_offset]; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_PREV_NODE_ID func_device_load_listranking_test_prev_node_id

ListRankingLink func_device_load_listranking_test_node_links(uint node_id)
{ 
    ListRankingLink link; 

    link.prev_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_PREV_NODE_ID(node_id); 
    link.next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id); 

    return link; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS func_device_load_listranking_test_node_links
#endif


#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING)
void func_device_store_listranking_test_next_node_id(uint node_id, uint next_node_id)
{ 
    uint node_offset = node_id * 2; 
    ssbo_list_ranking_links_[node_offset + 1] = next_node_id; 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_NEXT_NODE_ID func_device_store_listranking_test_next_node_id

void func_device_store_listranking_test_prev_node_id(uint node_id, uint prev_node_id)
{ 
    uint node_offset = node_id * 2; 
    ssbo_list_ranking_links_[node_offset] = prev_node_id; 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_PREV_NODE_ID func_device_store_listranking_test_prev_node_id
#endif



/* anchor/spliced-to-node pointers */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)
void func_device_store_anchor_to_node_id(uint anchor_id, uint node_id)
{
    ssbo_list_ranking_anchor_to_node_out_[anchor_id] = node_id; 
}
#define FUNC_DEVICE_STORE_PER_ANCHOR_NODEID func_device_store_anchor_to_node_id

void func_device_store_spliced_node_id(uint spliced_id, uint node_id)
{
    ssbo_list_ranking_spliced_node_id_out_[spliced_id] = node_id; 
}
#define FUNC_DEVICE_STORE_PER_SPLICED_NODEID func_device_store_spliced_node_id
#endif



/* node-to-anchor pointers */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)
void func_device_store_node_to_anchor_id(uint node_id, uint anchor_id)
{
    ssbo_list_ranking_node_to_anchor_out_[node_id] = anchor_id; 
}
#define FUNC_DEVICE_STORE_PER_NODE_ANCHORID func_device_store_node_to_anchor_id
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)
uint func_device_load_node_to_anchor_id(uint node_id)
{
    return ssbo_list_ranking_node_to_anchor_in_[node_id]; 
}
#define FUNC_DEVICE_LOAD_PER_NODE_ANCHORID func_device_load_node_to_anchor_id
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined (_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)

uint func_load_splicing_thread_node_id(uint thread_idx, uint splicing_iter = 1u)
{ /* node_id acts as a pointer to access node buffers */
    uint node_id = (splicing_iter == 0u) ? thread_idx
        : ssbo_list_ranking_anchor_to_node_in_[thread_idx]; 
    return node_id; 
}
#define FUNC_GET_NODE_ID_FOR_ANCHOR func_load_splicing_thread_node_id
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING)
uint func_device_load_spliced_node_id(uint spliced_id)
{
    return ssbo_list_ranking_spliced_node_id_in_[spliced_id]; 
}
#define FUNC_GET_NODE_ID_FOR_SPLICED func_device_load_spliced_node_id
#endif


/* node tags */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING)
void func_device_store_listranking_test_tag(uint node_id, uint tag)
{
    ssbo_list_ranking_tags_out_[node_id] = tag; 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_TAG func_device_store_listranking_test_tag
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING)
uint func_device_load_listranking_test_tag(uint node_id)
{ 
    return ssbo_list_ranking_tags_in_[node_id]; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG func_device_load_listranking_test_tag
#endif


/* node ranks */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)
uint func_device_load_node_rank(uint node_id)
{
    return ssbo_list_ranking_ranks_[node_id]; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK func_device_load_node_rank

void func_device_store_node_rank(uint node_id, uint rank)
{
    ssbo_list_ranking_ranks_[node_id] = rank; 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK func_device_store_node_rank

uint func_device_load_node_rank_with_head_tail_flags(uint node_id, out bool is_head, out bool is_tail)
{
    uint data = ssbo_list_ranking_ranks_[node_id]; 
    uint rank = data & 0x3fffffffu; 
    is_head = (((data >> 31) & 1u) != 0u);
    is_tail = (((data >> 30) & 1u) != 0u); 
    return rank; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK_WITH_HEAD_TAIL_BITS func_device_load_node_rank_with_head_tail_flags

void func_device_store_node_rank_with_head_tail_flags(uint node_id, uint rank, bool is_head, bool is_tail)
{
    ssbo_list_ranking_ranks_[node_id] = ((uint(is_head) << 31) | (uint(is_tail) << 30) | (rank & 0x3fffffffu)); 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_RANK_WITH_HEAD_TAIL_BITS func_device_store_node_rank_with_head_tail_flags
#endif



/* Per-anchor thread validation for all kernels after the compaction */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)

bool validate_thread_per_anchor(uint mapped_anchor_id, uint num_anchors)
{
    return mapped_anchor_id < num_anchors; 
}
#define FUNC_DEVICE_VALIDATE_THREAD_PER_ANCHOR validate_thread_per_anchor

/* Padded to odd number of iterations. *
 * Must be deterministic about the ping-pong buffer(s)  */
uint get_num_jump_iters(uint num_anchors)
{
    uint num_iters = uint(ceil(log2(float(num_anchors))) + .0001f);
    num_iters = ((num_iters % 2u) == 0u) ? (num_iters + 1u) : (num_iters); 
    return num_iters; 
}
#define FUNC_GET_NUM_JUMPS get_num_jump_iters

#endif



/* Packed Jumping Info */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES)
JumpingInfo func_device_load_per_anchor_jumping_info(uint anchor_id)
{
    uvec2 jump_data_packed = uvec2(
        ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[anchor_id*2u ], 
        ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[anchor_id*2u+1u]
    );  
     
    return DecodePointerJumpingInfo(jump_data_packed); 
}
#define FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO func_device_load_per_anchor_jumping_info

void func_device_store_per_anchor_jumping_pointer(uint anchor_id, JumpingInfo ji)
{
    uvec2 jump_data_packed = EncodePointerJumpingInfo(ji);
    ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[anchor_id*2u]    = jump_data_packed.x;
    ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[anchor_id*2u+1u] = jump_data_packed.y;  
}
#define FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO func_device_store_per_anchor_jumping_pointer
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)
bool func_is_loop_breaking_pass()
{
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
    return (pc_listranking_ranking_pass_with_broken_loops_ == 0); 
#else
    return false; 
#endif
}
#define IS_LOOP_BREAKING_PASS func_is_loop_breaking_pass 

bool func_is_loop_ranking_pass()
{
#if defined(_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD)
    return (pc_listranking_ranking_pass_with_broken_loops_ > 0); 
#else
    return false; 
#endif
}
#define IS_LOOP_RANKING_PASS func_is_loop_ranking_pass 

JumpingInfo func_device_init_per_anchor_jumping_info(uint anchor_id, uint splicing_iter)
{
    uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter);  

    JumpingInfo ji; 
    uint prev_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_PREV_NODE_ID(node_id); 
    uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id); 
    ji.jump_next_anchor_id = FUNC_DEVICE_LOAD_PER_NODE_ANCHORID(next_node_id); 

    ji.is_list_tail = (next_node_id == node_id); 
    ji.is_list_head = (prev_node_id == node_id); 

    if (IS_LOOP_BREAKING_PASS())
    {
        ji.data = node_id; /* unique for each anchor */
    }
    else
    {
        ji.data = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(node_id); 
    }



    return ji; 
}
#define FUNC_DEVICE_INIT_PER_ANCHOR_LIST_JUMPING_INFO func_device_init_per_anchor_jumping_info

JumpingInfo func_device_update_anchor_jumping_info(
    JumpingInfo ji, JumpingInfo ji_next, out bool jumped_to_end
){
    JumpingInfo ji_updated = ji; /* .is_list_head and .is_list_tail are not updated */ 
    /* Only update if the next jump has not reached the end of the list */
    jumped_to_end = ji_next.jump_next_anchor_id == ji.jump_next_anchor_id; 
    if (false == jumped_to_end)
    { 
        if (IS_LOOP_BREAKING_PASS())
            ji_updated.data = max(ji_updated.data, ji_next.data); /* find node with max code as head */ 
        else
            ji_updated.data += ji_next.data; /* accumulate rank */

        ji_updated.jump_next_anchor_id = ji_next.jump_next_anchor_id; 
    }

    return ji_updated;  
}
#define FUNC_DEVICE_UPDATE_ANCHOR_LIST_JUMPING_INFO func_device_update_anchor_jumping_info
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)
uint func_allocate_space_for_list(uint list_len)
{
    uint list_start = atomicAdd(ssbo_list_ranking_addressing_counters_[0], list_len); 
    return list_start;  
}
#define FUNC_DEVICE_ALLOC_LIST_ADDR func_allocate_space_for_list

void func_device_broadcast_list_topology(uint broadcast_anchor_id, uint list_len, uint list_start_addr)
{
    ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[broadcast_anchor_id * 2]     = list_start_addr; 
    ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[broadcast_anchor_id * 2 + 1] = list_len;  
}
#define FUNC_DEVICE_BROADCAST_LIST_TOPOLOGY func_device_broadcast_list_topology
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)
void func_device_retrieve_list_topology(uint broadcast_anchor_id, out uint list_len, out uint list_start_addr)
{
    list_start_addr = ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[broadcast_anchor_id * 2]; 
    list_len        = ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[broadcast_anchor_id * 2 + 1];  
}
#define FUNC_DEVICE_RETRIEVE_LIST_TOPOLOGY func_device_retrieve_list_topology

void func_device_store_list_topology(uint node_id, uint list_len, uint list_start_addr)
{
    ssbo_list_ranking_serialized_topo_[node_id * 2]     = list_start_addr; 
    ssbo_list_ranking_serialized_topo_[node_id * 2 + 1] = list_len;  
}
#define FUNC_DEVICE_STORE_LIST_TOPOLOGY func_device_store_list_topology 
#endif




#endif 



