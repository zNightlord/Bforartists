#ifndef NPR_STROKEGEN_LIST_RANKING_INPUTS_LIB_H
#define NPR_STROKEGEN_LIST_RANKING_INPUTS_LIB_H


struct ListRankingLink
{
    uint prev_node_id; 
    uint next_node_id; 
};

struct ContractionInfo
{
    uint id; 
    bool is_anchor; 
    
    /* list-end can be head or tail, 
     * depending on the jumping direction. 
     * by default is_list_end == is_list_tail */ 
    bool is_list_end; 
};

struct JumpingInfo
{
    uint jump_next_anchor_id; 
    uint data; 
    /* bool is_list_head; */ 
    bool is_list_tail; 
};
struct JumpingListInfo
{
    uint list_len; 
    uint list_addr; 
}; 
/* Note: "Specifically, structs cannot be used as input/output variables." 
 * See https://www.khronos.org/opengl/wiki/Data_Type_(GLSL)#Structs
*/

/* TODO: Instead of this port you fucking bit encoding functions from MPipeline */
#define BITMASK_LIST_RANKING_CONTRACTION_INFO_NODE_ID 0x3fffffffu
uint EncodeContractionInfo(bool is_anchor, bool is_list_end, uint anchor_id, uint next_node_id)
{
    uint encoded = 0u; 
    if (is_anchor)
    {
        encoded = (anchor_id | (1u << 31u) | (uint(is_list_end) << 30u)); 
    }
    else
    {
        encoded = next_node_id & BITMASK_LIST_RANKING_CONTRACTION_INFO_NODE_ID;
    }
    
    return encoded; 
}

ContractionInfo DecodeContractionInfo(uint encoded)
{
    bool is_anchor = ((encoded >> 31u) == 1u);
    if (is_anchor)
    {
        bool is_list_end = (((encoded >> 30u) & 1u) == 1u);
        return ContractionInfo(encoded & BITMASK_LIST_RANKING_CONTRACTION_INFO_NODE_ID, true, is_list_end); 
    }
    else
    {
        return ContractionInfo(encoded & BITMASK_LIST_RANKING_CONTRACTION_INFO_NODE_ID, false, false);
    }

    return ContractionInfo(0u, false, false); 
}

uvec2 EncodePointerJumpingInfo(JumpingInfo ji)
{
    uint packed_x = /* (uint(ji.is_list_head) |  */(ji.jump_next_anchor_id << 1)/* ) */;
    uint packed_y = (uint(ji.is_list_tail) | (ji.data << 1)); 
    return uvec2(packed_x, packed_y); 
}
JumpingInfo DecodePointerJumpingInfo(uvec2 encoded)
{
    JumpingInfo ji; 
    /* ji.is_list_head = (1u == (encoded.x & 1u)); */
    ji.jump_next_anchor_id = (encoded.x >> 1);
    ji.is_list_tail = (1u == (encoded.y & 1u));  
    ji.data = (encoded.y >> 1); 

    return ji; 
}

/* Global Configs */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES)

#define NUM_NODES_TOTAL ((ubo_list_ranking_splicing_.num_nodes))
#define SUB_BUFF_SIZE ((ubo_list_ranking_splicing_.subbuff_size))

#endif



/* Links */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) 
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


#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES)
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
    ssbo_list_ranking_spliced_node_id_[spliced_id] = node_id; 
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

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)
uint func_device_load_node_to_anchor_id(uint node_id)
{
    return ssbo_list_ranking_node_to_anchor_in_[node_id]; 
}
#define FUNC_DEVICE_LOAD_PER_NODE_ANCHORID func_device_load_node_to_anchor_id
#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined (_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) 

uint func_load_splicing_thread_node_id(uint thread_idx, uint splicing_iter)
{ /* node_id acts as a pointer to access node buffers */
    uint node_id = (splicing_iter == 0u) ? thread_idx
        : ssbo_list_ranking_anchor_to_node_in_[thread_idx]; 
    return node_id; 
}
#define FUNC_GET_NODE_ID_FOR_ANCHOR func_load_splicing_thread_node_id
#endif

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES)
uint func_device_load_spliced_node_id(uint spliced_id)
{
    return ssbo_list_ranking_spliced_node_id_[spliced_id]; 
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
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) 
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
#endif



/* Per-anchor thread validation for all kernels after the compaction */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)

bool validate_thread_per_anchor(uint mapped_anchor_id, uint num_anchors)
{
    return mapped_anchor_id < num_anchors; 
}
#define FUNC_DEVICE_VALIDATE_THREAD_PER_ANCHOR validate_thread_per_anchor

/* Remember we use the 0th iter to init data, 
 * actual jump iter starts form #1 */
uint get_num_jump_iters(uint num_anchors)
{
    uint num_iters = uint(ceil(log2(float(num_anchors))) + .0001f);
    return num_iters; 
}
#define FUNC_GET_NUM_JUMPS get_num_jump_iters

#endif




#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)
JumpingInfo func_device_init_per_anchor_jumping_info(uint anchor_id, uint splicing_iter)
{
    uint node_id = FUNC_GET_NODE_ID_FOR_ANCHOR(anchor_id, splicing_iter);  

    JumpingInfo ji; 
    uint next_node_id = FUNC_DEVICE_LOAD_LISTRANKING_NODE_NEXT_NODE_ID(node_id); 
    ji.jump_next_anchor_id = FUNC_DEVICE_LOAD_PER_NODE_ANCHORID(next_node_id); 
    ji.data = FUNC_DEVICE_LOAD_LISTRANKING_NODE_RANK(node_id); 

    ji.is_list_tail = (next_node_id == node_id); 

    return ji; 
}
#define FUNC_DEVICE_INIT_PER_ANCHOR_LIST_JUMPING_INFO func_device_init_per_anchor_jumping_info


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


JumpingInfo func_device_update_anchor_jumping_info(
    JumpingInfo ji, JumpingInfo ji_next, out bool jumped_to_end
){
    JumpingInfo ji_updated = ji; /* .is_list_head and .is_list_tail are not updated */ 
    /* Only update if the next jump has not reached the end of the list */
    jumped_to_end = ji_next.jump_next_anchor_id == ji.jump_next_anchor_id; 
    if (false == jumped_to_end)
    { 
        ji_updated.data += ji_next.data; 
        ji_updated.jump_next_anchor_id = ji_next.jump_next_anchor_id; 
    }

    return ji_updated; 
}
#define FUNC_DEVICE_UPDATE_ANCHOR_LIST_JUMPING_INFO func_device_update_anchor_jumping_info

#endif




#endif 



