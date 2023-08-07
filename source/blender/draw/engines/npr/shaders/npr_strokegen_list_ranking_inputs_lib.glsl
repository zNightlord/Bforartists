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
    bool is_list_head; 
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




#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_UPDATE_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_RANKING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING) 
/* buffers:
 ssbo_list_ranking_links_;
 ssbo_list_ranking_anchor_links_; 
 ssbo_list_ranking_output_;
 ubo_list_ranking_tagging_; */
 
ListRankingLink func_device_load_listranking_test_node_links(uint node_id)
{ 
    ListRankingLink link; 
    uint node_offset = node_id * 2; 
    link.prev_node_id = ssbo_list_ranking_links_[node_offset];
    link.next_node_id = ssbo_list_ranking_links_[node_offset + 1];

    return link; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS func_device_load_listranking_test_node_links

#endif




#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_UPDATE_ANCHORS) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)
/* buffers:
 ssbo_list_ranking_links_;
 ssbo_list_ranking_anchor_links_; 
 ssbo_list_ranking_output_;
 ubo_list_ranking_tagging_; */

uint func_device_load_listranking_test_tag(uint node_id, uint iter, uint subbuff_size)
{ 
    if (iter == 0u) /* init tag with node id */
        return node_id;

    uint subbuff_beg = subbuff_size * ((iter + 1) % 2u); /* output from last iter */
    uint node_offset = subbuff_beg + node_id; 
    return ssbo_list_ranking_output_[node_offset]; 
}
#define FUNC_DEVICE_LOAD_LISTRANKING_NODE_TAG func_device_load_listranking_test_tag

void func_device_store_listranking_test_tag(uint tag, uint node_id, uint iter, uint subbuff_size)
{ /* reuse rank buffer */
    uint subbuff_beg = subbuff_size * (iter % 2u); /* output to next iter */
    uint node_offset = subbuff_beg + node_id; 
    ssbo_list_ranking_output_[node_offset] = tag; 
}
#define FUNC_DEVICE_STORE_LISTRANKING_NODE_TAG func_device_store_listranking_test_tag

#endif





#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS)

void func_device_store_anchor_to_node_id(uint anchor_id, uint node_id)
{
    ssbo_list_ranking_anchor_to_node_[anchor_id] = node_id; 
}
#define FUNC_DEVICE_STORE_PER_ANCHOR_NODEID func_device_store_anchor_to_node_id

void func_device_cache_contraction_info(
    uint node_id, uint next_node_id, bool is_anchor, uint anchor_id)
{
    bool is_list_end = (next_node_id == node_id); 
    ssbo_list_ranking_node_to_anchor_[node_id] = EncodeContractionInfo(
        is_anchor, is_list_end, anchor_id, next_node_id
    ); 
}
#define FUNC_DEVICE_STORE_PER_NODE_CONTRACTION_INFO func_device_cache_contraction_info

#endif


/* Per-anchor thread validation for all kernels after the compaction */
#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_RANKING) || defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)

bool validate_thread_per_anchor(uint mapped_anchor_id, uint num_anchors)
{
    return mapped_anchor_id < num_anchors; 
}
#define FUNC_DEVICE_VALIDATE_THREAD_PER_ANCHOR validate_thread_per_anchor

uint get_num_anchors()
{
    return ssbo_list_ranking_atomics_.counter_anchors; 
}
#define FUNC_GET_NUM_ANCHORS get_num_anchors

/* Remember we use the 0th iter to init data, 
 * actual jump iter starts form #1 */
uint get_num_jump_iters()
{
    uint num_anchors = FUNC_GET_NUM_ANCHORS(); 
    uint num_iters = uint(ceil(log2(float(num_anchors))) + .0001f);

    return num_iters; 
}
#define FUNC_GET_NUM_JUMPS get_num_jump_iters


uint func_device_load_node_id_for_anchor(uint anchor_id)
{
    uint node_id = ssbo_list_ranking_anchor_to_node_[anchor_id]; 
    return node_id; 
}
#define FUNC_DEVICE_LOAD_PER_ANCHOR_NODEID(anchor_id) func_device_load_node_id_for_anchor(anchor_id)

#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_RANKING)

void func_device_load_contraction_info(uint node_id, inout ContractionInfo contraction_info)
{
    contraction_info = DecodeContractionInfo(ssbo_list_ranking_node_to_anchor_[node_id]);  
}
#define FUNC_DEVICE_LOAD_PER_NODE_CONTRACTION_INFO(node_id, contraction_info) \
     func_device_load_contraction_info(node_id, contraction_info)

void func_device_store_node_to_anchor(uint node_id, uint anchor_id)
{
    ssbo_list_ranking_node_to_anchor_[node_id] = anchor_id;  
}
#define FUNC_DEVICE_STORE_PER_NODE_TO_ANCHOR_ID func_device_store_node_to_anchor

void func_device_store_next_anchor_id(uint anchor_id, uint next_anchor_id)
{
    ssbo_list_ranking_anchor_to_next_anchor_[anchor_id] = next_anchor_id; 
}
#define FUNC_DEVICE_STORE_PER_ANCHOR_NEXT_ANCHORID(anchor_id, next_anchor_id) \
    func_device_store_next_anchor_id(anchor_id, next_anchor_id)

void func_device_store_per_node_sublist_rank(uint node_id, uint sublist_rank)
{
    ssbo_list_ranking_output_[node_id] = sublist_rank; 
}
#define FUNC_DEVICE_STORE_PER_NODE_SUBLIST_RANK(node_id, sublist_rank) func_device_store_per_node_sublist_rank(node_id, sublist_rank)

#endif




/*
buffers to add:
ssbo_list_ranking_per_anchor_sublist_jumping_info_[]; 
pc_listranking_jumping_iter_; 
*/

#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING)

uint func_device_load_next_anchor_id(uint anchor_id)
{
    return ssbo_list_ranking_anchor_to_next_anchor_[anchor_id]; 
}
#define FUNC_DEVICE_LOAD_PER_ANCHOR_NEXT_ANCHORID func_device_load_next_anchor_id

uint func_device_load_sublist_length(uint anchor_node_id)
{
    return ssbo_list_ranking_output_[anchor_node_id]; 
}
#define FUNC_DEVICE_LOAD_SUBLIST_LENGTH func_device_load_sublist_length

JumpingInfo func_device_init_per_anchor_jumping_info(uint anchor_id, uint node_id)
{
    JumpingInfo ji; 
    ji.jump_next_anchor_id = FUNC_DEVICE_LOAD_PER_ANCHOR_NEXT_ANCHORID(anchor_id); 
    ji.data = FUNC_DEVICE_LOAD_SUBLIST_LENGTH(node_id); 

    ListRankingLink link = FUNC_DEVICE_LOAD_LISTRANKING_NODE_LINKS(node_id); 
    ji.is_list_head = (link.prev_node_id == node_id); 
    ji.is_list_tail = (link.next_node_id == node_id); 

    return ji; 
}
#define FUNC_DEVICE_INIT_PER_ANCHOR_LIST_JUMPING_INFO func_device_init_per_anchor_jumping_info


JumpingInfo func_device_load_per_anchor_jumping_info(uint anchor_id, uint jump_iter, uint subbuff_size)
{
    uint subbuff_beg = subbuff_size * ((jump_iter + 1u) % 2u); /* output from last iter */
    uint ld_offset = (subbuff_beg + anchor_id);
    uvec2 jump_data_packed = uvec2(
        ssbo_list_ranking_per_anchor_sublist_jumping_info_[ld_offset*2u ], 
        ssbo_list_ranking_per_anchor_sublist_jumping_info_[ld_offset*2u + 1u]
    );  
     
    return DecodePointerJumpingInfo(jump_data_packed); 
}
#define FUNC_DEVICE_LOAD_PER_ANCHOR_LIST_JUMPING_INFO func_device_load_per_anchor_jumping_info


void func_device_store_per_anchor_jumping_pointer(uint anchor_id, uint jump_iter, JumpingInfo ji, uint subbuff_size)
{
    uint subbuff_beg = subbuff_size * (jump_iter % 2u); /* output to next iter */
    uint st_offset = (subbuff_beg + anchor_id);

    uvec2 jump_data_packed = EncodePointerJumpingInfo(ji);
    ssbo_list_ranking_per_anchor_sublist_jumping_info_[st_offset*2u]      = jump_data_packed.x;
    ssbo_list_ranking_per_anchor_sublist_jumping_info_[st_offset*2u + 1u] = jump_data_packed.y;  
}
#define FUNC_DEVICE_STORE_PER_ANCHOR_LIST_JUMPING_INFO func_device_store_per_anchor_jumping_pointer

JumpingInfo func_device_update_anchor_jumping_info(
    JumpingInfo ji, JumpingInfo ji_next, out bool jumped_to_end
){
    JumpingInfo ji_updated = ji; /* .is_list_head and .is_list_tail are not updated */ 
    /* Only update if the next jump has not reached the end of the lsit */
    jumped_to_end = ji_next.jump_next_anchor_id == ji.jump_next_anchor_id; 
    if (false == jumped_to_end)
    { 
        ji_updated.data += ji_next.data; 
        ji_updated.jump_next_anchor_id = ji_next.jump_next_anchor_id; 
    }

    return ji_updated; 
}
#define FUNC_DEVICE_UPDATE_ANCHOR_LIST_JUMPING_INFO func_device_update_anchor_jumping_info

void func_device_store_per_anchor_list_len_and_addr(
    uint tail_anchor_id, JumpingListInfo list_info, uint subbuff_size
){
    uint subbuff_beg = subbuff_size * 2u; /* only output at the last iter */ 
    uint st_offset = (subbuff_beg + tail_anchor_id * 2u);

    ssbo_list_ranking_per_anchor_sublist_jumping_info_[st_offset]      = list_info.list_len;
    ssbo_list_ranking_per_anchor_sublist_jumping_info_[st_offset + 1u] = list_info.list_addr;   
}
#define FUNC_DEVICE_STORE_TAIL_ANCHOR_LIST_INFO func_device_store_per_anchor_list_len_and_addr

#endif



#if defined(_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_RANKING)


#endif




#endif 



