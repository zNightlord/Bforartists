#pragma BLENDER_REQUIRE(npr_strokegen_segloopconv1d_lib.glsl)


/* input buffers:
 * ssbo_in_segloopconv1d_data_
 * ssbo_out_segloopconv1d_data_
 * ssbo_segloopconv1d_patch_table_
 * ubo_segloopconv1d_
 * -----------------------------------------------
*/

#define T_To_Uint(x) x
#define Uint_To_T(x) x

/* Generate segment topology for testing the right half of patch table.
* - that is, we let the seg tail cross the thread group border. */
void get_test_seg_topo(uint gl_invoc_id, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
{
    #define SECTION_LEN 1731u
    
    uint section_id =  gl_invoc_id / SECTION_LEN;
    uint sec_head_id = section_id * SECTION_LEN;
    uint sec_tail_id = (sec_head_id + SECTION_LEN) - 1u;

    seg_len = max(1, wang_hash(section_id) % SECTION_LEN);
    uint id_seg_lc = (gl_invoc_id - sec_head_id) / seg_len;

    seg_head_id = id_seg_lc * seg_len + sec_head_id;
    seg_tail_id = (seg_head_id + seg_len) - 1u;

    if (sec_tail_id < seg_tail_id)
        seg_tail_id = sec_tail_id;
    
    uint last_valid_elem_id = uint(ubo_segloopconv1d_.num_conv_items) - 1u; 
    if (seg_head_id <= last_valid_elem_id && last_valid_elem_id <= seg_tail_id)
        seg_tail_id = last_valid_elem_id; 

    seg_len = seg_tail_id - seg_head_id + 1u;

    #undef SECTION_LEN
}


#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    
    uint seg_head_id, seg_tail_id, seg_len;
    get_test_seg_topo(idx, seg_head_id, seg_tail_id, seg_len);
    
    /* Setup segment topo */
    bool hf = (idx == seg_head_id);
    bool tf = (idx == seg_tail_id);
        
    /* Build & Output conv patch table */
    _FUNC_COMPUTE_PATCH_TABLE(
        blockIdx, groupIdx,
        hf, tf, seg_head_id, seg_tail_id, seg_len 
    );
    
    /* Output debug info */
    uint addr_dbg_st = idx << 2;
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_head_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_tail_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_len);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = ((hf ? 1 : (tf ? 2 : 0)));
}

#endif



#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;

    uint seg_head_id, seg_tail_id, seg_len;
    get_test_seg_topo(idx, seg_head_id, seg_tail_id, seg_len);

    uint convData; 
    _FUNC_SETUP_SEGLOOP1DCONV(
        blockIdx, groupIdx, 
        /*out*/ convData
    );
    ssbo_in_segloopconv1d_data_[idx] = convData;

    
    
    for (uint d = 1; d <= MAX_CONV_RADIUS; ++d)
    {
        uint neighData = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
            d, blockIdx, groupIdx,
            seg_len, seg_head_id
        );
        
        convData += neighData; 
    }

    for (uint d = 1; d <= MAX_CONV_RADIUS; ++d)
    {
        uint neighData = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
            d, blockIdx, groupIdx,
            seg_len, seg_head_id
        );
 
        convData += neighData;
    }
    ssbo_out_segloopconv1d_data_[idx] = convData;
}

#endif