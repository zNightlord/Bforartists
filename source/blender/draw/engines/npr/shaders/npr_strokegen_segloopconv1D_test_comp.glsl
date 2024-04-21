#pragma BLENDER_REQUIRE(npr_strokegen_segloopconv1d_lib.glsl)







#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    
    bool seg_is_loop; 
    uint seg_head_id, seg_tail_id, seg_len;
    FUNC_GET_LOOP_TOPOLOGY(idx, seg_is_loop, seg_head_id, seg_tail_id, seg_len);
    
    /* Setup segment topo */
    bool hf = (idx == seg_head_id);
    bool tf = (idx == seg_tail_id);
        
    /* Build & Output conv patch table */
    _FUNC_COMPUTE_PATCH_TABLE(
        blockIdx, groupIdx,
        hf, tf, seg_head_id, seg_tail_id, seg_len 
    );
    
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
    /* Output debug info */
    uint addr_dbg_st = idx << 2;
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_head_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_tail_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_len);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = ((hf ? 1 : (tf ? 2 : 0)));
#endif
}

#endif



#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;

    bool seg_is_loop; 
    uint seg_head_id, seg_tail_id, seg_len;
    FUNC_GET_LOOP_TOPOLOGY(idx, seg_is_loop, seg_head_id, seg_tail_id, seg_len);

    DATA_TYPE_LOOPCONV1D orig_data; 
    _FUNC_SETUP_SEGLOOP1DCONV(
        blockIdx, groupIdx, 
        /*out*/ orig_data
    );
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
    ssbo_in_segloopconv1d_data_[idx] = floatBitsToUint(orig_data); // input is natively calculated rather than loaded from buffer
#endif

    T_CONV_TEMP_DATA conv_temp_data; 
    for (uint d = 1; d <= MAX_CONV_RADIUS; ++d)
    {
        DATA_TYPE_LOOPCONV1D neigh_data = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
            d, blockIdx, groupIdx,
            seg_len, seg_head_id
        );

        FUNC_CONVOLUTION(
            /*mov_left*/true, d, seg_is_loop, seg_head_id, seg_tail_id, idx/*item_id*/, 
            neigh_data, orig_data, /*inout*/conv_temp_data
        ); 
    }

    for (uint d = 1; d <= MAX_CONV_RADIUS; ++d)
    {
        DATA_TYPE_LOOPCONV1D neigh_data = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
            d, blockIdx, groupIdx,
            seg_len, seg_head_id
        );
 
        FUNC_CONVOLUTION(
            /*mov_left*/false, d, seg_is_loop, seg_head_id, seg_tail_id, idx, 
            neigh_data, orig_data, /*inout*/conv_temp_data); 
    }

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
        ssbo_out_segloopconv1d_data_[idx] = floatBitsToUint(conv_temp_data.val);
    #endif
}

#endif