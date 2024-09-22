#pragma BLENDER_REQUIRE(npr_strokegen_segloopconv1D_inputs_advanced_lib.glsl)


 
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


    /* Customizable ---------------------------------------------------------- 
     * Use this function to fetch the segment topology */
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
    /* Default style of convolution. Visits each neighbor once. 
     * Sufficient for most cases  */
#if !defined(SEGLOOPCONV1D_USE_ADVANCED_INPUT)
    for (uint d = 1; d <= MAX_CONV_RADIUS; ++d)
    {
        DATA_TYPE_LOOPCONV1D neigh_data = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
            d, blockIdx, groupIdx,
            seg_len, seg_head_id
        );

        /* Customizable ------------------------------------------------------------- */
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

        /* Customizable ------------------------------------------------------------- */
        FUNC_CONVOLUTION(
            /*mov_left*/false, d, seg_is_loop, seg_head_id, seg_tail_id, idx, 
            neigh_data, orig_data, /*inout*/conv_temp_data); 
    }
#else
    /* Advanced convolution. Use this when sophisticated logic is involved  */
    /* Customizable --------- */
    FUNC_CONVOLUTION_ADVANCED(idx, blockIdx, groupIdx, seg_is_loop, seg_head_id, seg_tail_id, seg_len, /*out*/conv_temp_data); 
#endif


    /* Post-process && Output convolution results. */
    /* Customizable ------------------------------------------------------------- */
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
        ssbo_out_segloopconv1d_data_[idx] = floatBitsToUint(conv_temp_data.val);
    #endif
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__SEG_DENOISING)
        ContourFlags efs_ori = decode_contour_flags(orig_data);
        
        uint num_snakes_left = uint(MAX_CONV_RADIUS); 
        uint num_snakes_right = uint(MAX_CONV_RADIUS);
        // if (!seg_is_loop) {
        //     num_snakes_left = min(MAX_CONV_RADIUS, idx - seg_head_id/*edge loop head, not the seg head*/);
        //     num_snakes_right = min(MAX_CONV_RADIUS, seg_tail_id - idx);
        // }
        uint num_pstv_cusps_left = (conv_temp_data.num_pstv_cusps_left);
        uint num_pstv_cusps_right = uint(conv_temp_data.num_pstv_cusps_right);

        // supress small noises 
        uint num_ptsv_cusps_total = num_pstv_cusps_left + num_pstv_cusps_right + (efs_ori.cusp_func_pstv ? 1u : 0u);
        uint num_snakes_total     = num_snakes_left     + num_snakes_right     + 1u;
        bool stochastic_pstv_cusp  = (num_snakes_total - 2u <= num_ptsv_cusps_total);
        bool stochastic_ngtv_cusp  =  (num_ptsv_cusps_total <= 2u);  
        set_contour_flags_dbg_flag_0(
            ((!efs_ori.cusp_func_pstv) && stochastic_pstv_cusp) 
            || (efs_ori.cusp_func_pstv && stochastic_ngtv_cusp), 
            efs_ori
        );  
        if (stochastic_pstv_cusp || stochastic_ngtv_cusp)
        {
            set_contour_seg_head(false, efs_ori);
            bool cusp_pstv_fixed = stochastic_pstv_cusp ? true : false; 
            // set_contour_cusp_flags(cusp_pstv_fixed, efs_ori); 
        }

        if (idx < get_num_items())
            ssbo_out_segloopconv1d_data_[idx] = encode_contour_flags(efs_ori);
    #endif
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0)
        float curvature = conv_temp_data.corner_curv; 
        if (idx < get_num_items())
            store_ssbo_contour_2d_sample_geometry__corner_curvature(idx, curvature, get_num_items()); 
    #endif
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_1)
        ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(idx); 
        if (idx < get_num_items())
        {
            cf.is_curv_minima = conv_temp_data.is_local_minima; 
            cf.seg_head = conv_temp_data.is_local_minima || cf.seg_head_contour || cf.seg_head_clipped; 
            cf.is_corner = conv_temp_data.is_local_maxima;
            store_ssbo_contour_2d_sample_topology__flags(idx, cf); 
        }
    #endif
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CALC_TANGENT_CURVATURE)
        if (idx < get_num_items())
            store_ssbo_contour_2d_sample_geometry__tangent(idx, conv_temp_data.tangent, get_num_items());             
    #endif
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2D_SAMPLE_VISIBILITY_DENOISING)
        if (idx < get_num_items())
        {
            ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(idx); 
            cf.occluded_filtered = conv_temp_data.filtered_occlusion; 
            store_ssbo_contour_2d_sample_topology__flags(idx, cf); 
        }
    #endif
    /* -------------------------------------------------------------------------- */
}

#endif