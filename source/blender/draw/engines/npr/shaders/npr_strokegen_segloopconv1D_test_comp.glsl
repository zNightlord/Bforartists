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
        CuspSegmentDenoiseData csdd_ori = decode_cusp_segment_denoise_data(orig_data); 
        ContourFlags cfs_ori = csdd_ori.cf;

        uint num_snakes_left = uint(MAX_CONV_RADIUS); 
        uint num_snakes_right = uint(MAX_CONV_RADIUS);
        bool is_non_looped_curve_head = (!seg_is_loop) && seg_head_id == idx; 
        // if (!seg_is_loop) {
        //     num_snakes_left = min(MAX_CONV_RADIUS, idx - seg_head_id/*edge loop head, not the seg head*/);
        //     num_snakes_right = min(MAX_CONV_RADIUS, seg_tail_id - idx);
        // }
        uint num_pstv_cusps_left = (conv_temp_data.num_pstv_cusps_left);
        uint num_pstv_cusps_right = uint(conv_temp_data.num_pstv_cusps_right);

        // ---------------------------------------------
        // Supress very short cusp segments 
        // ---------------------------------------------
        uint num_ptsv_cusps_total = num_pstv_cusps_left + num_pstv_cusps_right + (cfs_ori.cusp_func_pstv ? 1u : 0u);
        uint num_snakes_total     = num_snakes_left     + num_snakes_right     + 1u;

        float pstv_seg_len_total = conv_temp_data.pstv_seg_len_left + conv_temp_data.pstv_seg_len_right;
        float ngtv_seg_len_total = conv_temp_data.ngtv_seg_len_left + conv_temp_data.ngtv_seg_len_right;
        float seg_len_total = pstv_seg_len_total + ngtv_seg_len_total;
		
        const float cusp_len_threshold = 5e-3f; // TODO: adaptive to mesh size
        bool stochastic_pstv_seg = (pstv_seg_len_total / seg_len_total >= 0.95f) && ngtv_seg_len_total < cusp_len_threshold;
		bool stochastic_ngtv_seg = (ngtv_seg_len_total / seg_len_total >= 0.95f) && pstv_seg_len_total < cusp_len_threshold;
         
        if (cfs_ori.seg_head && (!is_non_looped_curve_head))
        {    
            if (stochastic_pstv_seg || stochastic_ngtv_seg)
            {  
                ContourCurveTopo cct = initialize_contour_curve_topo(cfs_ori.looped_curve, seg_head_id, seg_len); 
                uint cusp_seg_len = ssbo_contour_snake_seg_len_[idx]; 
                uint next_cusp_seg_head_id = move_contour_id_along_loop(cct, idx, float(cusp_seg_len));
                uint prev_cusp_seg_tail_id = move_contour_id_along_loop(cct, idx, -1.0f);
                uint prev_cusp_seg_len = ssbo_contour_snake_seg_len_[prev_cusp_seg_tail_id];
                uint prev_cusp_seg_head_id = move_contour_id_along_loop(cct, prev_cusp_seg_tail_id, -float(prev_cusp_seg_len - 1u));

                bool major_cusp_seg; // if this segment has the the most common cusp function sign 
                if (stochastic_pstv_seg) major_cusp_seg = cfs_ori.cusp_func_pstv; 
                if (stochastic_ngtv_seg) major_cusp_seg = !cfs_ori.cusp_func_pstv; 
				
                { // Long-range segment merging
                    uint merge_cusp_seg_head_id; 
                    if (major_cusp_seg) 
                    { // for the major cusp segment, we merge its prev noisy segment to the prev-prev major segment
                        uint prev_cusp_seg_tail_id = 
                            move_contour_id_along_loop(cct, idx, -1.0f);
                        uint prev_cusp_seg_len = ssbo_contour_snake_seg_len_[prev_cusp_seg_tail_id];
                        uint prev_cusp_seg_head_id = 
                            move_contour_id_along_loop(cct, prev_cusp_seg_tail_id, -float(prev_cusp_seg_len - 1u));

                        merge_cusp_seg_head_id = prev_cusp_seg_head_id; 
                    } else { 
                        // for the noise cusp segment, we merge it with its next long segment
                        uint next_cusp_seg_head_id = move_contour_id_along_loop(cct, idx, float(cusp_seg_len));

                        merge_cusp_seg_head_id = next_cusp_seg_head_id; 
                    }

					ContourFlags cfs_merged_seg_head; { // have to fetch cotour flags from input buffer to avoid racing condition
						CuspSegmentDenoiseData csdd;
						uvec4 conv_input_enc = uvec4(0u);
						Load4(ssbo_in_segloopconv1d_data_, merge_cusp_seg_head_id, conv_input_enc); 
						csdd = decode_cusp_segment_denoise_data(conv_input_enc);

						cfs_merged_seg_head = csdd.cf; 
					}
					set_contour_seg_head(false, cfs_merged_seg_head); 
					
					if (idx < get_num_items())
						ssbo_out_segloopconv1d_data_[merge_cusp_seg_head_id] = encode_contour_flags(cfs_merged_seg_head);
				}
 
				{ // Merge with the prev sub-segment
					// Merge with prev segment
					set_contour_seg_head(false, cfs_ori);
                    if (idx < get_num_items())
                       ssbo_out_segloopconv1d_data_[idx] = encode_contour_flags(cfs_ori);
				}
            }
        }



        #if defined(DEBUG)
            uint offset = get_num_items(); 
            uvec4 dbg_data = uvec4(
                floatBitsToUint(seg_len_total), // float(num_ptsv_cusps_total)), // pstv_seg_len_total), 
                floatBitsToUint(float(num_ptsv_cusps_total)), // pstv_seg_len_total), 
                floatBitsToUint(pstv_seg_len_total / seg_len_total), 
                floatBitsToUint(ngtv_seg_len_total / seg_len_total)
            ); 
            ssbo_out_segloopconv1d_data_[offset + idx*4+0] = dbg_data[0]; 
            ssbo_out_segloopconv1d_data_[offset + idx*4+1] = dbg_data[1]; 
            ssbo_out_segloopconv1d_data_[offset + idx*4+2] = dbg_data[2]; 
            ssbo_out_segloopconv1d_data_[offset + idx*4+3] = dbg_data[3]; 
        #endif

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
            cf.seg_head = 
				conv_temp_data.is_local_minima || conv_temp_data.is_local_maxima
				|| seg_head_id/*curve start, not "seg"*/ == idx// || cf.seg_head_contour 
				|| cf.seg_head_clipped; 
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