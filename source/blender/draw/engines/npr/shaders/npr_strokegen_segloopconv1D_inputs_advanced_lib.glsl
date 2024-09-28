
#pragma BLENDER_REQUIRE(npr_strokegen_segloopconv1D_lib.glsl)



#ifndef BNPR_STROKEGEN_SEGLOOPCONV1D_INPUT_ADVANCED_LIB
#define BNPR_STROKEGEN_SEGLOOPCONV1D_INPUT_ADVANCED_LIB


#if defined(SEGLOOPCONV1D_USE_ADVANCED_INPUT)

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)
        /* Advanced Convolution Function --------------------------------------------------- */
        void FUNC_CONVOLUTION_ADVANCED(
            uint idx, uint blockIdx, uint groupIdx,
            bool seg_is_loop, uint seg_head_id, uint seg_tail_id, uint seg_len,
            out T_CONV_TEMP_DATA conv_temp_data
        )
        {
        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
            DATA_TYPE_LOOPCONV1D orig_data = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                0, blockIdx, groupIdx,
                seg_len, seg_head_id
            );

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
        #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST


        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_PROCESSING)
            // "Robust corner detection using altitude to chord ratio accumulation"
            uint sample_id = idx;
            uint num_samples = ssbo_segloopconv1d_info_.num_conv_items;

            // TODO: we might need to avoid cusp-based contour segmentation for 2d samples
            uint sub_seg_rank = load_ssbo_contour_2d_sample_topology__seg_rank(sample_id, num_samples);
            uint sub_seg_len = load_ssbo_contour_2d_sample_topology__seg_len(sample_id, num_samples);
            ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id);
            bool is_sub_seg_loop = is_2d_sample_curve_looped(
                cf.looped_curve && cf.no_segmentation_on_contour_curve, 
                cf.curve_clipped, sub_seg_len == seg_len
            );

            bool short_seg = sub_seg_len < 42u; // (see the paper)

            #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0)
                conv_temp_data.corner_curv = .0f;
                vec2 pt = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                    0, blockIdx, groupIdx,
                    seg_len, seg_head_id
                );

                const uint L = 8; // 16u; // (in the paper they used 16 but it was too much for me)
                float C = .0f;
                uint filter_radius = L - 1u;
                // if (seg_len < L) continue; // no corner for short curves
                if (L <= sub_seg_len && (!short_seg))
                    for (uint d = 0; d < filter_radius; ++d)
                    {
                        bool valid_pt = true;
                        {
                            uint step_left = (filter_radius - d);
                            if (sub_seg_rank < step_left && !is_sub_seg_loop) valid_pt = false;
                            else
                            {
                                uint rank_chord_end = sub_seg_rank - step_left + L - 1u;
                                if (sub_seg_len <= rank_chord_end && !is_sub_seg_loop) valid_pt = false;
                            }
                        }
                        if (valid_pt)
                        {
                            vec2 pt_chord_beg = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                                filter_radius - d, blockIdx, groupIdx,
                                seg_len, seg_head_id
                            );

                            vec2 pt_chord_end = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
                                d + 1, blockIdx, groupIdx,
                                seg_len, seg_head_id
                            );

                            vec2 chord = pt_chord_end - pt_chord_beg;
                            vec2 a = pt - pt_chord_beg;
                            float dist = abs(chord.x * a.y - chord.y * a.x) / max(1e-10f, length(chord));
                            dist /= length(chord);

                            C += dist;
                        }
                    }

                conv_temp_data.corner_curv = C;

                if (sample_id < num_samples)
                {
                    vec4 dbg_col = vec4(C.xxx, 1.0f);
                    vec2 dbg_pix = pt;
                    // imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
                }
            #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0
            #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_1)
                conv_temp_data.is_local_maxima = false;
                conv_temp_data.is_local_minima = false;
                if (short_seg) return; // no corner for short sub-segs
                conv_temp_data.is_local_maxima = true;
                conv_temp_data.is_local_minima = true;

                float orig_curv = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                    0, blockIdx, groupIdx,
                    seg_len, seg_head_id
                );

                #define LOCAL_MAXIMA_RADIUS 2
                uint steps_left = is_sub_seg_loop ? LOCAL_MAXIMA_RADIUS : min(LOCAL_MAXIMA_RADIUS, sub_seg_rank);
                uint steps_right = is_sub_seg_loop ? LOCAL_MAXIMA_RADIUS : min(LOCAL_MAXIMA_RADIUS, sub_seg_len - sub_seg_rank - 1u);

                for (uint d = 1; d < steps_left; ++d)
                {
                    float neigh_curv = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                        d, blockIdx, groupIdx,
                        seg_len, seg_head_id
                    );

                    if (orig_curv < neigh_curv)
                        conv_temp_data.is_local_maxima = false;
                    else
                        conv_temp_data.is_local_minima = false;
                }
                for (uint d = 1; d < steps_right; ++d)
                {
                    float neigh_curv = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
                        d, blockIdx, groupIdx,
                        seg_len, seg_head_id
                    );

                    if (orig_curv < neigh_curv)
                        conv_temp_data.is_local_maxima = false;
                    else
                        conv_temp_data.is_local_minima = false;
                }

                if (orig_curv < 2.0f) // planar sections of the curve
                    conv_temp_data.is_local_maxima = false;
                if (2.0f <= orig_curv) // sharp point of the curve
                   conv_temp_data.is_local_minima = false;
				if (orig_curv < 0.3f)
					conv_temp_data.is_local_minima = true;

                if (sample_id < num_samples)
                {
                    vec4 dbg_col = vec4(1.0f);
                    if (conv_temp_data.is_local_maxima) dbg_col = vec4(1, 0, 0, 1);
                    if (conv_temp_data.is_local_minima) dbg_col = vec4(0, 1, 0, 1);
                    vec2 dbg_pix = pcs_screen_size_ * load_ssbo_contour_2d_sample_geometry__position(sample_id);
                    // imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
                }
            #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_1
            #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CALC_TANGENT_CURVATURE)
                vec2 pt = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                    0, blockIdx, groupIdx, seg_len, seg_head_id
                );

                const uint fit_radius = 8u;
                uint wls_window_sz = 2u * fit_radius + 1u;

                Curve2DWLSFitParams wls_params = curve_2d_wls_fit_init_params(); 
                float wls_l_l = .0f;
                float wls_l_r = .0f;
                for (uint d = 1; d <= fit_radius; ++d)
                {
                    vec2 pt_lt = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                        d, 
                        blockIdx, groupIdx, seg_len, seg_head_id
                    );
                    bool reach_left_seg_end = (sub_seg_rank < d && !is_sub_seg_loop); 

                    vec2 pt_rt = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
                        d, 
                        blockIdx, groupIdx, seg_len, seg_head_id
                    );
                    bool reach_right_seg_end = (sub_seg_len <= sub_seg_rank + d && !is_sub_seg_loop); 

                    // if (reach_left_seg_end || reach_right_seg_end) break; 

                    curve_2d_wls_fit_iter(wls_window_sz, true,  pt, pt_lt, /*inout*/wls_l_l, wls_params);
                    curve_2d_wls_fit_iter(wls_window_sz, false, pt, pt_rt, /*inout*/wls_l_r, wls_params);
                }

                vec2 tangent; float curv; 
                curve_2d_wls_fit_solve(wls_params, /*out*/tangent, curv);

                conv_temp_data.tangent = tangent; 

                if (sample_id < num_samples)
                {
                    vec4 dbg_col = vec4(tangent.xy * .5f + .5f, abs(curv), 1.0f);
                    vec2 dbg_pix = pt;
                    // imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
                }
            #endif
            #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2D_SAMPLE_VISIBILITY_DENOISING)
                conv_temp_data.num_occ_samples_left = 0u; 
                conv_temp_data.num_occ_samples_right = 0u; 
                uint occ = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                    0, blockIdx, groupIdx, seg_len, seg_head_id
                ); // == cf.occluded

                const uint filter_radius = 8u;
                bool lt_reach_clipped_seg_end = false;
                float num_lt_samples = 0; 
                bool rt_reach_clipped_seg_end = false;  
                float num_rt_samples = 0; 
                for (uint d = 1; d <= filter_radius; ++d)
                {
                    uint occ_lt = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                        d, blockIdx, groupIdx, seg_len, seg_head_id
                    );
                    convolution_iteration(true,  occ_lt, /*inout*/conv_temp_data, lt_reach_clipped_seg_end, num_lt_samples); 

                    uint occ_rt = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
                        d, blockIdx, groupIdx, seg_len, seg_head_id
                    );
                    convolution_iteration(false, occ_rt, /*inout*/conv_temp_data, rt_reach_clipped_seg_end, num_rt_samples); 
                }

                float num_occ_samples_left = float(conv_temp_data.num_occ_samples_left);
                float num_occ_samples_right = float(conv_temp_data.num_occ_samples_right);
                conv_temp_data.filtered_occlusion = 
                    (num_occ_samples_left + num_occ_samples_right) > .3f * (num_lt_samples + num_rt_samples); 

                if (sample_id < num_samples)
                {
                    vec4 dbg_col = vec4(
                        conv_temp_data.num_occ_samples_left, 
                        conv_temp_data.num_occ_samples_right, 
                        !conv_temp_data.filtered_occlusion, 
                        1.0f
                    );
                    vec2 dbg_pix = pcs_screen_size_ * load_ssbo_contour_2d_sample_geometry__position(sample_id);
                    // imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
                }
            #endif
        #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_PROCESSING
        }



    #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION



#endif // SEGLOOPCONV1D_USE_ADVANCED_INPUT


#endif
