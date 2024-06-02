
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

        
        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION)
            uint sample_id = idx; 
            uint num_samples = ssbo_segloopconv1d_info_.num_conv_items; 

            // TODO: we might need to avoid cusp-based contour segmentation for 2d samples
            uint sub_seg_rank = load_ssbo_contour_2d_sample_topology__seg_rank(sample_id, num_samples); 
            uint sub_seg_len = load_ssbo_contour_2d_sample_topology__seg_len(sample_id, num_samples); 
            ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id); 
            bool is_sub_seg_loop = cf.looped_curve && (!cf.curve_clipped) && (sub_seg_len == seg_len); 

            conv_temp_data.corner_scores = vec3(.0f); 
            vec2 pt = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                0, blockIdx, groupIdx,
                seg_len, seg_head_id
            ); 

            const uvec3 chord_lens = uvec3(10u, 20u, 30u); 
            vec3 C = vec3(.0f, .0f, .0f); 

            for (uint ic = 0u; ic < 3u; ++ic)
            {
                C[ic] = .0f; 
                
                uint chord_len = chord_lens[ic]; 
                uint filter_radius = chord_len - 1u;
                // if (seg_len < chord_len) continue; // no corner for short curves
                if (sub_seg_len < chord_len) continue; // no corner for short sub-segs

                for (uint d = 0; d < filter_radius; ++d)
                {
                    bool valid_pt = true; 
                    { 
                        uint step_left = (filter_radius - d); 
                        if (sub_seg_rank < step_left && !is_sub_seg_loop) valid_pt = false; 
                        else
                        {
                            uint rank_chord_end = sub_seg_rank - step_left + chord_len - 1u; 
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
                        vec2 a_proj = chord * (dot(a, chord) / dot(chord, chord)); 
                        float dist = sqrt(max(.0f, dot(a, a) - dot(a_proj, a_proj))); 
                        
                        C[ic] += dist; 
                    }
                }
            }
            
            float c_max = max(C[0], max(C[1], C[2])); 
            C /= max(c_max, 1e-10f); 
            conv_temp_data.corner_scores = C; 

            if (sample_id < num_samples)
            {
                vec4 dbg_col = vec4(1.0f); 
                vec2 dbg_pix = pt; 
                imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
            }
        #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION
        }



    #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION



#endif // SEGLOOPCONV1D_USE_ADVANCED_INPUT


#endif
