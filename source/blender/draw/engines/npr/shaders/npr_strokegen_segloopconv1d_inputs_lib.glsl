
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)


#ifndef BNPR_SEGLOOPCONV1D_LIB_INPUTS_INCLUDED
#define BNPR_SEGLOOPCONV1D_LIB_INPUTS_INCLUDED



uint get_num_items()
{
#if (_KERNEL_MULTI_COMPILE__SEGLOOPCONV1D_INFO_SSBO == 1)
    return uint(ssbo_segloopconv1d_info_.num_conv_items); 
#else
    return uint(ubo_segloopconv1d_.num_conv_items); 
#endif
}



/* Fetch topology --------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
    void FUNC_GET_LOOP_TOPOLOGY(uint elemId, out bool seg_is_loop, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
    {
        /* Generate segment topology for testing the right half of patch table.
            * - that is, we let the seg tail cross the thread group border. */
        #define SECTION_LEN 1731u
        seg_is_loop = true; 

        uint section_id =  elemId / SECTION_LEN;
        uint sec_head_id = section_id * SECTION_LEN;
        uint sec_tail_id = (sec_head_id + SECTION_LEN) - 1u;

        seg_len = max(1, wang_hash(section_id) % SECTION_LEN);
        uint id_seg_lc = (elemId - sec_head_id) / seg_len;

        seg_head_id = id_seg_lc * seg_len + sec_head_id;
        seg_tail_id = (seg_head_id + seg_len) - 1u;

        if (sec_tail_id < seg_tail_id)
            seg_tail_id = sec_tail_id;
        
        uint last_valid_elem_id = get_num_items() - 1u; 
        if (seg_head_id <= last_valid_elem_id && last_valid_elem_id <= seg_tail_id)
            seg_tail_id = last_valid_elem_id; 

        seg_len = seg_tail_id - seg_head_id + 1u;

        #undef SECTION_LEN
    }
#endif
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__SEG_DENOISING)
    void FUNC_GET_LOOP_TOPOLOGY(uint elem_id, out bool seg_is_loop, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
    {
        uint rank = ssbo_contour_snake_rank_[elem_id];
        seg_head_id = elem_id - rank; 

        uint len = ssbo_contour_snake_list_len_[elem_id];
        seg_len = len; 

        seg_tail_id = seg_head_id + len - 1u; 


        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)
            ContourFlags cfs = decode_contour_flags(ssbo_contour_snake_flags_[elem_id]);
        #endif
        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)
            uvec4 conv_input_enc = uvec4(0u);
            Load4(ssbo_in_segloopconv1d_data_, elem_id, conv_input_enc); 
            CuspSegmentDenoiseData csdd = decode_cusp_segment_denoise_data(conv_input_enc);

            ContourFlags cfs = csdd.cf; 
        #endif
        seg_is_loop = cfs.looped_curve; 
    }
#endif
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_PROCESSING)
    void FUNC_GET_LOOP_TOPOLOGY(uint elem_id, out bool seg_is_loop, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
    {
    	uint num_samples = ssbo_segloopconv1d_info_.num_conv_items; 

        ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(elem_id);
        ContourCurveTopo cct = load_contour_2d_sample_curve_topo(elem_id, cf, num_samples);

        seg_is_loop = cct.looped_curve;
        seg_head_id = cct.head_id;
        seg_tail_id = cct.tail_id;
        seg_len     = cct.len;
    }
#endif



/* Build patch table --------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)
void FUNC_DEVICE_STORE_LOOPCONV1D_PATCH_ID(uint stAddr, uint val)
{
    ssbo_segloopconv1d_patch_table_[stAddr] = val; // should be shared in most cases
}
#endif



/* Convolution --------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)
uint FUNC_DEVICE_LOAD_LOOPCONV1D_PATCH_ID(uint ldAddr)
{
    return ssbo_segloopconv1d_patch_table_[ldAddr]; // should be shared in most cases
}

bool should_init_conv_data(bool mov_left, uint mov_step) { return mov_left && (mov_step == 1u); }

#define T_CV DATA_TYPE_LOOPCONV1D
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
        struct T_CONV_TEMP_DATA
        {
            T_CV val; 
        }; 
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId)
        { /* Store is directly implemented in the shader code for more variety */
            return float(wang_hash(elemId * 73u) % 711u);
        }
        void FUNC_CONVOLUTION_LOOPCONV1D(
            bool mov_left, uint mov_step, bool seg_is_loop, uint seg_head_id, uint seg_tail_id, uint item_id, 
            T_CV neighbor_data, T_CV orig_data, inout T_CONV_TEMP_DATA conv_data)
        {
            if (should_init_conv_data(mov_left, mov_step)) conv_data.val = orig_data; // initialization
            conv_data.val += neighbor_data; // simple add
        }
    #endif
 
    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__SEG_DENOISING)
        struct T_CONV_TEMP_DATA
        {
            uint num_pstv_cusps_left;
            float pstv_seg_len_left; 
            float ngtv_seg_len_left; 
            uint num_pstv_cusps_right;  
            float pstv_seg_len_right;
            float ngtv_seg_len_right; 

            bool is_seg_head;

            vec3 vpos_prev;  
        }; 
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId) 
        {
            uvec4 encoded_data; 
            Load4(ssbo_in_segloopconv1d_data_, elemId, encoded_data); 
            return encoded_data; 
        }
        
        #define SEG_DENOISE_PROBE_DISTANCE 0.08f

        void FUNC_CONVOLUTION_LOOPCONV1D(
            bool mov_left, uint mov_step, bool seg_is_loop, uint seg_head_id, uint seg_tail_id, uint item_id, 
            T_CV neighbor_data, T_CV orig_data, inout T_CONV_TEMP_DATA conv_data)
        {
            CuspSegmentDenoiseData csdd_neigh = decode_cusp_segment_denoise_data(neighbor_data);
            if (1u == mov_step) 
            { // init/reset convolution data
                CuspSegmentDenoiseData csdd_ori = decode_cusp_segment_denoise_data(orig_data);
                
                conv_data.is_seg_head = csdd_ori.cf.seg_head; 
                if (mov_left) {
                    conv_data.num_pstv_cusps_left = 0u;
                    conv_data.pstv_seg_len_left = 0.0f; 
                    conv_data.ngtv_seg_len_left = 0.0f;
                } else {
                    conv_data.num_pstv_cusps_right = 0u;
                    conv_data.pstv_seg_len_right = 0.0f;
                    conv_data.ngtv_seg_len_right = 0.0f; 
                }
                
                conv_data.vpos_prev = csdd_ori.vpos_ws; 
            }

            bool moved_across_end = // exit when moving across the end of a non-looped curve
                (mov_left && (seg_head_id + mov_step > item_id)) 
                || (!mov_left && (seg_tail_id < item_id + mov_step));
            bool skip_when_not_looped = (!seg_is_loop) && moved_across_end; 
            if (skip_when_not_looped) return; 

            float probed_distance = mov_left ? // constrain the Euclidean radius of the convolution window
                conv_data.pstv_seg_len_left + conv_data.ngtv_seg_len_left 
                : conv_data.pstv_seg_len_right + conv_data.ngtv_seg_len_right;
            if (probed_distance > SEG_DENOISE_PROBE_DISTANCE) return;

            float edge_len = length(csdd_neigh.vpos_ws - conv_data.vpos_prev); 
            if (mov_left)
            {
                if (csdd_neigh.cf.cusp_func_pstv) {
                    conv_data.num_pstv_cusps_left += 1;
                    conv_data.pstv_seg_len_left += edge_len; 
                } else {
                    conv_data.ngtv_seg_len_left += edge_len;  
                }
            }
            else
            {
                if (csdd_neigh.cf.cusp_func_pstv) {
                    conv_data.num_pstv_cusps_right += 1;
                    conv_data.pstv_seg_len_right += edge_len; 
                } else {
                    conv_data.ngtv_seg_len_right += edge_len;  
                }
            }
            conv_data.vpos_prev = csdd_neigh.vpos_ws; 
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0)
        struct T_CONV_TEMP_DATA
        {
            float corner_curv; 
        }; 
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId) 
        {
			vec2 pix_coord = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(elemId); 
            return pix_coord;
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_1)
        struct T_CONV_TEMP_DATA
        {
            bool is_local_maxima; 
            bool is_local_minima; 
        };
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId) 
        {
            uint num_samples = get_num_items(); 
			float corner_curv = load_ssbo_contour_2d_sample_geometry__corner_curvature(elemId, num_samples); 
            return corner_curv;
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CALC_TANGENT_CURVATURE)
        struct T_CONV_TEMP_DATA
        {
            vec2 tangent;
        };
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId) 
        {
			vec2 pix_coord = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(elemId); 
            return pix_coord;
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2D_SAMPLE_VISIBILITY_DENOISING)
        struct T_CONV_TEMP_DATA
        {
            uint num_occ_samples_left;
            uint num_occ_samples_right;  
            bool filtered_occlusion; 
        };
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId) 
        {
			ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(elemId); 
            uint seg_key = (cf.occluded ? 1u : 0u); 
        #define SEG_KEY_AT_2D_SAMPLE_CLIPPED_END 0xffffffffu 
            if (cf.seg_head_clipped || cf.seg_tail_clipped) seg_key = SEG_KEY_AT_2D_SAMPLE_CLIPPED_END;
            return seg_key; 
        }
        void convolution_iteration(bool traverse_lt, T_CV seg_key, 
            inout T_CONV_TEMP_DATA conv_temp_data, 
			inout bool reach_clipped_seg_end, inout float conv_counter
        ){
            if (seg_key == SEG_KEY_AT_2D_SAMPLE_CLIPPED_END || reach_clipped_seg_end) 
            { // reached clipped point, do nothing & return
                reach_clipped_seg_end = true; 
                return; 
            }
            conv_counter++; 

            if (seg_key == 1u)
                if (traverse_lt) conv_temp_data.num_occ_samples_left++; 
			    else conv_temp_data.num_occ_samples_right++; 
        }
    #endif
#undef T_CV
#endif


#endif


