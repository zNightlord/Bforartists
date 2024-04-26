
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)


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

        uint len = ssbo_contour_edge_list_len_[elem_id];
        seg_len = len; 

        seg_tail_id = seg_head_id + len - 1u; 

        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)
            ContourFlags cfs = decode_contour_flags(ssbo_contour_edge_flags_[elem_id]);
        #endif
        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)
            ContourFlags cfs = decode_contour_flags(ssbo_in_segloopconv1d_data_[elem_id]);
        #endif
        seg_is_loop = cfs.looped_curve; 
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
            T_CV neighbor_data, T_CV ori_data, inout T_CONV_TEMP_DATA conv_data)
        {
            if (should_init_conv_data(mov_left, mov_step)) conv_data.val = ori_data; // initialization
            conv_data.val += neighbor_data; // simple add
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__SEG_DENOISING)
        struct T_CONV_TEMP_DATA
        {
            uint num_neg_cusp_funcs; 
        }; 
        T_CV FUNC_DEVICE_LOAD_LOOPCONV1D_DATA(uint elemId)
        {
            return ssbo_in_segloopconv1d_data_[elemId];
        }
        void FUNC_CONVOLUTION_LOOPCONV1D(
            bool mov_left, uint mov_step, bool seg_is_loop, uint seg_head_id, uint seg_tail_id, uint item_id, 
            T_CV neighbor_data, T_CV ori_data, inout T_CONV_TEMP_DATA conv_data)
        {
            ContourFlags efs_neigh = decode_contour_flags(neighbor_data);
            ContourFlags efs_ori = decode_contour_flags(ori_data);
            
            conv_data.num_neg_cusp_funcs = 0u; 
            // efs_conv.seg_head
        }
    #endif
#undef T_CV
#endif


#endif


