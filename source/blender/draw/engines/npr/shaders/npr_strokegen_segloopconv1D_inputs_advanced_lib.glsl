
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

        
        #if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST)
            vec2 orig_data = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                0, blockIdx, groupIdx,
                seg_len, seg_head_id
            ); 

            for (uint d = 0; d < MAX_CONV_RADIUS; ++d)
            {
                vec2 arc_beg_pos = _FUNC_LOAD_CONV_DATA_LDS_LEFT(
                    MAX_CONV_RADIUS - d, blockIdx, groupIdx,
                    seg_len, seg_head_id 
                );

                vec2 arc_end_pos = _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
                    d + 1, blockIdx, groupIdx,
                    seg_len, seg_head_id
                ); 
                 
            }
        #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST
        }



    #endif // _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION



#endif // SEGLOOPCONV1D_USE_ADVANCED_INPUT


#endif
