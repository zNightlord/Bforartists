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

#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)

/* Generate segment topology for testing the right half of patch table.
* - that is, we let the seg tail cross the thread group border. */
void get_seg_topo_for_validate_patch_last_tail(uint gl_invoc_id, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
{ 
#define SECTION_LEN 231u
    uint section_id =  gl_invoc_id / SECTION_LEN;
    uint sec_head_id = section_id * SECTION_LEN;
    uint sec_tail_id = (sec_head_id + SECTION_LEN) - 1u;
    
    seg_len = wang_hash(section_id) % SECTION_LEN;
    uint id_seg_lc = (gl_invoc_id - sec_head_id) / seg_len;

    seg_head_id = id_seg_lc * seg_len + sec_head_id;
    seg_tail_id = (seg_head_id + seg_len) - 1u;
    
    if (sec_tail_id < seg_tail_id)
    {
        seg_tail_id = sec_tail_id;
        seg_len = seg_tail_id - seg_head_id + 1u;
    }
    
#undef SECTION_LEN
}


void main()
{
    const uint groupIdx = gl_LocalInvocationID.x;
    const uint idx = gl_GlobalInvocationID.x;
    const uint blockIdx = gl_WorkGroupID.x;
    
    uint seg_head_id, seg_tail_id, seg_len;
    get_seg_topo_for_validate_patch_last_tail(idx, seg_head_id, seg_tail_id, seg_len);
    
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
/**#define STORE_DBG_DATA_UINT(val) \
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (val); \
#endif*/
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_head_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_tail_id);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = (seg_len);
    ssbo_debug_segloopconv1d_data_[addr_dbg_st++] = ((hf ? 1 : (tf ? 2 : 0)));
}
#endif
