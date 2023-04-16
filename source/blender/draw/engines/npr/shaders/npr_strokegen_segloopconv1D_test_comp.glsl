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

/* Generate segment topology for testing. */
void get_test_seg_topo(uint gl_invoc_id, out uint seg_head_id, out uint seg_tail_id, out uint seg_len)
{
#define SECTION_LEN 2173u
    uint section_id =  gl_invoc_id / SECTION_LEN;
    uint sec_head_id = section_id * SECTION_LEN;
    uint sec_tail_id = (sec_head_id + SECTION_LEN) - 1u;
    
    uint id_seg_lc = (gl_invoc_id - section_id) / SECTION_LEN;
    seg_len = wang_hash(section_id) % SECTION_LEN;
    
    seg_head_id = id_seg_lc * seg_len + sec_head_id;
    seg_tail_id = min(sec_tail_id, (seg_head_id + seg_len) - 1u);
    
#undef SECTION_LEN
}


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
        
    _FUNC_COMPUTE_PATCH_TABLE(
        blockIdx, groupIdx,
        hf, tf, seg_head_id, seg_tail_id, seg_len 
    );
}
#endif
