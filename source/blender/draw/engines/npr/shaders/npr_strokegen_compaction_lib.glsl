#ifndef BNPR_COMPACTION__INCLUDED
#define BNPR_COMPACTION__INCLUDED

/* Macro expansion, for details, see
/* ---------------------------------------
/* https://stackoverflow.com/questions/1489932/how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg */
#define CAT(x, y) CAT_(x, y)
#define CAT_(x, y) x ## y

/* Parallel Compaction --------------------------------------------------------------------------------------------- */
/* !!! Must invoke this for EVERY thread in the work group !!! 
 * not matter if your thread maps to "invalid" element(s) */
/* Code gen Input: */
/* #define COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN --- wether use the default compaction code generator */
/* #define GLOBAL_COUNTER XXX --- Global Compaction Counter, should be cleared to 0*/
/* #define CP_TAG XXX --- Tag, defined in shader info for default code gen */
#define DECL_LDS_DIGIT_PER_LANE(tag) \
    shared uint CAT(LDS_digit_per_lane_, tag)[32u]; \

#define DECL_LDS_OFFSET_PER_LANE_SLOT(tag) \
    shared uint CAT(LDS_offset_per_lane_slot_, tag)[32u]; \

#define DECL_LDS_HIST_BLK(tag) \
    shared uint CAT(LDS_hist_blk_, tag); \

#define DECL_LDS_SCAN_BLOCK_OFFSET(tag) \
    shared uint CAT(LDS_scan_block_offset_, tag); \


#define DECL_COMPACTION_FUNC(tag, GLOBAL_ATOMIC_COUNTER) \
uint CAT(compact_, tag)(bool val, uint groupIdx, uint multiplier = 1u)                             \
{                                                                                                  \
	const uint wave_id = groupIdx >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */ \
    const uint lane_id = groupIdx % 32u;                                                           \
    const uint num_waves = gl_WorkGroupSize.x >> 5u;                                               \
                                                                                                   \
	/* Clear LDS counters */                                                                       \
    if (wave_id == 0u)                                                                             \
    {                                                                                              \
        if (lane_id == 0u) LDS_HIST_BLK = 0u;                                                      \
        LDS_DIGIT_PER_LANE[lane_id] = 0u;                                                          \
    }                                                                                              \
    barrier();                                                                                     \
	/*  w0    w1    w2    w3         LDS_DIGIT_PER_LANE                                            \
	 * 00:1  04:1  08:0  12:1  l0    0000                                                          \
	 * 01:0  05:0  09:0  13:0  l1    0000                                                          \
	 * 02:1  06:1  10:0  14:0  l2    0000                                                          \
	 * 03:1  07:0  11:1  15:0  l3    0000                                                          \
	 * ---------------------------                                                                 \
	 * LDS_HIST_BLK                                                                                \
	 *  0                                                                                          \
	*/                                                                                             \
                                                                                                   \
    /* Mark 1/0 at bit #wave_id */                                                                 \
    uint compact_bitval = (val ? 1u : 0u);                                                         \
    uint lds_compact_input = compact_bitval << wave_id;                                            \
    atomicOr(LDS_DIGIT_PER_LANE[lane_id], lds_compact_input);                                      \
    barrier();                                                                                     \
	/*  w0    w1    w2    w3         LDS_DIGIT_PER_LANE                                            \
	 * 00:1  04:1  08:0  12:1  l0    1011                                                          \
	 * 01:0  05:0  09:0  13:0  l1    0000                                                          \
	 * 02:1  06:1  10:0  14:0  l2    0011                                                          \
	 * 03:1  07:0  11:1  15:0  l3    0101                                                          \
	*/                                                                                             \
                                                                                                   \
    /* Prefix sum on lane sums */                                                                  \
    uint lane_digit = LDS_DIGIT_PER_LANE[lane_id];                                                 \
    uint wave_mask = (~(0xffffffffu << wave_id));                                                  \
	/*  w0    w1    w2    w3                                                                       \
	 * 0000  0001  0011  0111        wave_mask                                                     \
	 *                                                                                             \
	 *  w0    w1    w2    w3         lane_digit                                                    \
	 * 00:1  04:1  08:0  12:1  l0    1011                                                          \
	 * 01:0  05:0  09:0  13:0  l1    0000                                                          \
	 * 02:1  06:1  10:0  14:0  l2    0011                                                          \
	 * 03:1  07:0  11:1  15:0  l3    0101                                                          \
	*/                                                                                             \
    uint lane_digit_masked = lane_digit & wave_mask;                                               \
    uint num_1_bits_low = bitCount(lane_digit_masked);                                             \
    uint lane_offset = num_1_bits_low;                                                             \
    barrier();                                                                                     \
                                                                                                   \
    if (wave_id == 0u)                                                                             \
    {                                                                                              \
        uint lane_digit_sum = bitCount(lane_digit); /* digit sum */                                \
        LDS_OFFSET_PER_LANE_SLOT[lane_id] = atomicAdd(LDS_HIST_BLK, lane_digit_sum);               \
    }                                                                                              \
    barrier();                                                                                     \
                                                                                                   \
    memoryBarrier();                                                                               \
    /* Add block sum to global counter. */                                                         \
    if (groupIdx == gl_WorkGroupSize.x - 1u)                                                       \
    {                                                                                              \
        LDS_SCAN_BLOCK_OFFSET = atomicAdd(                                                         \
            (GLOBAL_ATOMIC_COUNTER),                                                               \
            LDS_HIST_BLK * (multiplier)                                                            \
        );                                                                                         \
    }                                                                                              \
    memoryBarrier();                                                                               \
    barrier();                                                                                     \
                                                                                                   \
    /* Compute final offset */                                                                     \
    uint local_offset = LDS_OFFSET_PER_LANE_SLOT[lane_id] + lane_offset;                           \
    local_offset *= (multiplier);                                                                  \
    uint blk_offset = LDS_SCAN_BLOCK_OFFSET;                                                       \
                                                                                                   \
    uint scanres = local_offset + blk_offset;                                                      \
	return scanres;                                                                                \
}                                                                                                  \


/* ----------------------------------------------------------------------------
 * Default code generation, 
 * for those kernels only use one simple compaction op 
*/
#if !defined(COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN) 

DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, GLOBAL_COUNTER)

    #undef LDS_HIST_BLK
    #undef LDS_OFFSET_PER_LANE_SLOT
    #undef LDS_DIGIT_PER_LANE
#undef CP_TAG

#endif








/* /////////////////////////////////////////////////////////////////////////
 * Custom code generation 
*/

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0)
// nothing special here, just to demonstrate how to use the custom code gen
#define CP_TAG normal_line
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_dbg_vnor_lines)
#undef CP_TAG
#endif



#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1) || defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
#define CP_TAG general_dbg_lines
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines)
#undef CP_TAG
#endif




#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT)
#define CP_TAG contour_edge
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_contour_edges)
    #undef LDS_HIST_BLK
    #undef LDS_OFFSET_PER_LANE_SLOT
    #undef LDS_DIGIT_PER_LANE
#undef CP_TAG

#define CP_TAG dbg_edge
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_dbg_edge_lines)
    #undef LDS_HIST_BLK
    #undef LDS_OFFSET_PER_LANE_SLOT
    #undef LDS_DIGIT_PER_LANE
#undef CP_TAG

#define CP_TAG draw_faces
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_draw_faces)
    #undef LDS_HIST_BLK
    #undef LDS_OFFSET_PER_LANE_SLOT
    #undef LDS_DIGIT_PER_LANE
#undef CP_TAG

#endif /* _KERNEL_MULTICOMPILE__GEOM_EXTRACT */



#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_1)
#define CP_TAG visibility_contour_split
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, ssbo_bnpr_mesh_pool_counters_.num_contour_edges)
#undef CP_TAG
#endif



#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__COMPACT)
#define CP_TAG temporal_contour_record
    DECL_LDS_DIGIT_PER_LANE(CP_TAG)
    #define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_OFFSET_PER_LANE_SLOT(CP_TAG)
    #define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_digit_per_lane_, CP_TAG)
    DECL_LDS_HIST_BLK(CP_TAG)
    #define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
    DECL_LDS_SCAN_BLOCK_OFFSET(CP_TAG)
    #define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG) 

    DECL_COMPACTION_FUNC(CP_TAG, (ssbo_temporal_record_counters_[(pc_obj_id_ * MAX_TEMPORAL_FRAMES) + (pc_frame_id_ % MAX_TEMPORAL_FRAMES)]))
#undef CP_TAG
#endif



#endif




