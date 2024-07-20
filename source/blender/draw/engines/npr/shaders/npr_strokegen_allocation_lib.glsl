

#ifndef BNPR_ALLOCATION__INCLUDED
#define BNPR_ALLOCATION__INCLUDED


/* Parallel Allocation --------------------------------------------------------------------------------------------- */
#define DECL_LDS_COUNTERS_PER_LANE(tag) \
    shared uint CAT(LDS_counters_per_lane_, tag)[32u]; \

#define DECL_LDS_COUNTERS_PER_LANE_SUM(tag) \
    shared uint CAT(LDS_counters_per_lane_sum_, tag)[32u]; \

#define DECL_LDS_ALLOC_BLK_COUNTER(tag) \
    shared uint CAT(LDS_alloc_blk_counter_, tag); \

#define DECL_LDS_ALLOC_BLOCK_OFFSET(tag) \
    shared uint CAT(LDS_alloc_blk_offset_, tag); \

#define DECL_ALLOCATION_FUNC(tag, GLOBAL_ATOMIC_COUNTER) \
uint CAT(alloc_, tag)(uint groupIdx, uint num_items)                                               \
{ /* allocation */                                                                                 \
    const uint wave_id = groupIdx >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */ \
    const uint lane_id = groupIdx % 32u;                                                           \
    const uint num_waves = gl_WorkGroupSize.x >> 5u;                                               \
                                                                                                   \
    /* Clear LDS counters */                                                                       \
    if (wave_id == 0u)                                                                             \
    {                                                                                              \
        LDS_COUNTERS_PER_LANE[lane_id] = 0u;                                                       \
        LDS_COUNTERS_PER_LANE_SUM[lane_id] = 0u;                                                   \
        if (lane_id == 0u)                                                                         \
            LDS_ALLOC_BLK_COUNTER = 0u;                                                            \
    }                                                                                              \
    barrier();                                                                                     \
                                                                                                   \
    uint lane_offset = atomicAdd(LDS_COUNTERS_PER_LANE[lane_id], num_items);                       \
    barrier();                                                                                     \
                                                                                                   \
    if (wave_id == 0u)                                                                             \
        LDS_COUNTERS_PER_LANE_SUM[lane_id] =                                                       \
            atomicAdd(LDS_ALLOC_BLK_COUNTER, LDS_COUNTERS_PER_LANE[lane_id]);                      \
    barrier();                                                                                     \
    memoryBarrier();                                                                               \
                                                                                                   \
    uint lane_sum = LDS_COUNTERS_PER_LANE_SUM[lane_id];                                            \
    uint group_offset = lane_offset + lane_sum;                                                    \
                                                                                                   \
    /* Add block sum to global counter. */                                                         \
    if (groupIdx == gl_WorkGroupSize.x - 1u)                                                       \
    {                                                                                              \
        LDS_ALLOC_BLOCK_OFFSET = atomicAdd(                                                        \
            (GLOBAL_ATOMIC_COUNTER),                                                               \
            LDS_ALLOC_BLK_COUNTER                                                                  \
        );                                                                                         \
    }                                                                                              \
    memoryBarrier();                                                                               \
    barrier();                                                                                     \
                                                                                                   \
    uint global_offset = group_offset + LDS_ALLOC_BLOCK_OFFSET;                                    \
                                                                                                   \
    return global_offset;                                                                          \
}                                                                                                  \






#if defined(_KERNEL_MULTICOMPILE__FACE_SPLIT_WORK_GEN)
#define AC_TAG split_face

DECL_LDS_COUNTERS_PER_LANE(AC_TAG)
#define LDS_COUNTERS_PER_LANE CAT(LDS_counters_per_lane_, AC_TAG)
DECL_LDS_COUNTERS_PER_LANE_SUM(AC_TAG)
#define LDS_COUNTERS_PER_LANE_SUM CAT(LDS_counters_per_lane_sum_, AC_TAG)
DECL_LDS_ALLOC_BLK_COUNTER(AC_TAG)
#define LDS_ALLOC_BLK_COUNTER CAT(LDS_alloc_blk_counter_, AC_TAG)
DECL_LDS_ALLOC_BLOCK_OFFSET(AC_TAG)
#define LDS_ALLOC_BLOCK_OFFSET CAT(LDS_alloc_blk_offset_, AC_TAG)

DECL_ALLOCATION_FUNC(AC_TAG, ssbo_face_split_counters_[pcs_split_iter_].num_split_faces)
#undef AC_TAG
#endif



#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT)
#define AC_TAG draw_face

DECL_LDS_COUNTERS_PER_LANE(AC_TAG)
#define LDS_COUNTERS_PER_LANE CAT(LDS_counters_per_lane_, AC_TAG)
DECL_LDS_COUNTERS_PER_LANE_SUM(AC_TAG)
#define LDS_COUNTERS_PER_LANE_SUM CAT(LDS_counters_per_lane_sum_, AC_TAG)
DECL_LDS_ALLOC_BLK_COUNTER(AC_TAG)
#define LDS_ALLOC_BLK_COUNTER CAT(LDS_alloc_blk_counter_, AC_TAG)
DECL_LDS_ALLOC_BLOCK_OFFSET(AC_TAG)
#define LDS_ALLOC_BLOCK_OFFSET CAT(LDS_alloc_blk_offset_, AC_TAG)

DECL_ALLOCATION_FUNC(AC_TAG, ssbo_bnpr_mesh_pool_counters_.num_draw_faces)
#undef AC_TAG
#endif



#if defined(_KERNEL_MULTICOMPILE__EXTRACT_MESH_CONTOUR_DATA)
#define AC_TAG raster_frags

DECL_LDS_COUNTERS_PER_LANE(AC_TAG)
#define LDS_COUNTERS_PER_LANE CAT(LDS_counters_per_lane_, AC_TAG)
DECL_LDS_COUNTERS_PER_LANE_SUM(AC_TAG)
#define LDS_COUNTERS_PER_LANE_SUM CAT(LDS_counters_per_lane_sum_, AC_TAG)
DECL_LDS_ALLOC_BLK_COUNTER(AC_TAG)
#define LDS_ALLOC_BLK_COUNTER CAT(LDS_alloc_blk_counter_, AC_TAG)
DECL_LDS_ALLOC_BLOCK_OFFSET(AC_TAG)
#define LDS_ALLOC_BLOCK_OFFSET CAT(LDS_alloc_blk_offset_, AC_TAG)

DECL_ALLOCATION_FUNC(AC_TAG, ssbo_bnpr_mesh_pool_counters_.num_frags)
#undef AC_TAG
#endif



#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
#define AC_TAG general_dbg_lines

DECL_LDS_COUNTERS_PER_LANE(AC_TAG)
#define LDS_COUNTERS_PER_LANE CAT(LDS_counters_per_lane_, AC_TAG)
DECL_LDS_COUNTERS_PER_LANE_SUM(AC_TAG)
#define LDS_COUNTERS_PER_LANE_SUM CAT(LDS_counters_per_lane_sum_, AC_TAG)
DECL_LDS_ALLOC_BLK_COUNTER(AC_TAG)
#define LDS_ALLOC_BLK_COUNTER CAT(LDS_alloc_blk_counter_, AC_TAG)
DECL_LDS_ALLOC_BLOCK_OFFSET(AC_TAG)
#define LDS_ALLOC_BLOCK_OFFSET CAT(LDS_alloc_blk_offset_, AC_TAG)

DECL_ALLOCATION_FUNC(AC_TAG, ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines)
#undef AC_TAG
#endif



#endif




