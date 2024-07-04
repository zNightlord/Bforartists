/**
 * Shared structures, enums & defines between C++ and GLSL.
 * Can also include some math functions but they need to be simple enough to be valid in both
 * language.
 */

#ifndef USE_GPU_SHADER_CREATE_INFO
#  pragma once

#include "BLI_math_bits.h"
#  include "BLI_memory_utils.hh"
#  include "DRW_gpu_wrapper.hh"

#  include "draw_manager.hh"
#  include "draw_pass.hh"

#  include "bnpr_defines.hh"

#  include "GPU_shader_shared.hh"

namespace blender::npr::strokegen {

using namespace draw;

#endif


  /* ----------------------------------------------------- */
  /** \name Noise Functions
   * \{ */
  static inline uint wang_hash(uint seed)
  {
    seed = (seed ^ 61) ^ (seed >> 16);
    seed *= 9;
    seed = seed ^ (seed >> 4);
    seed *= 0x27d4eb2d;
    seed = seed ^ (seed >> 15);
    return seed;
  }
  /** \} */


  /* ----------------------------------------------------- */
  /** \name Thread Group Counting
   * \{ */
  static inline uint compute_num_threads(
    uint numWorkItems,
    uint numItemsPerThread = 1u, uint numThreadsPerItem = 1u
  )
  {
    uint numThreads = numWorkItems;

    if (numItemsPerThread != 1u)
      numThreads = ((numWorkItems + numItemsPerThread - 1u) / numItemsPerThread);
    else if (numThreadsPerItem != 1u)
      numThreads = numThreadsPerItem * numWorkItems;

    return numThreads;
  }

  static inline uint compute_num_groups(
    uint numWorkItems, uint groupSize,
    uint numItemsPerThread = 1u, uint numThreadsPerItem = 1u
  )
  {
    const uint numThreads = compute_num_threads(
      numWorkItems, numItemsPerThread,
      numThreadsPerItem
    );

#if !defined(GPU_SHADER)
    return math::max(1u, ((numThreads + groupSize - 1) / groupSize));
#else
      return max(1u, ((numThreads + groupSize - 1) / groupSize));
#endif
  }
  /** \} */



  /* ----------------------------------------------------- */
  /** \name GPU Scan Testing
   * \{ */
  static inline uint tree_seg_scan_encode_upsweep_hfs(uint hf_partialSum, uint hf_orig)
  {
    return ((hf_orig << 1) | hf_partialSum);
  }

  static inline void tree_seg_scan_decode_upsweep_hfs(
#if !defined(GPU_SHADER)
    uint hfs_encoded, uint& out_hf_orig, uint& out_hf_partialSum
#else
    uint hfs_encoded, out uint out_hf_orig, out uint out_hf_partialSum
#endif
  )
  {
    out_hf_partialSum = (hfs_encoded & 1);
    hfs_encoded >>= 1;
    out_hf_orig = (hfs_encoded & 1);
  }

  static inline uint tree_seg_scan_decode_upsweep_hfs_get_origHF(
#if !defined(GPU_SHADER)
    uint hfs_encoded
#else
    uint hfs_encoded
#endif
  )
  {
    hfs_encoded >>= 1;
    return (hfs_encoded & 1);
  }

  static inline uint tree_seg_scan_decode_upsweep_hfs_get_sumHF(
#if !defined(GPU_SHADER)
    uint hfs_encoded
#else
    uint hfs_encoded
#endif
  )
  {
    return (hfs_encoded & 1);
  }


  struct UBData_TreeScan
  {
    uint num_scan_items;
    uint num_valid_scan_threads;
    uint num_thread_groups;
    uint dummy;
  };
  BLI_STATIC_ASSERT_ALIGN(UBData_TreeScan, 16)


  /* Segmented scan: pack hf and data together */
#ifndef GPU_SHADER
#  define uvec2 uint2
#  define uvec3 uint3
#  define uvec4 uint4
#  define vec2 float2
#  define vec3 float3
#  define vec4 float4
#endif
    
  static inline uint segscan_encode_hf_with_val(
    uint val, uint hfs/* can have 2 flags for per-block aggregates; this only happens for segscan */)
  {
    uint encoded = (val << 2u) | (hfs & 0x00000003u);
    return encoded;
  }
  static inline uint segscan_decode_hf(uint encoded)
  {
    return encoded & 0x3u; 
  }
  static inline uint segscan_decode_val(uint encoded)
  {
    return (encoded >> 2u);
  }

  // Scan data type, encoder/decoders are pre-defined here
  // use these to spawn different shader infos
  struct SSBOData_SegScanType_uint { // scan on uint
    uint val;
    uint hf;
  };
  static inline uint segscan_uint_hf_encode(SSBOData_SegScanType_uint data)
  {
    uint encoded = segscan_encode_hf_with_val(data.val, data.hf); 
    return encoded;
  }
  static inline SSBOData_SegScanType_uint segscan_uint_hf_decode(uint encoded)
  {
    SSBOData_SegScanType_uint data;
    data.val = segscan_decode_val(encoded);
    data.hf = segscan_decode_hf(encoded);
    return data;
  }

  struct SSBOData_SegScanType_uvec2 { // scan on uvec2
    uvec2 val;
    uint hf;
  };
  static inline uvec2 segscan_uvec2_hf_encode(SSBOData_SegScanType_uvec2 data)
  {
    uvec2 encoded;
    encoded.x = segscan_encode_hf_with_val(data.val.x, data.hf);
    encoded.y = data.val.y;
    return encoded;
  }
  static inline SSBOData_SegScanType_uvec2 segscan_uvec2_hf_decode(uvec2 encoded)
  {
    SSBOData_SegScanType_uvec2 data;
    data.val.x = segscan_decode_val(encoded.x);
    data.hf = segscan_decode_hf(encoded.x);
    data.val.y = encoded.y;
    return data;
  }

  struct SSBOData_SegScanType_uvec3 { // scan on uvec3
    uvec3 val;
    uint hf;
  };
  static inline uvec3 segscan_uvec3_hf_encode(SSBOData_SegScanType_uvec3 data)
  {
    uvec3 encoded;
    encoded.x = segscan_encode_hf_with_val(data.val.x, data.hf);
    encoded.y = data.val.y;
    encoded.z = data.val.z;
    return encoded;
  }
  static inline SSBOData_SegScanType_uvec3 segscan_uvec3_hf_decode(uvec3 encoded)
  {
    SSBOData_SegScanType_uvec3 data;
    data.val.x = segscan_decode_val(encoded.x);
    data.hf = segscan_decode_hf(encoded.x);
    data.val.y = encoded.y;
    data.val.z = encoded.z;
    return data;
  }

  struct SSBOData_SegScanType_uvec4 {  // scan on uvec4
    uvec4 val;
    uint hf;
  };
  static inline uvec4 segscan_uvec4_hf_encode(SSBOData_SegScanType_uvec4 data)
  {
    uvec4 encoded;
    encoded.x = segscan_encode_hf_with_val(data.val.x, data.hf);
    encoded.y = data.val.y;
    encoded.z = data.val.z;
    encoded.w = data.val.w; 
    return encoded;
  }
  static inline SSBOData_SegScanType_uvec4 segscan_uvec4_hf_decode(uvec4 encoded)
  {
    SSBOData_SegScanType_uvec4 data;
    data.val.x = segscan_decode_val(encoded.x);
    data.hf = segscan_decode_hf(encoded.x);
    data.val.y = encoded.y;
    data.val.z = encoded.z;
    data.val.w = encoded.w; 
    return data;
  }

  struct SSBOData_SegScanType_float {  // scan on float
    float val;
    uint hf;
  };
  static inline uvec2 segscan_float_hf_encode(SSBOData_SegScanType_float data)
  {
    uvec2 encoded;
#ifdef GPU_SHADER
    encoded.x = floatBitsToUint(data.val);
#else
    encoded.x = float_as_uint(data.val);
#endif

    encoded.y = data.hf;

    return encoded;
  }
  static inline SSBOData_SegScanType_float segscan_float_hf_decode(uvec2 encoded)
  {
    SSBOData_SegScanType_float data;
    #ifdef GPU_SHADER
      data.val = uintBitsToFloat(encoded.x);
    #else
      data.val = uint_as_float(encoded.x);
    #endif

    data.hf = encoded.y;

    return data;
  }

  struct SSBOData_SegScanType_vec3 {  // scan on vec2/vec3
    vec3 val;
    uint hf;
  };
  static inline uvec4 segscan_vec3_hf_encode(SSBOData_SegScanType_vec3 data)
  {
    uvec4 encoded;
#ifdef GPU_SHADER
    encoded.xyz = floatBitsToUint(data.val);
#else
    encoded.x = float_as_uint(data.val.x);
    encoded.y = float_as_uint(data.val.y);
    encoded.z = float_as_uint(data.val.z);
#endif

    encoded.w = data.hf;

    return encoded;
  }
  static inline SSBOData_SegScanType_vec3 segscan_vec3_hf_decode(uvec4 encoded)
  {
    SSBOData_SegScanType_vec3 data;
#ifdef GPU_SHADER
    data.val = uintBitsToFloat(encoded.xyz);
#else
    data.val.x = uint_as_float(encoded.x);
    data.val.y = uint_as_float(encoded.y);
    data.val.z = uint_as_float(encoded.z);
#endif

    data.hf = encoded.w;

    return data;
  }
#ifndef GPU_SHADER
  #undef uvec2
  #undef uvec3
  #undef uvec4

#  define SSBOData_SegScanTestEncoded uint4
#  define SSBOData_SegScanTestDecodeFunc segscan_uvec4_hf_decode
#  define SSBOData_SegScanTest SSBOData_SegScanType_uvec4
#  define SSBOData_SegScanZeroValue ((SSBOData_SegScanTest{uint4(0u, 0u, 0u, 0u), 1u}))
#endif
  /** \} */



  /* ----------------------------------------------------- */
  /** \name GPU 1D Segmented Looped Convolution Testing
   * \{ */
  struct UBData_SegLoopConv1D
  {
    uint num_thread_groups;
    uint num_conv_items;
    uint dummy2;
    uint dummy3;
  };
  BLI_STATIC_ASSERT_ALIGN(UBData_SegLoopConv1D, 16);
  #ifndef GPU_SHADER
    #define NPR_SEGLOOPCONV1D_TEST_DATA_TYPE float
  #endif
#if defined(GPU_SHADER) && defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)

#endif
  /** \} */



  /* ----------------------------------------------------- */
  /** \name GPU List Ranking Testing
   * \{ */
  struct UBData_ListRanking
  {
    uint num_thread_groups;
    uint num_nodes;
    uint subbuff_size;
    uint num_tagging_iters;
    uint dummy0;
    uint dummy1;
    uint dummy2;
    uint dummy4;
  };
  BLI_STATIC_ASSERT_ALIGN(UBData_ListRanking, 16);

  struct SSBOData_ListRankingInputs {
    uint num_nodes;
    uint dummy0;
    uint dummy1;
    uint dummy2;
  };


#define BNPR_LIST_RANKING_MAX_SUBLIST_LEN 5u
  static inline uint ComputeTaggingIters(uint numNodes)
  {
    uint iters = 0;
    uint l = numNodes; /* max #nodes between anchors after tagging for #iters */
    while (l > (BNPR_LIST_RANKING_MAX_SUBLIST_LEN >> 1u))
    /* apply tagging iters to make sub-list size <= BNPR_LIST_RANKING_MAX_SUBLIST_LEN */
    {
      /* http://www.graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2 */
      bool power_of_2 = (l != 0u) && (0u == (l & (l - 1u)));
#ifndef GPU_SHADER
      l = (31 - bitscan_reverse_uint(l)) + 1; /* (highest_set_bit) + 1 == #bits for a tag */
#else
      l = findMSB(l) < 0 ? 0u : (findMSB(l) + 1);
#endif
      if (!power_of_2) iters++; /* ceiling. */
    }
    return iters;
  }
  /** \} */



  /* -------------------------------------------------------------------- */
  /** \name Geometry Extraction from GPUBatch(es)
   * \{ */
  struct SSBOData_StrokeGenMeshPoolCounters
  {
    uint num_verts;
    uint num_edges;
    uint num_faces;
    uint num_contour_edges; // 4
    uint num_contour_verts;
    uint num_filtered_edges;
    uint num_filtered_verts;
    uint num_dbg_vnor_lines;  // 4
    uint num_dbg_vpdir_lines;  
    uint num_dbg_edge_lines;  
    uint num_draw_faces;  
    uint num_frags;  // 4
    uint num_2d_samples;
    uint dummy0;
    uint dummy1;
    uint dummy2; // 4 
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenMeshPoolCounters, 16)

  static inline uint get_prim_base_addr_ibo(bool is_ibo_fmt_u16, uint prim_id, uint num_verts_per_prim)
  {
    uint prim_vert_beg = prim_id * num_verts_per_prim;
    return (is_ibo_fmt_u16 ? (prim_vert_beg / 2u) : (prim_vert_beg));
  }

  struct SSBOData_StrokeGenEdgeSplitCounters {
    uint num_split_edges_pass_1;
    uint num_split_edges;
    uint dummy_0; 
    uint dummy_1;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenEdgeSplitCounters, 16)

  struct SSBOData_StrokeGenEdgeCollapseCounters {
    uint num_collapsed_edges_pass_1;
    uint dummy;
    uint dummy_0;
    uint dummy_1;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenEdgeCollapseCounters, 16)

  struct SSBOData_StrokeGenEdgeFlipCounters {
    uint num_flip_edges_pass_1;
    uint dummy;
    uint dummy_0;
    uint dummy_1;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenEdgeFlipCounters, 16)

  struct SSBOData_StrokeGenFaceSplitCounters {
    uint num_split_faces;
    uint dummy_0;
    uint dummy_1;
    uint dummy_2;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenFaceSplitCounters, 16)

  struct SSBOData_StrokeGenDynamicMeshCounters {
    uint num_verts;
    uint num_edges;
    uint dummy_0;
    uint dummy_1; 
  }; 
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenDynamicMeshCounters, 16)



  /* Addressing */
#define MAX_NUM_CONTOUR_EDGES ((6000000)) /* ssbo buffer cannot be larger than 8192*8192*/
  /* ! Reuse the whole buffer ! */
  /*x2*/
  static inline uint mesh_pool_addr__edgedir(uint contour_edge_id)
  { 
    uint base_addr = 0; 
    return base_addr + contour_edge_id * 2u;
  } 
  /*x4*/
  static inline uint mesh_pool_addr__zwhclip(uint contour_edge_id)
  {
    uint base_addr = mesh_pool_addr__edgedir(MAX_NUM_CONTOUR_EDGES);
    return base_addr + contour_edge_id * 4;
  }
  /*x4(+1debug)*/
  static inline uint mesh_pool_addr__edgeuv(uint contour_edge_id)
  {
    uint base_addr = mesh_pool_addr__zwhclip(MAX_NUM_CONTOUR_EDGES);
    return base_addr + contour_edge_id * 5/* 4 */;
  }
  /** } */

  
  /* Addressing - topology for filtered mesh edge/verts 
   * - data is temporary and only used in mesh filtering passes
   * - hence we reuse this buffer
  */
  static inline uint ssbo_edge_to_edges_addr__edge_to_selected_edge(uint edge_id, uint num_all_edges)
  { 
    return num_all_edges * 4u + edge_id;
  }
  static inline uint ssbo_vert_to_edge_list_header_addr__vert_to_selected_vert(
    uint vert_id, uint num_all_verts
  ){
    return ((num_all_verts + 3u) / 4u) * 4u + vert_id; 
  }
  /** } */


  /* Addressing - tex2d_contour_pix_marks_ */
#ifndef GPU_SHADER
#  define CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE (int2(8u, 4u))
#else
#  define CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE (uvec2(8u, 4u))
#endif
  /** } */

#ifdef GPU_SHADER
  #define GROUP_SIZE_BNPR_SCAN_SWEEP 1024u
  #define GROUP_SIZE_BNPR_SCAN_AGGRG 1024u  // for recursive-scan should be as large as possible
#endif


#ifdef __cplusplus

using UBO_ViewMatrices = draw::UniformBuffer<ViewMatrices>;
using SSBO_IndirectDrawArgs = draw::StorageBuffer<DrawCommand, true>;
using SSBO_IndirectDispatchArgs = draw::StorageBuffer<DispatchCommand>;


// Persistent Mesh Buffers ----------------
// about 64MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerEdge = draw::StorageArrayBuffer<T, MAX_NUM_EDGES_PER_BATCH * Stride, true>;
// about 32MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerSelectedEdge = draw::StorageArrayBuffer<T, MAX_NUM_EDGES_PER_BATCH * Stride / 2, true>;
// about 42MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerFace = draw::StorageArrayBuffer<T, MAX_NUM_EDGES_PER_BATCH * Stride, true>;
// about 32MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerVert = draw::StorageArrayBuffer<T, MAX_NUM_VERTS_PER_BATCH * Stride, true>;
// about 16MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerSelectedVert = draw::StorageArrayBuffer<T, MAX_NUM_VERTS_PER_BATCH * Stride / 2, true>;
// about 1MB for single stride
template<typename T, size_t Stride>
using SSBO_StrokeGenMeshBufPerContour = draw::StorageArrayBuffer<T, MAX_NUM_CONTOUR_EDGES_PER_BATCH * Stride, true>;


// Reused Buffers --------------------------
using SSBO_StrokeGenReusedLarge = StorageArrayBuffer<uint, 2048 * 2048 * 16, true>;
using SSBO_StrokeGenReusedMedium = StorageArrayBuffer<uint, 2048 * 2048 * 8, true>;
using SSBO_StrokeGenReusedSmall = StorageArrayBuffer<uint, 2048 * 2048 * 4, true>;
using SSBO_StrokeGenReusedTiny = StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_StrokeGenReusedMinimum = StorageArrayBuffer<uint, 2048 * 2048, true>;

using SSBO_StrokeGenMeshPoolCounters = StorageBuffer<SSBOData_StrokeGenMeshPoolCounters>;
using SSBO_StrokeGenEdgeSplitCounters = StorageArrayBuffer<SSBOData_StrokeGenEdgeSplitCounters, MAX_CONSEQ_EDGE_SPLITS + 1, true>;
using SSBO_StrokeGenEdgeCollapseCounters = StorageArrayBuffer<SSBOData_StrokeGenEdgeCollapseCounters, MAX_CONSEQ_EDGE_COLLAPSES + 1, true>;
using SSBO_StrokeGenEdgeFlipCounters = StorageArrayBuffer<SSBOData_StrokeGenEdgeFlipCounters, MAX_CONSEQ_EDGE_FLIPS + 1, true>;
using SSBO_StrokeGenFaceSplitCounters = StorageArrayBuffer<SSBOData_StrokeGenFaceSplitCounters, MAX_CONSEQ_FACE_SPLITS + 1, true>;
using SSBO_StrokeGenDynamicMeshCounters = StorageBuffer<SSBOData_StrokeGenDynamicMeshCounters>; 
using SSBO_StrokeGenTemporalRecordCounters =
  StorageArrayBuffer<uint, 4u * ((MAX_TEMPORAL_RECOREDED_FRAMES + 3u) / 4u), true>; 

// Buffers for testing parallel primitives --------------------
using SSBO_BnprScanData = StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_BnprScanAggregates = StorageArrayBuffer<uint, 512 * 16, true>;
using UBO_BnprTreeScan = UniformBuffer<UBData_TreeScan>;
using SSBO_BnprTreeScan = StorageBuffer<UBData_TreeScan>;

using SSBO_SegLoopConv1DData = StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_SegLoopConvDebugData = StorageArrayBuffer<uint, 2048 * 2048 * 8, true>;
using SSBO_SegLoopConvPatchTable = StorageArrayBuffer<uint, 4096 * 64, true>;
using UBO_SegLoopConv1D = UniformBuffer<UBData_SegLoopConv1D>;
using SSBO_SegLoopConv1DInfo = StorageBuffer<UBData_SegLoopConv1D>; 

using SSBO_ListRankingLinksStagingBuf = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST * 2, false/* We need to init this from CPU */>;
using SSBO_ListRankingLinks = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST * 2, true>;
using SSBO_ListRankingTags = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingRanks = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingSerializedTopo = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST * 2, true>; /* list-len&start */
using SSBO_ListRankingAnchorJumpingInfo = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST * 8, true>; /* also cache per-node data*/

/* #anchors is around #nodes/4 */
using SSBO_ListRankingAnchorToNode = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingAnchorToNextAnchor = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingSplicedNodeId = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingNodeToAnchor = draw::StorageArrayBuffer<uint, NUM_ITEMS_BNPR_LIST_RANK_TEST, true>;
using SSBO_ListRankingCounters = draw::StorageArrayBuffer<uint, BNPR_LIST_RANK_ANCHOR_COUNTER_BUFFER_SIZE>;
using SSBO_ListRankingAllocationCounters = draw::StorageArrayBuffer<uint, 4>;
using UBO_ListRanking = draw::UniformBuffer<UBData_ListRanking>;
using SSBO_ListRankingInputs = draw::StorageBuffer<SSBOData_ListRankingInputs>;

}

// namespace blender::bnpr
#endif
