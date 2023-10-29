/* SPDX-License-Identifier: GPL-2.0-or-later */

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

#  include "GPU_shader_shared.h"

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

  struct SSBOData_SegScanTest
  {
#ifndef GPU_SHADER
    uint3 val;
#else
    uvec3 val;
#endif
    uint hf;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_SegScanTest, 16)
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

#define NPR_SEGLOOPCONV1D_TEST_DATA_TYPE uint
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
    uint num_contour_edges;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenMeshPoolCounters, 16)

  static inline uint get_prim_base_addr_ibo(bool is_ibo_fmt_u16, uint prim_id, uint num_verts_per_prim)
  {
    uint prim_vert_beg = prim_id * num_verts_per_prim;
    return (is_ibo_fmt_u16 ? (prim_vert_beg / 2u) : (prim_vert_beg));
  }

  /* The posnor vbo of a Surface GPUBatch
   * See "PosNorLoop" in mesh_extractors\extract_mesh_vbo_pos_nor.cc
   */
  struct SSBO_Data_PosNor {
    float x;
    float y;
    float z;
    uint packedNormal; 
  };
  BLI_STATIC_ASSERT_ALIGN(SSBO_Data_PosNor, 16)


  /* Addressing - buf_strokegen_mesh_pool */
  /*x6*/
  static inline uint mesh_pool_addr__wpos(uint contour_edge_id) 
  {
    return contour_edge_id * 6u;
  }
  /*x2*/
  static inline uint mesh_pool_addr__edgedir(uint contour_edge_id, uint num_contour_edges)
  {
    uint base_addr = mesh_pool_addr__wpos(num_contour_edges);
    return base_addr + contour_edge_id * 2u;
  }
  /*x4*/
  static inline uint mesh_pool_addr__zwhclip(uint contour_edge_id, uint num_contour_edges)
  {
    uint base_addr = mesh_pool_addr__edgedir(num_contour_edges, num_contour_edges);
    return base_addr + contour_edge_id * 4;
  }
  /*x4*/
  static inline uint mesh_pool_addr__edgeuv(uint contour_edge_id, uint num_contour_edges)
  {
    uint base_addr = mesh_pool_addr__zwhclip(num_contour_edges, num_contour_edges);
    return base_addr + contour_edge_id * 4;
  }
  /** } */

  /* Addressing - tex2d_contour_pix_marks_ */
#ifndef GPU_SHADER
#  define CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE (int2(8u, 4u))
#else
#  define CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE (uvec2(8u, 4u))
#endif
  /** } */


#ifdef __cplusplus

// Template to set buffer size in compile time
using UBO_ViewMatrices = draw::UniformBuffer<ViewMatrices>;
using SSBO_IndirectDrawArgs = draw::StorageBuffer<DrawCommand, true>;
using SSBO_IndirectDispatchArgs = draw::StorageBuffer<DispatchCommand>;

using SSBO_StrokeGenTest = draw::StorageArrayBuffer<uint, 4096 * 4, true>; 

using SSBO_StrokeGenMeshVertHashTable = draw::StorageArrayBuffer<uint, MAX_VERT_HASH_TABLE_SIZE, true>;
using SSBO_StrokeGenMeshLarge = draw::StorageArrayBuffer<uint, 2048 * 2048 * 16, true>;
using SSBO_StrokeGenMeshSmall = draw::StorageArrayBuffer<uint, 2048 * 2048 * 4, true>;
using SSBO_StrokeGenMeshTiny = draw::StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_StrokeGenMeshLarge_Float = draw::StorageArrayBuffer<float, 2048 * 2048 * 12, true>;
using SSBO_StrokeGenMeshSmall_Float = draw::StorageArrayBuffer<float, 2048 * 2048 * 4, true>;

using SSBO_StrokeGenMeshPoolCounters = draw::StorageBuffer<SSBOData_StrokeGenMeshPoolCounters>;

using SSBO_BnprScanData = draw::StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_BnprScanAggregates = draw::StorageArrayBuffer<uint, 512 * 16, true>;
using UBO_BnprTreeScan = draw::UniformBuffer<UBData_TreeScan>;

using SSBO_SegLoopConv1DData = draw::StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_SegLoopConvDebugData = draw::StorageArrayBuffer<uint, 2048 * 2048 * 8, true>;
using SSBO_SegLoopConvPatchTable = draw::StorageArrayBuffer<uint, 4096 * 64, true>;
using UBO_SegLoopConv1D = draw::UniformBuffer<UBData_SegLoopConv1D>;

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

}

// namespace blender::bnpr
#endif
