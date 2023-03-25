/* SPDX-License-Identifier: GPL-2.0-or-later */

/**
 * Shared structures, enums & defines between C++ and GLSL.
 * Can also include some math functions but they need to be simple enough to be valid in both
 * language.
 */

#ifndef USE_GPU_SHADER_CREATE_INFO
#  pragma once

#  include "BLI_memory_utils.hh"
#  include "DRW_gpu_wrapper.hh"

#  include "draw_manager.hh"
#  include "draw_pass.hh"

#  include "bnpr_defines.hh"

#  include "GPU_shader_shared.h"

namespace blender::bnpr {

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


  
  /* -------------------------------------------------------------------- */
  /** \name Geometry Extraction from GPUBatch(es)
   * \{ */
  struct SSBOData_StrokeGenMeshPoolArgs
  {
    uint num_verts;
    uint num_edges;
    uint num_faces;
    uint dummy;
  };
  BLI_STATIC_ASSERT_ALIGN(SSBOData_StrokeGenMeshPoolArgs, 16)

  static inline uint get_prim_base_addr_ibo(bool is_ibo_fmt_u16, uint prim_id, uint num_verts_per_prim)
  {
    uint prim_vert_beg = prim_id * num_verts_per_prim;
    return (is_ibo_fmt_u16 ? (prim_vert_beg / 2u) : (prim_vert_beg));
  }

  static inline uint get_vbo_addr(bool is_ibo_fmt_u16, uint ibo_data, uint prim_id, uint num_verts_per_prim)
  {
    uint prim_vert_beg = prim_id * num_verts_per_prim;
    
    uint ibo_data_16h = (ibo_data >> 16u);
    uint ibo_data_16l = (ibo_data & 0xFFFFu);
    return (is_ibo_fmt_u16
              ? ( /* TODO: not sure about which 16 bits to use really */ 
                ((prim_vert_beg & 1u) == 0u) ? ibo_data_16l : ibo_data_16h)
              : ibo_data
            );
  }
  /** } */





#ifdef __cplusplus

// Template to set buffer size in compile time
using SSBO_StrokeGenTest = draw::StorageArrayBuffer<uint, 4096 * 4, true>;

using SSBO_StrokeGenMeshPool = draw::StorageArrayBuffer<uint, 2048 * 2048 * 4, true>;
using SSBO_StrokeGenMeshPoolArgs = draw::StorageBuffer<SSBOData_StrokeGenMeshPoolArgs>;
  
using SSBO_BnprScanData = draw::StorageArrayBuffer<uint, 2048 * 2048 * 2, true>;
using SSBO_BnprScanAggregates = draw::StorageArrayBuffer<uint, 512 * 16, true>;
using UBO_BnprTreeScan = draw::UniformBuffer<UBData_TreeScan>;

}

// namespace blender::bnpr
#endif
