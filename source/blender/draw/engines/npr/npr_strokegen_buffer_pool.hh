/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "bnpr_shader_shared.hh"

#include <set>

namespace blender::npr::strokegen {
class Instance;

class GPUBufferPoolModule {
  friend class StrokeGenPassModule;
  friend class StrokegenMeshRasterPass; 

 private:
  /** Instance */
  Instance &instance;

  /** GPU Buffers */
  UBO_ViewMatrices ubo_view_matrices_;
  UBO_ViewMatrices ubo_view_matrices_cache_; // for visualizing contour curves in different view


  SSBO_StrokeGenMeshPoolCounters ssbo_bnpr_mesh_pool_counters_;
  SSBO_StrokeGenMeshPoolCounters ssbo_bnpr_mesh_pool_counters_prev_; // keep counters from last mesh extraction iter

  SSBO_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_[MAX_REMESHING_ITERS][2]; 
  SSBO_StrokeGenEdgeSplitCounters ssbo_edge_split_counters_;

  SSBO_IndirectDrawArgs ssbo_bnpr_mesh_pool_draw_args_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_filtered_edge_; 
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_filtered_vert_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_split_edge_; 
  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_edge_dispatch_args_; 

  // Persistent Buffers -------------------------------------------------------------------
  SSBO_StrokeGenMeshBufPerVert<float, 3> ssbo_vbo_full_;                // 96MB    Total
  SSBO_StrokeGenMeshBufPerVert<uint, 1> ssbo_vert_to_edge_list_header_; // 32MB    128MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 4> ssbo_edge_to_vert_;             // 256MB   384MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 4> ssbo_edge_to_edges_;            // 256MB   640MB
  SSBO_StrokeGenMeshBufPerVert<uint, 1> ssbo_vert_flags_;               // 32MB    672MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 1> ssbo_edge_flags_;               // 64MB    736MB
  SSBO_StrokeGenMeshBufPerContour<uint, 2> ssbo_contour_to_contour_;    // 
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_edge_rank_;     // 
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_edge_list_len_; // 
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_edge_list_head_;// 

  // Reusable Large Buffers ---------------------------------------------------------------
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_0_;                    // 256MB  Total  
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_1_;                   // 128MB  384MB
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_2_;                   // 128MB  512MB
  // finally contains contour edge data
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_3_;                    // 256MB  768MB
  // temporally used for hashing
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_4_;                    // 256MB 1024MB
  SSBO_StrokeGenReusedSmall ssbo_mesh_buffer_reuse_5_;                    // 64MB  1088MB
  SSBO_StrokeGenReusedSmall ssbo_mesh_buffer_reuse_6_;                    // 64MB  1152MB


  // Reused Buffer Scheme for Basic Meshing ------------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vert_spatial_map_headers_()
  {
    return ssbo_mesh_buffer_reuse_4_; 
  }
  inline GPUStorageBuf *reused_ssbo_vert_merged_id_()
  {
    return ssbo_mesh_buffer_reuse_5_; 
  }


  // Reused Buffer Scheme for Edge Splitting --------------------------------------------
  inline GPUStorageBuf *reused_ssbo_per_edge_split_info_()
  {
    return ssbo_mesh_buffer_reuse_1_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_split_edge_info_()
  {
    return ssbo_mesh_buffer_reuse_2_; 
  }


  // Reused Buffer Scheme for Mesh Filtering -----------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vert_spatial_map_payloads_()
  {
    return ssbo_mesh_buffer_reuse_0_;  
  }
  inline GPUStorageBuf *reused_ssbo_edge_spatial_map_payloads_()
  {
    return ssbo_mesh_buffer_reuse_0_; 
  }
  inline void reused_ssbo_wedge_flooding_pointers_(uint iter, GPUStorageBuf *&buf_in, GPUStorageBuf*& buf_out) const
  {
    if (iter % 2u == 0u) {
      buf_in = ssbo_mesh_buffer_reuse_6_;
      buf_out = ssbo_mesh_buffer_reuse_0_; 
    }else {
      buf_in = ssbo_mesh_buffer_reuse_0_;
      buf_out = ssbo_mesh_buffer_reuse_6_; 
    }
  }
  inline GPUStorageBuf *reused_ssbo_filtered_edge_to_edge_()
  {
    return ssbo_contour_to_contour_; 
  }
  inline GPUStorageBuf *reused_ssbo_filtered_vert_to_vert_()
  {
    return ssbo_contour_edge_rank_; 
  }
  inline GPUStorageBuf *reused_ssbo_edge_quadric_data()
  { // TODO: we can actually share one buffer with vert_quadric data if this also goes iterative
    return ssbo_mesh_buffer_reuse_5_; 
  }
  inline void reused_ssbo_vert_quadric_data_(int iter, GPUStorageBuf*& buf_in, GPUStorageBuf*& buf_out)
  {
    buf_in = iter % 2 == 0 ? ssbo_mesh_buffer_reuse_3_ : ssbo_mesh_buffer_reuse_0_;
    buf_out = iter % 2 == 0 ? ssbo_mesh_buffer_reuse_0_ : ssbo_mesh_buffer_reuse_3_; 
  }
  inline void reused_ssbo_filtered_normal_edge_(int iter, GPUStorageBuf*& buf_in, GPUStorageBuf*& buf_out)
  {
    buf_in = iter % 2 == 0 ? ssbo_mesh_buffer_reuse_4_ : ssbo_mesh_buffer_reuse_2_;
    buf_out = iter % 2 == 0 ? ssbo_mesh_buffer_reuse_2_ : ssbo_mesh_buffer_reuse_4_; 
  }
  inline GPUStorageBuf *reused_ssbo_filtered_normal_vert_()
  {
    return ssbo_mesh_buffer_reuse_1_; 
  }

  // Reused Buffers for per-batch Contour Processing
  inline GPUStorageBuf *reused_ssbo_bnpr_mesh_pool_()
  {
    return ssbo_mesh_buffer_reuse_3_; 
  }
  inline GPUStorageBuf *reused_ssbo_edge_to_contour_()
  {
    return ssbo_mesh_buffer_reuse_6_; 
  }



  SSBO_BnprScanData       ssbo_in_scan_data_;
  SSBO_BnprScanData       ssbo_out_scan_data_;
  SSBO_BnprScanAggregates ssbo_scan_block_sum_;
  UBO_BnprTreeScan        ubo_bnpr_tree_scan_infos_;

  SSBO_SegLoopConv1DData     ssbo_in_segloopconv1d_data_;
  SSBO_SegLoopConv1DData     ssbo_out_segloopconv1d_data_;
  SSBO_SegLoopConvDebugData  ssbo_debug_segloopconv1d_data_;
  SSBO_SegLoopConvPatchTable ssbo_segloopconv1d_patch_table_;
  UBO_SegLoopConv1D          ubo_segloopconv1d_;

  SSBO_ListRankingCounters            ssbo_list_ranking_anchor_counters_; // x1 slot per splicing iteration
  SSBO_ListRankingCounters            ssbo_list_ranking_splice_counters_; // x1 slot per re-link iteration
  SSBO_ListRankingAllocationCounters  ssbo_list_ranking_addressing_counters_; // allocate space for serialized lists
  SSBO_IndirectDispatchArgs           ssbo_list_ranking_indirect_dispatch_args_per_anchor[NUM_ITERS_BNPR_LIST_RANK_SPLICE+1];
  SSBO_IndirectDispatchArgs           ssbo_list_ranking_indirect_dispatch_args_per_spliced[NUM_ITERS_BNPR_LIST_RANK_RELINK+1];
  SSBO_ListRankingTags                ssbo_list_ranking_tags_[2];
  SSBO_ListRankingRanks               ssbo_list_ranking_ranks_;
  SSBO_ListRankingLinksStagingBuf     ssbo_list_ranking_links_staging_buf_;
  SSBO_ListRankingLinks               ssbo_list_ranking_links_;
  SSBO_ListRankingAnchorToNode        ssbo_list_ranking_anchor_to_node_[2];
  SSBO_ListRankingSplicedNodeId       ssbo_list_ranking_spliced_node_id_[NUM_ITERS_BNPR_LIST_RANK_RELINK];
  SSBO_ListRankingSerializedTopo      ssbo_list_ranking_serialized_topo_;

  SSBO_ListRankingAnchorToNextAnchor  ssbo_list_ranking_anchor_to_next_anchor_;
  SSBO_ListRankingAnchorJumpingInfo   ssbo_list_ranking_per_anchor_sublist_jumping_info_[2];
  SSBO_ListRankingNodeToAnchor        ssbo_list_ranking_node_to_anchor_;


  SSBO_ListRankingLinks ssbo_list_ranking_debug_;

  UBO_ListRanking                     ubo_list_ranking_splicing_;
  SSBO_ListRankingInputs              ssbo_list_ranking_inputs_; 

  // Note: Whenever you found weird shader error:
  // Check following things:
  // *) if the data type of your buffer is correct, this often happens with a custom data type
  // *) the undefined variable "_ubo/_ssbo_xxxx_"
  //    Check the if there is any commented code at the end of the shader file.
  //    "// XXX" will fuck up the #pragma .... at the beginning of your code
  // *) Check the included XXXlib.glsl, if a #endif is at the end without blank lines followed,
  //    it fucks up with your #pragma
  // *) Check if you have included the shared .hh file
  //    with sth. like .typedef_source("bnpr_shader_shared.hh")
  // *) Check if you have multiple instances of a XXX.glsl included in the shaderinfo,
  //    for example
  // GPU_SHADER_CREATE_INFO(bnpr_scan_uint_add)
  //   .typedef_source("bnpr_shader_shared.hh")
  //   . ... ... ...;
  // GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_compact_anchors)
  //   .typedef_source("bnpr_shader_shared.hh")
  //   . ... ... ...
  //   .additional_info("bnpr_scan_uint_add")
  //   . ... ... ...;
  // *) Check if you uses any weird keyword as variable names, such as "packed"
  //
  //
  // *) When using parallel primitive lib (scan, loopconv, etc)
  //    with multiple kernels in one shader file
  //    and you dont need the primitive in a certain kernel,
  //    use define guard to exclude the scan library like this:
  //    (shaderinfo) .define("BNPR_STROKEGEN_SCAN_NO_SUBGROUP_CODEGEN_LIB", "1")
  //    Also, remember to check if there are any other confliction with the primitive library
  // 
  //
  // For "syntax error, unexpected $undefined at token "<undefined>""
  // *) Check if your multi-line macro has space after the '\'
  // *) Check any other garbage characters (something other than a printable ASCII character, a space, a tab, or a newline).
  //    see https://stackoverflow.com/questions/10877386/opengl-shader-compilation-errors-unexpected-undefined-at-token-undefined
  //

  /** Mapped Buffers */
  /** Temp data for testing list ranking  */
  Vector<int> listranking_test_nodes_prev_next;
  std::vector<int> listranking_test_nodes_rank;
  std::vector<int> listranking_test_nodes_head;
  std::vector<int> listranking_test_nodes_tail;
  std::vector<int> listranking_test_nodes_list_len;
  std::set<int> listranking_test_head_nodes;

  bool listranking_test_data_uploaded;
  bool listranking_test_data_validated;



 public:
  GPUBufferPoolModule(Instance &inst) :
    instance(inst)
    , listranking_test_data_uploaded(false)
    , listranking_test_data_validated(false)
  {
  }
  ~GPUBufferPoolModule()
  {
  }

  void on_begin_sync(const DRWView* drw_view, bool upload_list_ranking_test_data, bool update_view_matrices_for_contour_extraction = true);
  void sync_object(Object *ob);
  void end_sync();


  /* -------------------------------------------------------------------- */
  /** \name List Ranking Validation
   * \{ */
  void build_list_ranking_testing_data(bool test_loop_lists = false);
  /** \} */
};
}  // namespace blender::npr::strokegen
