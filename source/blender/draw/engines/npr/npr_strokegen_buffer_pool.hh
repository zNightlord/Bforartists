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

  SSBO_StrokeGenTest ssbo_bnpr_test_;

  SSBO_StrokeGenMeshLarge     ssbo_bnpr_mesh_pool_; // contains contour edge data
  SSBO_StrokeGenMeshPoolCounters ssbo_bnpr_mesh_pool_counters_;
  SSBO_IndirectDrawArgs ssbo_bnpr_mesh_pool_draw_args_;
  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_edge_dispatch_args_; 

  // TODO: these buffers are huge. consider re-use them for each mesh. 
  SSBO_StrokeGenMeshGiant_Float ssbo_vbo_full_;
  SSBO_StrokeGenMeshLarge ssbo_edge_to_vert_;
  SSBO_StrokeGenMeshLarge ssbo_edge_to_edges_;
  SSBO_StrokeGenMeshMedium ssbo_edge_to_contour_;
  SSBO_StrokeGenMeshTiny ssbo_contour_to_contour_;
  SSBO_StrokeGenMeshMinimum ssbo_contour_edge_rank_;
  SSBO_StrokeGenMeshMinimum ssbo_contour_edge_list_len_;
  SSBO_StrokeGenMeshMinimum ssbo_contour_edge_list_head_; 
  SSBO_StrokeGenMeshVertHashTable ssbo_vert_spatial_map_headers_;
  SSBO_StrokeGenMeshLarge ssbo_mesh_buffer_reuse_0_;
  SSBO_StrokeGenMeshLarge ssbo_vert_merged_id_;

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

  /** CPU Buffers */
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
