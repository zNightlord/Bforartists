/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "bnpr_shader_shared.hh"

#include <set>

namespace blender::gpu {
class StorageBuf;
}

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
  float rot_angle_ubo_view_matrices_cache_;
  UBO_ViewMatrices ubo_view_matrices_cache_; // for visualizing contour curves in different view
  UBO_ViewMatrices ubo_view_matrices_cache_2_; // for visualizing contour curves in different view


  SSBO_StrokeGenMeshPoolCounters ssbo_bnpr_mesh_pool_counters_;
  SSBO_StrokeGenMeshPoolCounters ssbo_bnpr_mesh_pool_counters_prev_; // keep counters from last mesh extraction iter

  SSBO_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_[2];

  inline GPUStorageBuf *ssbo_dyn_mesh_counters_in_()
  {
    return ssbo_dyn_mesh_counters_[0];
  }
  inline GPUStorageBuf *ssbo_dyn_mesh_counters_out_() { return ssbo_dyn_mesh_counters_[1]; }
  SSBO_StrokeGenEdgeSplitCounters ssbo_edge_split_counters_;
  SSBO_StrokeGenEdgeCollapseCounters ssbo_edge_collapse_counters_;
  SSBO_StrokeGenEdgeFlipCounters ssbo_edge_flip_counters_;
  SSBO_StrokeGenFaceSplitCounters ssbo_face_split_counters_;

  SSBO_StrokeGenTemporalRecordCounters ssbo_temporal_record_counters_; 


  // Dispatch Args --------------------------------------------------------
  SSBO_IndirectDrawArgs ssbo_bnpr_contour_mesh_draw_args_;
  SSBO_IndirectDrawArgs ssbo_bnpr_2d_sample_draw_args_;
  SSBO_IndirectDrawArgs ssbo_bnpr_vert_debug_draw_args_;

  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_filtered_edge_; 
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_filtered_vert_;

  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_split_edge_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_collapsed_edge_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_flip_edge_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_split_face_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_remeshed_edges_;
  SSBO_IndirectDispatchArgs ssbo_indirect_dispatch_args_per_remeshed_verts_;

  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_edge_dispatch_args_;
  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_vert_dispatch_args_;
  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_frag_dispatch_args_;
  SSBO_IndirectDispatchArgs ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_;

  SSBO_IndirectDispatchArgs ssbo_bnpr_temporal_record_dispatch_args_; 



  // Persistent Buffers ------------------------------------------------------------------------
  // TODO: Consider compress some of them when the data became static
  SSBO_StrokeGenMeshBufPerVert<float, 3> ssbo_vbo_full_;                        // 96MB    Total
  SSBO_StrokeGenMeshBufPerVert<uint, 1> ssbo_vert_to_edge_list_header_;         // 32MB    128MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 4> ssbo_edge_to_vert_;                     // 256MB   384MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 4> ssbo_edge_to_edges_;                    // 256MB   640MB
  SSBO_StrokeGenMeshBufPerVert<uint, 1> ssbo_vert_flags_;                       // 32MB    672MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 1> ssbo_edge_flags_;                       // 64MB
  SSBO_StrokeGenMeshBufPerVert<float, 3> ssbo_vnor_;                            // 96MB
  SSBO_StrokeGenMeshBufPerVert<uint, 2> ssbo_vcurv_max_;                        // 32MB
  SSBO_StrokeGenMeshBufPerEdge<uint, 1> ssbo_subd_edge_tree_node_up_; 

  SSBO_StrokeGenMeshBufPerSelectedEdge<uint, 1> ssbo_selected_edge_to_edge_;    // 32MB    
  SSBO_StrokeGenMeshBufPerSelectedVert<uint, 1> ssbo_selected_vert_to_vert_;    // 16MB    
  SSBO_StrokeGenMeshBufPerEdge<uint, 6> ssbo_dbg_lines_;                        // 256MB

  SSBO_StrokeGenMeshBufPerContour<uint, 2> ssbo_contour_to_contour_;       //
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_rank_;       // 
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_list_len_;   // 
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_list_head_;  //
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_seg_rank_;   //
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_seg_len_;    //
  SSBO_StrokeGenMeshBufPerContour<uint, 6> ssbo_contour_snake_vpos_;       //
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_flags_;      //


  SSBO_StrokeGenMeshBufPerContour<uint, 32> ssbo_contour_temporal_records_[MAX_TEMPORAL_FRAMES];
  SSBO_StrokeGenMeshBufPerContour<uint, 1> ssbo_contour_snake_to_temporal_record_;  

  // Reusable Large Buffers -------------------------------------------------------------
  // note: don't ssbo_mesh_buffer_reuse_0_
  // from interpo contour tessellation till the end of append_subpass_setup_contour_edge_data
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_0_;                    // 256MB  Total  
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_1_;                   // 128MB  384MB
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_2_;                   // 128MB  512MB
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_3_;                    // 256MB  768MB
  SSBO_StrokeGenReusedLarge ssbo_mesh_buffer_reuse_4_;                    // 256MB 1024MB
  SSBO_StrokeGenReusedSmall ssbo_mesh_buffer_reuse_5_;                    // 64MB  1088MB
  SSBO_StrokeGenReusedSmall ssbo_mesh_buffer_reuse_6_;                    // 64MB  1152MB
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_7_;                   // 128MB 1280MB
  // note: don't ssbo_mesh_buffer_reuse_8_
  // from interpo contour tessellation till the end of append_subpass_setup_contour_edge_data
  SSBO_StrokeGenReusedMedium ssbo_mesh_buffer_reuse_8_;                   // 128MB 1280MB
  // Free after subpass_serialize_contour_edges 
  SSBO_StrokeGenMeshBufPerContour<uint, 8> ssbo_contour_edge_transfer_data_;  // 64MB
  // Free after subpass_visibility_split_contour_edges
  SSBO_StrokeGenMeshBufPerContour<uint, 5> ssbo_contour_raster_data_;         // 40MB
  // Notes: DO NOT use reuse_4 when remeshing, it holds per-vertex remesh len


  // Reused Buffer Scheme for Temporal Coherence   -----------------------------------------------
  // Long life-time buffers  -------------------------------
  inline GPUStorageBuf *reused_ssbo_subd_edge_tree_node_dw_()
  { // [append_subpasses_loop_subdiv, interpolated contour tessellation)
    return ssbo_mesh_buffer_reuse_3_;
  }
  inline GPUStorageBuf *reused_ssbo_contour_vert_to_old_edge_()
  { // [interpolated contour tessellation, append_subpass_setup_contour_edge_data]
    return ssbo_mesh_buffer_reuse_0_;
  }
  inline GPUStorageBuf *reused_ssbo_edge_to_temporal_record_()
  { // [interpolated contour tessellation, append_subpass_setup_contour_edge_data]
    return ssbo_mesh_buffer_reuse_8_;
  }



  // Reused Buffer Scheme for Basic Meshing ------------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vert_spatial_map_headers_()
  {
    return ssbo_mesh_buffer_reuse_4_; 
  }
  inline GPUStorageBuf *reused_ssbo_vert_merged_id_()
  {
    return ssbo_mesh_buffer_reuse_5_; 
  }
  inline GPUStorageBuf *reused_ssbo_vert_spatial_map_payloads_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }
  inline GPUStorageBuf *reused_ssbo_edge_spatial_map_payloads_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }
  inline void reused_ssbo_wedge_flooding_pointers_(uint iter,
                                                   GPUStorageBuf *&buf_in,
                                                   GPUStorageBuf *&buf_out) const
  {
    if (iter % 2u == 0u) {
      buf_in = ssbo_mesh_buffer_reuse_6_;
      buf_out = ssbo_mesh_buffer_reuse_0_;
    }
    else {
      buf_in = ssbo_mesh_buffer_reuse_0_;
      buf_out = ssbo_mesh_buffer_reuse_6_;
    }
  }


  // Reused Buffer Scheme throughout Edge Split, Collapse, Flip and Relocation ------
  inline GPUStorageBuf *reused_ssbo_vtx_remesh_len_()
  {
    return ssbo_mesh_buffer_reuse_4_; 
  }


  // Reused Buffer Scheme throughout Sqrt-3/Loop Subdiv ------
  // Don't reuse buffers reused by edge split/flip
  inline GPUStorageBuf *reused_ssbo_vpos_subd_()
  {
    return ssbo_mesh_buffer_reuse_8_;
  }
  inline GPUStorageBuf *reused_ssbo_epos_subd_()
  {
    return ssbo_mesh_buffer_reuse_7_; 
  }

  // Reused Buffer Scheme throughout Loop Subdiv ------
  // Don't reuse buffers reused by edge split
  inline GPUStorageBuf* reused_ssbo_subd_edge_vert_to_old_edge_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }


  // Reused Buffer Scheme for Edge Split --------------------------------------------
  inline GPUStorageBuf *reused_ssbo_per_edge_split_info_()
  {
    return ssbo_mesh_buffer_reuse_1_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_split_edge_info_()
  {
    return ssbo_mesh_buffer_reuse_2_; 
  }
  

  // Reused Buffer Scheme for Edge Collapse --------------------------------------------
  inline GPUStorageBuf *reused_ssbo_per_edge_collapse_info_in_(int step)
  {
    return (step % 2 == 0) ? ssbo_mesh_buffer_reuse_0_ : ssbo_mesh_buffer_reuse_1_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_edge_collapse_info_out_(int step)
  {
    return (step % 2 == 0) ? ssbo_mesh_buffer_reuse_1_ : ssbo_mesh_buffer_reuse_0_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_collapse_edge_info_()
  {
    return ssbo_mesh_buffer_reuse_2_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_vert_collapse_wedge_id_()
  {
    return ssbo_mesh_buffer_reuse_5_; 
  }



  // Reused Buffer Scheme for Edge Flip --------------------------------------------
  inline GPUStorageBuf *reused_ssbo_per_edge_flip_info_()
  {
    return ssbo_mesh_buffer_reuse_1_; 
  }
  inline GPUStorageBuf *reused_ssbo_per_flip_edge_info_()
  {
    return ssbo_mesh_buffer_reuse_2_; 
  }
  inline GPUStorageBuf *reused_ssbo_vertex_edge_flip_info_()
  {
    return ssbo_mesh_buffer_reuse_5_; 
  }

  // Reused Buffer Scheme for Face Split --------------------------------------------
  inline GPUStorageBuf *reused_ssbo_per_face_split_info_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }


  

  // Reused Buffer Scheme for Vertex Position Filtering -----------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vpos_temp_()
  {
    return ssbo_mesh_buffer_reuse_0_; 
  }


  // Reused Buffer Scheme for Quadric-based Filtering -----------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vert_quadric_data_in_(int step)
  {
    return step % 2u == 0u ? ssbo_mesh_buffer_reuse_3_ : ssbo_mesh_buffer_reuse_0_;
  }
  inline GPUStorageBuf *reused_ssbo_vert_quadric_data_out_(int step)
  {
    return step % 2u == 0u ? ssbo_mesh_buffer_reuse_0_ : ssbo_mesh_buffer_reuse_3_;
  }

  // Reused Buffer Scheme for Curvature Smoothing -----------------------------------------------
  inline GPUStorageBuf *reused_ssbo_vcurv_max_temp_()
  {
    return ssbo_mesh_buffer_reuse_1_;
  }




  // Reused buffers for Contour Processing -------------------------------------------------------
  // lifetime [append_subpass_extract_contour_edges, append_subpass_setup_contour_edge_data]
  inline GPUStorageBuf *reused_ssbo_edge_to_contour_()
  {
    return ssbo_mesh_buffer_reuse_6_; 
  }
  inline GPUStorageBuf *reused_ssbo_contour_temp_data_()
  {
    return ssbo_mesh_buffer_reuse_4_; 
  }

  // lifetime [append_subpass_extract_contour_edges, append_pass_remeshed_surface_depth_drawcall]
  inline GPUStorageBuf *reused_ssbo_face_to_vert_draw_depth_()
  {
    return ssbo_mesh_buffer_reuse_2_; 
  }

  inline GPUStorageBuf *reused_ssbo_frag_to_contour_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }
  inline GPUStorageBuf *reused_ssbo_frag_raster_data_()
  {
    return ssbo_mesh_buffer_reuse_3_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_contour_fragment_idmapping_()
  { // only used within the segscan passes
    return ssbo_mesh_buffer_reuse_8_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_contour_fragment_idmapping_()
  { 
    return reused_ssbo_frag_to_contour_();
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_contour_visibility_split_0_()
  {  // only used within the segscan passes
    return ssbo_mesh_buffer_reuse_8_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_contour_visibility_split_0_()
  { 
    return ssbo_mesh_buffer_reuse_2_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_contour_visibility_split_1_()
  {  // only used within the segscan passes
    return ssbo_mesh_buffer_reuse_1_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_contour_visibility_split_1_()
  {
    return ssbo_mesh_buffer_reuse_4_; 
  }
  inline GPUStorageBuf *reused_ssbo_contour_visibility_split_info_()
  {  
    return ssbo_mesh_buffer_reuse_8_;
  }
  inline GPUStorageBuf *reused_ssbo_frag_seg_head_to_visibility_split_contour_()
  {
    return ssbo_mesh_buffer_reuse_1_;
  }
  

  // lifetime [append_subpass_list_ranking, append_subpass_serialize_contour_edges]
  inline GPUStorageBuf *reused_ssbo_contour_edge_rank_()
  {
    return ssbo_mesh_buffer_reuse_5_;
  }
  inline GPUStorageBuf *reused_ssbo_contour_edge_list_len_()
  {
    return ssbo_mesh_buffer_reuse_2_;
  }
  inline GPUStorageBuf *reused_ssbo_contour_edge_list_head_info_()
  {
    return ssbo_mesh_buffer_reuse_8_;
  }

  // lifetime [append_subpass_serialize_contour_edges, append_subpass_contour_segmentation)
  inline GPUStorageBuf *reused_ssbo_in_segloopconv1d_data_contour_seg_denoise()
  {
    return ssbo_mesh_buffer_reuse_4_;
  }

  // lifetime (append_subpass_contour_segmentation, append_subpass_calc_contour_edges_draw_data)
  inline GPUStorageBuf *reused_ssbo_contour_to_start_sample_()
  {
    return ssbo_mesh_buffer_reuse_8_;
  }
  inline GPUStorageBuf *reused_ssbo_2d_sample_to_contour_()
  {
    return ssbo_mesh_buffer_reuse_4_;
  }
  inline GPUStorageBuf *reused_ssbo_contour_arc_len_param_()
  {
    return ssbo_mesh_buffer_reuse_2_;
  }

  // life time (append_subpass_contour_segmentation, strokegen_contour_2d_resample_eval_position]
  inline GPUStorageBuf *reused_ssbo_contour_2d_resample_raster_data_()
  {
    return ssbo_mesh_buffer_reuse_0_;
  }

  // shorter-lifetime buffers
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_()
  { // lifetime within segscan
    return ssbo_mesh_buffer_reuse_1_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_2d_resampler_accumulate_curvlen_()
  {
    return reused_ssbo_contour_arc_len_param_(); 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_2d_resampler_alloc_samples_()
  { // lifetime within scan, racing with above 2 buffers
    return ssbo_mesh_buffer_reuse_7_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_2d_resampler_alloc_samples_()
  {
    return reused_ssbo_contour_to_start_sample_();
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_2d_resample_contour_idmapping_()
  { // lifetime within scan, racing with above 2 buffers
    return ssbo_mesh_buffer_reuse_3_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_2d_resample_contour_idmapping_()
  {
    return reused_ssbo_2d_sample_to_contour_();
  }

  inline GPUStorageBuf *reused_ssbo_contour_2d_sample_geometry_()
  {
    return ssbo_contour_edge_transfer_data_; 
  }
  inline GPUStorageBuf *reused_ssbo_contour_2d_sample_topology_()
  {
    return ssbo_contour_raster_data_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_2d_sample_segmentation_0_()
  {
    return ssbo_mesh_buffer_reuse_0_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_2d_sample_segmentation_1_()
  {
    return ssbo_mesh_buffer_reuse_3_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_2d_sample_segmentation_0_()
  {
    return ssbo_mesh_buffer_reuse_7_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_2d_sample_segmentation_1_()
  {
    return ssbo_mesh_buffer_reuse_1_;
  }
  inline GPUStorageBuf *reused_ssbo_stroke_mesh_pool_()
  {
    return ssbo_mesh_buffer_reuse_0_; 
  }


  // lifetime [append_subpass_calc_contour_edges_draw_data, append_draw_contour_subpass]
  inline GPUStorageBuf *reused_ssbo_bnpr_mesh_pool_()
  {
    return ssbo_mesh_buffer_reuse_3_;
  }








  // Scan Working Buffers -----------------------------------
  SSBO_BnprScanData         ssbo_in_scan_data_;
  SSBO_BnprScanData         ssbo_out_scan_data_;
  SSBO_BnprScanAggregates   ssbo_scan_block_sum_;
  UBO_BnprTreeScan          ubo_bnpr_tree_scan_infos_;
  SSBO_BnprTreeScan         ssbo_tree_scan_infos_[32]; 
  SSBO_IndirectDispatchArgs ssbo_scan_dispatch_args_;
  inline GPUStorageBuf *reused_ssbo_tree_scan_infos_contour_fragment_idmapping_()
  {
    return ssbo_tree_scan_infos_[0];
  }

  inline GPUStorageBuf *reused_ssbo_tree_scan_infos_contour_visibility_split_()
  {
    return ssbo_tree_scan_infos_[0];
  }

  inline GPUStorageBuf *reused_ssbo_tree_scan_infos_contour_segmentation_()
  {
    return ssbo_tree_scan_infos_[0];
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_contour_segmentation_step_0()
  {
    return ssbo_mesh_buffer_reuse_5_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_contour_segmentation_step_0()
  {
    return ssbo_mesh_buffer_reuse_7_;
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_input_contour_segmentation_step_1()
  {
    return ssbo_mesh_buffer_reuse_6_; 
  }
  inline GPUStorageBuf *reused_ssbo_tree_scan_output_contour_segmentation_step_1()
  {
    return ssbo_mesh_buffer_reuse_8_;
  }

  inline GPUStorageBuf *reused_ssbo_tree_scan_infos_2d_resampler_()
  {
    return ssbo_tree_scan_infos_[0];
  }


  // Segmented Loop Convolution Working Buffers -----------------------------------
  SSBO_SegLoopConv1DData     ssbo_in_segloopconv1d_data_;
  SSBO_SegLoopConv1DData     ssbo_out_segloopconv1d_data_;
  SSBO_SegLoopConvDebugData  ssbo_debug_segloopconv1d_data_;
  SSBO_SegLoopConvPatchTable ssbo_segloopconv1d_patch_table_;
  UBO_SegLoopConv1D          ubo_segloopconv1d_;
  SSBO_SegLoopConv1DInfo     ssbo_segloopconv1d_info_;
  SSBO_IndirectDispatchArgs  ssbo_segloopconv1d_dispatch_args_;

  // List Ranking Working Buffers ------------------------------------------------------------------------------------
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
    , rot_angle_ubo_view_matrices_cache_(.0f)
  {
  }
  ~GPUBufferPoolModule()
  {
  }


  void on_begin_sync(
    const DRWView* drw_view, bool upload_list_ranking_test_data,
    bool update_view_matrices_for_contour_extraction = true,
    bool update_view_matrices_for_dbg_view_rotation = false
  );
  void sync_object(Object *ob);
  void end_sync();


  /* -------------------------------------------------------------------- */
  /** \name List Ranking Validation
   * \{ */
  void build_list_ranking_testing_data(bool test_loop_lists = false);
  /** \} */


  /* -------------------------------------------------------------------- */
  /** \name View Debugging
   * \{ */
  static void get_dbg_view_rot_angle(float& rot_ang);
  void should_dbg_rotate_view_matrices_cache(float& rot_ang, bool& dbg_rotate_view_matrix);
  /** \} */
};
}  // namespace blender::npr::strokegen
