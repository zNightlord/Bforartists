/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Shader module that manage shader libraries, deferred compilation,
 * and static shader usage.
 */

#include "gpu_shader_create_info.hh"

#include "npr_strokegen_shader.hh"

namespace blender::npr::strokegen {

/* -------------------------------------------------------------------- */
/** \name Module
 *
 * \{ */

StrokeGenShaderModule *StrokeGenShaderModule::g_shader_module = nullptr;

StrokeGenShaderModule *StrokeGenShaderModule::module_get()
{
  if (g_shader_module == nullptr) {
    /* TODO(@fclem) thread-safety. */
    g_shader_module = new StrokeGenShaderModule();
  }
  return g_shader_module;
}

void StrokeGenShaderModule::module_free()
{
  if (g_shader_module != nullptr) {
    /* TODO(@fclem) thread-safety. */
    delete g_shader_module;
    g_shader_module = nullptr;
  }
}

StrokeGenShaderModule::StrokeGenShaderModule()
{
  for (GPUShader *&shader : shaders_) {
    shader = nullptr;
  }

#ifdef DEBUG
  /* Ensure all shader are described. */
  for (auto i : IndexRange(MAX_SHADER_TYPE)) {
    const char *name = static_shader_create_info_name_get(eShaderType(i));
    if (name == nullptr) {
      std::cerr << "bnpr: Missing case for eShaderType(" << i
                << ") in static_shader_create_info_name_get().";
      BLI_assert(0);
    }
    const GPUShaderCreateInfo *create_info = GPU_shader_create_info_get(name);
    BLI_assert_msg(create_info != nullptr, "bnpr: Missing create info for static shader.");
  }
#endif
}

StrokeGenShaderModule::~StrokeGenShaderModule()
{
  for (GPUShader *&shader : shaders_) {
    DRW_SHADER_FREE_SAFE(shader);
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Static shaders
 *
 * \{ */

const char *StrokeGenShaderModule::static_shader_create_info_name_get(eShaderType shader_type)
{ // returns the name of ShaderCreateInfo
  switch (shader_type) {
    case DEPTH:
      return "basic_depth_mesh";
    case POINTCLOUD_DEPTH:
      return "basic_depth_pointcloud";
    case CURVES_DEPTH:
      return "basic_depth_curves";
    case DEPTH_CONSERVATIVE:
      return "basic_depth_mesh_conservative";
    case POINTCLOUD_DEPTH_CONSERVATIVE:
      return "basic_depth_pointcloud_conservative";
    // case COMPUTE_TEST:
    //   return "bnpr_strokegen_test_xxx";

    case GPU_MESHING_BOOSTRAP:
      return "bnpr_geom_extract_boostrap"; 
    case CONTOUR_GEOM_EXTRACT:
      return "bnpr_geom_extract";
    case CONTOUR_GEOM_EXTRACT_IBO_16BIT:
      return "bnpr_geom_extract_ibo_16bits";
    case FILL_DRAW_ARGS_CONTOUR_EDGES:
      return "bnpr_geom_fill_draw_args_contour_edges";
    case FILL_DRAW_ARGS_REMESHED_SURFACE_DEPTH:
      return "bnpr_geom_fill_draw_args_remeshed_surface_depth"; 
    case FILL_DISPATCH_ARGS_CONTOUR_EDGES:
      return "strokegen_fill_dispatch_args_per_contour_edge";
    case FILL_DISPATCH_ARGS_CONTOUR_VERTS:
      return "strokegen_fill_dispatch_args_per_contour_vert";
    case FILL_DISPATCH_ARGS_CONTOUR_FRAGS:
      return "strokegen_fill_dispatch_args_per_contour_fragment";
    case FILL_DISPATCH_ARGS_CONTOUR_2D_SAMPLES:
      return "strokegen_fill_dispatch_args_per_contour_2d_sample"; 
    case EXTRACT_CURRENT_MESH_CONTOUR_DATA:
      return "bnpr_geom_extract_mesh_contour_data";


    case CLEAR_FRAG_TO_CONTOUR_IDMAPPING:
      return "strokegen_clear_frag_to_contour_idmapping";
    case PREP_FRAG_TO_CONTOUR_IDMAPPING:
      return "strokegen_prep_segscan_frag_to_contour_idmapping";
    case FINISH_FRAG_TO_CONTOUR_IDMAPPING:
      return "strokegen_finish_segscan_frag_to_contour_idmapping"; 
    case CONTOUR_FRAG_VISIBILITY_TEST:
      return "strokegen_contour_frag_visibility_test";
    case CONTOUR_FRAG_SETUP_VISIBILITY_SEGMENTATION:
      return "strokegen_contour_frag_setup_visibility_segmentation"; 
    case CONTOUR_VISIBILITY_SPLIT_STEP_0:
      return "strokegen_contour_visibility_split_step_0";
    case CONTOUR_VISIBILITY_SPLIT_STEP_1:
      return "strokegen_contour_visibility_split_step_1";
    case CONTOUR_VISIBILITY_SPLIT_STEP_2:
      return "strokegen_contour_visibility_split_step_2";
    case CONTOUR_VISIBILITY_SPLIT_STEP_3:
      return "strokegen_contour_visibility_split_step_3";

    case FILL_CONTOUR_EDGE_RANKING_INPUTS:
      return "strokegen_fill_cotour_edge_ranking_inputs"; 

    case SERIALIZE_RANKED_CONTOUR_EDGES:
      return "strokegen_serialize_contour_edges_pass_0";
    case SETUP_CONTOUR_SEGMENTATION:
      return "strokegen_setup_contour_segmentation";
    case FINISH_CONTOUR_SEGMENTATION:
      return "strokegen_finish_contour_segmentation";
    case PREP_CUSP_SEGMENTATION:
      return "strokegen_prep_contour_cusp_segmentation";

    case RECORD_TEMPORAL_CONTOUR_DATA:
      return "strokegen_record_temporal_contour_data"; 

    case CALC_CONTOUR_EDGES_DRAW_DATA:
      return "strokegen_calc_contour_edges_draw_data";

    case INDIRECT_DRAW_CONTOUR_MESH_DEPTH:
      return "bnpr_geom_draw_remeshed_surface_depth"; 
    case INDIRECT_DRAW_CONTOUR_EDGES:
      return "bnpr_geom_draw_contour_edges";

    case CONTOUR_ARCLEN_PARAMETERIZATION:
      return "strokegen_contour_2d_resample_arclen_parameterization";
    case CONTOUR_ALLOC_2D_SAMPLES:
      return "strokegen_contour_2d_resample_alloc_samples";
    case CONTOUR_ALLOC_2D_SAMPLES_FINISH:
      return "strokegen_contour_2d_resample_alloc_samples_finish";
    case CLEAR_2D_SAMPLE_TO_CONTOUR_IDMAPPING:
        return "strokegen_contour_2d_resample_idmapping_clear_buffer"; 
    case PREP_2D_SAMPLE_TO_CONTOUR_IDMAPPING:
      return "strokegen_contour_2d_resample_idmapping_setup_segscan";
    case FINISH_2D_SAMPLE_TO_CONTOUR_IDMAPPING:
      return "strokegen_contour_2d_resample_idmapping_finish";
    case CONTOUR_2D_SAMPLES_EVAL_POSITION:
      return "strokegen_contour_2d_resample_eval_position";
    case CONTOUR_2D_SAMPLES_EVAL_TOPO_STEP_0:
        return "strokegen_contour_2d_resample_eval_topo_step_0"; 
    case CONTOUR_2D_SAMPLES_SEGMENTATION_PREP_SEGTAILS:
      return "strokegen_contour_2d_resample_segmentation_prep_seg_tails";
    case CONTOUR_2D_SAMPLES_EVAL_TOPO_SETUP_SEGSCAN:
      return "strokegen_contour_2d_resample_segmentation_setup_segscan";
    case CONTOUR_2D_SAMPLES_EVAL_TOPO_FINISH_SEGSCAN:
      return "strokegen_contour_2d_resample_eval_topo_finish";
    case CONTOUR_2D_SAMPLES_REJECT_FAKE_CORNERS:
      return "strokegen_contour_2d_resample_eval_topo_remove_fake_corners";
    case CONTOUR_2D_SAMPLES_CALC_RENDER_DATA:
      return "strokegen_calc_contour_2d_curve_render_data";
    case FILL_DRAW_ARGS_CONTOUR_2D_SAMPLES:
      return "bnpr_geom_fill_draw_args_contour_2d_samples"; 
    case INDIRECT_DRAW_CONTOUR_2D_SAMPLES:
      return "bnpr_geom_draw_contour_2d_samples"; 

    case CONTOUR_PIXEL_COMPRESS:
      return "bnpr_compress_contour_pixels";
    case CONTOUR_PIXEL_COMPRESS_DBG:
      return "bnpr_compress_contour_pixels_dbg";

    case MESH_COLLECT_VBO:
      return "bnpr_geom_extract_collect_verts"; 
    case MESH_VERT_MERGE_INIT:
      return "bnpr_meshing_merge_verts_bootstrap";
    case MESH_VERT_MERGE_HASH:
      return "bnpr_meshing_merge_verts_spatial_hashing";
    case MESH_VERT_MERGE_REINDEX:
      return "bnpr_meshing_merge_verts_deduplicate";

    case MESH_COLLECT_EDGE_ADJ_IBO:
      return "bnpr_geom_extract_collect_edges";
    case MESH_COLLECT_EDGE_ADJ_IBO_16BIT:
      return "bnpr_geom_extract_collect_edges_16bits"; 
    case MESH_EDGE_ADJACENCY_INIT:
      return "bnpr_meshing_merge_edges_bootstrap"; 
    case MESH_EDGE_ADJACENCY_HASH:
      return "bnpr_meshing_merge_edges_hashing";
    case MESH_EDGE_ADJACENCY_HASH_FINISH:
      return "bnpr_meshing_merge_edges_hashing_finish"; 
    case MESH_EDGE_ADJACENCY_FILL:
      return "bnpr_meshing_merge_edges_fill";

    case MESH_WEDGE_FLOODING_ITER:
      return "bnpr_meshing_wedge_flooding_iter";
    case MESH_WEDGE_FLOODING_FIRST_ITER_INIT:
      return "bnpr_meshing_wedge_flooding_iter_init"; 
    case MESH_WEDGE_FLOODING_LAST_ITER_OUTPUT_FLAGS:
      return "bnpr_meshing_wedge_flooding_last_iter_output_flags"; 
    case MESH_WEDGE_FLOODING_LAST_ITER_COMPACTION:
      return "bnpr_meshing_wedge_flooding_last_iter_compaction";

    case MESH_WEDGE_MARK_EDGES_ON_SELECTION_BORDER:
      return "bnpr_meshing_edge_selection_mark_selection_border"; 

    case MESH_SELECT_VERTS_FROM_SELECTED_EDGES:
      return "strokegen_select_verts_from_selected_edges"; 
    case MESH_EXPAND_VERTS_FROM_SELECTED_EDGES:
      return "strokegen_expand_verts_from_selected_edges";
    case MESH_COMPACT_SELECTED_VERTS:
      return "strokegen_compact_selected_verts";  

    case FILL_DISPATCH_ARGS_FILTERED_EDGES:
      return "bnpr_meshing_fill_dispatch_args_per_filtered_edge";
    case FILL_DISPATCH_ARGS_FILTERED_VERTS:
      return "bnpr_meshing_fill_dispatch_args_per_filtered_vert";

    case MESH_FILTER_VNOR_FILTERING:
      return "bnpr_meshing_surf_filtering_vnor_filtering";
    case MESH_FILTER_VPOS_FILTERING:
      return "bnpr_meshing_surf_filtering_vpos_filtering"; 

    case MESH_FILTER_VQUADRIC_DIFFUSION:
      return "bnpr_meshing_surf_filtering_vquadric_diffusion"; 
    case MESH_FILTER_QUADRIC_VPOS_FILTERING:
      return "bnpr_meshing_surf_filtering_quadric_vpos_filtering";
    case MESH_FILTER_SUBDIV_VPOS_SMOOTH:
      return "bnpr_meshing_surf_filtering_subdiv_vpos_smoothing"; 
    case MESH_FILTER_SUBDIV_VPOS_SMOOTH_FINISH:
      return "bnpr_meshing_surf_filtering_subdiv_vpos_smoothing_finish";
    case MESH_FILTER_LOOP_SUBDIV_EDGE_POINTS:
      return "bnpr_meshing_surf_filtering_subdiv_edge_points"; 
    case MESH_FILTER_VCURV_SMOOTHING:
      return "bnpr_meshing_surf_filtering_vcurv_smoothing";
    case MESH_FILTER_VCURV_SMOOTHING_OUTPUT_REMESH_LEN:
      return "bnpr_meshing_surf_filtering_vcurv_smoothing_output_remesh_len"; 

    case FILL_DISPATCH_ARGS_REMESHED_EDGES:
      return "strokegen_remeshing_fill_dispatch_args_per_remeshed_edge";
    case FILL_DISPATCH_ARGS_REMESHED_VERTS:
      return "strokegen_remeshing_fill_dispatch_args_per_remeshed_vert"; 

    case MESH_OP_SPLIT_EDGE_INIT:
      return "bnpr_meshing_edge_split_init"; 
    case MESH_OP_SPLIT_EDGE_COMPACT:
      return "bnpr_meshing_edge_split_compact";
    case FILL_DISPATCH_ARGS_SPLIT_EDGES:
      return "strokegen_remeshing_fill_dispatch_args_per_split_edge";
    case MESH_OP_SPLIT_EDGE_EXCLUDE_BORDER:
      return "bnpr_meshing_edge_split_exclude_border"; 
    case MESH_OP_SPLIT_EDGE_RESOLVE_CONFLICT:
      return "bnpr_meshing_edge_split_resolve_conflict";
    case MESH_OP_SPLIT_EDGE_EXECUTE:
      return "bnpr_meshing_edge_split_execute"; 

    case MESH_OP_COLLAPSE_EDGE_INIT: 
      return "bnpr_meshing_edge_collapse_init"; 
    case MESH_OP_COLLAPSE_EDGE_COMPACT: 
      return "bnpr_meshing_edge_collapse_compact"; 
    case FILL_DISPATCH_ARGS_COLLAPSE_EDGES: 
      return "strokegen_remeshing_fill_dispatch_args_per_collapsed_edge"; 
    case MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_0: 
      return "bnpr_meshing_edge_collapse_resolve_conflict_0"; 
    case MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_1: 
      return "bnpr_meshing_edge_collapse_resolve_conflict_1"; 
    case MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_2: 
      return "bnpr_meshing_edge_collapse_resolve_conflict_2"; 
    case MESH_OP_COLLAPSE_EDGE_EXECUTE: 
      return "bnpr_meshing_edge_collapse_execute"; 

    case MESH_OP_FLIP_EDGE_INIT: 
      return "bnpr_meshing_edge_flip_init"; 
    case MESH_OP_FLIP_EDGE_COMPACT: 
      return "bnpr_meshing_edge_flip_compact"; 
    case FILL_DISPATCH_ARGS_FLIP_EDGES:
      return "strokegen_remeshing_fill_dispatch_args_per_flip_edge"; 
    case MESH_OP_FLIP_EDGE_VALIDATE: 
      return "bnpr_meshing_edge_flip_validate"; 
    case MESH_OP_FLIP_EDGE_RESOLVE_CONFLICT: 
      return "bnpr_meshing_edge_flip_resolve_conflict"; 
    case MESH_OP_FLIP_EDGE_EXECUTE: 
      return "bnpr_meshing_edge_flip_execute";

    case MESH_OP_SPLIT_FACE_INIT:
      return "bnpr_meshing_face_split_init"; 
    case MESH_OP_SPLIT_FACE_GENERATE_WORKS:
      return "bnpr_meshing_face_split_work_generation";
    case FILL_DISPATCH_ARGS_SPLIT_FACES:
      return "strokegen_remeshing_fill_dispatch_args_per_split_face"; 
    case MESH_OP_SPLIT_FACE_EXECUTE:
      return "bnpr_meshing_face_split_execute"; 

    case MESH_ANALYSE_VERT_NORMAL:
      return "bnpr_geom_analysis_order_0_vert_normal";
    case MESH_ANALYSE_VERT_NORMAL_VORONOIAREA:
      return "bnpr_geom_analysis_order_0_vert_normal_voroarea";
    case MESH_ANALYSE_VERT_NORMAL_TOPOFLAGS:
      return "bnpr_geom_analysis_order_0_vert_normal_topoflags";
    case MESH_ANALYSE_VERT_NORMAL_VORONOIAREA_TOPOFLAGS:
      return "bnpr_geom_analysis_order_0_vert_normal_voroarea_topoflags"; 
    case MESH_ANALYSE_EDGE_FEATURES:
      return "bnpr_geom_analysis_feature_edges"; 

    case MESH_ANALYSE_VERT_CONTOUR_GRADIENT:
      return "bnpr_geom_analysis_order_1_contour_grad_"; 

    case FILL_DRAW_ARGS_DBG_LINES:
      return "strokegen_remeshing_fill_draw_args_dbg_lines"; 
    case INDIRECT_DRAW_DBG_LINES:
      return "bnpr_geom_draw_debug_lines"; 

    case MESH_ANALYSE_VERT_CURV_PASS_0:
      return "bnpr_geom_analysis_order_1_vert_curv_pass_0";
    case MESH_ANALYSE_VERT_CURV_PASS_1_RUSINKIEWICZ:
      return "bnpr_geom_analysis_order_1_main_curvature_rusinkiewicz"; 
    case MESH_ANALYSE_VERT_CURV_PASS_1_JACQUES:
      return "bnpr_geom_analysis_order_1_main_curvature_jacques";

    case MESH_LOOP_SUBD_TREE_BUILD_NODES_UPWARDS_FOR_FACE_EDGES:
      return "strokegen_loop_subdiv_tree_build_nodes_upwards_for_face_edges"; 
    case MESH_LOOP_SUBD_TREE_INIT_NODES_DOWNWARDS:
      return "strokegen_loop_subdiv_tree_build_nodes_downwards_init_"; 
    case MESH_LOOP_SUBD_TREE_BUILD_NODES_DOWNWARDS:
      return "strokegen_loop_subdiv_tree_build_nodes_downwards_calc_";

    case MESH_COMPACT_NEW_TEMPORAL_RECORDS:
      return "strokegen_compact_new_temporal_contour_records"; 
    case MESH_CALCULATE_NEW_TEMPORAL_RECORDS_FIRST_ITER:
      return "strokegen_calculate_new_temporal_contour_records_first_iter";
    case MESH_CALCULATE_NEW_TEMPORAL_RECORDS_TRACE_ITER:
      return "strokegen_calculate_new_temporal_contour_records_trace_iter";
    case MESH_CALCULATE_NEW_TEMPORAL_RECORDS_LAST_ITER:
      return "strokegen_calculate_new_temporal_contour_records_last_iter"; 
    case MESH_RECONSTRUCT_OLD_TEMPORAL_RECORDS:
      return "strokegen_reconstruct_temporal_contour_records"; 
    case FILL_DISPATCH_ARGS_TEMPORAL_RECORDS:
      return "strokegen_fill_dispatch_args_per_temporal_record"; 

    case SCAN_TEST_AGGREGATE:
      return "bnpr_scan_test_aggregate";
    case SCAN_TEST_UPSWEEP:
      return "bnpr_scan_test_upsweep";
    case SCAN_TEST_DWSWEEP:
      return "bnpr_scan_test_dwsweep";

    case SCAN_UINT_ADD_UPSWEEP:
      return "bnpr_scan_uint_add_upsweep"; 
    case SCAN_UINT_ADD_AGGREGATE:
        return "bnpr_scan_uint_add_aggregate";
    case SCAN_UINT_ADD_DWSWEEP:
        return "bnpr_scan_uint_add_dwsweep";

    case SEGSCAN_TEST_UPSWEEP:
      return "bnpr_segscan_test_upsweep";
    case SEGSCAN_TEST_AGGREGATE:
      return "bnpr_segscan_test_aggregate";
    case SEGSCAN_TEST_DWSWEEP:
      return "bnpr_segscan_test_dwsweep";

    case SEGSCAN_UINT_ADD_UPSWEEP:
      return "bnpr_segscan_uint_add_upsweep";
    case SEGSCAN_UINT_ADD_AGGREGATE:
      return "bnpr_segscan_uint_add_aggregate";
    case SEGSCAN_UINT_ADD_DWSWEEP:
      return "bnpr_segscan_uint_add_dwsweep";

    case SEGSCAN_FLOAT_ADD_UPSWEEP:
      return "bnpr_segscan_float_add_upsweep";
    case SEGSCAN_FLOAT_ADD_AGGREGATE:
      return "bnpr_segscan_float_add_aggregate";
    case SEGSCAN_FLOAT_ADD_DWSWEEP:
      return "bnpr_segscan_float_add_dwsweep";

    case SEGSCAN_UINT_MIN_UPSWEEP:
      return "bnpr_segscan_uint_min_upsweep";
    case SEGSCAN_UINT_MIN_AGGREGATE:
      return "bnpr_segscan_uint_min_aggregate";
    case SEGSCAN_UINT_MIN_DWSWEEP:
      return "bnpr_segscan_uint_min_dwsweep"; 

    case SCAN_FILL_DISPTACH_ARGS:
      return "strokegen_scan_fill_dispatch_args";

    case CONV1D_FILL_DISPATCH_ARGS:
      return "strokegen_segloopconv1d_fill_dispatch_args";
    case CONV1D_TEST_BUILD_PATCH:
      return "strokegen_segloopconv1D_test_build_patch";
    case CONV1D_TEST_CONVOLUTION:
      return "strokegen_segloopconv1D_test_convolution";
    case CONV1D_SEG_DENOISE_BUILD_PATCH:
      return "strokegen_segloopconv1D_seg_denoising_build_patch";
    case CONV1D_SEG_DENOISE_CONVOLUTION:
      return "strokegen_segloopconv1D_seg_denoising_convolution";
    case CONV1D_2D_SAMPLE_BUILD_PATCH:
        return "strokegen_segloopconv1D_samples_2d_build_patch"; 
    case CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_0:
      return "strokegen_segloopconv1D_corner_detection_step_0_convolution";
    case CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_1:
      return "strokegen_segloopconv1D_corner_detection_step_1_convolution";
    case CONV1D_2D_SAMPLE_CALC_TANGENT_CURVATURE:
      return "strokegen_segloopconv1D_calc_2d_sample_tangent_curv_convolution";
    case CONV1D_2D_SAMPLE_DENOISE_VISIBILITY:
      return "strokegen_segloopconv1D_denoise_2d_sample_visibility_convolution";

    case LISTRANKING_INIT_ANCHORS:
      return "strokegen_list_ranking_test_tagging";
    case LISTRANKING_COMPACT_ANCHORS:
      return "strokegen_list_ranking_test_compact_anchors";
    case LISTRANKING_SPLICE_OUT_NODES:
      return "strokegen_list_ranking_test_splice_out_nodes";
    case LISTRANKING_FILL_DISPATCH_ARGS:
      return "strokegen_list_ranking_test_fill_dispatch_args";
    case LISTRANKING_SUBLIST_POINTER_JUMPING:
     return "strokegen_list_ranking_test_sublist_pointer_jumping";
    case LISTRANKING_LOOPED_POINTER_JUMPING:
     return "strokegen_list_ranking_test_looped_pointer_jumping";
    case LISTRANKING_MARK_LOOP_HEAD_TAIL:
      return "strokegen_list_ranking_test_mark_loop_head_tail";
    case LISTRANKING_RELINKING:
     return "strokegen_list_ranking_test_relinking";
    case LISTRANKING_LOOPED_RELINKING:
     return "strokegen_list_ranking_test_looped_relinking";
    case LISTRANKING_SETUP_INPUTS:
     return "strokegen_list_ranking_setup_input_data";
    case LISTRANKING_OUTPUT_DATA_PASS_0:
      return "strokegen_list_ranking_test_output_pass_0";
    case LISTRANKING_OUTPUT_DATA_PASS_1:
      return "strokegen_list_ranking_test_output_pass_1";

    /* To avoid compiler warning about missing case. */
    case MAX_SHADER_TYPE:
      return "";
  }
  return "";
}

GPUShader *StrokeGenShaderModule::static_shader_get(eShaderType shader_type)
{
  if (shaders_[shader_type] == nullptr) {
    const char *shader_name = static_shader_create_info_name_get(shader_type);

    shaders_[shader_type] = GPU_shader_create_from_info_name(shader_name);

    if (shaders_[shader_type] == nullptr) {
      fprintf(stderr, "bnpr: error: Could not compile static shader \"%s\"\n", shader_name);
    }
    BLI_assert(shaders_[shader_type] != nullptr);
  }
  return shaders_[shader_type];
}


/** \} */

}  // namespace blender::bnpr
