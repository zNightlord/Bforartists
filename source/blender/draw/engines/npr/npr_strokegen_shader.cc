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
    case COMPUTE_TEST:
      return "bnpr_strokegen_test_xxx";

    case GPU_MESHING_BOOSTRAP:
      return "bnpr_geom_extract_boostrap"; 
    case CONTOUR_GEOM_EXTRACT:
      return "bnpr_geom_extract";
    case CONTOUR_GEOM_EXTRACT_IBO_16BIT:
      return "bnpr_geom_extract_ibo_16bits";
    case FILL_DRAW_ARGS_CONTOUR_EDGES:
      return "bnpr_geom_fill_draw_args_contour_edges";
    case FILL_DISPATCH_ARGS_CONTOUR_EDGES:
      return "strokegen_fill_dispatch_args_per_contour_edge";
    case COMPUTE_CONTOUR_EDGE_RASTER_DATA:
      return "bnpr_geom_extract_calc_contour_edge_raster_data"; 
    case INDIRECT_DRAW_CONTOUR_EDGES:
      return "bnpr_geom_draw_contour_edges";


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
    case MESH_WEDGE_FLOODING_SELECT_VERTS:
      return "bnpr_meshing_compact_filtered_verts";

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

    case MESH_FILTERING_EDGE_NORMAL:
      return "bnpr_meshing_mesh_filtering_edge_normal_"; 
    case MESH_FILTERING_EDGE_QUADRIC:
      return "bnpr_meshing_mesh_filtering_edge_quadric_";
    case MESH_FILTERING_VERT_QUADRIC_INIT:
      return "bnpr_meshing_mesh_filtering_init_vert_quadric_"; 
    case MESH_FILTERING_VERT_QUADRIC_DIFFUSION:
      return "bnpr_meshing_mesh_filtering_diffuse_vert_quadric_"; 
    case MESH_FILTERING_MOVE_VERTS:
      return "bnpr_meshing_mesh_filtering_move_verts_";


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

    case MESH_ANALYSE_VERT_NORMAL:
      return "bnpr_geom_analysis_order_0_vert_normal";
    case MESH_ANALYSE_VERT_NORMAL_VORONOIAREA:
      return "bnpr_geom_analysis_order_0_vert_all"; 

    case FILL_DRAW_ARGS_DBG_LINES:
      return "strokegen_remeshing_fill_draw_args_dbg_lines"; 
    case INDIRECT_DRAW_DBG_LINES:
      return "bnpr_geom_draw_debug_lines"; 

    case MESH_ANALYSE_VERT_CURV_PASS_0:
      return "bnpr_geom_analysis_order_1_vert_curv_pass_0";
    case MESH_ANALYSE_VERT_CURV_PASS_1:
      return "bnpr_geom_analysis_order_1_main_curvature"; 


    case SCAN_TEST_AGGREGATE:
      return "bnpr_scan_test_aggregate";
    case SCAN_TEST_UPSWEEP:
      return "bnpr_scan_test_upsweep";
    case SCAN_TEST_DWSWEEP:
      return "bnpr_scan_test_dwsweep";
    case SEGSCAN_TEST_UPSWEEP:
      return "bnpr_segscan_test_upsweep";
    case SEGSCAN_TEST_AGGREGATE:
      return "bnpr_segscan_test_aggregate";
    case SEGSCAN_TEST_DWSWEEP:
      return "bnpr_segscan_test_dwsweep";
    case CONV1D_TEST_BUILD_PATCH:
      return "strokegen_segloopconv1D_test_build_patch";
    case CONV1D_TEST_CONVOLUTION:
      return "strokegen_segloopconv1D_test_convolution";
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
    case LISTRANKING_OUTPUT_DATA:
      return "strokegen_list_ranking_test_output"; 

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
