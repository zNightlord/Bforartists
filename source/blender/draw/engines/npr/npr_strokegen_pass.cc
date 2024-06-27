/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 */

#include "npr_strokegen_pass.hh"

#include "DEG_depsgraph_query.hh"
#include "gpu_storage_buffer_private.hh"
#include "NOD_geometry_exec.hh"

#include <iomanip>

namespace blender::npr::strokegen
{
  using namespace blender;

  PassSimple& StrokeGenPassModule::get_compute_pass(eType passType, int pass_id)
  {
    switch (passType) {
      case SCAN_TEST:
        return pass_scan_test;
      case SEGSCAN_TEST:
        return pass_segscan_test;
      case SEGLOOPCONV_TEST:
        return pass_conv1d_test;
      case LIST_RANKING_TEST:
        return pass_listranking_test;

      case GEOM_EXTRACTION:
        return pass_extract_geom_arr[pass_id];
      case CONTOUR_PROCESS:
        return pass_process_contours; 
      case COMPRESS_CONTOUR_PIXELS:
        return pass_compress_contour_pixels; 
    }
    return pass_comp_test;
  }

  PassMain &StrokeGenPassModule::get_render_pass(eType passType, int pass_id)
  {
    switch (passType) {
      case INDIRECT_DRAW_CONTOUR_EDGES:
        return pass_draw_contour_edges;
      case INDIRECT_DRAW_CONTOUR_2D_SAMPLES:
        return pass_draw_contour_2d_samples; 
      case INDIRECT_DRAW_REMESHED_DEPTH:
        return pass_draw_remeshed_surface_depth_[pass_id]; 
      case INDIRECT_DRAW_DBG_VNOR:
        return pass_draw_debug_lines_; 
    }
    return pass_draw_contour_edges; 
  }


  void StrokeGenPassModule::prepare_validation_passes(int frame_counter)
  {
    ScanSettings scan_test_settings;
    scan_test_settings.is_validation_shader = true;
    scan_test_settings.frame_counter = frame_counter;
    scan_test_settings.use_indirect_dispatch = false;
    scan_test_settings.ssbo_scan_infos_ = nullptr;
    scan_test_settings.ssbo_in_scan_data_ = buffers_.ssbo_in_scan_data_;
    scan_test_settings.ssbo_out_scan_data_ = buffers_.ssbo_out_scan_data_;
    scan_test_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
    scan_test_settings.shader_upsweep = eShaderType::SCAN_TEST_UPSWEEP;
    scan_test_settings.shader_aggregate = eShaderType::SCAN_TEST_AGGREGATE;
    scan_test_settings.shader_dwsweep = eShaderType::SCAN_TEST_DWSWEEP; 
    append_subpass_scan(scan_test_settings, pass_scan_test);

    ScanSettings segscan_test_settings;
    segscan_test_settings.is_validation_shader = true; 
    segscan_test_settings.frame_counter = frame_counter; 
    segscan_test_settings.use_indirect_dispatch = false;
    segscan_test_settings.ssbo_scan_infos_ = nullptr; 
    segscan_test_settings.ssbo_in_scan_data_ = buffers_.ssbo_in_scan_data_;
    segscan_test_settings.ssbo_out_scan_data_ = buffers_.ssbo_out_scan_data_;
    segscan_test_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
    segscan_test_settings.shader_upsweep = eShaderType::SEGSCAN_TEST_UPSWEEP;
    segscan_test_settings.shader_aggregate = eShaderType::SEGSCAN_TEST_AGGREGATE;
    segscan_test_settings.shader_dwsweep = eShaderType::SEGSCAN_TEST_DWSWEEP;
    append_subpass_segscan(segscan_test_settings, pass_segscan_test);

    SegLoopConv1DSettings segloopconv1d_test_settings;
    segloopconv1d_test_settings.is_validation_shader = true;
    segloopconv1d_test_settings.use_indirect_dispatch = false;
    segloopconv1d_test_settings.ssbo_segloopconv1d_info_ = buffers_.ssbo_segloopconv1d_info_;
    segloopconv1d_test_settings.ssbo_in_segloopconv1d_data_ = buffers_.ssbo_in_segloopconv1d_data_;
    segloopconv1d_test_settings.ssbo_out_segloopconv1d_data_ = buffers_.ssbo_out_segloopconv1d_data_;
    segloopconv1d_test_settings.ssbo_segloopconv1d_patch_table_ = buffers_.ssbo_segloopconv1d_patch_table_;
    segloopconv1d_test_settings.shader_build_patch_table = CONV1D_TEST_BUILD_PATCH;
    segloopconv1d_test_settings.lazy_dispatch = false;
    segloopconv1d_test_settings.shader_convolution = CONV1D_TEST_CONVOLUTION; 
    append_subpass_segloopconv1d(segloopconv1d_test_settings, pass_conv1d_test);

    append_subpass_list_ranking(ListRankingPassUsage::TestListRanking, pass_listranking_test, true);
  }

  void StrokeGenPassModule::on_begin_sync(int frame_counter)
  {
    init_mesh_extraction_passes();
    pass_draw_contour_edges.init_pass(shaders_, textures_, StrokegenMeshRasterPass::DRAW_CONTOUR_EDGES);
    pass_draw_contour_2d_samples.init_pass(shaders_, textures_, StrokegenMeshRasterPass::DRAW_CONTOUR_2D_SAMPLES);

    init_surface_depth_passes();

    pass_draw_debug_lines_.init_pass(shaders_, textures_, StrokegenMeshRasterPass::DBG_LINES); 

    prepare_validation_passes(frame_counter);

    num_total_mesh_tris = num_total_mesh_verts = num_total_mesh_edges = 0; // TODO: these should go to the UBO

    // fetch ui inputs
    const DRWContextState *draw_ctx = DRW_context_state_get();
    const Scene *scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);

    meshing_params.max_num_remesh_dbg_iters = 100000; // scene_eval->npr.npr_test_val_0;

    int val_1 = scene_eval->npr.npr_test_val_1; 
    surf_dbg_ctx.dbg_lines = (0 < val_1);
    surf_dbg_ctx.dbg_vert_normal = (1 == val_1);
    surf_dbg_ctx.dbg_vert_curv = (2 == val_1);
    surf_dbg_ctx.dbg_edges = (3 == val_1);
    if (val_1 == 4)
      surf_dbg_ctx.dbg_vert_normal = surf_dbg_ctx.dbg_vert_curv = true;
    if (val_1 == 5)
      surf_dbg_ctx.dbg_vert_normal = surf_dbg_ctx.dbg_edges = true;
    if (val_1 = 6)
      surf_dbg_ctx.dbg_vert_curv = surf_dbg_ctx.dbg_edges = true; 

    meshing_params.edge_visualize_mode = -1;
    if (surf_dbg_ctx.dbg_edges) {
      meshing_params.edge_visualize_mode = (int)(scene_eval->npr.npr_test_val_2 + 1e-10f); 
    }

    meshing_params.seconds_sync_view_mat = (int)(scene_eval->npr.npr_test_val_3 + 1e-10f);
    meshing_params.num_vtx_smooth_iters = (int)(scene_eval->npr.npr_test_val_4 + 1e-10f);

    meshing_params.contour_mode = (int)(scene_eval->npr.npr_test_val_5 + 1e-10f);
    meshing_params.visualize_contour_edges = 0 < meshing_params.contour_mode;
    meshing_params.iters_test_sqrt_subdiv = (int)(scene_eval->npr.npr_test_val_6 + 1e-10f); 

    pass_draw_contour_edges.draw_settings.draw_hidden_lines = scene_eval->npr.npr_test_val_7 > .5f;
    pass_draw_contour_2d_samples.draw_settings.draw_hidden_lines = scene_eval->npr.npr_test_val_7 > .5f; 
    pass_draw_debug_lines_.draw_settings.draw_hidden_lines = scene_eval->npr.npr_test_val_7 > .5f;

    meshing_params.num_quadric_diffusion_iters = scene_eval->npr.npr_test_val_8;
    meshing_params.quadric_deviation = scene_eval->npr.npr_test_val_9;
    meshing_params.position_regularization_scale = scene_eval->npr.npr_test_val_10;
    
    meshing_params.num_edge_flooding_iters = (int)(scene_eval->npr.npr_test_val_11 + 1e-10f);

    meshing_params.remeshing_targ_edge_len = scene_eval->npr.npr_test_val_12;
    meshing_params.remeshing_split_iters = (int)(scene_eval->npr.npr_test_val_13 + 1e-10f); 
    meshing_params.remeshing_collapse_iters = (int)(scene_eval->npr.npr_test_val_14 + 1e-10f);
    meshing_params.remeshing_flip_iters = (int)(scene_eval->npr.npr_test_val_15 + 1e-10f);
    meshing_params.remeshing_iters = (int)(scene_eval->npr.npr_test_val_16 + 1e-10f);
    meshing_params.remeshing_delaunay_flip_iters = (int)(scene_eval->npr.npr_test_val_17 + 1e-10f);

    surf_dbg_ctx.dbg_line_length = scene_eval->npr.npr_test_val_18;
    meshing_params.subdiv_type = (int)(scene_eval->npr.npr_test_val_19 + 1e-10f);
    meshing_params.subdiv_use_crease = .0f < scene_eval->npr.npr_test_val_20; 

    meshing_params.denoise_cusp_segmentation = .0f < scene_eval->npr.npr_test_val_21; 
    meshing_params.visibility_thresh = scene_eval->npr.npr_test_val_22; 
  }

  void StrokeGenPassModule::on_end_sync()
  {
    /* Debug draw */
    rebuild_pass_dbg_geom_drawcall(surf_dbg_ctx);

    /* Post processing after iterated over all meshes. */
    rebuild_pass_process_contours();
    if (meshing_params.visualize_contour_edges) { // I need to hide contour edges for better debug visuals, for now
      rebuild_pass_contour_edge_drawcall();
      rebuild_pass_contour_2d_samples_drawcall(); 
      rebuild_pass_compress_contour_pixels();  
    }
  }


  /**
   * \brief In each render loop, re-initialize the compute pass for geometry extraction.
   */
  void StrokeGenPassModule::init_mesh_extraction_passes()
  {
    curr_mesh_id_extract_geom = -1; 

    for (StrokegenMeshComputePass& pass_extract_geom : pass_extract_geom_arr) {
      pass_extract_geom.init(); 
    }
    boostrap_before_extract_first_batch = true;

    surf_dbg_ctx.dbg_line_length = 1.0f;
  }

  void StrokeGenPassModule::init_surface_depth_passes()
  {
    curr_mesh_id_surf_depth = -1;
    for (StrokegenMeshRasterPass &surf_depth_pass : pass_draw_remeshed_surface_depth_) {
      surf_depth_pass.init_pass(
          shaders_, textures_, StrokegenMeshRasterPass::REMESHED_SURFACE_DEPTH);
    }
  }

  void StrokeGenPassModule::GetSurfaceAnalysisContext_InitPass(SurfaceAnalysisContext &surf_analysis_ctx) const
  {
    surf_analysis_ctx.set_calc_feature_edges(true, false); 

    surf_analysis_ctx.order_0_only_selected = false;
    surf_analysis_ctx.set_calc_vert_normal(true, false);
    surf_analysis_ctx.ssbo_vnor_ = buffers_.ssbo_vnor_;
    surf_analysis_ctx.set_calc_vert_voronoi_area(true);
    surf_analysis_ctx.ssbo_varea_ = buffers_.ssbo_mesh_buffer_reuse_5_;
    surf_analysis_ctx.set_calc_vert_topo_flags(true); 

    surf_analysis_ctx.order_1_only_selected = false;
    surf_analysis_ctx.set_calc_vert_curvature(true,
                                              SurfaceAnalysisContext::CurvatureEstimator::Jacques,
                                              false, false);
    surf_analysis_ctx.ssbo_edge_vtensors_ = buffers_.ssbo_mesh_buffer_reuse_7_;
    surf_analysis_ctx.ssbo_vcurv_tensor_ = buffers_.ssbo_mesh_buffer_reuse_1_;
    surf_analysis_ctx.ssbo_vcurv_pdirs_k1k2_ = buffers_.ssbo_mesh_buffer_reuse_2_;

  }

  void StrokeGenPassModule::GetSurfaceAnalysisContext_CurvatureForAdaptiveRemeshing(
      SurfaceAnalysisContext &surf_analysis_ctx) const
  {
    surf_analysis_ctx.set_calc_feature_edges(true, true); 

    surf_analysis_ctx.order_0_only_selected = false;
    surf_analysis_ctx.set_calc_vert_normal(true, false);
    surf_analysis_ctx.ssbo_vnor_ = buffers_.ssbo_vnor_;
    surf_analysis_ctx.set_calc_vert_voronoi_area(true);
    surf_analysis_ctx.ssbo_varea_ = buffers_.ssbo_mesh_buffer_reuse_5_;
    surf_analysis_ctx.set_calc_vert_topo_flags(true); 

    surf_analysis_ctx.order_1_only_selected = true;
    surf_analysis_ctx.set_calc_vert_curvature(true,
                                              SurfaceAnalysisContext::CurvatureEstimator::Jacques,
                                              false, false);
    surf_analysis_ctx.ssbo_edge_vtensors_ = buffers_.ssbo_mesh_buffer_reuse_7_;
    surf_analysis_ctx.ssbo_vcurv_tensor_ = buffers_.ssbo_mesh_buffer_reuse_1_;
    surf_analysis_ctx.ssbo_vcurv_pdirs_k1k2_ = buffers_.ssbo_mesh_buffer_reuse_2_;
  }

  void StrokeGenPassModule::GetSurfaceAnalysisContext_VertexRelocationPass(
      SurfaceAnalysisContext &surf_analysis_ctx) const
  {
    surf_analysis_ctx.order_0_only_selected = false;
    surf_analysis_ctx.set_calc_vert_normal(true, false);
    surf_analysis_ctx.ssbo_vnor_ = buffers_.ssbo_vnor_;
    surf_analysis_ctx.set_calc_vert_voronoi_area(false);
    surf_analysis_ctx.ssbo_varea_ = buffers_.ssbo_mesh_buffer_reuse_5_;
    surf_analysis_ctx.set_calc_vert_topo_flags(false); 

    surf_analysis_ctx.order_1_only_selected = false;
    surf_analysis_ctx.set_calc_vert_curvature(false,
                                              SurfaceAnalysisContext::CurvatureEstimator::Jacques,
                                  false, false);
    surf_analysis_ctx.ssbo_edge_vtensors_ = buffers_.ssbo_mesh_buffer_reuse_7_;
    surf_analysis_ctx.ssbo_vcurv_tensor_ = buffers_.ssbo_mesh_buffer_reuse_1_;
    surf_analysis_ctx.ssbo_vcurv_pdirs_k1k2_ = buffers_.ssbo_mesh_buffer_reuse_2_;
  }

  void StrokeGenPassModule::GetSurfaceAnalysisContext_ContourInsertionPass(SurfaceAnalysisContext &surf_analysis_ctx) const
  {
    surf_analysis_ctx.order_0_only_selected = false;
    surf_analysis_ctx.set_calc_vert_normal(true, true);
    surf_analysis_ctx.ssbo_vnor_ = buffers_.ssbo_vnor_;
    surf_analysis_ctx.set_calc_vert_voronoi_area(false);
    surf_analysis_ctx.ssbo_varea_ = buffers_.ssbo_mesh_buffer_reuse_5_;
    surf_analysis_ctx.set_calc_vert_topo_flags(true); 

    surf_analysis_ctx.order_1_only_selected = false;
    surf_analysis_ctx.set_calc_vert_curvature(false, SurfaceAnalysisContext::CurvatureEstimator::Jacques, false, false);
    surf_analysis_ctx.ssbo_edge_vtensors_ = buffers_.ssbo_mesh_buffer_reuse_7_;
    surf_analysis_ctx.ssbo_vcurv_tensor_ = buffers_.ssbo_mesh_buffer_reuse_1_;
    surf_analysis_ctx.ssbo_vcurv_pdirs_k1k2_ = buffers_.ssbo_mesh_buffer_reuse_2_;
  }

  void StrokeGenPassModule::GetSurfaceAnalysisContext_CuspDetectionPass(
      SurfaceAnalysisContext &surf_analysis_ctx) const
  {
    surf_analysis_ctx.order_0_only_selected = false;
    surf_analysis_ctx.set_calc_vert_normal(true, false);
    surf_analysis_ctx.ssbo_vnor_ = buffers_.ssbo_vnor_;
    surf_analysis_ctx.set_calc_vert_voronoi_area(false);
    surf_analysis_ctx.ssbo_varea_ = buffers_.ssbo_mesh_buffer_reuse_5_;
    surf_analysis_ctx.set_calc_vert_topo_flags(false); 

    surf_analysis_ctx.order_1_only_selected = true;
    surf_analysis_ctx.set_calc_vert_curvature(true,
                                              SurfaceAnalysisContext::CurvatureEstimator::Jacques,
                                              false,
                                              true);
    surf_analysis_ctx.ssbo_edge_vtensors_ = buffers_.ssbo_mesh_buffer_reuse_7_;
    surf_analysis_ctx.ssbo_vcurv_tensor_ = buffers_.ssbo_mesh_buffer_reuse_1_;
    surf_analysis_ctx.ssbo_vcurv_pdirs_k1k2_ = buffers_.ssbo_mesh_buffer_reuse_2_;
  }

  void StrokeGenPassModule::append_subpasses_estimate_curvature_for_adaptive_remeshing(
    ResourceHandle &rsc_handle, int num_edges, int num_verts, bool output_dbg_lines)
  {
    const bool debug_remesh_edge_len = false; 

    SurfaceAnalysisContext surf_analysis_ctx;
    GetSurfaceAnalysisContext_CurvatureForAdaptiveRemeshing(surf_analysis_ctx);
    auto surf_dbg_ctx_cpy = surf_dbg_ctx;
    surf_dbg_ctx_cpy.dbg_vert_curv = debug_remesh_edge_len ? false : output_dbg_lines; 
    append_subpass_surf_geom_analysis(
        rsc_handle, num_verts, num_edges, surf_analysis_ctx, surf_dbg_ctx_cpy);

    // Note: in fact, this smoothing is not very necessary.
    // we don't really have that much of nan/inf curvatures
    // but I still do this to ensure the curvature is robust
    int num_steps_vcurv_smooth = 2 * 1 /*must be even*/; 
    for (int step_vcurv_smooth = 0; step_vcurv_smooth < num_steps_vcurv_smooth; ++step_vcurv_smooth) {
      append_subpass_vertex_curv_smoothing(
          step_vcurv_smooth, num_verts, num_edges,
          false, step_vcurv_smooth == num_steps_vcurv_smooth - 1
        ); 
    }

    if (debug_remesh_edge_len && output_dbg_lines)
    {
      SurfaceAnalysisContext surf_ctx_cpy = surf_analysis_ctx;
      surf_ctx_cpy.ssbo_vcurv_pdirs_k1k2_ = buffers_.reused_ssbo_vtx_remesh_len_();
      surf_ctx_cpy.output_curvature_tensors = false; 
      // surf_ctx_cpy.output_vertex_facing_flag = true; 

      surf_dbg_ctx_cpy = surf_dbg_ctx;
      surf_dbg_ctx_cpy.dbg_vert_curv = true;
      
      append_subpass_surf_geom_analysis(
          rsc_handle, num_verts, num_edges, surf_ctx_cpy, surf_dbg_ctx_cpy);  
    }
  }

  void StrokeGenPassModule::append_subpasses_sqrt_subdiv(int num_edges, int num_verts)
  {
    for (int iter_subdiv = 0; iter_subdiv < meshing_params.iters_test_sqrt_subdiv; ++iter_subdiv)
    {
      append_subpass_vertex_relocation(Sqrt3SubdivSmoothCache, num_edges, num_verts, 0, true);

      append_subpass_split_faces(0, num_edges, num_verts);
      // update elem counters after face split
      append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
      append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);
      
      append_subpass_flip_edges(SqrtSubdiv, 0, num_edges, num_verts);

      append_subpass_vertex_relocation(Sqrt3SubdivSmoothApply, num_edges, num_verts, 0, true);
    }
  }

  void StrokeGenPassModule::append_subpasses_loop_subdiv(int num_edges, int num_verts)
  {
    for (int iter_subdiv = 0; iter_subdiv < meshing_params.iters_test_sqrt_subdiv; ++iter_subdiv)
    {
      append_subpass_vertex_relocation(LoopSubdivSmoothCache, num_edges, num_verts, 0, true);

      /* try to split each edge once, but not theoretically guaranteed,
       * some edge may not be split and hence fuck up the local subdivision */
      for (int iter_split = 0; iter_split < 6; ++iter_split)
      {
        append_subpass_split_edges(LoopSubdivSplit, iter_split, num_edges, num_verts);
        // update elem counters after split
        append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
        append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);
      }  

      append_subpass_flip_edges(LoopSubdivFlip, 0, num_edges, num_verts);

      append_subpass_vertex_relocation(LoopSubdivSmoothApply, num_edges, num_verts, 0, true);
    }
  }

  /**
   * \brief Add a subpass for extracting geometry from given GPUBatch.
   * \param ob Mesh Object
   * \param gpu_batch_line_adj Mesh geometry stored in GPUBatch, ib stored with line adjacency info.
   */
  void StrokeGenPassModule::append_per_mesh_pass(
      Object* ob,
      gpu::Batch* gpu_batch_line_adj,
      gpu::Batch* gpu_batch_surf,
      ResourceHandle& rsc_handle,
      const DRWView* drw_view
  ){
    curr_mesh_id_extract_geom++; 

    if (boostrap_before_extract_first_batch) {
      auto &sub = pass_extract_geom().sub("bootstrap meshing passes");
      sub.shader_set(shaders_.static_shader_get(GPU_MESHING_BOOSTRAP));

      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      sub.dispatch(int3(1, 1, 1));

      curr_mesh_id_extract_geom = -1; // reset counter

      return; // bootstrapping done. will re-enter this func for actual work
    }

    pass_extract_geom().debug_name = "extract mesh pass"; 

    int batch_resource_index = (int)rsc_handle.resource_index();

    /* see example in "GPU_batch_draw_parameter_get" */
    gpu::Batch *edge_batch = static_cast<gpu::Batch *>(gpu_batch_line_adj);
    if (edge_batch == nullptr || edge_batch->elem == nullptr)
      fprintf(stderr, "StrokeGen Error: empty mesh when rebuilding pass 'extract_mesh_geom'");
    int num_edges = edge_batch->elem_()->index_len_get() / 4;
    // 4 indices per primitive, basically 4 verts around the edge

    gpu::IndexBuf *ib = edge_batch->elem_();
    if (ib->index_len_get() == 0)
      return;
    /* Hack to get ibo format */
    gpu::GPUIndexBufType ib_type = // the actual "index_type_" is protected
        (ib->size_get() / ib->index_len_get() == sizeof(uint32_t)) ?
          gpu::GPU_INDEX_U32 :
          gpu::GPU_INDEX_U16;

    gpu::Batch *batch_surf = static_cast<gpu::Batch *>(gpu_batch_surf);
    gpu::IndexBuf *ibo_surf = batch_surf->elem_();
    /* Layout of multi vbos, see
     * "DRW_batch_requested(cache->batch.surface, GPU_PRIM_TRIS)" in
     * "DRW_mesh_batch_cache_create_requested"
     *
     * For the pos vbo of a Surface GPUBatch
     * See mesh_extractors\extract_mesh_vbo_pos.cc
    */
    gpu::VertBuf *vbo_surf_pos = batch_surf->verts_(1);
    int num_tris = ibo_surf->index_len_get() / 3;
    int num_verts = vbo_surf_pos->vertex_len;



    // Copy VBO/IBO ----------------------------------------------------------
    // Note: this also works
    /* GPU_storagebuf_copy_sub_from_vertbuf(
     *  buffers_.ssbo_vbo_full_, gpu_batch_surf->verts[0], 4 * num_total_mesh_verts, 0, 4 * num_verts
     * ); */
    append_subpass_cpy_vbo(gpu_batch_surf, batch_resource_index, num_verts);
    append_subpass_cpy_line_adj_ibo(gpu_batch_line_adj, ib_type, num_edges);

    // Basic Meshing -----------------------------------------------------------
    append_subpass_meshing_merge_verts(num_verts);
    append_subpass_meshing_wedge_adjacency(num_edges, num_verts);


    // Mesh Selection -----------------------------------------------------------
    EdgeFloodingOptions flooding_options;
    flooding_options.compact_edges = true;
    flooding_options.output_selected_to_edge = true;
    flooding_options.output_edge_to_selected = false;
    append_subpass_select_remeshed_edges(num_edges, num_verts, flooding_options);
    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);

    append_subpass_mark_selection_border_edges(num_edges, num_verts); 

    SelectVertsFromEdgesContext vtx_sel_ctx_flooding;
    vtx_sel_ctx_flooding.active_selection_slots = {0, -1, -1, -1}; // select flooded verts at slot#0
    vtx_sel_ctx_flooding.expand_selection = true;
    vtx_sel_ctx_flooding.input_selection_slots_for_expansion = {0, -1, -1, -1};
    vtx_sel_ctx_flooding.output_selection_slots_for_expansion = { 1, -1, -1, -1};
    vtx_sel_ctx_flooding.compact_verts = true;
    vtx_sel_ctx_flooding.selection_slots_for_compation = {0, 1, -1, -1}; 
    vtx_sel_ctx_flooding.compact_all_slots_selected = false;

    SelectVertsFromEdgesContext vtx_sel_ctx_vnor = vtx_sel_ctx_flooding;
    append_subpass_select_verts_from_selected_edges(vtx_sel_ctx_vnor, num_edges, num_verts); 

    // Init vert curvatures -------------------------------------------------------------------
    SurfaceAnalysisContext surf_analysis_ctx_init;
    GetSurfaceAnalysisContext_InitPass(surf_analysis_ctx_init);
    auto surf_dbg_ctx_cpy = surf_dbg_ctx;
    surf_dbg_ctx_cpy.dbg_vert_normal = false;
    surf_dbg_ctx_cpy.dbg_vert_curv = false;
    append_subpass_surf_geom_analysis(
        rsc_handle, num_verts, num_edges, surf_analysis_ctx_init, surf_dbg_ctx_cpy);



    // GPU Remesher -------------------------------------------------------------------------
    int dbg_step = 0;
    bool step_dbg_remesh = true;
    auto should_remesh_when_dbg = [&]() {
      return ((!step_dbg_remesh) || (dbg_step < meshing_params.max_num_remesh_dbg_iters));
    }; 

    for (int iter_remesh = 0; iter_remesh < meshing_params.remeshing_iters; ++iter_remesh)
    {
      int num_edge_split_iters = meshing_params.remeshing_split_iters;
        for (int iter_edge_split = 0; iter_edge_split < num_edge_split_iters; ++iter_edge_split) {
          if (should_remesh_when_dbg()) {
            append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
            append_subpass_split_edges(LongEdge, iter_edge_split, num_edges, num_verts);
            dbg_step++;   
          }
        }
      // update elem counters after split
      append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true); 
      append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);


      EdgeFlipOptiGoal opti_goal = Delaunay;
      int num_edge_flip_iters = meshing_params.remeshing_delaunay_flip_iters;
        for (int iter_edge_flip = 0; iter_edge_flip < num_edge_flip_iters; ++iter_edge_flip) {
          if (should_remesh_when_dbg()) {
            append_subpass_flip_edges(opti_goal, iter_edge_flip, num_edges, num_verts);
            dbg_step++; 
          }
        }

      append_subpasses_estimate_curvature_for_adaptive_remeshing(
          rsc_handle,
          num_edges,
          num_verts,
          false /*iter_remesh == meshing_params.remeshing_iters - 1*/);

      int num_edge_collapse_iters = meshing_params.remeshing_collapse_iters;
        for (int iter_edge_collapse = 0; iter_edge_collapse < num_edge_collapse_iters;
             ++iter_edge_collapse) {
          if (should_remesh_when_dbg()) {
              append_subpass_collapse_edges(iter_remesh, iter_edge_collapse, num_edges, num_verts);
              dbg_step++;   
          }
        }

      opti_goal = Valence; 
      num_edge_flip_iters = meshing_params.remeshing_flip_iters;
        for (int iter_edge_flip = 0; iter_edge_flip < num_edge_flip_iters; ++iter_edge_flip) {
          if (should_remesh_when_dbg()) {
              append_subpass_flip_edges(opti_goal, iter_edge_flip, num_edges, num_verts);
              dbg_step++; 
          }
        }


      { // Vertex Relocation
        VertexRelocationMode relocation_mode = TangentialSmoothing; // QuadricFiltering; 
        // vertex normal & curvature required
        SurfaceAnalysisContext surf_analysis_ctx_remesh;
        GetSurfaceAnalysisContext_VertexRelocationPass(surf_analysis_ctx_remesh);
        auto surf_dbg_ctx_cpy = surf_dbg_ctx;
        surf_dbg_ctx_cpy.dbg_vert_normal = false;
        surf_dbg_ctx_cpy.dbg_vert_curv = false;
        append_subpass_surf_geom_analysis(
            rsc_handle, num_verts, num_edges, surf_analysis_ctx_remesh, surf_dbg_ctx_cpy);
        append_subpass_vertex_relocation(
          relocation_mode,
          num_edges,
          num_verts,
          meshing_params.num_quadric_diffusion_iters,
          true
        );
      }
    }

    if (meshing_params.contour_mode != ContourType::Raw) {
      if (meshing_params.subdiv_type == 0)
        append_subpasses_sqrt_subdiv(num_edges, num_verts);
      else
        append_subpasses_loop_subdiv(num_edges, num_verts);
    }

    // test interpolated contour tessellation
    if (meshing_params.contour_mode != ContourType::Raw) {
      SurfaceAnalysisContext surf_analysis_ctx_contour;
      GetSurfaceAnalysisContext_ContourInsertionPass(surf_analysis_ctx_contour); 
      auto surf_dbg_ctx_cpy = surf_dbg_ctx;
      surf_dbg_ctx_cpy.dbg_vert_normal = false;
      surf_dbg_ctx_cpy.dbg_vert_curv = false;
      append_subpass_surf_geom_analysis(
          rsc_handle, num_verts, num_edges, surf_analysis_ctx_contour, surf_dbg_ctx_cpy);
      
      for (int iter_contour_insertion = 0; iter_contour_insertion < 6; ++iter_contour_insertion) {
        append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
        append_subpass_split_edges(InterpContour, iter_contour_insertion, num_edges, num_verts);
        // update elem counters after split
        append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
        append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);
      }
    }
    
    // test cusp detection
    {
      SurfaceAnalysisContext surf_analysis_ctx_contour;
      GetSurfaceAnalysisContext_CuspDetectionPass(surf_analysis_ctx_contour);
      auto surf_dbg_ctx_cpy = surf_dbg_ctx;
      surf_dbg_ctx_cpy.dbg_vert_normal = surf_dbg_ctx.dbg_vert_normal;
      surf_dbg_ctx_cpy.dbg_vert_curv = surf_dbg_ctx.dbg_vert_curv;
      append_subpass_surf_geom_analysis(rsc_handle, num_verts, num_edges, surf_analysis_ctx_contour, surf_dbg_ctx_cpy);
    }



    // Contour Processing --------------------------------------------------------
    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, false);
    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);
    append_subpass_extract_contour_edges(
        gpu_batch_line_adj,
        rsc_handle,
        edge_batch,
        num_edges,
        ib_type,
        meshing_params.edge_visualize_mode,
        meshing_params.contour_mode
      );

    append_subpass_fill_dispatch_args_contour_edges(pass_extract_geom(), false);
    append_subpass_setup_contour_edge_data();

    // Note: this should be  appening at the end
    num_total_mesh_tris += num_tris;
    num_total_mesh_verts += num_verts;
    num_total_mesh_edges += num_edges; 
  }


  void StrokeGenPassModule::bind_rsc_for_bnpr_meshing_surf_filtering_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf>& sub,
      int num_verts, int num_edges,
      int& out_num_ssbo)
  {
    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_dyn_mesh_counters_out_());
    sub.bind_ssbo(2, buffers_.ssbo_edge_to_edges_);
    sub.bind_ssbo(3, buffers_.ssbo_edge_to_vert_);
    sub.bind_ssbo(4, buffers_.ssbo_vert_to_edge_list_header_);
    sub.bind_ssbo(5, buffers_.ssbo_selected_edge_to_edge_);
    sub.bind_ssbo(6, buffers_.ssbo_selected_vert_to_vert_);
    sub.bind_ssbo(7, buffers_.ssbo_vbo_full_);
    sub.bind_ssbo(8, buffers_.ssbo_vnor_);
    sub.bind_ssbo(9, buffers_.ssbo_edge_flags_);
    sub.bind_ssbo(10, buffers_.ssbo_vert_flags_);
    out_num_ssbo = 11; 

    sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);

    sub.push_constant("pcs_edge_count_", num_edges);
    sub.push_constant("pcs_vert_count_", num_verts);
  }

  // TODO: Store quadrics for selected verts in a compacted manner. 
  void StrokeGenPassModule::append_subpass_vertex_relocation(
      VertexRelocationMode mode,
      int num_edges,
      int num_verts,
      int num_vnor_filter_iters,
      bool only_selected_verts)
  {
    int ssbo_offset = 0; 
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      bind_rsc_for_bnpr_meshing_surf_filtering_(
          sub, num_verts, num_edges, ssbo_offset); 
    };

    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, only_selected_verts);


    if (mode == TangentialSmoothing) {
      int step_vpos_filter = 0;
      int num_smooth_steps = ((meshing_params.num_vtx_smooth_iters + 1) / 2) * 2; 
      for (; step_vpos_filter < num_smooth_steps; ++step_vpos_filter)
      {
        {
          auto &sub = pass_extract_geom().sub("bnpr_meshing_surf_filtering_vpos_filtering");
          sub.shader_set(shaders_.static_shader_get(MESH_FILTER_VPOS_FILTERING));

          bind_src(sub);
          sub.bind_ssbo(ssbo_offset, buffers_.reused_ssbo_vpos_temp_());
          sub.bind_ssbo(ssbo_offset+1, buffers_.reused_ssbo_vtx_remesh_len_()); 
          sub.push_constant("pcs_vpos_filtering_iter_", step_vpos_filter);
          sub.push_constant("pcs_num_vpos_filtering_iters_", meshing_params.num_vtx_smooth_iters); 

          sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
          sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
        }
      }
    }


    if (mode == QuadricFiltering)
    { // Quadric Filtering
      int step_vpos_filter = 0;
      for (; step_vpos_filter < meshing_params.num_vtx_smooth_iters; ++step_vpos_filter) {
        {
          auto bind_src_vq = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
                                 int vq_filter_step) {
            sub.bind_ssbo(ssbo_offset, buffers_.reused_ssbo_vert_quadric_data_in_(vq_filter_step));
            sub.bind_ssbo(ssbo_offset+1, buffers_.reused_ssbo_vert_quadric_data_out_(vq_filter_step)); 
            sub.push_constant("pcs_vq_filtering_iter_", vq_filter_step);
            sub.push_constant("pcs_filtered_quadric_type_", (int)0);
            sub.push_constant("pcs_quadric_deviation_", meshing_params.quadric_deviation);
            sub.push_constant("pcs_geodist_deviation_", 1.0f);
            sub.push_constant("pcs_position_regularization_scale_", meshing_params.position_regularization_scale);
          };

          int num_vq_filter_iters = num_vnor_filter_iters;
          int step_vq_filter = 0;
          for (; step_vq_filter < num_vq_filter_iters; ++step_vq_filter) {
            auto &sub = pass_extract_geom().sub("bnpr_meshing_surf_filtering_vquadric_diffusion");
            sub.shader_set(shaders_.static_shader_get(MESH_FILTER_VQUADRIC_DIFFUSION));

            bind_src(sub);
            bind_src_vq(sub, step_vq_filter);

            sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
            sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
          }

          {
            auto &sub = pass_extract_geom().sub(
                "bnpr_meshing_surf_filtering_quadric_vpos_filtering");
            sub.shader_set(shaders_.static_shader_get(MESH_FILTER_QUADRIC_VPOS_FILTERING));

            bind_src(sub);
            bind_src_vq(sub, step_vq_filter);

            sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
            sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
          }
        }
      }
    }


    /* Vertex interpolations for Subdivision */
    const int SUBDIV_VPOS_FILTER_TYPE__SQRT3 = 0; // match shader macro
    const int SUBDIV_VPOS_FILTER_TYPE__LOOP = 1;  // match shader macro
    int subdiv_mode = SUBDIV_VPOS_FILTER_TYPE__SQRT3;
    if (mode == LoopSubdivSmoothCache || mode == LoopSubdivSmoothApply)
      subdiv_mode = SUBDIV_VPOS_FILTER_TYPE__LOOP;

    auto bind_subd_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> & sub, int subdType)
    {
      sub.bind_ssbo(ssbo_offset, buffers_.reused_ssbo_vpos_subd_());
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_epos_subd_());
      sub.push_constant("pcs_subdiv_type_", subdType);
      sub.push_constant("pcs_loop_subd_enable_crease_", meshing_params.subdiv_use_crease ? 1 : 0); 
    };

    // Smooth old vertices
    if (mode == Sqrt3SubdivSmoothCache || mode == LoopSubdivSmoothCache) 
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_surf_filtering_subdiv_vpos_smoothing");
      sub.shader_set(shaders_.static_shader_get(MESH_FILTER_SUBDIV_VPOS_SMOOTH));

      bind_src(sub);
      bind_subd_src(sub, subdiv_mode);

      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    if (mode == Sqrt3SubdivSmoothApply || mode == LoopSubdivSmoothApply) {
      auto &sub = pass_extract_geom().sub(
          "bnpr_meshing_surf_filtering_subdiv_vpos_smoothing_finish");
      sub.shader_set(shaders_.static_shader_get(MESH_FILTER_SUBDIV_VPOS_SMOOTH_FINISH));

      bind_src(sub);
      bind_subd_src(sub, subdiv_mode); 
    
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }

    /* For loop subd, we also need to cache new points on old edges */
    if (mode == LoopSubdivSmoothCache)
    {
      append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true); 

      {
        auto &sub = pass_extract_geom().sub(
           "bnpr_meshing_surf_filtering_subdiv_edge_points"
        );
        sub.shader_set(shaders_.static_shader_get(MESH_FILTER_LOOP_SUBDIV_EDGE_POINTS));

        bind_src(sub);
        bind_subd_src(sub, subdiv_mode);

        sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
        sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
      }
    }
  }

  void StrokeGenPassModule::append_subpass_vertex_curv_smoothing(
    int iter_smooth, int num_verts, int num_edges, bool only_selected_verts, bool output_remesh_len)
  {
    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, only_selected_verts);

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_surf_filtering_vcurv_smoothing");
      if (false == output_remesh_len)
        sub.shader_set(shaders_.static_shader_get(MESH_FILTER_VCURV_SMOOTHING));
      else
        sub.shader_set(shaders_.static_shader_get(MESH_FILTER_VCURV_SMOOTHING_OUTPUT_REMESH_LEN)); 

      int ssbo_offset = 0; 
      bind_rsc_for_bnpr_meshing_surf_filtering_(sub, num_verts, num_edges, ssbo_offset);
      sub.bind_ssbo(ssbo_offset, buffers_.ssbo_vcurv_max_);
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_vcurv_max_temp_());
      sub.bind_ssbo(ssbo_offset + 2, buffers_.reused_ssbo_vtx_remesh_len_());
      sub.push_constant("pcs_vcurv_smooth_iter_", iter_smooth); 

      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_mark_selection_border_edges(int num_edges, int num_verts)
  {
    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(1, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(2, buffers_.ssbo_selected_edge_to_edge_);
      sub.bind_ssbo(3, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(4, buffers_.ssbo_dyn_mesh_counters_out_());
      sub.bind_ssbo(5, buffers_.ssbo_edge_flags_);
      
      sub.push_constant("pcs_vert_count_", num_verts);
      sub.push_constant("pcs_edge_count_", num_edges);
    };

    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true); 
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_selection_mark_selection_border");
      sub.shader_set(shaders_.static_shader_get(MESH_WEDGE_MARK_EDGES_ON_SELECTION_BORDER)); 
      bind_rsc(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
  }


  void StrokeGenPassModule::append_subpass_select_verts_from_selected_edges(
      SelectVertsFromEdgesContext ctx, int num_edges, int num_verts)
  {
    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub, int4 active_sel_slots) {
      sub.bind_ssbo(0, buffers_.ssbo_edge_to_edges_);            
      sub.bind_ssbo(1, buffers_.ssbo_edge_to_vert_);             
      sub.bind_ssbo(2, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(3, buffers_.ssbo_edge_flags_);               
      sub.bind_ssbo(4, buffers_.ssbo_vert_flags_);               
      sub.bind_ssbo(5, buffers_.ssbo_selected_edge_to_edge_);
      sub.bind_ssbo(6, buffers_.ssbo_selected_vert_to_vert_);
      sub.bind_ssbo(7, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(8, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(9, buffers_.ssbo_dyn_mesh_counters_out_());

      sub.push_constant("pcs_vert_count_", num_verts);
      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_vertex_selection_slots_", active_sel_slots); 
    }; 

    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, true);
    {
      auto &sub = pass_extract_geom().sub("strokegen_select_verts_from_selected_edges");
      sub.shader_set(shaders_.static_shader_get(MESH_SELECT_VERTS_FROM_SELECTED_EDGES));
      bind_rsc(sub, ctx.active_selection_slots);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_); 
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }

    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, false);
    if (ctx.expand_selection) {
      auto &sub = pass_extract_geom().sub("strokegen_expand_verts_from_selected_edges");
      sub.shader_set(shaders_.static_shader_get(MESH_EXPAND_VERTS_FROM_SELECTED_EDGES));
      bind_rsc(sub, ctx.input_selection_slots_for_expansion);
      sub.push_constant("pcs_vertex_selection_slots_out_",
                        ctx.output_selection_slots_for_expansion); 
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
    
    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, false);
    if (ctx.compact_verts) {
      auto &sub = pass_extract_geom().sub("strokegen_compact_selected_verts");
      sub.shader_set(shaders_.static_shader_get(MESH_COMPACT_SELECTED_VERTS));

      bind_rsc(sub, ctx.selection_slots_for_compation);
      sub.push_constant("pcs_vertex_select_all_slots_", ctx.compact_all_slots_selected ? 1 : 0); 

      sub.dispatch(/*buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_*/int3(compute_num_groups(num_verts, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT), 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
  }

  void StrokeGenPassModule::append_subpass_cpy_vbo(gpu::Batch *gpu_batch_surf, int batch_resource_index, int num_verts)
  {
    auto &sub = pass_extract_geom().sub("bnpr_geom_extract_collect_verts");
      
    sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_COLLECT_VBO)); 
      
    sub.bind_ssbo(0, &(gpu_batch_surf->elem));  // TODO: investigate whether double pointer is necessary
    sub.bind_ssbo(1, &(gpu_batch_surf->verts[1])); 
    sub.bind_ssbo(2, buffers_.ssbo_vbo_full_);
    sub.bind_ssbo(3, DRW_manager_get()->matrix_buf.current());
    sub.bind_ssbo(4, buffers_.ssbo_bnpr_mesh_pool_counters_prev_); 
    sub.bind_ssbo(5, buffers_.ssbo_bnpr_mesh_pool_counters_); 
      
    sub.bind_ubo(0, buffers_.ubo_view_matrices_); 
    sub.push_constant("pcs_rsc_handle_", batch_resource_index); 
    sub.push_constant("pcs_meshbatch_num_verts_", num_verts); 
      
    int num_groups = compute_num_groups(num_verts, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT); 
    sub.dispatch(int3(num_groups, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_cpy_line_adj_ibo(
      gpu::Batch *gpu_batch_line_adj,
      gpu::GPUIndexBufType ib_type,
      int num_added_edges)
  {
    auto &sub = pass_extract_geom().sub("bnpr_geom_extract_collect_edges");

    if (ib_type == gpu::GPU_INDEX_U16)
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_COLLECT_EDGE_ADJ_IBO_16BIT));
    else
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_COLLECT_EDGE_ADJ_IBO));

    sub.bind_ssbo(0, &(gpu_batch_line_adj->elem));
    sub.bind_ssbo(1, buffers_.ssbo_edge_to_vert_);

    sub.push_constant("pcs_edge_count_", num_added_edges);

    int num_groups = compute_num_groups(num_added_edges, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
    sub.dispatch(int3(num_groups, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_meshing_merge_verts(int num_verts_in, bool debug)
  {
    const int hashmap_size = std::min(MAX_GPU_HASH_TABLE_SIZE, num_verts_in * 16);

    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub)
    {
      sub.bind_ssbo(0, buffers_.reused_ssbo_vert_spatial_map_headers_());
      sub.bind_ssbo(1, buffers_.reused_ssbo_vert_spatial_map_payloads_());
      sub.bind_ssbo(2, buffers_.reused_ssbo_vert_merged_id_());
      sub.bind_ssbo(3, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(4, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(5, buffers_.ssbo_vert_flags_); 
      sub.push_constant("pcs_hash_map_size_", hashmap_size);
      sub.push_constant("pcs_vert_count_", (int)num_verts_in);
    }; 

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_merge_verts_bootstrap");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_INIT));

      bind_rsc(sub);
      /* TODO: use GL buffer copy function rather than this */
      int num_groups = compute_num_groups(hashmap_size, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_merge_verts_spatial_hashing");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_HASH));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_verts_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_merge_verts_deduplicate");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_REINDEX));

      bind_rsc(sub); 

      int num_groups = compute_num_groups(num_verts_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }

  void StrokeGenPassModule::append_subpass_meshing_wedge_adjacency(
    int num_edges_in, int num_verts_in, bool debug
  )
  {
    const int hashmap_size = std::min(MAX_GPU_HASH_TABLE_SIZE, num_edges_in * 16);

    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub)
    {
      GPUStorageBuf *reused_ssbo_wedge_flooding_pointers_in_ = nullptr;
      GPUStorageBuf *reused_ssbo_wedge_flooding_pointers_out_ = nullptr;
      buffers_.reused_ssbo_wedge_flooding_pointers_(
          1/*iter=="-1"*/, reused_ssbo_wedge_flooding_pointers_in_, reused_ssbo_wedge_flooding_pointers_out_
      ); 

      sub.bind_ssbo(0, buffers_.reused_ssbo_vert_merged_id_());
      sub.bind_ssbo(1, buffers_.reused_ssbo_vert_spatial_map_headers_());
      sub.bind_ssbo(2, buffers_.reused_ssbo_edge_spatial_map_payloads_());
      sub.bind_ssbo(3, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(4, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(5, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(6, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(7, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(8, buffers_.ssbo_edge_flags_);
      sub.bind_ssbo(9, buffers_.ssbo_dyn_mesh_counters_in_());
      sub.bind_ssbo(10, buffers_.ssbo_dyn_mesh_counters_out_());
      sub.bind_ssbo(11, buffers_.ssbo_edge_split_counters_);
      sub.bind_ssbo(12, buffers_.ssbo_bnpr_vert_debug_draw_args_); 
      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_); 
      sub.push_constant("pcs_hash_map_size_", hashmap_size);
      sub.push_constant("pcs_edge_count_", num_edges_in);
      sub.push_constant("pcs_vert_count_", num_verts_in); 
    };

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_adj_bootstrap");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_INIT));

      bind_rsc(sub);
      /* TODO: use GL buffer copy function rather than this */
      int num_groups = compute_num_groups(hashmap_size, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_adj_hashing");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_HASH));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_edges_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_adj_finish");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_FILL));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_edges_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }

  void StrokeGenPassModule::append_subpass_select_remeshed_edges(int num_edges, int num_verts, EdgeFloodingOptions options)
  {
    auto bind_flooding_rsc = [&](
        draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
        int flooding_iter = 0
      ) {
      GPUStorageBuf *reused_ssbo_wedge_flooding_pointers_in_ = nullptr;
      GPUStorageBuf *reused_ssbo_wedge_flooding_pointers_out_ = nullptr;
      buffers_.reused_ssbo_wedge_flooding_pointers_(flooding_iter,
                                                    reused_ssbo_wedge_flooding_pointers_in_,
                                                    reused_ssbo_wedge_flooding_pointers_out_);

      sub.bind_ssbo(0, reused_ssbo_wedge_flooding_pointers_in_);
      sub.bind_ssbo(1, reused_ssbo_wedge_flooding_pointers_out_);
      sub.bind_ssbo(2, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(3, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(4, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(5, buffers_.ssbo_selected_edge_to_edge_);
      sub.bind_ssbo(6, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(7, buffers_.ssbo_dyn_mesh_counters_out_());
      sub.bind_ssbo(8, buffers_.ssbo_edge_flags_);
      sub.bind_ssbo(9, buffers_.ssbo_vbo_full_);

      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_); 

      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_vert_count_", num_verts); 
      sub.push_constant("pcs_output_selected_edge_to_edge_",
                        (options.compact_edges && options.output_selected_to_edge) ? 1 : 0);
      sub.push_constant("pcs_output_edge_to_selected_edge_",
                        (options.compact_edges && options.output_edge_to_selected) ? 1 : 0);
    };

    const int num_flooding_iters = meshing_params.num_edge_flooding_iters + 1/* +1 init pass */;
    // Must be even number !!!
    int flooding_iter = 0;
    for (; flooding_iter < num_flooding_iters; ++flooding_iter) {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_wedge_flooding");

      if (flooding_iter == 0) {
        // first iter, setup selection seeds
        sub.shader_set(shaders_.static_shader_get(MESH_WEDGE_FLOODING_FIRST_ITER_INIT)); 
      }
      else if (flooding_iter < num_flooding_iters - 1)
        sub.shader_set(shaders_.static_shader_get(MESH_WEDGE_FLOODING_ITER));
      else {
        // last iter, output results
        if (options.compact_edges)
          sub.shader_set(shaders_.static_shader_get(MESH_WEDGE_FLOODING_LAST_ITER_COMPACTION));
        else
          sub.shader_set(shaders_.static_shader_get(MESH_WEDGE_FLOODING_LAST_ITER_OUTPUT_FLAGS));
      }

      bind_flooding_rsc(sub, flooding_iter);

      int num_groups = compute_num_groups(num_edges, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }

  void StrokeGenPassModule::append_subpass_split_edges(EdgeSplitMode mode, int iter_split, int num_edges, int num_verts)
  {
    auto bind_src = [&](
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub, float remesh_edge_len,
        GPUStorageBuf* ssbo_14_override = nullptr
      ) {
        sub.bind_ssbo(0, buffers_.ssbo_dyn_mesh_counters_in_()); // in
        sub.bind_ssbo(1, buffers_.ssbo_dyn_mesh_counters_out_()); // out
        sub.bind_ssbo(2, buffers_.ssbo_edge_split_counters_);
        sub.bind_ssbo(3, buffers_.ssbo_vbo_full_);
        sub.bind_ssbo(4, buffers_.ssbo_edge_to_vert_);
        sub.bind_ssbo(5, buffers_.ssbo_edge_to_edges_);
        sub.bind_ssbo(6, buffers_.ssbo_vert_to_edge_list_header_);
        sub.bind_ssbo(7, buffers_.ssbo_vert_flags_);
        sub.bind_ssbo(8, buffers_.ssbo_edge_flags_);
        sub.bind_ssbo(9, buffers_.reused_ssbo_per_edge_split_info_()); 
        sub.bind_ssbo(10, buffers_.reused_ssbo_per_split_edge_info_());
        sub.bind_ssbo(11, buffers_.ssbo_selected_edge_to_edge_);
        sub.bind_ssbo(12, buffers_.ssbo_bnpr_mesh_pool_counters_);
        sub.bind_ssbo(13, buffers_.reused_ssbo_vtx_remesh_len_());
        if (ssbo_14_override == nullptr)
          sub.bind_ssbo(14, buffers_.reused_ssbo_epos_subd_());
        else
          sub.bind_ssbo(14, ssbo_14_override /*override the name ssbo_epos_subd_ in shader using macro*/); 
        sub.bind_ssbo(15, buffers_.ssbo_vnor_);

        sub.push_constant("pcs_split_iter_", iter_split); 
        sub.push_constant("pcs_remesh_edge_len_", remesh_edge_len); 
        sub.push_constant("pcs_edge_count_", num_edges); 
        sub.push_constant("pcs_vert_count_", num_verts);

        // Test contour insertion
        sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);
        sub.push_constant("pcs_split_mode_", (int)mode); 
        // ----------------------
    };

    float remesh_len_scaled = meshing_params.remeshing_targ_edge_len / 100.0f; 
    if (iter_split == 0u) {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_split_init");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_SPLIT_EDGE_INIT));
      bind_src(sub, remesh_len_scaled);
      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_split_compact");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_SPLIT_EDGE_COMPACT));
      bind_src(sub, remesh_len_scaled);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_split_edge"); 
      sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_SPLIT_EDGES));
      sub.bind_ssbo(0, buffers_.ssbo_edge_split_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_indirect_dispatch_args_per_split_edge_);
      sub.push_constant("pcs_edge_split_dispatch_group_size_", (int)(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT));
      sub.push_constant("pcs_split_iter_", iter_split);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_split_exclude_border");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_SPLIT_EDGE_EXCLUDE_BORDER));
      bind_src(sub, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_split_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_split_resolve_conflict");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_SPLIT_EDGE_RESOLVE_CONFLICT));
      bind_src(sub, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_split_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_split_execute");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_SPLIT_EDGE_EXECUTE));
      bind_src(sub, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_split_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_collapse_edges(int iter_remesh, int iter_collapse, int num_edges, int num_verts)
  {
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
                        int pingpong_id, float remesh_edge_len, GPUStorageBuf* ssbo_12_override = nullptr) {
      sub.bind_ssbo(0, buffers_.ssbo_dyn_mesh_counters_in_());   // in
      sub.bind_ssbo(1, buffers_.ssbo_dyn_mesh_counters_out_());    // out
      sub.bind_ssbo(2, buffers_.ssbo_edge_collapse_counters_);
      sub.bind_ssbo(3, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(4, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(5, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(6, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(7, buffers_.ssbo_vert_flags_);
      sub.bind_ssbo(8, buffers_.ssbo_edge_flags_);
      sub.bind_ssbo(9, buffers_.reused_ssbo_per_edge_collapse_info_in_(pingpong_id));
      sub.bind_ssbo(10, buffers_.reused_ssbo_per_edge_collapse_info_out_(pingpong_id));
      sub.bind_ssbo(11, buffers_.reused_ssbo_per_collapse_edge_info_());
      if (ssbo_12_override != nullptr)
        sub.bind_ssbo(12, ssbo_12_override);  // overwrite, we don't have slot for that
      else
        sub.bind_ssbo(12, buffers_.reused_ssbo_per_vert_collapse_wedge_id_());
      sub.bind_ssbo(13, buffers_.ssbo_selected_edge_to_edge_);
      sub.bind_ssbo(14, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.push_constant("pcs_collapse_iter_", iter_collapse);
      sub.push_constant("pcs_remesh_edge_len_", remesh_edge_len);
      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_vert_count_", num_verts);
    };

    float remesh_len_scaled = meshing_params.remeshing_targ_edge_len / 100.0f;
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_init");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_INIT));
      bind_src(sub, 0, remesh_len_scaled);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_compact");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_COMPACT));
      bind_src(sub, 0, remesh_len_scaled,
        buffers_.reused_ssbo_vtx_remesh_len_() // override ssbo slot
      );
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_collapsed_edge");
      sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_COLLAPSE_EDGES));
      sub.bind_ssbo(0, buffers_.ssbo_edge_collapse_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_indirect_dispatch_args_per_collapsed_edge_);
      sub.push_constant("pcs_edge_collapse_dispatch_group_size_", (int)(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT));
      sub.push_constant("pcs_collapse_iter_", iter_collapse);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_resolve_conflict_0");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_0));
      bind_src(sub, 1, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_collapsed_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_resolve_conflict_1");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_1));
      bind_src(sub, 2, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_collapsed_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_resolve_conflict_2");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_RESOLVE_CONFLICT_2));
      bind_src(sub, 3, meshing_params.remeshing_targ_edge_len);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_collapsed_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_collapse_execute");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_OP_COLLAPSE_EDGE_EXECUTE));
      bind_src(sub, 4, meshing_params.remeshing_targ_edge_len,
        buffers_.reused_ssbo_vtx_remesh_len_() // override ssbo slot
      );
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_collapsed_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_flip_edges(StrokeGenPassModule::EdgeFlipOptiGoal opti_goal,
                                                      int iter_flip,
                                                      int num_edges,
                                                      int num_verts)
  {
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_dyn_mesh_counters_in_());  // in
      sub.bind_ssbo(1, buffers_.ssbo_dyn_mesh_counters_out_());  // out
      sub.bind_ssbo(2, buffers_.ssbo_edge_flip_counters_);
      sub.bind_ssbo(3, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(4, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(5, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(6, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(7, buffers_.ssbo_edge_flags_);
      sub.bind_ssbo(8, buffers_.reused_ssbo_per_edge_flip_info_());
      sub.bind_ssbo(9, buffers_.reused_ssbo_per_flip_edge_info_());
      sub.bind_ssbo(10, buffers_.reused_ssbo_vertex_edge_flip_info_());
      sub.bind_ssbo(11, buffers_.ssbo_selected_edge_to_edge_);
      sub.bind_ssbo(12, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(13, buffers_.ssbo_vert_flags_); 
      sub.push_constant("pcs_flip_opti_goal_type_", (int)opti_goal); 
      sub.push_constant("pcs_flip_iter_", iter_flip);
      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_vert_count_", num_verts);
    };

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_flip_init");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_FLIP_EDGE_INIT));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE); 
    }

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_flip_compact");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_FLIP_EDGE_COMPACT));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }

    {
      auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_flip_edge");
      sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_FLIP_EDGES));

      sub.bind_ssbo(0, buffers_.ssbo_edge_flip_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_indirect_dispatch_args_per_flip_edge_);
      sub.push_constant("pcs_edge_flip_dispatch_group_size_", (int)(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)); 
      sub.push_constant("pcs_flip_iter_", iter_flip);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_flip_validate");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_FLIP_EDGE_VALIDATE));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_flip_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_flip_resolve_conflict");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_FLIP_EDGE_RESOLVE_CONFLICT));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_flip_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_edge_flip_execute");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_FLIP_EDGE_EXECUTE));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_flip_edge_);
      sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_split_faces(int iter_split, int num_edges, int num_verts)
  {
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0,  buffers_.ssbo_bnpr_mesh_pool_counters_); 
      sub.bind_ssbo(1,  buffers_.ssbo_selected_edge_to_edge_); 
      sub.bind_ssbo(2,  buffers_.ssbo_dyn_mesh_counters_in_()); 
      sub.bind_ssbo(3,  buffers_.ssbo_dyn_mesh_counters_out_()); 
      sub.bind_ssbo(4,  buffers_.ssbo_face_split_counters_); 
      sub.bind_ssbo(5,  buffers_.ssbo_vbo_full_); 
      sub.bind_ssbo(6,  buffers_.ssbo_edge_to_vert_); 
      sub.bind_ssbo(7,  buffers_.ssbo_edge_to_edges_); 
      sub.bind_ssbo(8,  buffers_.ssbo_vert_to_edge_list_header_); 
      sub.bind_ssbo(9,  buffers_.ssbo_vert_flags_); 
      sub.bind_ssbo(10, buffers_.ssbo_edge_flags_); 
      sub.bind_ssbo(11, buffers_.reused_ssbo_per_face_split_info_()); 
      sub.push_constant("pcs_split_iter_", iter_split);
      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_vert_count_", num_verts);
    };

    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_face_split_init");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_SPLIT_FACE_INIT)); 
      bind_src(sub);
      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_face_split_work_generation");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_SPLIT_FACE_GENERATE_WORKS));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
    {
      auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_split_face");
      sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_SPLIT_FACES));
      
      sub.bind_ssbo(0, buffers_.ssbo_face_split_counters_); 
      sub.bind_ssbo(1, buffers_.ssbo_indirect_dispatch_args_per_split_face_);
      sub.push_constant("pcs_face_split_dispatch_group_size_", (int)(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT));
      sub.push_constant("pcs_split_iter_", iter_split);
      
      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
    {
      auto &sub = pass_extract_geom().sub("bnpr_meshing_face_split_execute");
      sub.shader_set(shaders_.static_shader_get(MESH_OP_SPLIT_FACE_EXECUTE));
      bind_src(sub);
      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_split_face_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_remeshed_edges_(int num_static_edges, bool only_selected_edges)
  {
    auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_remeshed_edge");
    sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_REMESHED_EDGES));

    sub.bind_ssbo(0, buffers_.ssbo_dyn_mesh_counters_out_());
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(2, buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_); 

    sub.push_constant("pcs_edge_count_", num_static_edges);
    sub.push_constant("pcs_only_selected_elems_", only_selected_edges ? 1 : 0); 
    sub.push_constant("pcs_remeshed_edges_dispatch_group_size_", (int)GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE); 
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_remeshed_verts_(
      int num_static_verts, bool only_selected_elems_)
  {
    auto &sub = pass_extract_geom().sub("strokegen_remeshing_fill_dispatch_args_per_remeshed_vert");
    sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_REMESHED_VERTS));
    sub.bind_ssbo(0, buffers_.ssbo_dyn_mesh_counters_out_());
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(2, buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
    sub.push_constant("pcs_vert_count_", num_static_verts);
    sub.push_constant("pcs_only_selected_elems_", only_selected_elems_ ? 1 : 0); 
    sub.push_constant("pcs_remeshed_verts_dispatch_group_size_", (int)GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE); 
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_contour_edges(PassSimple& pass, bool all_contour_edges)
  { // fill dispatch args for contour edges
    auto &sub = pass.sub("strokegen_fill_dispatch_args_per_contour_edge");
    sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DISPATCH_ARGS_CONTOUR_EDGES));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
    sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
    int kernel_size = GROUP_SIZE_STROKEGEN_GEOM_EXTRACT;
    sub.push_constant("pc_per_contour_edge_dispatch_group_size_", kernel_size);
    sub.push_constant("pc_dispatch_for_all_edges_", all_contour_edges ? 1 : 0);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_contour_verts(PassSimple &pass)
  {  // fill dispatch args for contour edges
    auto &sub = pass.sub("strokegen_fill_dispatch_args_per_contour_vert");
    sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DISPATCH_ARGS_CONTOUR_VERTS));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
    int kernel_size = GROUP_SIZE_STROKEGEN_GEOM_EXTRACT;
    sub.push_constant("pc_per_contour_vert_dispatch_group_size_", kernel_size);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_contour_frags(PassSimple& pass, bool all_contour_frags)
  {
    auto &sub = pass.sub("strokegen_fill_dispatch_args_per_contour_fragment");
    sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DISPATCH_ARGS_CONTOUR_FRAGS));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
    sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
    int kernel_size = GROUP_SIZE_STROKEGEN_GEOM_EXTRACT;
    sub.push_constant("pc_per_contour_frag_dispatch_group_size_", kernel_size);
    sub.push_constant("pc_dispatch_for_all_frags_", all_contour_frags ? 1 : 0);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_contour_2d_samples(PassSimple &pass)
  {
    auto &sub = pass.sub("strokegen_fill_dispatch_args_per_contour_2d_sample");
    sub.shader_set(shaders_.static_shader_get(FILL_DISPATCH_ARGS_CONTOUR_2D_SAMPLES));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
    int kernel_size = GROUP_SIZE_STROKEGEN_GEOM_EXTRACT;
    sub.push_constant("pc_per_contour_2d_sample_dispatch_group_size_", kernel_size);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND); 
  }

  void StrokeGenPassModule::append_subpass_surf_geom_analysis(
      ResourceHandle& rsc_handle, int num_verts, int num_edges, const SurfaceAnalysisContext & ctx,
      const SurfaceDebugContext & dbg_options
  ){
    int ssbo_offset_base = 0;
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub, bool only_selected_verts, bool output_dbg_lines) {
      sub.bind_ssbo(0,  buffers_.ssbo_dyn_mesh_counters_out_()); 
      sub.bind_ssbo(1,  buffers_.ssbo_bnpr_mesh_pool_counters_); 
      sub.bind_ssbo(2,  buffers_.ssbo_selected_edge_to_edge_); 
      sub.bind_ssbo(3,  buffers_.ssbo_selected_vert_to_vert_); 
      sub.bind_ssbo(4,  buffers_.ssbo_edge_to_vert_); 
      sub.bind_ssbo(5,  buffers_.ssbo_edge_to_edges_); 
      sub.bind_ssbo(6,  buffers_.ssbo_vert_to_edge_list_header_); 
      sub.bind_ssbo(7,  buffers_.ssbo_vert_flags_); 
      sub.bind_ssbo(8, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(9, buffers_.ssbo_dbg_lines_); 
      ssbo_offset_base = 10;

      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_); 

      sub.push_constant("pcs_vert_count_", num_verts);
      sub.push_constant("pcs_edge_count_", num_edges);
      sub.push_constant("pcs_rsc_handle", (int)rsc_handle.resource_index());
      int dbg_lines = output_dbg_lines ? 1 : 0; 
      sub.push_constant("pcs_output_dbg_geom_", dbg_lines);
      sub.push_constant("pcs_dbg_geom_scale_", dbg_options.dbg_line_length);
      sub.push_constant("pcs_only_selected_verts_", only_selected_verts ? 1 : 0); 
    };

    // Calculate Feature Edges
    if (ctx.calc_feature_edges) {
      append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, ctx.only_selected_edges);
      {
        auto &sub = pass_extract_geom().sub("bnpr_geom_analysis_feature_edges");
        sub.shader_set(shaders_.static_shader_get(MESH_ANALYSE_EDGE_FEATURES));

        bind_src(sub, false, false);
        sub.bind_ssbo(ssbo_offset_base, buffers_.ssbo_edge_flags_);
        sub.push_constant("pcs_only_selected_edges_", ctx.only_selected_edges);

        sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }
    }


    // Calculate Order-0 Vertex Attributes
    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, ctx.order_0_only_selected); 
    if (ctx.calc_vert_normal)
    {
      auto &sub = pass_extract_geom().sub("bnpr_geom_analysis_order_0_vert_normal");
      eShaderType shaderType =
          ctx.calc_vert_voronoi_area && ctx.calc_vert_topo_flags ?
              MESH_ANALYSE_VERT_NORMAL_VORONOIAREA_TOPOFLAGS :
          ctx.calc_vert_voronoi_area ?
              MESH_ANALYSE_VERT_NORMAL_VORONOIAREA :
          ctx.calc_vert_topo_flags ?
              MESH_ANALYSE_VERT_NORMAL_TOPOFLAGS :
              MESH_ANALYSE_VERT_NORMAL; 
      sub.shader_set(shaders_.static_shader_get(shaderType));


      bind_src(sub, ctx.order_0_only_selected, dbg_options.dbg_vert_normal);
      sub.bind_ssbo(ssbo_offset_base, ctx.ssbo_vnor_);
      if (shaderType == MESH_ANALYSE_VERT_NORMAL_VORONOIAREA_TOPOFLAGS) {
        sub.bind_ssbo(ssbo_offset_base + 1, ctx.ssbo_varea_);
        sub.bind_ssbo(ssbo_offset_base + 2, buffers_.ssbo_edge_flags_);
      }
      else if (shaderType == MESH_ANALYSE_VERT_NORMAL_VORONOIAREA)
        sub.bind_ssbo(ssbo_offset_base + 1, ctx.ssbo_varea_);
      else if (shaderType == MESH_ANALYSE_VERT_NORMAL_TOPOFLAGS)
        sub.bind_ssbo(ssbo_offset_base + 1, buffers_.ssbo_edge_flags_);

      sub.push_constant("pcs_output_vertex_facing_flag_", ctx.output_vertex_facing_flag); 

      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }

    // Calculate Order-1 Vertex Attributes
    append_subpass_fill_dispatch_args_remeshed_edges_(num_edges, ctx.order_1_only_selected); 
    append_subpass_fill_dispatch_args_remeshed_verts_(num_verts, ctx.order_1_only_selected); 
    int ssbo_offset_base_1 = 0; 
    auto bind_src_order_1 = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
                                bool output_dbg_lines = false) {
          bind_src(sub, ctx.order_1_only_selected, output_dbg_lines); 
          sub.bind_ssbo(ssbo_offset_base + 0, ctx.ssbo_vnor_);
          sub.bind_ssbo(ssbo_offset_base + 1, ctx.ssbo_varea_);
          sub.bind_ssbo(ssbo_offset_base + 2, buffers_.ssbo_edge_flags_);
          ssbo_offset_base_1 = ssbo_offset_base + 3;
        };
    if (ctx.calc_vert_curvature)
    {
      if (ctx.curvature_estimator == SurfaceAnalysisContext::Rusinkiewicz) {
        auto &sub = pass_extract_geom().sub("bnpr_geom_analysis_order_1_vert_curv_pass_0");
        sub.shader_set(shaders_.static_shader_get(MESH_ANALYSE_VERT_CURV_PASS_0));

        bind_src_order_1(sub); 
        sub.bind_ssbo(ssbo_offset_base_1 + 0, ctx.ssbo_edge_vtensors_);

        sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
      }
      {
        auto &sub = pass_extract_geom().sub("bnpr_geom_analysis_order_1_main_curvature");
        sub.shader_set(shaders_.static_shader_get(
          ctx.curvature_estimator == SurfaceAnalysisContext::Rusinkiewicz ?
          MESH_ANALYSE_VERT_CURV_PASS_1_RUSINKIEWICZ : MESH_ANALYSE_VERT_CURV_PASS_1_JACQUES)
        );

        bind_src_order_1(sub, dbg_options.dbg_vert_curv);
        sub.bind_ssbo(ssbo_offset_base_1 + 0, ctx.ssbo_edge_vtensors_);
        sub.bind_ssbo(ssbo_offset_base_1 + 1, ctx.ssbo_vcurv_pdirs_k1k2_);
        sub.bind_ssbo(ssbo_offset_base_1 + 2, buffers_.ssbo_vcurv_max_);
        sub.push_constant("pcs_output_curv_tensors_", ctx.output_curvature_tensors);
        sub.push_constant("pcs_output_maxcurv_with_cusp_function_", ctx.output_maxcurv_with_cusp_function ? 1 : 0); 

        sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_verts_);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
      }
    }

  }

  void StrokeGenPassModule::rebuild_pass_dbg_geom_drawcall(SurfaceDebugContext dbg_ctx)
  {
    Vector<int> dbgLineTypes;
    if (dbg_ctx.dbg_vert_normal)
      dbgLineTypes.append(SurfaceDebugContext::DbgLineType::vnor); 
    if (dbg_ctx.dbg_vert_curv)
      dbgLineTypes.append(SurfaceDebugContext::DbgLineType::vcurv);
    if (dbg_ctx.dbg_edges)
      dbgLineTypes.append(SurfaceDebugContext::DbgLineType::edges);

    for (int dbg_line_type : dbgLineTypes) {
      {
        auto &sub = pass_draw_debug_lines_.sub("fill_draw_args_debug_lines_");

        sub.shader_set(shaders_.static_shader_get(FILL_DRAW_ARGS_DBG_LINES));

        sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
        sub.bind_ssbo(1, buffers_.ssbo_bnpr_vert_debug_draw_args_);
        sub.push_constant("pcs_line_type_", dbg_line_type);

        sub.dispatch(int3(1, 1, 1));
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }

      pass_draw_debug_lines_.append_draw_dbg_lines_subpass(shaders_, buffers_, dbg_line_type);
    }
  }


  void StrokeGenPassModule::append_subpass_extract_contour_edges(
      gpu::Batch* gpu_batch_line_adj,
      ResourceHandle& rsc_handle,
      gpu::Batch* edge_batch,
      int num_edges,
      gpu::GPUIndexBufType ib_type,
      int edge_visualize_mode, int contour_visualize_mode
  )
  {
    auto bind_rsc = [&](ResourceHandle &rsc_handle,
                        gpu::Batch *edge_batch,
                        gpu::GPUIndexBufType ib_type,
                        int edge_visualize_mode,
                        draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub)
    {
      sub.bind_ssbo(0, buffers_.ssbo_edge_to_vert_ /*&(gpu_batch_line_adj->elem)*/);
      sub.bind_ssbo(1, buffers_.ssbo_vbo_full_ /*&(gpu_batch_line_adj->verts[0])*/);
      sub.bind_ssbo(2, buffers_.reused_ssbo_contour_temp_data_());
      sub.bind_ssbo(3, DRW_manager_get()->matrix_buf.current());
      sub.bind_ssbo(4, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(5, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
      sub.bind_ssbo(6, buffers_.reused_ssbo_edge_to_contour_());
      sub.bind_ssbo(7, buffers_.ssbo_dyn_mesh_counters_out_());
      sub.bind_ssbo(8, buffers_.ssbo_edge_flags_);
      sub.bind_ssbo(9, buffers_.ssbo_edge_to_edges_); 
      sub.bind_ssbo(10, buffers_.reused_ssbo_face_to_vert_draw_depth_());
      // for debugging
      sub.bind_ssbo(11, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(12, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(13, buffers_.ssbo_dbg_lines_);
      sub.bind_ssbo(14, buffers_.ssbo_vert_flags_);
      // --------------
      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);
      sub.push_constant("pcs_ib_fmt_u16", ib_type == gpu::GPU_INDEX_U16 ? 1 : 0);
      sub.push_constant("pcs_num_edges", (int)edge_batch->elem_()->index_len_get() / 4);
      sub.push_constant(
          // not used for now, but I suspect this could be interfering with complex meshes
          "pcs_num_ib_offset",
          (int)edge_batch->elem_()->index_start_get() +
              (int)edge_batch->elem_()->index_base_get());
      sub.push_constant("pcs_rsc_handle", (int)rsc_handle.resource_index());
      sub.push_constant("pcs_edge_visualize_mode_", edge_visualize_mode);
      sub.push_constant("pcs_chain_interpo_contour_", contour_visualize_mode == 2 ? 0 : 1);
      float2 fb_res = textures_.get_contour_raster_screen_res();
      sub.push_constant("pcs_screen_size_", fb_res);
      sub.push_constant("pcs_dbg_geom_scale_", surf_dbg_ctx.dbg_line_length); 
    };

    {
      auto &sub = pass_extract_geom().sub(
         boostrap_before_extract_first_batch ? "bnpr_geom_extract_boostrap" : "bnpr_geom_extract");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_GEOM_EXTRACT));

      bind_rsc(rsc_handle, edge_batch, ib_type, edge_visualize_mode, sub);

      sub.dispatch(buffers_.ssbo_indirect_dispatch_args_per_remeshed_edges_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_setup_contour_edge_data()
  {
    {
      auto &sub = pass_extract_geom().sub("bnpr_geom_extract_mesh_contour_data");
      sub.shader_set(shaders_.static_shader_get(eShaderType::EXTRACT_CURRENT_MESH_CONTOUR_DATA));

      sub.bind_ssbo(0, buffers_.reused_ssbo_contour_temp_data_());
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
      sub.bind_ssbo(3, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(4, buffers_.reused_ssbo_edge_to_contour_());
      sub.bind_ssbo(5, buffers_.ssbo_contour_to_contour_);
      sub.bind_ssbo(6, buffers_.ssbo_vert_to_edge_list_header_);
      sub.bind_ssbo(7, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(8, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(9, buffers_.ssbo_list_ranking_inputs_);
      sub.bind_ssbo(10, buffers_.ssbo_contour_edge_transfer_data_);
      sub.bind_ssbo(11, buffers_.ssbo_vcurv_max_);
      sub.bind_ssbo(12, buffers_.ssbo_contour_raster_data_); 
      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_); 
      float2 fb_res = textures_.get_contour_raster_screen_res(); 
      sub.push_constant("pcs_screen_size_", fb_res); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_serialize_contour_edges()
  {
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.reused_ssbo_contour_edge_rank_());
      sub.bind_ssbo(1, buffers_.reused_ssbo_contour_edge_list_len_());
      sub.bind_ssbo(2, buffers_.reused_ssbo_contour_edge_list_head_info_());
      sub.bind_ssbo(3, buffers_.ssbo_contour_edge_transfer_data_);
      sub.bind_ssbo(4, buffers_.ssbo_contour_snake_rank_);
      sub.bind_ssbo(5, buffers_.ssbo_contour_snake_list_len_);
      sub.bind_ssbo(6, buffers_.ssbo_contour_snake_list_head_);
      sub.bind_ssbo(7, buffers_.ssbo_contour_snake_vpos_);
      sub.bind_ssbo(8, buffers_.ssbo_contour_snake_flags_);
      sub.bind_ssbo(9, buffers_.ssbo_contour_to_contour_);
      sub.bind_ssbo(10, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(11, buffers_.ssbo_segloopconv1d_info_);
      sub.bind_ssbo(12, buffers_.reused_ssbo_in_segloopconv1d_data_contour_seg_denoise());
      sub.bind_ssbo(13, buffers_.ssbo_list_ranking_addressing_counters_); 
    };
    
    {
      auto &sub = pass_process_contours.sub("strokegen_serialize_contour_edges_pass_0");
      sub.shader_set(shaders_.static_shader_get(SERIALIZE_RANKED_CONTOUR_EDGES));

      bind_src(sub); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
  }

  void StrokeGenPassModule::append_subpass_contour_segmentation()
  {
    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_contour_snake_rank_);
      sub.bind_ssbo(1, buffers_.ssbo_contour_snake_list_len_);
      sub.bind_ssbo(2, buffers_.ssbo_contour_snake_list_head_);
      sub.bind_ssbo(3, buffers_.ssbo_contour_snake_flags_);
      sub.bind_ssbo(4, buffers_.reused_ssbo_tree_scan_input_contour_segmentation_step_0());
      sub.bind_ssbo(5, buffers_.reused_ssbo_tree_scan_input_contour_segmentation_step_1());
      sub.bind_ssbo(6, buffers_.reused_ssbo_tree_scan_output_contour_segmentation_step_0());
      sub.bind_ssbo(7, buffers_.reused_ssbo_tree_scan_output_contour_segmentation_step_1());
      sub.bind_ssbo(8, buffers_.ssbo_contour_snake_seg_rank_);
      sub.bind_ssbo(9, buffers_.ssbo_contour_snake_seg_len_);
      sub.bind_ssbo(10, buffers_.reused_ssbo_tree_scan_infos_contour_segmentation_());
      sub.bind_ssbo(11, buffers_.ssbo_bnpr_mesh_pool_counters_);
    }; 

    {
      auto &sub = pass_process_contours.sub("strokegen_setup_contour_segmentation");
      sub.shader_set(shaders_.static_shader_get(SETUP_CONTOUR_SEGMENTATION));

      bind_rsc(sub); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }

    {
      ScanSettings scan_settings;
      scan_settings.is_validation_shader = false;
      scan_settings.frame_counter = 0; 
      scan_settings.use_indirect_dispatch = true;
      scan_settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_contour_segmentation_();
      scan_settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_contour_segmentation_step_0();
      scan_settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_contour_segmentation_step_0();
      scan_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      scan_settings.shader_upsweep = SEGSCAN_UINT_ADD_UPSWEEP;
      scan_settings.shader_aggregate = SEGSCAN_UINT_ADD_AGGREGATE;
      scan_settings.shader_dwsweep= SEGSCAN_UINT_ADD_DWSWEEP;

      append_subpass_segscan(scan_settings, pass_process_contours); 
    }
    {
      ScanSettings scan_settings;
      scan_settings.is_validation_shader = false;
      scan_settings.frame_counter = 0; 
      scan_settings.use_indirect_dispatch = true;
      scan_settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_contour_segmentation_();
      scan_settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_contour_segmentation_step_1();
      scan_settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_contour_segmentation_step_1();
      scan_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      scan_settings.shader_upsweep = SEGSCAN_UINT_ADD_UPSWEEP;
      scan_settings.shader_aggregate = SEGSCAN_UINT_ADD_AGGREGATE;
      scan_settings.shader_dwsweep= SEGSCAN_UINT_ADD_DWSWEEP;

      append_subpass_segscan(scan_settings, pass_process_contours); 
    }

    {
      auto &sub = pass_process_contours.sub("strokegen_finish_contour_segmentation");
      sub.shader_set(shaders_.static_shader_get(FINISH_CONTOUR_SEGMENTATION));

      bind_rsc(sub);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::bind_rsc_for_contour_2d_resample_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
    float2 screen_res, float pcs_sample_rate, int& out_ssbo_offset
  ){
    sub.bind_ssbo(0, buffers_.reused_ssbo_contour_2d_resample_raster_data_());
    sub.bind_ssbo(1, buffers_.reused_ssbo_contour_2d_sample_geometry_());
    sub.bind_ssbo(2, buffers_.reused_ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_());
    sub.bind_ssbo(3, buffers_.reused_ssbo_contour_arc_len_param_());
    sub.bind_ssbo(4, buffers_.reused_ssbo_tree_scan_input_2d_resampler_alloc_samples_());
    sub.bind_ssbo(5, buffers_.reused_ssbo_contour_to_start_sample_());
    sub.bind_ssbo(6, buffers_.reused_ssbo_tree_scan_input_2d_resample_contour_idmapping_());
    sub.bind_ssbo(7, buffers_.reused_ssbo_2d_sample_to_contour_());
    sub.bind_ssbo(8, buffers_.ssbo_contour_snake_rank_);
    sub.bind_ssbo(9, buffers_.ssbo_contour_snake_list_len_);
    sub.bind_ssbo(10, buffers_.ssbo_contour_snake_list_head_);
    sub.bind_ssbo(11, buffers_.ssbo_contour_snake_flags_);
    sub.bind_ssbo(12, buffers_.ssbo_contour_snake_vpos_);
    sub.bind_ssbo(13, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(14, buffers_.reused_ssbo_tree_scan_infos_2d_resampler_());
    out_ssbo_offset = 15;

    sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);

    sub.push_constant("pcs_screen_size_", screen_res);
    sub.push_constant("pcs_sample_rate_", pcs_sample_rate); 
  }

  void StrokeGenPassModule::bind_rsc_for_contour_2d_sample_evaluation_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub,
      float2 screen_res, float pcs_sample_rate, int &out_ssbo_offset)
  {
    sub.bind_ssbo(0, buffers_.reused_ssbo_contour_2d_sample_geometry_());
    sub.bind_ssbo(1, buffers_.reused_ssbo_contour_2d_sample_topology_());
    sub.bind_ssbo(2, buffers_.reused_ssbo_contour_arc_len_param_());
    sub.bind_ssbo(3, buffers_.reused_ssbo_contour_to_start_sample_());
    sub.bind_ssbo(4, buffers_.reused_ssbo_2d_sample_to_contour_());
    sub.bind_ssbo(5, buffers_.ssbo_contour_snake_flags_); 
    sub.bind_ssbo(6, buffers_.ssbo_contour_snake_rank_); 
    sub.bind_ssbo(7, buffers_.ssbo_contour_snake_list_head_); 
    sub.bind_ssbo(8, buffers_.ssbo_contour_snake_list_len_); 
    sub.bind_ssbo(9, buffers_.ssbo_contour_snake_seg_rank_); 
    sub.bind_ssbo(10, buffers_.ssbo_contour_snake_seg_len_); 
    sub.bind_ssbo(11, buffers_.ssbo_bnpr_mesh_pool_counters_); 
    out_ssbo_offset = 12;

    sub.bind_image(0, textures_.tex2d_contour_dbg_); 

    sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);

    sub.push_constant("pcs_screen_size_", screen_res);
    sub.push_constant("pcs_sample_rate_", pcs_sample_rate);
  }

  void StrokeGenPassModule::append_subpass_contour_arclen_parameterization(float2 screen_res, float sample_rate)
  {
    int ssbo_offset = 0; 
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_arclen_parameterization");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_ARCLEN_PARAMETERIZATION)); 

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);
      sub.bind_image(0, textures_.tex2d_contour_dbg_); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE); 
    }
    {
      ScanSettings settings;
      settings.use_indirect_dispatch = true;
      settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_2d_resampler_();
      settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_();
      settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_2d_resampler_accumulate_curvlen_();
      settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      settings.shader_upsweep = SEGSCAN_FLOAT_ADD_UPSWEEP;
      settings.shader_aggregate = SEGSCAN_FLOAT_ADD_AGGREGATE;
      settings.shader_dwsweep = SEGSCAN_FLOAT_ADD_DWSWEEP; 
      settings.is_validation_shader = false;
      settings.frame_counter = 0;

      append_subpass_segscan(settings, pass_process_contours); 
    }
  }


  void StrokeGenPassModule::append_subpass_2d_sample_segmentation(
      float2 screen_res, float sample_rate, bool is_segmentation_by_curve_pass
  ) {
    int ssbo_offset_ = 0;
    { // Mark seg tails
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_segmentation_prep_seg_tails");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_SEGMENTATION_PREP_SEGTAILS));

      bind_rsc_for_contour_2d_sample_evaluation_(sub, screen_res, sample_rate, ssbo_offset_);
      sub.push_constant("pcs_segment_by_seg_", is_segmentation_by_curve_pass ? 0 : 1);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // Setup Segmentation
      auto &sub = pass_process_contours.sub(
          "strokegen_contour_2d_resample_segmentation_setup_segscan"
          );
      sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_EVAL_TOPO_SETUP_SEGSCAN));

      bind_rsc_for_contour_2d_sample_evaluation_(sub, screen_res, sample_rate, ssbo_offset_);
      sub.bind_ssbo(ssbo_offset_ + 0,
                    buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_0_());
      sub.bind_ssbo(ssbo_offset_ + 1,
                    buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_1_());
      sub.bind_ssbo(ssbo_offset_ + 2, buffers_.reused_ssbo_tree_scan_infos_2d_resampler_()); 
      sub.push_constant("pcs_segment_by_seg_", is_segmentation_by_curve_pass ? 0 : 1); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }


    ScanSettings settings;
    settings.use_indirect_dispatch = true;
    settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_2d_resampler_();
    settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
    settings.shader_upsweep = SEGSCAN_UINT_ADD_UPSWEEP;
    settings.shader_aggregate = SEGSCAN_UINT_ADD_AGGREGATE;
    settings.shader_dwsweep = SEGSCAN_UINT_ADD_DWSWEEP;

    settings.is_validation_shader = false;
    settings.frame_counter = 0;
    {
      settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_0_();
      settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_2d_sample_segmentation_0_();
      append_subpass_segscan(settings, pass_process_contours); 
    }
    {
      settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_1_();
      settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_2d_sample_segmentation_1_();
      append_subpass_segscan(settings, pass_process_contours); 
    }

    { // Finish Segmentation
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_eval_topo_finish");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_EVAL_TOPO_FINISH_SEGSCAN));

      sub.bind_ssbo(0, buffers_.reused_ssbo_contour_2d_sample_geometry_()); 
      sub.bind_ssbo(1, buffers_.reused_ssbo_contour_2d_sample_topology_()); 
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_counters_); 
      sub.bind_ssbo(3, buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_0_()); 
      sub.bind_ssbo(4, buffers_.reused_ssbo_tree_scan_output_2d_sample_segmentation_0_()); 
      sub.bind_ssbo(5, buffers_.reused_ssbo_tree_scan_input_2d_sample_segmentation_1_()); 
      sub.bind_ssbo(6, buffers_.reused_ssbo_tree_scan_output_2d_sample_segmentation_1_()); 

      sub.bind_image(0, textures_.tex2d_contour_dbg_); 

      sub.push_constant("pcs_segment_by_seg_", is_segmentation_by_curve_pass ? 0 : 1);
      sub.push_constant("pcs_screen_size_", screen_res);
      sub.push_constant("pcs_sample_rate_", sample_rate);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }


  void StrokeGenPassModule::append_subpass_contour_generate_2d_samples(float2 screen_res, float sample_rate)
  {
    int ssbo_offset = 0; 

    // Convert 3D Contour to 2D Samples
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_alloc_samples");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_ALLOC_2D_SAMPLES));

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);
      sub.bind_image(0, textures_.tex2d_contour_dbg_);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    {
      ScanSettings settings;
      settings.use_indirect_dispatch = true;
      settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_2d_resampler_();
      settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_2d_resampler_alloc_samples_();
      settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_2d_resampler_alloc_samples_();
      settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      settings.shader_upsweep = SCAN_UINT_ADD_UPSWEEP;
      settings.shader_aggregate = SCAN_UINT_ADD_AGGREGATE;
      settings.shader_dwsweep = SCAN_UINT_ADD_DWSWEEP; 
      settings.is_validation_shader = false;
      settings.frame_counter = 0;

      append_subpass_scan(settings, pass_process_contours);
    }
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_alloc_samples_finish");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_ALLOC_2D_SAMPLES_FINISH));

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);
      sub.bind_image(0, textures_.tex2d_contour_dbg_);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }

    append_subpass_fill_dispatch_args_contour_2d_samples(pass_process_contours);
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_idmapping_clear_buffer");
      sub.shader_set(shaders_.static_shader_get(CLEAR_2D_SAMPLE_TO_CONTOUR_IDMAPPING));

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_idmapping_setup_segscan");
      sub.shader_set(shaders_.static_shader_get(PREP_2D_SAMPLE_TO_CONTOUR_IDMAPPING));

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    {
      ScanSettings settings;
      settings.use_indirect_dispatch = true;
      settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_2d_resampler_();
      settings.ssbo_in_scan_data_ =
          buffers_.reused_ssbo_tree_scan_input_2d_resample_contour_idmapping_();
      settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      settings.ssbo_out_scan_data_ =
          buffers_.reused_ssbo_tree_scan_output_2d_resample_contour_idmapping_();
      settings.shader_upsweep = SEGSCAN_UINT_MIN_UPSWEEP;
      settings.shader_aggregate = SEGSCAN_UINT_MIN_AGGREGATE;
      settings.shader_dwsweep = SEGSCAN_UINT_MIN_DWSWEEP;
      settings.is_validation_shader = false;
      settings.frame_counter = 0;

      append_subpass_segscan(settings, pass_process_contours);
    }
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_idmapping_finish");
      sub.shader_set(shaders_.static_shader_get(FINISH_2D_SAMPLE_TO_CONTOUR_IDMAPPING));

      bind_rsc_for_contour_2d_resample_(sub, screen_res, sample_rate, ssbo_offset);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }

    // Evaluate Geom & Topo
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_eval_position");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_EVAL_POSITION));

      bind_rsc_for_contour_2d_sample_evaluation_(sub, screen_res, sample_rate, ssbo_offset);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_contour_2d_resample_raster_data_()); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }


    // Segmentation
    for (int i = 0; i < 2; ++i)
    {
      int ssbo_offset_ = 0; 
      bool is_segmentation_by_curve_pass = (i == 0);
      { // Find seg heads
        auto &sub = pass_process_contours.sub("strokegen_contour_2d_resample_eval_topo_step_0");
        sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_EVAL_TOPO_STEP_0));

        bind_rsc_for_contour_2d_sample_evaluation_(sub, screen_res, sample_rate, ssbo_offset_);
        sub.bind_ssbo(ssbo_offset_ + 0, buffers_.ssbo_segloopconv1d_info_);
        sub.push_constant("pcs_segment_by_seg_", is_segmentation_by_curve_pass ? 0 : 1);

        sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }

      append_subpass_2d_sample_segmentation(screen_res, sample_rate, is_segmentation_by_curve_pass);
    }

    { // Find 2D Corners & Eval tangent and curvature
      SegLoopConv1DSettings settings;
      settings.use_indirect_dispatch = true;
      settings.ssbo_segloopconv1d_info_ = buffers_.ssbo_segloopconv1d_info_;
      settings.ssbo_in_segloopconv1d_data_ = buffers_.reused_ssbo_contour_2d_sample_geometry_();
      settings.ssbo_segloopconv1d_patch_table_ = buffers_.ssbo_segloopconv1d_patch_table_;
      settings.ssbo_out_segloopconv1d_data_ = buffers_.reused_ssbo_contour_2d_sample_topology_();
      settings.shader_build_patch_table = CONV1D_2D_SAMPLE_BUILD_PATCH;
      settings.shader_convolution = CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_0;
      settings.lazy_dispatch = false;
      settings.is_validation_shader = false;
      append_subpass_segloopconv1d(settings, pass_process_contours);

      settings.shader_convolution = CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_1;
      settings.lazy_dispatch = true;
      append_subpass_segloopconv1d(settings, pass_process_contours);


      append_subpass_2d_sample_segmentation(screen_res, sample_rate, false);
      { 
        auto& sub = pass_process_contours.sub("strokegen_contour_2d_resample_eval_topo_remove_fake_corners");
        sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_REJECT_FAKE_CORNERS));

        bind_rsc_for_contour_2d_sample_evaluation_(sub, screen_res, sample_rate, ssbo_offset);

        sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }
      append_subpass_2d_sample_segmentation(screen_res, sample_rate, false);


      settings.shader_convolution = CONV1D_2D_SAMPLE_CALC_TANGENT_CURVATURE;
      settings.lazy_dispatch = true;
      append_subpass_segloopconv1d(settings, pass_process_contours);

      settings.shader_convolution = CONV1D_2D_SAMPLE_DENOISE_VISIBILITY;
      settings.lazy_dispatch = true;
      append_subpass_segloopconv1d(settings, pass_process_contours);
    }


    { // Output 2D stroke mesh
      auto &sub = pass_process_contours.sub("strokegen_calc_contour_2d_curve_render_data");
      sub.shader_set(shaders_.static_shader_get(CONTOUR_2D_SAMPLES_CALC_RENDER_DATA));

      sub.bind_ssbo(0, buffers_.reused_ssbo_contour_2d_sample_geometry_());
      sub.bind_ssbo(1, buffers_.reused_ssbo_contour_2d_sample_topology_());
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(3, buffers_.reused_ssbo_stroke_mesh_pool_());
      sub.bind_image(0, textures_.tex2d_contour_dbg_);
      sub.push_constant("pcs_screen_size_", screen_res);
      sub.push_constant("pcs_stroke_width_", 5.0f);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }




  void StrokeGenPassModule::append_subpass_calc_contour_edges_draw_data()
  {
    int ssbo_offset = 0; 
    auto bind_src = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_contour_snake_rank_);
      sub.bind_ssbo(1, buffers_.ssbo_contour_snake_list_len_);
      sub.bind_ssbo(2, buffers_.ssbo_contour_snake_list_head_);
      sub.bind_ssbo(3, buffers_.ssbo_contour_snake_vpos_);
      sub.bind_ssbo(4, buffers_.ssbo_contour_snake_flags_);
      sub.bind_ssbo(5, buffers_.ssbo_bnpr_mesh_pool_counters_);
      ssbo_offset = 6; 

      sub.bind_ubo(0, buffers_.ubo_view_matrices_);
      sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res()); 
    };

    {
      auto &sub = pass_process_contours.sub("strokegen_calc_contour_edges_draw_data"); 
      sub.shader_set(shaders_.static_shader_get(CALC_CONTOUR_EDGES_DRAW_DATA));
      bind_src(sub); 
      sub.bind_ssbo(ssbo_offset, buffers_.reused_ssbo_bnpr_mesh_pool_());

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_vert_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }


  void StrokeGenPassModule::rebuild_pass_process_contours()
  {
    pass_process_contours.init();

    // Split contour edges based on rasterized visibility
    append_subpass_contour_edges_soft_rasterization();
    append_subpass_visibility_split_contour_edges(); 

    // List ranking to generate curves
    append_subpass_fill_contour_list_ranking_inputs(); 
    append_subpass_fill_dispatch_args_contour_edges(pass_process_contours, true); 
    append_subpass_list_ranking(ContourEdgeLinking, pass_process_contours, true);
    // Reorder contour verts & generate curves
    append_subpass_fill_dispatch_args_contour_edges(pass_process_contours, true);
    append_subpass_serialize_contour_edges();

    // Segmentation based on cusp/(TODO: visibility)
    append_subpass_fill_dispatch_args_contour_verts(pass_process_contours); 
    if (meshing_params.denoise_cusp_segmentation){
      SegLoopConv1DSettings conv1d_settings;
      conv1d_settings.is_validation_shader = false;
      conv1d_settings.use_indirect_dispatch = true;
      conv1d_settings.lazy_dispatch = false;
      conv1d_settings.ssbo_segloopconv1d_info_ = buffers_.ssbo_segloopconv1d_info_;
      conv1d_settings.ssbo_segloopconv1d_patch_table_ = buffers_.ssbo_segloopconv1d_patch_table_;
      conv1d_settings.ssbo_in_segloopconv1d_data_ = buffers_.reused_ssbo_in_segloopconv1d_data_contour_seg_denoise();
      conv1d_settings.ssbo_out_segloopconv1d_data_ = buffers_.ssbo_contour_snake_flags_;
      conv1d_settings.shader_build_patch_table = CONV1D_SEG_DENOISE_BUILD_PATCH;
      conv1d_settings.shader_convolution = CONV1D_SEG_DENOISE_CONVOLUTION;
    
      append_subpass_segloopconv1d(conv1d_settings, pass_process_contours); 
    }
    append_subpass_contour_segmentation();

    // 2D resampling
    float sample_rate = 4.0f; 
    append_subpass_contour_arclen_parameterization(textures_.get_contour_raster_screen_res(), sample_rate); 
    append_subpass_contour_generate_2d_samples(textures_.get_contour_raster_screen_res(), sample_rate);

    append_subpass_calc_contour_edges_draw_data();
  }

  void StrokeGenPassModule::rebuild_pass_contour_edge_drawcall()
  {
    auto &sub = pass_draw_contour_edges.sub("fill_draw_args_contour_edges");

    sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DRAW_ARGS_CONTOUR_EDGES));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_bnpr_contour_mesh_draw_args_);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);

    pass_draw_contour_edges.append_draw_contour_subpass(shaders_, buffers_, textures_);
  }

  void StrokeGenPassModule::rebuild_pass_contour_2d_samples_drawcall()
  {
    StrokegenMeshRasterPass &pass = pass_draw_contour_2d_samples;

    {
      auto &sub = pass.sub("bnpr_geom_fill_draw_args_contour_2d_samples");
      sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DRAW_ARGS_CONTOUR_2D_SAMPLES)); 
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_2d_sample_draw_args_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    pass_draw_contour_2d_samples.append_draw_contour_2d_subpass(shaders_, buffers_, textures_); 
  }

  void StrokeGenPassModule::append_pass_remeshed_surface_depth_drawcall()
  {
    StrokegenMeshRasterPass &pass = pass_draw_remeshed_surface_depth_[++curr_mesh_id_surf_depth];

    {
      auto &sub = pass.sub("bnpr_geom_fill_draw_args_remeshed_surface_depth");

      sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DRAW_ARGS_REMESHED_SURFACE_DEPTH));

      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_contour_mesh_draw_args_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    pass.append_draw_remeshed_surface_depth_subpass(shaders_, buffers_, textures_);
  }

  void StrokeGenPassModule::append_subpass_contour_edges_soft_rasterization()
  {
    GPUSamplerState depthSamplerState;
    depthSamplerState.filtering = GPU_SAMPLER_FILTERING_DEFAULT;
    depthSamplerState.extend_x = depthSamplerState.extend_yz = GPU_SAMPLER_EXTEND_MODE_CLAMP_TO_BORDER; 

    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub)
    {
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.reused_ssbo_frag_to_contour_());
      sub.bind_ssbo(2, buffers_.ssbo_contour_raster_data_);
      sub.bind_ssbo(3, buffers_.reused_ssbo_frag_raster_data_());
      sub.bind_ssbo(4, buffers_.reused_ssbo_tree_scan_infos_contour_fragment_idmapping_()); 
      sub.bind_ssbo(5, buffers_.reused_ssbo_tree_scan_input_contour_fragment_idmapping_());
      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_); 
      sub.bind_image(0, textures_.tex2d_contour_dbg_);
      sub.bind_texture(0, textures_.tex_remeshed_surf_depth_, depthSamplerState); 

      sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res());
    };

    append_subpass_fill_dispatch_args_contour_frags(pass_process_contours, true);

    eGPUBarrier barrier_all = GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND |
                          GPU_BARRIER_SHADER_IMAGE_ACCESS | GPU_BARRIER_TEXTURE_FETCH; 

    { // Clear frag to contour id mapping, dispatch per-fragment
      auto &sub = pass_process_contours.sub("strokegen_clear_frag_to_contour_idmapping");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CLEAR_FRAG_TO_CONTOUR_IDMAPPING));

      bind_rsc(sub);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    append_subpass_fill_dispatch_args_contour_edges(pass_process_contours, true); 
    { // Seed at the segment heads
      auto &sub = pass_process_contours.sub("strokegen_prep_segscan_frag_to_contour_idmapping");
      sub.shader_set(shaders_.static_shader_get(eShaderType::PREP_FRAG_TO_CONTOUR_IDMAPPING));

      bind_rsc(sub);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND); 
    }

    { // Segscan to broadcast seeds to all fragments 
      ScanSettings scan_settings;
      scan_settings.is_validation_shader = false;
      scan_settings.frame_counter = 0; 
      scan_settings.use_indirect_dispatch = true;
      scan_settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_contour_fragment_idmapping_();
      scan_settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_contour_fragment_idmapping_();
      scan_settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_contour_fragment_idmapping_();
      scan_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      scan_settings.shader_upsweep = SEGSCAN_UINT_MIN_UPSWEEP;
      scan_settings.shader_aggregate = SEGSCAN_UINT_MIN_AGGREGATE;
      scan_settings.shader_dwsweep= SEGSCAN_UINT_MIN_DWSWEEP;

      append_subpass_segscan(scan_settings, pass_process_contours); 
    }

    { // Finish segscan
      auto &sub = pass_process_contours.sub("strokegen_finish_segscan_frag_to_contour_idmapping");
      sub.shader_set(shaders_.static_shader_get(eShaderType::FINISH_FRAG_TO_CONTOUR_IDMAPPING));

      bind_rsc(sub);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND); 
    }

    { // Visibility test
      auto &sub = pass_process_contours.sub("strokegen_contour_frag_visibility_test");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_FRAG_VISIBILITY_TEST));

      bind_rsc(sub);
      sub.push_constant("pcs_visibility_thresh_", meshing_params.visibility_thresh); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND |
                  GPU_BARRIER_SHADER_IMAGE_ACCESS | GPU_BARRIER_TEXTURE_FETCH); 
    }
  }

  void StrokeGenPassModule::append_subpass_visibility_split_contour_edges()
  {
    int ssbo_offset = 0; 
    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.reused_ssbo_frag_to_contour_());
      sub.bind_ssbo(2, buffers_.ssbo_contour_raster_data_);
      sub.bind_ssbo(3, buffers_.reused_ssbo_frag_raster_data_());
      sub.bind_ssbo(4, buffers_.reused_ssbo_tree_scan_infos_contour_visibility_split_());
      ssbo_offset = 5; 

      sub.bind_image(0, textures_.tex2d_contour_dbg_);
      sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);
    };

    { 
      auto &sub = pass_process_contours.sub("strokegen_contour_frag_setup_visibility_segmentation");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_FRAG_SETUP_VISIBILITY_SEGMENTATION));

      bind_rsc(sub);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_tree_scan_input_contour_visibility_split_0_());
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_tree_scan_input_contour_visibility_split_1_()); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    {
      ScanSettings segscan_test_settings;
      segscan_test_settings.is_validation_shader = false;
      segscan_test_settings.frame_counter = 0;
      segscan_test_settings.use_indirect_dispatch = true;
      segscan_test_settings.ssbo_scan_infos_ = buffers_.reused_ssbo_tree_scan_infos_contour_visibility_split_();
      segscan_test_settings.ssbo_in_scan_data_ = buffers_.reused_ssbo_tree_scan_input_contour_visibility_split_0_();
      segscan_test_settings.ssbo_out_scan_data_ = buffers_.reused_ssbo_tree_scan_output_contour_visibility_split_0_();
      segscan_test_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      segscan_test_settings.shader_upsweep = eShaderType::SEGSCAN_UINT_ADD_UPSWEEP;
      segscan_test_settings.shader_aggregate = eShaderType::SEGSCAN_UINT_ADD_AGGREGATE;
      segscan_test_settings.shader_dwsweep = eShaderType::SEGSCAN_UINT_ADD_DWSWEEP;

      append_subpass_segscan(segscan_test_settings, pass_process_contours);
    }

    {
      ScanSettings segscan_test_settings;
      segscan_test_settings.is_validation_shader = false;
      segscan_test_settings.frame_counter = 0;
      segscan_test_settings.use_indirect_dispatch = true;
      segscan_test_settings.ssbo_scan_infos_ =
          buffers_.reused_ssbo_tree_scan_infos_contour_visibility_split_();
      segscan_test_settings.ssbo_in_scan_data_ =
          buffers_.reused_ssbo_tree_scan_input_contour_visibility_split_1_();
      segscan_test_settings.ssbo_out_scan_data_ =
          buffers_.reused_ssbo_tree_scan_output_contour_visibility_split_1_();
      segscan_test_settings.ssbo_scan_block_sum_ = buffers_.ssbo_scan_block_sum_;
      segscan_test_settings.shader_upsweep = eShaderType::SEGSCAN_UINT_ADD_UPSWEEP;
      segscan_test_settings.shader_aggregate = eShaderType::SEGSCAN_UINT_ADD_AGGREGATE;
      segscan_test_settings.shader_dwsweep = eShaderType::SEGSCAN_UINT_ADD_DWSWEEP;

      append_subpass_segscan(segscan_test_settings, pass_process_contours);
    }

    {
      auto &sub = pass_process_contours.sub("strokegen_contour_visibility_split_step_0");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_VISIBILITY_SPLIT_STEP_0));

      bind_rsc(sub);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_contour_visibility_split_info_());

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    {
      auto &sub = pass_process_contours.sub("strokegen_contour_visibility_split_step_1");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_VISIBILITY_SPLIT_STEP_1));

      bind_rsc(sub);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_contour_visibility_split_info_());
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_frag_seg_head_to_visibility_split_contour_());
      sub.bind_ssbo(ssbo_offset + 2, buffers_.reused_ssbo_tree_scan_output_contour_visibility_split_0_());
      sub.bind_ssbo(ssbo_offset + 3, buffers_.reused_ssbo_tree_scan_output_contour_visibility_split_1_());
      sub.bind_ssbo(ssbo_offset + 4, buffers_.ssbo_contour_to_contour_);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_frag_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND |
                  GPU_BARRIER_SHADER_IMAGE_ACCESS | GPU_BARRIER_TEXTURE_FETCH);
    }
    // step 1 pass has generated new contour edges
    append_subpass_fill_dispatch_args_contour_edges(pass_process_contours, true);

    {
      auto &sub = pass_process_contours.sub("strokegen_contour_visibility_split_step_2");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_VISIBILITY_SPLIT_STEP_2));

      bind_rsc(sub);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_contour_visibility_split_info_());
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_frag_seg_head_to_visibility_split_contour_());
      sub.bind_ssbo(ssbo_offset + 2, buffers_.ssbo_contour_to_contour_);
      sub.bind_ssbo(ssbo_offset + 3, buffers_.ssbo_contour_edge_transfer_data_);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
    {
      auto &sub = pass_process_contours.sub("strokegen_contour_visibility_split_step_3");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_VISIBILITY_SPLIT_STEP_3));

      bind_rsc(sub);
      sub.bind_ssbo(ssbo_offset + 0, buffers_.reused_ssbo_contour_visibility_split_info_());
      sub.bind_ssbo(ssbo_offset + 1, buffers_.reused_ssbo_frag_seg_head_to_visibility_split_contour_());
      sub.bind_ssbo(ssbo_offset + 2, buffers_.ssbo_contour_to_contour_);
      sub.bind_ssbo(ssbo_offset + 3, buffers_.ssbo_contour_edge_transfer_data_);

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
  }

  void StrokeGenPassModule::append_subpass_fill_contour_list_ranking_inputs()
  {
    auto &sub = pass_process_contours.sub("strokegen_fill_cotour_edge_ranking_inputs");
    sub.shader_set(shaders_.static_shader_get(FILL_CONTOUR_EDGE_RANKING_INPUTS));

    sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_list_ranking_inputs_);

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND); 
  }

  void StrokeGenPassModule::rebuild_pass_compress_contour_pixels(bool debug)
  {
    pass_compress_contour_pixels.init();
    {
      auto &sub = pass_compress_contour_pixels.sub("pass_compress_contour_pixels");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_PIXEL_COMPRESS));

      sub.bind_image(0, textures_.tex_contour_raster);
      sub.bind_image(1, textures_.tex2d_contour_pix_marks_);
      sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res());

      int2 compressed_tex_res = int2(
          GPU_texture_width(textures_.tex2d_contour_pix_marks_),
          GPU_texture_height(textures_.tex2d_contour_pix_marks_)
      );
      sub.dispatch(int3(compressed_tex_res.x, compressed_tex_res.y, 1));
      sub.barrier(eGPUBarrier::GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS); 
    }
    if (debug)
    {
      auto &sub = pass_compress_contour_pixels.sub("pass_compress_contour_pixels_debug");
      sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_PIXEL_COMPRESS_DBG));

      sub.bind_image(0, textures_.tex_contour_raster);
      sub.bind_image(1, textures_.tex2d_contour_pix_marks_);
      sub.bind_image(2, textures_.tex2d_contour_dbg_);
      sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res());

      int2 screen_res = int2(GPU_texture_width(textures_.tex_contour_raster),
                             GPU_texture_height(textures_.tex_contour_raster));
      sub.dispatch(int3(screen_res.x, screen_res.y, 1));
      sub.barrier(eGPUBarrier::GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }


  void StrokeGenPassModule::append_subpass_scan(StrokeGenPassModule::ScanSettings scan_settings, PassSimple& pass)
  {
    if (scan_settings.is_validation_shader) pass.init();

    // prep tree-scan args
    if (scan_settings.use_indirect_dispatch) {
      auto &sub = pass.sub("strokegen_scan_fill_dispatch_args");
      sub.shader_set(shaders_.static_shader_get(SCAN_FILL_DISPTACH_ARGS));

      sub.bind_ssbo(0, scan_settings.ssbo_scan_infos_);
      sub.bind_ssbo(1, buffers_.ssbo_scan_dispatch_args_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    { // upsweep for tree-scan
      auto &sub = pass.sub("strokegen_scan_test_upsweep");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_upsweep));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, scan_settings.ssbo_in_scan_data_);
      sub.bind_ssbo(1, scan_settings.ssbo_out_scan_data_);
      sub.bind_ssbo(2, scan_settings.ssbo_scan_block_sum_);
      if (scan_settings.use_indirect_dispatch)
        sub.bind_ssbo(3, scan_settings.ssbo_scan_infos_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      if (scan_settings.use_indirect_dispatch)
        sub.dispatch(buffers_.ssbo_scan_dispatch_args_);
      else
        sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));

      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
    {  // reduction for tree-scan
      auto &sub = pass.sub("strokegen_scan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_aggregate));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, scan_settings.ssbo_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    {  // down sweep for tree-scan
      auto &sub = pass.sub("strokegen_scan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_dwsweep));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, scan_settings.ssbo_out_scan_data_);
      sub.bind_ssbo(1, scan_settings.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      if (scan_settings.use_indirect_dispatch)
        sub.dispatch(buffers_.ssbo_scan_dispatch_args_);
      else
        sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));

      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  // TODO: we need to have single-pass version for this
  void StrokeGenPassModule::append_subpass_segscan(
    ScanSettings scan_settings, 
    PassSimple &pass)
  {
    if (scan_settings.is_validation_shader) pass.init();

    // prep tree-scan args
    if (scan_settings.use_indirect_dispatch)
    {
      auto& sub = pass.sub("strokegen_scan_fill_dispatch_args");
      sub.shader_set(shaders_.static_shader_get(SCAN_FILL_DISPTACH_ARGS));

      sub.bind_ssbo(0, scan_settings.ssbo_scan_infos_);
      sub.bind_ssbo(1, buffers_.ssbo_scan_dispatch_args_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND); 
    }

    { // upsweep for tree-scan
      auto& sub = pass.sub("strokegen_segscan_test_upsweep");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_upsweep));

      sub.bind_ssbo(0, scan_settings.ssbo_in_scan_data_);
      sub.bind_ssbo(1, scan_settings.ssbo_out_scan_data_);
      sub.bind_ssbo(2, scan_settings.ssbo_scan_block_sum_);
      if (scan_settings.use_indirect_dispatch)
        sub.bind_ssbo(3, scan_settings.ssbo_scan_infos_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      if (scan_settings.is_validation_shader)
        sub.push_constant("pcs_segscan_test_random_seed_", scan_settings.frame_counter);  

      if (!scan_settings.use_indirect_dispatch)
        sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      else {
        sub.dispatch(buffers_.ssbo_scan_dispatch_args_);
      }

      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
    { // reduction for tree-scan
      auto& sub = pass.sub("strokegen_segscan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_aggregate));

      sub.bind_ssbo(0, scan_settings.ssbo_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // down sweep for tree-scan
      auto& sub = pass.sub("strokegen_segscan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(scan_settings.shader_dwsweep));

      sub.bind_ssbo(0, scan_settings.ssbo_out_scan_data_);
      sub.bind_ssbo(1, scan_settings.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      if (!scan_settings.use_indirect_dispatch)
        sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      else{
        sub.dispatch(buffers_.ssbo_scan_dispatch_args_);
      }

      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::append_subpass_segloopconv1d(
      SegLoopConv1DSettings settings, PassSimple& pass)
  {
    if (settings.is_validation_shader)
      pass.init();

    if (settings.use_indirect_dispatch) {
      auto &sub = pass.sub("strokegen_segloopconv1d_fill_dispatch_args");
      sub.shader_set(shaders_.static_shader_get(CONV1D_FILL_DISPATCH_ARGS));

      sub.bind_ssbo(0, settings.ssbo_segloopconv1d_info_);
      sub.bind_ssbo(1, buffers_.ssbo_segloopconv1d_dispatch_args_);
      sub.push_constant("pc_segloopconv1d_dispatch_group_size_", (int)GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }

    if (!settings.lazy_dispatch)
    {
      auto &sub = pass.sub("strokegen_segloopconv1D_build_patch");
      sub.shader_set(shaders_.static_shader_get(settings.shader_build_patch_table));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, settings.ssbo_segloopconv1d_patch_table_);
      sub.bind_ssbo(1, settings.ssbo_segloopconv1d_info_);
      if (settings.is_validation_shader)
        sub.bind_ssbo(2, buffers_.ssbo_debug_segloopconv1d_data_);
      else if (settings.shader_build_patch_table == CONV1D_SEG_DENOISE_BUILD_PATCH) {
        sub.bind_ssbo(2, buffers_.ssbo_contour_snake_rank_);
        sub.bind_ssbo(3, buffers_.ssbo_contour_snake_list_len_);
        sub.bind_ssbo(4, buffers_.ssbo_contour_snake_flags_);
      }
      else if (settings.shader_build_patch_table == CONV1D_2D_SAMPLE_BUILD_PATCH) {
        sub.bind_ssbo(2, buffers_.reused_ssbo_contour_2d_sample_topology_());
        sub.bind_ssbo(3, buffers_.reused_ssbo_contour_2d_sample_geometry_());
        sub.bind_image(0, textures_.tex2d_contour_dbg_); 
        sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res()); 
      }
      sub.bind_ubo(0, buffers_.ubo_segloopconv1d_);

      if (settings.use_indirect_dispatch)
        sub.dispatch(buffers_.ssbo_segloopconv1d_dispatch_args_);
      else
        sub.dispatch(int3(buffers_.ubo_segloopconv1d_.num_thread_groups, 1, 1));

      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS/*for dbg texture*/);
    }
    {
      auto &sub = pass.sub("strokegen_segloopconv1D_convolution");
      sub.shader_set(shaders_.static_shader_get(settings.shader_convolution));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, settings.ssbo_segloopconv1d_patch_table_);
      sub.bind_ssbo(1, settings.ssbo_in_segloopconv1d_data_);
      sub.bind_ssbo(2, settings.ssbo_out_segloopconv1d_data_);
      sub.bind_ssbo(3, settings.ssbo_segloopconv1d_info_); 
      if (settings.is_validation_shader)
        sub.bind_ssbo(4, buffers_.ssbo_debug_segloopconv1d_data_);
      else if (settings.shader_convolution == CONV1D_SEG_DENOISE_CONVOLUTION) {
        sub.bind_ssbo(4, buffers_.ssbo_contour_snake_rank_);
        sub.bind_ssbo(5, buffers_.ssbo_contour_snake_list_len_);
      }
      else if (settings.shader_convolution == CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_0
                || settings.shader_convolution == CONV1D_2D_SAMPLE_CORNER_CONVOLUTION_STEP_1
                || settings.shader_convolution == CONV1D_2D_SAMPLE_CALC_TANGENT_CURVATURE
                || settings.shader_convolution == CONV1D_2D_SAMPLE_DENOISE_VISIBILITY) {
        sub.bind_ssbo(4, buffers_.reused_ssbo_contour_2d_sample_topology_());
        sub.bind_ssbo(5, buffers_.reused_ssbo_contour_2d_sample_geometry_());
        sub.bind_image(0, textures_.tex2d_contour_dbg_); 
        sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res()); 
      }
      sub.bind_ubo(0, buffers_.ubo_segloopconv1d_);

      if (settings.use_indirect_dispatch)
        sub.dispatch(buffers_.ssbo_segloopconv1d_dispatch_args_);
      else
        sub.dispatch(int3(buffers_.ubo_segloopconv1d_.num_thread_groups, 1, 1));

      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS /*for dbg texture*/);
    }
  }

  void StrokeGenPassModule::rebuild_pass_list_ranking_pointer_jumping(
      PassSimple &pass_listranking,
      int num_splice_iters,
      const int num_jump_iters, int jumping_info_offset = 0,
      bool loop_breaking_pass = false, bool loop_ranking_pass = false
    ){
    // Spliced Pointer Jumping
    for (int jump_iter = 0; jump_iter < num_jump_iters; jump_iter++)
    {
      auto& sub = pass_listranking.sub("strokegen_list_ranking_test_sublist_pointer_jumping");
      if (!(loop_breaking_pass || loop_ranking_pass))
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_SUBLIST_POINTER_JUMPING));
      else // we are dealing with looped list(s), and we need one more jumping pass to determine head&tail nodes
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_LOOPED_POINTER_JUMPING));

      int jumping_info_subbuff_in = jump_iter + jumping_info_offset;
      if (loop_ranking_pass && jump_iter > 0)
        jumping_info_subbuff_in += 1; // we do nothing in the 0th jumping pass hence did not swap the buffers
      sub.bind_ssbo(0, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[(jumping_info_subbuff_in)%2]);
      sub.bind_ssbo(1, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[(jumping_info_subbuff_in+1)%2]);
      sub.bind_ssbo(2, buffers_.ssbo_list_ranking_node_to_anchor_);
      sub.bind_ssbo(3, buffers_.ssbo_list_ranking_anchor_to_node_[(num_splice_iters)%2]); // note: double check this
      sub.bind_ssbo(4, buffers_.ssbo_list_ranking_links_);
      sub.bind_ssbo(5, buffers_.ssbo_list_ranking_ranks_);
      sub.bind_ssbo(6, buffers_.ssbo_list_ranking_anchor_counters_);
      sub.bind_ssbo(7, buffers_.ssbo_list_ranking_serialized_topo_);
      sub.bind_ssbo(8, buffers_.ssbo_list_ranking_inputs_);
      sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
      sub.push_constant("pc_listranking_splice_iter_", num_splice_iters);
      sub.push_constant("pc_listranking_jumping_iter_", jump_iter);
      int is_looped_ranking_pass = loop_ranking_pass ? 1 : 0;
      sub.push_constant("pc_listranking_ranking_pass_with_broken_loops_", is_looped_ranking_pass);

      // note: double check this
      sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[num_splice_iters]);

      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
  }

  // Note: we must do this within 15 ms(0.64ms x 23iters from willey's algo.), otherwise it for nothing.
  void StrokeGenPassModule::append_subpass_list_ranking(ListRankingPassUsage passType, 
                                                      PassSimple &pass_listranking,
                                                      bool looped_pass_list_ranking
  )
  {
    // Build render passes
    bool custom_pass = (passType != ListRankingPassUsage::TestListRanking); 

    int num_splice_iters = 3;
    for (int splice_iter = 0; splice_iter < num_splice_iters; ++splice_iter)
    {
      rebuild_pass_list_ranking_fill_args(pass_listranking, true, false, splice_iter, (int)GROUP_SIZE_BNPR_LIST_RANK_TEST, custom_pass);

      if (splice_iter == 0)
      {
        { // copy from staging buffer,
          auto &sub = pass_listranking.sub("strokegen_list_ranking_test_upload_cpu_data");
          sub.shader_set(shaders_.static_shader_get(LISTRANKING_SETUP_INPUTS));
          if (!custom_pass)
            sub.bind_ssbo(0, buffers_.ssbo_list_ranking_links_staging_buf_);
          else if (passType == ListRankingPassUsage::ContourEdgeLinking)
            sub.bind_ssbo(0, buffers_.ssbo_contour_to_contour_); 
          sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
          sub.bind_ssbo(2, buffers_.ssbo_list_ranking_inputs_); 
          sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
          sub.push_constant("pc_listranking_custom_", custom_pass ? 1 : 0); 

          sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[0]);
          sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
        }
      }

      int tag_iter;
      int num_tag_iters = splice_iter == 0 ? 3 : splice_iter == 1 ? 2 : 1;
      for (tag_iter = 0; tag_iter < num_tag_iters; ++tag_iter)
      {
        {
          auto &sub = pass_listranking.sub("strokegen_list_ranking_test_tagging");
          sub.shader_set(shaders_.static_shader_get(LISTRANKING_INIT_ANCHORS));
          sub.bind_ssbo(0, buffers_.ssbo_list_ranking_tags_[tag_iter % 2]);
          sub.bind_ssbo(1, buffers_.ssbo_list_ranking_tags_[(tag_iter + 1) % 2]);
          sub.bind_ssbo(2, buffers_.ssbo_list_ranking_anchor_to_node_[splice_iter % 2]);
          sub.bind_ssbo(3, buffers_.ssbo_list_ranking_anchor_to_node_[(splice_iter + 1) % 2]);
          sub.bind_ssbo(4, buffers_.ssbo_list_ranking_links_);
          sub.bind_ssbo(5, buffers_.ssbo_list_ranking_anchor_counters_);
          sub.bind_ssbo(6, buffers_.ssbo_list_ranking_splice_counters_);
          sub.bind_ssbo(7, buffers_.ssbo_list_ranking_ranks_);
          sub.bind_ssbo(8, buffers_.ssbo_list_ranking_addressing_counters_);
          sub.bind_ssbo(9, buffers_.ssbo_list_ranking_inputs_); 
          sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
          sub.push_constant("pc_listranking_splice_iter_", splice_iter);
          sub.push_constant("pc_listranking_tagging_iter_", tag_iter);

          sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splice_iter]);
          sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
        }
      }

      {
        auto& sub = pass_listranking.sub("strokegen_list_ranking_test_compaction");
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_COMPACT_ANCHORS));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_tags_[tag_iter % 2]);
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_anchor_to_node_[splice_iter % 2]); // pingpong
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_anchor_to_node_[(splice_iter + 1) % 2]);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_spliced_node_id_[splice_iter]);
        sub.bind_ssbo(6, buffers_.ssbo_list_ranking_anchor_to_next_anchor_);
        sub.bind_ssbo(7, buffers_.ssbo_list_ranking_node_to_anchor_);
        sub.bind_ssbo(8, buffers_.ssbo_list_ranking_anchor_counters_);
        sub.bind_ssbo(9, buffers_.ssbo_list_ranking_splice_counters_);
        sub.bind_ssbo(10, buffers_.ssbo_list_ranking_inputs_);
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", splice_iter);
        sub.push_constant("pc_num_splice_iters_", num_splice_iters);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splice_iter]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }

      int splice_disptatch_args_slot = splice_iter + 1;
      rebuild_pass_list_ranking_fill_args(pass_listranking, false, true, splice_disptatch_args_slot, (int)GROUP_SIZE_BNPR_LIST_RANK_TEST, custom_pass);
      {
        auto& sub = pass_listranking.sub("strokegen_list_ranking_test_splice_out_nodes");
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_SPLICE_OUT_NODES));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_tags_[tag_iter % 2]);
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_anchor_to_node_[splice_iter % 2]); // pingpong
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_anchor_to_node_[(splice_iter + 1) % 2]);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_spliced_node_id_[splice_iter]);
        sub.bind_ssbo(6, buffers_.ssbo_list_ranking_anchor_to_next_anchor_);
        sub.bind_ssbo(7, buffers_.ssbo_list_ranking_node_to_anchor_);
        sub.bind_ssbo(8, buffers_.ssbo_list_ranking_anchor_counters_);
        sub.bind_ssbo(9, buffers_.ssbo_list_ranking_splice_counters_);
        sub.bind_ssbo(10, buffers_.ssbo_list_ranking_inputs_);
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", splice_iter);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[splice_disptatch_args_slot]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }
    }
    rebuild_pass_list_ranking_fill_args(pass_listranking, true, false, num_splice_iters, GROUP_SIZE_BNPR_LIST_RANK_TEST, custom_pass);

    const int num_jump_iters = MAX_NUM_JUMPS_BNPR_LIST_RANK_TEST;
    if (!looped_pass_list_ranking){
      rebuild_pass_list_ranking_pointer_jumping(pass_listranking, num_splice_iters, num_jump_iters);
    }else
    {
      rebuild_pass_list_ranking_pointer_jumping(
          pass_listranking, num_splice_iters, num_jump_iters, 0, true, false
          );
      {
        auto& sub = pass_listranking.sub("strokegen_list_ranking_test_mark_loop_head_tail");
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_MARK_LOOP_HEAD_TAIL));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_anchor_to_node_[num_splice_iters%2]);
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[num_jump_iters%2]);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[(num_jump_iters+1)%2]);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_node_to_anchor_);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(6, buffers_.ssbo_list_ranking_anchor_counters_);
        sub.bind_ssbo(7, buffers_.ssbo_list_ranking_inputs_); 
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", num_splice_iters);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[num_splice_iters]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }
      rebuild_pass_list_ranking_pointer_jumping(
          pass_listranking, num_splice_iters, num_jump_iters, 1, false, true);
    }

    { // Relinking spliced nodes back to the list
      int num_relink_iters = num_splice_iters;
      int input_jump_cache_pingpong = 0;
      for (int relink_iter = 0; relink_iter < num_relink_iters; relink_iter++)
      {
        auto& sub = pass_listranking.sub("strokegen_list_ranking_test_relink");
        if (!looped_pass_list_ranking)
          sub.shader_set(shaders_.static_shader_get(LISTRANKING_RELINKING));
        else
          sub.shader_set(shaders_.static_shader_get(LISTRANKING_LOOPED_RELINKING));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_spliced_node_id_[num_relink_iters - 1 - relink_iter]);  // note: double check this
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_splice_counters_);
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[input_jump_cache_pingpong]);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_serialized_topo_);
        sub.bind_ssbo(6, buffers_.ssbo_list_ranking_inputs_); 
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_relink_iter_", relink_iter);
        sub.push_constant("pc_listranking_num_relink_iters_", num_relink_iters);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[num_relink_iters - relink_iter]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }
    }

    // Output list rank, len, addr to custom buffers
    if (passType != ListRankingPassUsage::TestListRanking)
    {
      for (int i = 0; i < 2; ++i)
      {
        auto &sub = pass_listranking.sub("strokegen_list_ranking_test_output");
        sub.shader_set(shaders_.static_shader_get(
          i == 0 ? LISTRANKING_OUTPUT_DATA_PASS_0 : LISTRANKING_OUTPUT_DATA_PASS_1
        ));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_serialized_topo_);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_inputs_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_addressing_counters_);
        if (passType == ListRankingPassUsage::ContourEdgeLinking) {
          sub.bind_ssbo(4, buffers_.reused_ssbo_contour_edge_rank_());
          sub.bind_ssbo(5, buffers_.reused_ssbo_contour_edge_list_len_());
          sub.bind_ssbo(6, buffers_.reused_ssbo_contour_edge_list_head_info_());
        }
        sub.push_constant("pcs_contour_edge_linking_output_pass_type_", passType); 

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[0]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
      }
    }
  }


  /**
   * \brief Add a sub pass for filling indirect dispatch args to the list ranking pass.
   * \param per_anchor true if kernel is dispatched x1 thread per anchor node.
   * \param per_spliced true if kernel is dispatched x1 thread per spliced node.
   * \param splicing_or_relinking_iter the iteration of splicing or relinking process.
   */
  void StrokeGenPassModule::rebuild_pass_list_ranking_fill_args(
      PassSimple& pass_listranking, bool per_anchor, bool per_spliced, int splicing_or_relinking_iter,
      int group_size_x, bool custom_pass
      ){
    auto& sub = pass_listranking.sub("strokegen_list_ranking_test_fill_dispatch_args");
    sub.shader_set(shaders_.static_shader_get(eShaderType::LISTRANKING_FILL_DISPATCH_ARGS));

    sub.bind_ssbo(0, buffers_.ssbo_list_ranking_anchor_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_list_ranking_splice_counters_);
    sub.bind_ssbo(2, buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splicing_or_relinking_iter]);
    sub.bind_ssbo(3, buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[splicing_or_relinking_iter]);
    sub.bind_ssbo(4, buffers_.ssbo_list_ranking_inputs_);
    sub.push_constant("pc_listranking_custom_", custom_pass ? 1 : 0); 
    int dispatch_granularity_tag = per_anchor ? 0 : (per_spliced ? 1 : 2);
    sub.push_constant("pc_listranking_indirect_arg_granularity_", dispatch_granularity_tag); /* 0 := anchors, 1:= spliced nodes */
    sub.push_constant("pc_listranking_counter_buffer_slot_id_", splicing_or_relinking_iter); /* which counter to use */
    sub.push_constant("pc_listranking_dispatch_group_size_", group_size_x); /* which counter to use */

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
  }

  void StrokeGenPassModule::print_list_ranking_nodes(
    int head_node_id,
    uint* computed_ranks, uint* computed_topo, uint* computed_links
  ) const
{
    int list_len_gt = buffers_.listranking_test_nodes_list_len[head_node_id];
    uint list_len_out = computed_topo[2 * head_node_id + 1];
    std::cout << std::left // output each list in a column
      << "length: " << list_len_gt << " | " << list_len_out << " | " << std::endl;
    int curr_node_id = head_node_id;

    for (int i = 0; i < list_len_gt; ++i)
    {
      std::stringstream out_str;
      {
        uint rank_gt = buffers_.listranking_test_nodes_rank[curr_node_id];

        uint rank_out = computed_ranks[curr_node_id];
        uint list_len_out = computed_topo[2 * curr_node_id + 1];
        uint list_addr_out = computed_topo[2 * curr_node_id];
        uint list_broadcast_anchor_id = computed_links[2 * curr_node_id + 1];
        // if (test_looped_pass_list_ranking)
        {
          rank_out = rank_out & 0x3fffffffu;
        }
        int rank_out_adjusted = (list_len_gt - rank_out);

        out_str << curr_node_id << "_r" << rank_out_adjusted << "_l" << list_len_out << "_b" << list_broadcast_anchor_id
          << "_R" << rank_gt;

        std::cout << std::left << std::setw(20) << out_str.str();
      }

      if ((i + 1) != list_len_gt)
        std::cout << std::left << std::setw(2) << "   ->  ";

      if (i % 4 == 1)
        std::cout << std::endl;

      curr_node_id = buffers_.listranking_test_nodes_prev_next[curr_node_id * 2 + 1];
    }
    std::cout << std::endl;
  }

  bool StrokeGenPassModule::validate_list_ranking() const
  {
    SSBO_ListRankingRanks& buf_final_ranks = buffers_.ssbo_list_ranking_ranks_;
    buf_final_ranks.read();
    uint* computed_ranks = reinterpret_cast<uint *>(buf_final_ranks.data());

    SSBO_ListRankingSerializedTopo& buf_final_topo = buffers_.ssbo_list_ranking_serialized_topo_;
    buf_final_topo.read();
    uint* computed_topo = reinterpret_cast<uint *>(buf_final_topo.data());

    SSBO_ListRankingLinks& buf_computed_links = buffers_.ssbo_list_ranking_links_;
    buf_computed_links.read(); // .next should be pointing to anchor id for broadcasting topo data
    uint* computed_links = reinterpret_cast<uint *>(buf_computed_links.data());

    // Print linked lists
    const uint num_nodes = buffers_.ubo_list_ranking_splicing_.num_nodes;
    if (num_nodes <= 64)
    {
      for (int head_node_id : buffers_.listranking_test_head_nodes)
      {
        print_list_ranking_nodes(head_node_id, computed_ranks, computed_topo, computed_links);
      }
    }


    if (!test_looped_pass_list_ranking)
      for (size_t i = 0; i < num_nodes; ++i)
      {
        uint rank_gt      = buffers_.listranking_test_nodes_rank[i];
        uint list_len_gt  = buffers_.listranking_test_nodes_list_len[i];
        uint rank_out     = computed_ranks[i];
        uint list_len_out = computed_topo[2*i+1];

        uint head_node_id_gt = buffers_.listranking_test_nodes_head[i]; 

        rank_out = rank_out & 0x3fffffffu;

        if (list_len_gt != list_len_out) {
          fprintf(stderr, "strokegen error: incorrect list length, list ranking test failed: ");
          print_list_ranking_nodes(head_node_id_gt, computed_ranks, computed_topo, computed_links);

          return false;
        }

        if (rank_gt != (list_len_gt - 1 - rank_out))
        {
          fprintf(stderr, "strokegen error: incorrect rank, list ranking test failed: ");
          print_list_ranking_nodes(head_node_id_gt, computed_ranks, computed_topo, computed_links);

          return false;
        }
      }
    else // looped lists
    {
      for (int head_node_id : buffers_.listranking_test_head_nodes)
      { // now that I understood the heads computed here can be different from GPU outputs
        // since in shader the head is computed from anchors with max node id,
        // and head-anchors may not have the max node id among ALL list nodes,
        // which was the criterion I pick head in the CPU code.
        uint list_len_gt = buffers_.listranking_test_nodes_list_len[head_node_id];
        uint prev_list_len_out = computed_topo[2 * head_node_id + 1];
        uint head_list_addr_out = computed_topo[2 * head_node_id];
        if (list_len_gt != prev_list_len_out)
          fprintf(stderr, "incorrect list len");

        uint prev_node_id = head_node_id;
        uint prev_rank_out = computed_ranks[prev_node_id];

        auto decode_rank = [](uint raw_rank, uint list_len) {
          // decode rank, pack with flags
          raw_rank = raw_rank & 0x3fffffffu;
          raw_rank = list_len - raw_rank;
          return raw_rank;
        };

        auto valid_rank = [](uint rank, uint rank_pre, uint list_len)
        {
          return (rank < list_len) && (0 <= rank) && (rank_pre < list_len) && (0 <= rank_pre)
            && (rank == ((rank_pre + 1 + list_len) % list_len));
        };

        for (int i = 0; (i + 1) < list_len_gt; ++i) {
          uint curr_node_id = buffers_.listranking_test_nodes_prev_next[prev_node_id * 2 + 1];

          bool fucked = false;

          uint curr_list_len_out = computed_topo[2 * curr_node_id + 1];
          if (curr_list_len_out != list_len_gt)
          {
            fprintf(stderr, "incorrect list len");
            fucked = true;
          }

          uint curr_list_addr_out = computed_topo[2 * curr_node_id];
          if (curr_list_addr_out != head_list_addr_out)
          {
            fprintf(stderr, "incorrect list addr");
            fucked = true;
          }

          uint curr_rank_out = computed_ranks[curr_node_id];
          curr_rank_out = (curr_rank_out & 0x3fffffffu); 
          if (curr_rank_out > list_len_gt || curr_rank_out > curr_list_len_out) {
            fprintf(stderr, "incorrect list addr");
            fucked = true;
          }

          uint curr_rank_out_adjusted = decode_rank(curr_rank_out, list_len_gt);
          uint prev_rank_out_adjusted = decode_rank(prev_rank_out, list_len_gt);
          if (false == valid_rank(curr_rank_out_adjusted,
                                  prev_rank_out_adjusted,
                                  list_len_gt)
          ) {
            fprintf(stderr, "incorrect list rank");
            fucked = true;
          }

          if (fucked) {
            print_list_ranking_nodes(head_node_id, computed_ranks, computed_topo, computed_links);
            break;
          }

          prev_node_id = curr_node_id;
          prev_rank_out = curr_rank_out;
        }
        std::cout << std::endl;
      }
    }


    return true;
  }
}
