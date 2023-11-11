/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 */

#include "npr_strokegen_pass.hh"

#include "gpu_batch_private.hh"

#include <iomanip>

namespace blender::npr::strokegen
{
  using namespace blender;

  PassSimple& StrokeGenPassModule::get_compute_pass(eType passType)
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
        return pass_extract_geom;
      case CONTOUR_PROCESS:
        return pass_process_contours; 
      case COMPRESS_CONTOUR_PIXELS:
        return pass_compress_contour_pixels; 
    }
    return pass_comp_test;
  }

  PassMain &StrokeGenPassModule::get_contour_edge_draw_pass(eType passType)
  {
    switch (passType) {
      case INDIRECT_DRAW_CONTOUR_EDGES:
        return pass_draw_contour_edges;     
    }
    return pass_draw_contour_edges; 
  }



  void StrokeGenPassModule::on_begin_sync()
  {
    init_per_mesh_pass();
    pass_draw_contour_edges.init_pass(shaders_, textures_); 

    rebuild_pass_scan_test();
    rebuild_pass_segscan_test();
    rebuild_pass_conv_test();
    append_subpass_list_ranking(ListRankingPassType::Test, pass_listranking_test, true);

    num_total_mesh_tris = num_total_mesh_verts = num_total_mesh_edges = 0; // TODO: these should go to the UBO
  }

  void StrokeGenPassModule::on_end_sync()
  {
    /* Post processing after iterated over all meshes. */
    rebuild_pass_process_contours(); 
    rebuild_pass_contour_edge_drawcall();
    rebuild_pass_compress_contour_pixels();
  }


  /**
   * \brief In each render loop, re-initialize the compute pass for geometry extraction.
   */
  void StrokeGenPassModule::init_per_mesh_pass()
  {
    pass_extract_geom.init();
    boostrap_before_extract_first_batch = true;
  }

  /**
   * \brief Add a subpass for extracting geometry from given GPUBatch.
   * \param ob Mesh Object
   * \param gpu_batch_line_adj Mesh geometry stored in GPUBatch, ib stored with line adjacency info.
   */
  void StrokeGenPassModule::append_per_mesh_pass(
      Object* ob,
      GPUBatch* gpu_batch_line_adj,
      GPUBatch* gpu_batch_surf, 
      ResourceHandle& rsc_handle,
      const DRWView* drw_view
      ){

    if (boostrap_before_extract_first_batch) {
      auto &sub = pass_extract_geom.sub("bootstrap meshing passes");
      sub.shader_set(shaders_.static_shader_get(GPU_MESHING_BOOSTRAP));

      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_prev_); 
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      sub.dispatch(int3(1, 1, 1));

      return; // bootstrapping done. will re-enter this func for actual work
    }


    int batch_resource_index = (int)rsc_handle.resource_index(); 

    /* see example in "GPU_batch_draw_parameter_get" */
    gpu::Batch* edge_batch = static_cast<gpu::Batch*>(gpu_batch_line_adj);
    if (edge_batch == nullptr || edge_batch->elem == nullptr)
      fprintf(stderr, "StrokeGen Error: empty mesh when rebuilding pass 'extract_mesh_geom'");
    int num_edges = edge_batch->elem_()->index_len_get() / 4; // 4 indices per primitive, basically 4 verts around the edge

    gpu::IndexBuf* ib = edge_batch->elem_();
    /* Hack to get ibo format */
    gpu::GPUIndexBufType ib_type = // the actual "index_type_" is protected
        (ib->size_get() / ib->index_len_get() == sizeof(uint32_t))
          ? gpu::GPU_INDEX_U32 : gpu::GPU_INDEX_U16;
    gpu::VertBuf *vbo_edge_adj = edge_batch->verts_(0);
    int num_verts_edge_vbo = vbo_edge_adj->vertex_len; 


    gpu::Batch *batch_surf = static_cast<gpu::Batch *>(gpu_batch_surf); 
    gpu::IndexBuf *ibo_surf = batch_surf->elem_();
    /* Layout of multi vbos, see
     * "DRW_batch_requested(cache->batch.surface, GPU_PRIM_TRIS)" in
     * "DRW_mesh_batch_cache_create_requested"
     *
     * For the data layout of "posnor" vbo, see
     * "extract_mesh_vbo_pos_nor.cc"
     */
    gpu::VertBuf *vbo_surf_posnor = batch_surf->verts_(1);
    GPUVertBuf *vbo_surf_posnor_gpuverbuf = gpu_batch_surf->verts[1];
    int num_tris = ibo_surf->index_len_get() / 3;
    int num_verts = vbo_surf_posnor->vertex_len; 

    /* Cache per-edge geometry data */
    // Note: this also works
    /*GPU_storagebuf_copy_sub_from_vertbuf(buffers_.ssbo_vbo_full_,
                                         gpu_batch_surf->verts[0],
                                         4 * num_total_mesh_verts,
                                         0,
                                         4 * num_verts);*/
    append_subpass_merge_vbo(gpu_batch_surf, batch_resource_index, num_verts);
    append_subpass_merge_line_adj_ibo(gpu_batch_line_adj, num_edges, ib_type);

    append_subpass_meshing_merge_verts(num_verts);
    append_subpass_meshing_edge_adjacency(num_edges); 

    append_subpass_extract_contour_edges(gpu_batch_line_adj, rsc_handle, edge_batch, num_edges, ib_type);

    append_subpass_fill_dispatch_args_contour_edges(pass_extract_geom, false);
    append_subpass_process_contour_edges();


    // Note: this should be  appening at the end
    num_total_mesh_tris += num_tris;
    num_total_mesh_verts += num_verts;
    num_total_mesh_edges += num_edges; 
  }

  void StrokeGenPassModule::append_subpass_merge_vbo(GPUBatch *gpu_batch_surf, int batch_resource_index, int num_verts)
  {
    auto &sub = pass_extract_geom.sub("merge batch vbo");
      
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
    sub.push_constant("pcs_vert_id_offset_", num_total_mesh_verts);
      
    int num_groups = compute_num_groups(num_verts, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT); 
    sub.dispatch(int3(num_groups, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_merge_line_adj_ibo(
      GPUBatch *gpu_batch_line_adj,
      int num_added_edges,
      gpu::GPUIndexBufType ib_type)
  {
    auto &sub = pass_extract_geom.sub("merge batch edge adj ibo");

    if (ib_type == gpu::GPU_INDEX_U16)
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_COLLECT_EDGE_ADJ_IBO_16BIT));
    else
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_COLLECT_EDGE_ADJ_IBO));

    sub.bind_ssbo(0, &(gpu_batch_line_adj->elem));
    sub.bind_ssbo(1, buffers_.ssbo_edge_to_vert_);

    sub.push_constant("pcs_edge_count_", num_added_edges);
    sub.push_constant("pcs_vert_id_offset_", num_total_mesh_verts);
    sub.push_constant("pcs_edge_id_offset_", num_total_mesh_edges);

    int num_groups = compute_num_groups(num_added_edges, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
    sub.dispatch(int3(num_groups, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::append_subpass_meshing_merge_verts(int num_verts_in, bool debug)
  {
    const int hashmap_size = std::min(MAX_GPU_HASH_TABLE_SIZE, num_verts_in * 16);

    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_vert_spatial_map_headers_);
      sub.bind_ssbo(1, buffers_.ssbo_mesh_buffer_reuse_0_);
      sub.bind_ssbo(2, buffers_.ssbo_vert_merged_id_);
      sub.bind_ssbo(3, buffers_.ssbo_vbo_full_);
      sub.bind_ssbo(4, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.push_constant("pcs_hash_map_size_", hashmap_size);
      sub.push_constant("pcs_vert_count_", (int)num_verts_in); 
    }; 

    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_merge_verts_bootstrap");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_INIT));

      bind_rsc(sub);
      /* TODO: use GL buffer copy function rather than this */
      int num_groups = compute_num_groups(hashmap_size, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_merge_verts_spatial_hashing");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_HASH));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_verts_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_merge_verts_deduplicate");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_VERT_MERGE_REINDEX));

      bind_rsc(sub); 

      int num_groups = compute_num_groups(num_verts_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }

  void StrokeGenPassModule::append_subpass_meshing_edge_adjacency(int num_edges_in, bool debug)
  {
    const int hashmap_size = std::min(MAX_GPU_HASH_TABLE_SIZE, num_edges_in * 16);

    auto bind_rsc = [&](draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub) {
      sub.bind_ssbo(0, buffers_.ssbo_vert_merged_id_);
      sub.bind_ssbo(1, buffers_.ssbo_vert_spatial_map_headers_);
      sub.bind_ssbo(2, buffers_.ssbo_mesh_buffer_reuse_0_);
      sub.bind_ssbo(3, buffers_.ssbo_edge_to_vert_);
      sub.bind_ssbo(4, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(5, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.push_constant("pcs_hash_map_size_", hashmap_size);
      sub.push_constant("pcs_edge_count_", (int)num_edges_in);
    };

    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_edge_adj_bootstrap");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_INIT));

      bind_rsc(sub);
      /* TODO: use GL buffer copy function rather than this */
      int num_groups = compute_num_groups(hashmap_size, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_edge_adj_hashing");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_HASH));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_edges_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
    {
      auto &sub = pass_extract_geom.sub("bnpr_meshing_edge_adj_finish");
      sub.shader_set(shaders_.static_shader_get(eShaderType::MESH_EDGE_ADJACENCY_FILL));

      bind_rsc(sub);

      int num_groups = compute_num_groups(num_edges_in, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT);
      sub.dispatch(int3(num_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }

  void StrokeGenPassModule::append_subpass_fill_dispatch_args_contour_edges(PassSimple& pass, bool all_contour_edges)
  { // fill dispatch args for contour edges
    {
      auto &sub = pass.sub("fill_dispatch_args_per_contour_edge");
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
  }


  void StrokeGenPassModule::append_subpass_extract_contour_edges(
      GPUBatch *gpu_batch_line_adj,
      ResourceHandle &rsc_handle,
      gpu::Batch *edge_batch,
      int num_edges,
      gpu::GPUIndexBufType ib_type)
  {
    auto &sub = pass_extract_geom.sub(
        boostrap_before_extract_first_batch ? "extract_geom_boostrap" : "Extract"
        );
    // if (ib_type == gpu::GPU_INDEX_U16)
    //   sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_GEOM_EXTRACT_IBO_16BIT));
    // else
    sub.shader_set(shaders_.static_shader_get(eShaderType::CONTOUR_GEOM_EXTRACT));

    sub.bind_ssbo(0, buffers_.ssbo_edge_to_vert_/*&(gpu_batch_line_adj->elem)*/);
    sub.bind_ssbo(1, buffers_.ssbo_vbo_full_ /*&(gpu_batch_line_adj->verts[0])*/);
    sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_);
    sub.bind_ssbo(3, DRW_manager_get()->matrix_buf.current());
    sub.bind_ssbo(4, buffers_.ssbo_bnpr_mesh_pool_counters_);
    sub.bind_ssbo(5, buffers_.ssbo_edge_to_contour_);
    sub.bind_ubo(0, buffers_.ubo_view_matrices_cache_);
    sub.push_constant("pcs_ib_fmt_u16", ib_type == gpu::GPU_INDEX_U16 ? 1 : 0);
    sub.push_constant("pcs_num_edges", (int)edge_batch->elem_()->index_len_get() / 4);
    sub.push_constant(
        // not used for now, but I suspect this could be interfering with complex meshes
        "pcs_num_ib_offset",
        (int)edge_batch->elem_()->index_start_get() + (int)edge_batch->elem_()->index_base_get()
        );
    sub.push_constant("pcs_rsc_handle", (int)rsc_handle.resource_index());
    sub.push_constant("pcs_vert_id_offset_", num_total_mesh_verts);
    sub.push_constant("pcs_edge_id_offset_", num_total_mesh_edges);

    sub.dispatch(int3(
            compute_num_groups(num_edges, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT, 1),
            1,
            1)
        );
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  void StrokeGenPassModule::rebuild_pass_process_contours()
  {
    pass_process_contours.init();
    append_subpass_fill_dispatch_args_contour_edges(pass_process_contours, true); 
    // append_subpass_process_contour_edges();
    append_subpass_list_ranking(
        StrokeGenPassModule::ListRankingPassType::ContourEdgeLinking,
        pass_process_contours, true
      );
  }

  void StrokeGenPassModule::append_subpass_process_contour_edges()
  {
    {
      auto &sub = pass_extract_geom.sub("calc contour edge raster data");
      sub.shader_set(shaders_.static_shader_get(eShaderType::COMPUTE_CONTOUR_EDGE_RASTER_DATA));

      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_counters_prev_);
      sub.bind_ssbo(3, buffers_.ssbo_edge_to_edges_);
      sub.bind_ssbo(4, buffers_.ssbo_edge_to_contour_);
      sub.bind_ssbo(5, buffers_.ssbo_contour_to_contour_);
      sub.bind_ssbo(6, buffers_.ssbo_list_ranking_inputs_); 
      sub.bind_ubo(0, buffers_.ubo_view_matrices_);
      float2 fb_res = textures_.get_contour_raster_screen_res(); 
      sub.push_constant("pcs_screen_size_", fb_res); 

      sub.dispatch(buffers_.ssbo_bnpr_mesh_contour_edge_dispatch_args_);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::rebuild_pass_contour_edge_drawcall()
  {
    {
      auto &sub = pass_draw_contour_edges.sub("fill_draw_args_contour_edges");

      sub.shader_set(shaders_.static_shader_get(eShaderType::FILL_DRAW_ARGS_CONTOUR_EDGES));

      sub.bind_ssbo(0, buffers_.ssbo_bnpr_mesh_pool_counters_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_mesh_pool_draw_args_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
    }
    
    pass_draw_contour_edges.append_draw_subpass(shaders_, buffers_, textures_);
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
      sub.bind_image(2, textures_.tex2d_contour_pix_marks_dbg_);
      sub.push_constant("pcs_screen_size_", textures_.get_contour_raster_screen_res());

      int2 screen_res = int2(GPU_texture_width(textures_.tex_contour_raster),
                             GPU_texture_height(textures_.tex_contour_raster));
      sub.dispatch(int3(screen_res.x, screen_res.y, 1));
      sub.barrier(eGPUBarrier::GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_SHADER_IMAGE_ACCESS);
    }
  }


  void StrokeGenPassModule::rebuild_pass_scan_test()
  {
    pass_scan_test.init();
    { // upsweep for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_upsweep");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_UPSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_in_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_out_scan_data_);
      sub.bind_ssbo(2, buffers_.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // reduction for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_AGGREGATE));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // down sweep for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_DWSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_out_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::rebuild_pass_segscan_test()
  {
    pass_segscan_test.init();
    { // upsweep for tree-scan
      auto& sub = pass_segscan_test.sub("strokegen_segscan_test_upsweep");
      sub.shader_set(shaders_.static_shader_get(SEGSCAN_TEST_UPSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_in_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_out_scan_data_);
      sub.bind_ssbo(2, buffers_.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // reduction for tree-scan
      auto& sub = pass_segscan_test.sub("strokegen_segscan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(SEGSCAN_TEST_AGGREGATE));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // down sweep for tree-scan
      auto& sub = pass_segscan_test.sub("strokegen_segscan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(SEGSCAN_TEST_DWSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_out_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  void StrokeGenPassModule::rebuild_pass_conv_test()
  {
    pass_conv1d_test.init();
    {
      auto &sub = pass_conv1d_test.sub("strokegen_segloopconv1D_test_build_patch");
      sub.shader_set(shaders_.static_shader_get(CONV1D_TEST_BUILD_PATCH));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_segloopconv1d_patch_table_);
      sub.bind_ssbo(1, buffers_.ssbo_debug_segloopconv1d_data_);
      sub.bind_ubo(0, buffers_.ubo_segloopconv1d_);

      sub.dispatch(int3(buffers_.ubo_segloopconv1d_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    {
      auto &sub = pass_conv1d_test.sub("strokegen_segloopconv1D_test_convolution");
      sub.shader_set(shaders_.static_shader_get(CONV1D_TEST_CONVOLUTION));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_segloopconv1d_patch_table_);
      sub.bind_ssbo(1, buffers_.ssbo_in_segloopconv1d_data_);
      sub.bind_ssbo(2, buffers_.ssbo_out_segloopconv1d_data_);
      sub.bind_ssbo(3, buffers_.ssbo_debug_segloopconv1d_data_);
      sub.bind_ubo(0, buffers_.ubo_segloopconv1d_);

      sub.dispatch(int3(buffers_.ubo_segloopconv1d_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
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
      sub.bind_ssbo(7, buffers_.ssbo_list_ranking_addressing_counters_);
      sub.bind_ssbo(8, buffers_.ssbo_list_ranking_serialized_topo_);
      sub.bind_ssbo(9, buffers_.ssbo_list_ranking_inputs_);
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
  void StrokeGenPassModule::append_subpass_list_ranking(ListRankingPassType passType, 
                                                      PassSimple &pass_listranking,
                                                      bool looped_pass_list_ranking
  )
  {
    // Build render passes
    bool custom_pass = (passType != ListRankingPassType::Test); 

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
          else if (passType == ListRankingPassType::ContourEdgeLinking)
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
    if (passType != ListRankingPassType::Test)
    { 
      auto& sub = pass_listranking.sub("strokegen_list_ranking_test_output");
      sub.shader_set(shaders_.static_shader_get(LISTRANKING_OUTPUT_DATA));

      sub.bind_ssbo(0, buffers_.ssbo_list_ranking_ranks_); 
      sub.bind_ssbo(1, buffers_.ssbo_list_ranking_serialized_topo_);
      sub.bind_ssbo(2, buffers_.ssbo_list_ranking_inputs_);
      if (passType == ListRankingPassType::ContourEdgeLinking)
      {
        sub.bind_ssbo(3, buffers_.ssbo_contour_edge_rank_);
        sub.bind_ssbo(4, buffers_.ssbo_contour_edge_list_len_);
        sub.bind_ssbo(5, buffers_.ssbo_contour_edge_list_head_);
      }

      sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[0]);
      sub.barrier(GPU_BARRIER_SHADER_STORAGE | GPU_BARRIER_COMMAND);
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
        if (test_looped_pass_list_ranking)
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

        if (rank_gt != (list_len_gt - 1 - rank_out))
        {
          fprintf(stderr, "strokegen error: incorrect rank, list ranking test failed: ");
          return false;
        }

        if (list_len_gt != list_len_out)
        {
          fprintf(stderr, "strokegen error: incorrect list length, list ranking test failed: ");
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
