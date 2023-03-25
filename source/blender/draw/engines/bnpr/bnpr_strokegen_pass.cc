/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 */

#include "bnpr_strokegen_pass.hh"

#include "gpu_batch_private.hh"

namespace blender::bnpr
{
  using namespace blender;

  PassSimple& StrokeGenPassModule::get_compute_pass(eType passType)
  {
    switch (passType) {
    case SCAN_TEST:
      return pass_scan_test;
    case SEGSCAN_TEST:
      return pass_segscan_test;
    case GEOM_EXTRACTION:
      return pass_extract_geom;
    }
    return pass_comp_test;
  }



  void StrokeGenPassModule::on_begin_sync()
  {
    rebuild_pass_scan_test();
    rebuild_pass_segscan_test();
    
    reset_pass_extract_mesh_geom();
  }



  /**
   * \brief In each render loop, re-initialize the compute pass for geometry extraction.
   */
  void StrokeGenPassModule::reset_pass_extract_mesh_geom()
  {
    pass_extract_geom.init();
  }
  
  /**
   * \brief Add a subpass for extracting geometry from given GPUBatch.
   * \param ob Mesh Object
   * \param gpu_batch_line_adj Mesh geometry stored in GPUBatch, ib stored with line adjacency info.
   */
  void StrokeGenPassModule::rebuild_sub_pass_extract_mesh_geom(Object* ob, GPUBatch* gpu_batch_line_adj)
  {
    const std::string pass_name = "extract_geom_";

    gpu::Batch* batch = static_cast<gpu::Batch*>(gpu_batch_line_adj); /* see example in "GPU_batch_draw_parameter_get" */ 
    if (batch == nullptr || batch->elem == nullptr)
      fprintf(stderr, "StrokeGen Error: empty mesh when rebuilding pass 'extract_mesh_geom'");
    int num_edges = batch->elem_()->index_len_get() / 4; // 4 indices per primitive, basically 4 verts around the edge

    gpu::IndexBuf* ib = batch->elem_();
    /* Hack to get ibo format */
    gpu::GPUIndexBufType ib_type = // the actual "index_type_" is protected
      (ib->size_get() / ib->index_len_get() == sizeof(uint32_t)) 
        ? gpu::GPU_INDEX_U32 : gpu::GPU_INDEX_U16;
    
    
    /* Cache per-edge geometry data */
    {
      auto& sub = pass_extract_geom.sub(pass_name.c_str());
      sub.shader_set(shaders_.static_shader_get(eShaderType::COMPUTE_GEOM_EXTRACT));
      sub.bind_ssbo(0, &(gpu_batch_line_adj->elem)); // TODO: investigate whether double pointer is necessary
      sub.bind_ssbo(1, &(gpu_batch_line_adj->verts[0])); // TODO: investigate what geom->vert[1~15] does
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_mesh_pool_);
      sub.push_constant("pcs_ib_fmt_u16", ib_type == gpu::GPU_INDEX_U16 ? 1 : 0);
      sub.push_constant("pcs_num_verts", (int)batch->elem_()->index_len_get());
      sub.push_constant(
        "pcs_num_ib_offset",
        (int)batch->elem_()->index_start_get() + (int)batch->elem_()->index_base_get()
      );

      sub.dispatch(int3(
        compute_num_groups(num_edges, GROUP_SIZE_STROKEGEN_GEOM_EXTRACT, 1),
        1, 1)
      );
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }

  

  void StrokeGenPassModule::rebuild_pass_scan_test()
  {
    pass_scan_test.init();
    { // upsweep for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_upsweep");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_UPSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_in_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_out_scan_data_);
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // reduction for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_AGGREGATE));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // down sweep for tree-scan
      auto& sub = pass_scan_test.sub("strokegen_scan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(eShaderType::SCAN_TEST_DWSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_out_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_scan_block_sum_);
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
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_in_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_out_scan_data_);
      sub.bind_ssbo(2, buffers_.ssbo_bnpr_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // reduction for tree-scan
      auto& sub = pass_segscan_test.sub("strokegen_segscan_test_aggregate");
      sub.shader_set(shaders_.static_shader_get(SEGSCAN_TEST_AGGREGATE));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_scan_block_sum_);

      sub.dispatch(int3(1, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
    { // down sweep for tree-scan
      auto& sub = pass_segscan_test.sub("strokegen_segscan_test_dwsweep");
      sub.shader_set(shaders_.static_shader_get(SEGSCAN_TEST_DWSWEEP));

      // Note: keep the same slot binding as in shader_create_info
      sub.bind_ssbo(0, buffers_.ssbo_bnpr_out_scan_data_);
      sub.bind_ssbo(1, buffers_.ssbo_bnpr_scan_block_sum_);
      sub.bind_ubo(0, buffers_.ubo_bnpr_tree_scan_infos_);

      sub.dispatch(int3(buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups, 1, 1));
      sub.barrier(GPU_BARRIER_SHADER_STORAGE);
    }
  }
}
