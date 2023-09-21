/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 */

#include "npr_strokegen_pass.hh"

#include "gpu_batch_private.hh"

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
      case GEOM_EXTRACTION:
        return pass_extract_geom;
      case LIST_RANKING_TEST:
        return pass_listranking_test;
    }
    return pass_comp_test;
  }



  void StrokeGenPassModule::on_begin_sync()
  {
    rebuild_pass_scan_test();
    rebuild_pass_segscan_test();
    rebuild_pass_conv_test();
    rebuild_pass_list_ranking();

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

  // Note: we must do this within 15 ms(0.64ms x 23iters from willey's algo.), otherwise it for nothing.
  void StrokeGenPassModule::rebuild_pass_list_ranking()
  {
    // Build render passes
    pass_listranking_test.init();

    int num_splice_iters = 3;
    for (int splice_iter = 0; splice_iter < num_splice_iters; ++splice_iter)
    {
      rebuild_pass_list_ranking_fill_args(true, false, splice_iter, (int)GROUP_SIZE_BNPR_LIST_RANK_TEST);

      if (splice_iter == 0)
      {
        { // copy from staging buffer,
          auto &sub = pass_listranking_test.sub("strokegen_list_ranking_test_upload_cpu_data");
          sub.shader_set(shaders_.static_shader_get(LISTRANKING_UPLOAD_CPU_DATA));
          sub.bind_ssbo(0, buffers_.ssbo_list_ranking_links_staging_buf_);
          sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
          sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);

          sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[0]);
          sub.barrier(GPU_BARRIER_SHADER_STORAGE);
        }
      }

      int tag_iter;
      int num_tag_iters = splice_iter == 0 ? 3 : splice_iter == 1 ? 2 : 1;
      for (tag_iter = 0; tag_iter < num_tag_iters; ++tag_iter)
      {
        {
          auto &sub = pass_listranking_test.sub("strokegen_list_ranking_test_tagging");
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
          sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
          sub.push_constant("pc_listranking_splice_iter_", splice_iter);
          sub.push_constant("pc_listranking_tagging_iter_", tag_iter);

          sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splice_iter]);
          sub.barrier(GPU_BARRIER_SHADER_STORAGE);
        }
      }

      {
        auto& sub = pass_listranking_test.sub("strokegen_list_ranking_test_compaction");
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
        sub.bind_ssbo(10, buffers_.ssbo_list_ranking_debug_);
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", splice_iter);
        sub.push_constant("pc_num_splice_iters_", num_splice_iters);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splice_iter]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }

      int splice_disptatch_args_slot = splice_iter + 1;
      rebuild_pass_list_ranking_fill_args(false, true, splice_disptatch_args_slot, (int)GROUP_SIZE_BNPR_LIST_RANK_TEST);
      {
        auto& sub = pass_listranking_test.sub("strokegen_list_ranking_test_splice_out_nodes");
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
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", splice_iter);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[splice_disptatch_args_slot]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }
    }
    rebuild_pass_list_ranking_fill_args(true, false, num_splice_iters, GROUP_SIZE_BNPR_LIST_RANK_TEST);
    const int num_jump_iters = MAX_NUM_JUMPS_BNPR_LIST_RANK_TEST;
    {  // Sublist Pointer Jumping
      for (size_t jump_iter = 0; jump_iter < num_jump_iters; jump_iter++)
      {
        auto& sub = pass_listranking_test.sub("strokegen_list_ranking_test_sublist_pointer_jumping");
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_SUBLIST_POINTER_JUMPING));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[jump_iter%2]);
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[(jump_iter+1)%2]);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_node_to_anchor_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_anchor_to_node_[(num_splice_iters)%2]); // note: double check this
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(6, buffers_.ssbo_list_ranking_anchor_counters_);
        sub.bind_ssbo(7, buffers_.ssbo_list_ranking_addressing_counters_);
        sub.bind_ssbo(8, buffers_.ssbo_list_ranking_serialized_topo_);
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_splice_iter_", num_splice_iters);
        sub.push_constant("pc_listranking_jumping_iter_", (int)jump_iter);
        sub.push_constant("pc_listranking_ranking_pass_with_broken_loops_", 0);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[num_splice_iters]); // note: double check this
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
      }
    }

    { // Relinking spliced nodes back to the list
      int num_relink_iters = num_splice_iters;
      for (int relink_iter = 0; relink_iter < num_relink_iters; relink_iter++)
      {
        auto& sub = pass_listranking_test.sub("strokegen_list_ranking_test_relink");
        sub.shader_set(shaders_.static_shader_get(LISTRANKING_RELINKING));

        sub.bind_ssbo(0, buffers_.ssbo_list_ranking_spliced_node_id_[num_relink_iters - 1 - relink_iter]);  // note: double check this
        sub.bind_ssbo(1, buffers_.ssbo_list_ranking_links_);
        sub.bind_ssbo(2, buffers_.ssbo_list_ranking_ranks_);
        sub.bind_ssbo(3, buffers_.ssbo_list_ranking_splice_counters_);
        sub.bind_ssbo(4, buffers_.ssbo_list_ranking_per_anchor_sublist_jumping_info_[(num_jump_iters)%2]);
        sub.bind_ssbo(5, buffers_.ssbo_list_ranking_serialized_topo_);
        sub.bind_ubo(0, buffers_.ubo_list_ranking_splicing_);
        sub.push_constant("pc_listranking_relink_iter_", relink_iter);
        sub.push_constant("pc_listranking_num_relink_iters_", num_relink_iters);

        sub.dispatch(buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[num_relink_iters - relink_iter]);
        sub.barrier(GPU_BARRIER_SHADER_STORAGE);
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
    bool per_anchor, bool per_spliced, int splicing_or_relinking_iter, int group_size_x
  ){
    auto& sub = pass_listranking_test.sub("strokegen_list_ranking_test_fill_dispatch_args");
    sub.shader_set(shaders_.static_shader_get(eShaderType::LISTRANKING_FILL_DISPATCH_ARGS));

    sub.bind_ssbo(0, buffers_.ssbo_list_ranking_anchor_counters_);
    sub.bind_ssbo(1, buffers_.ssbo_list_ranking_splice_counters_);
    sub.bind_ssbo(2, buffers_.ssbo_list_ranking_indirect_dispatch_args_per_anchor[splicing_or_relinking_iter]);
    sub.bind_ssbo(3, buffers_.ssbo_list_ranking_indirect_dispatch_args_per_spliced[splicing_or_relinking_iter]);
    int dispatch_granularity_tag = per_anchor ? 0 : (per_spliced ? 1 : 2);
    sub.push_constant("pc_listranking_indirect_arg_granularity_", dispatch_granularity_tag); /* 0 := anchors, 1:= spliced nodes */
    sub.push_constant("pc_listranking_counter_buffer_slot_id_", splicing_or_relinking_iter); /* which counter to use */
    sub.push_constant("pc_listranking_dispatch_group_size_", group_size_x); /* which counter to use */

    sub.dispatch(int3(1, 1, 1));
    sub.barrier(GPU_BARRIER_SHADER_STORAGE);
  }

  bool StrokeGenPassModule::validate_list_ranking() const
  {
    SSBO_ListRankingRanks& buf_final_ranks = buffers_.ssbo_list_ranking_ranks_;
    buf_final_ranks.read();
    uint* computed_ranks = reinterpret_cast<uint *>(buf_final_ranks.data());

    SSBO_ListRankingSerializedTopo& buf_final_topo = buffers_.ssbo_list_ranking_serialized_topo_;
    buf_final_topo.read();
    uint* computed_topo = reinterpret_cast<uint *>(buf_final_topo.data());

    SSBO_ListRankingLinks& buf_final_links = buffers_.ssbo_list_ranking_links_;
    buf_final_links.read();
    uint* computed_links = reinterpret_cast<uint *>(buf_final_links.data());

    const uint num_nodes = buffers_.ubo_list_ranking_splicing_.num_nodes;
    uint tail_out_0 = computed_links[2 * 0 + 1];
    for (size_t i = 0; i < num_nodes; ++i)
    {
      uint rank_gt      = buffers_.listranking_test_nodes_rank[i];
      uint list_len_gt  = buffers_.listranking_test_nodes_list_len[i];
      uint tail_gt      = buffers_.listranking_test_nodes_tail[i];
      uint rank_out     = computed_ranks[i];
      uint list_len_out = computed_topo[2*i+1];
      uint tail_out     = computed_links[2*i+1]; // in the end every next pointer links to the tail node

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

    return true;
  }
}
