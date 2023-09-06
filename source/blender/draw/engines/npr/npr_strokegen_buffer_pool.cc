/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_instance.hh"
#include "bnpr_defines.hh"

namespace blender::npr::strokegen
{
  void GPUBufferPoolModule::on_begin_sync()
  {
    UBO_BnprTreeScan& ubo_tree_scan = ubo_bnpr_tree_scan_infos_;
    {
      ubo_bnpr_tree_scan_infos_.num_scan_items = NUM_ITEMS_BNPR_SCAN_TEST;
      ubo_bnpr_tree_scan_infos_.num_valid_scan_threads = compute_num_threads(
        NUM_ITEMS_BNPR_SCAN_TEST,
        2u
      );
      ubo_bnpr_tree_scan_infos_.num_thread_groups = compute_num_groups(
        NUM_ITEMS_BNPR_SCAN_TEST,
        GROUP_SIZE_BNPR_SCAN_TEST_SWEEP,
        2u
      );
      ubo_bnpr_tree_scan_infos_.dummy = 0u;
    }
    ubo_bnpr_tree_scan_infos_.push_update();


    UBO_SegLoopConv1D& ubo_segloopconv1d = ubo_segloopconv1d_;
    {
      ubo_segloopconv1d.num_thread_groups = compute_num_groups(
        NUM_ITEMS_SEGLOOPCONV1D_TEST,
        GROUP_SIZE_SEGLOOPCONV1D_TEST
      );
      ubo_segloopconv1d.num_conv_items = NUM_ITEMS_SEGLOOPCONV1D_TEST;
    }
    ubo_segloopconv1d_.push_update();

    UBO_ListRanking& ubo_list_ranking_splicing__ = ubo_list_ranking_splicing_;
    {
      ubo_list_ranking_splicing__.num_thread_groups = compute_num_groups(
        NUM_ITEMS_BNPR_LIST_RANK_TEST,
        GROUP_SIZE_BNPR_LIST_RANK_TEST
      );
      ubo_list_ranking_splicing__.num_nodes = NUM_ITEMS_BNPR_LIST_RANK_TEST;
      ubo_list_ranking_splicing__.subbuff_size = 4 * ((ubo_list_ranking_splicing__.num_nodes + 3) / 4);
      ubo_list_ranking_splicing__.num_tagging_iters = 1u;
    }
    ubo_list_ranking_splicing_.push_update();


    {
      // ssbo_list_ranking_indirect_dispatch_args_per_anchor[0].clear_to_zero();
      // GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE | GPU_BARRIER_COMMAND);
      
      // DispatchCommand* data = ssbo_list_ranking_indirect_dispatch_args_per_anchor[0].data();
      // ssbo_list_ranking_indirect_dispatch_args_per_anchor[0].read();
      // uint dispatch_args[4] = {
      //   compute_num_groups(NUM_ITEMS_BNPR_LIST_RANK_TEST, GROUP_SIZE_BNPR_LIST_RANK_TEST),
      //   1u,
      //   1u,
      //   1u
      // };
      // memcpy(data, dispatch_args, 4 * sizeof(uint));
      // ssbo_list_ranking_indirect_dispatch_args_per_anchor[0].push_update();  // glBufferSubData
      
      // GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE | GPU_BARRIER_COMMAND);
    }

    if (false == listranking_test_data_uploaded)
    { // list ranking test: build nodes on CPU & upload to ssbo
      build_list_ranking_testing_data(); // build listranking_test_nodes_prev_next

      ssbo_list_ranking_links_staging_buf_.resize(NUM_ITEMS_BNPR_LIST_RANK_TEST * 2);
      ssbo_list_ranking_links_staging_buf_.clear_to_zero();
      GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE);

      uint* data = ssbo_list_ranking_links_staging_buf_.data();
      ssbo_list_ranking_links_staging_buf_.read();
      memcpy(data, listranking_test_nodes_prev_next.data(), ssbo_list_ranking_links_staging_buf_.size() * sizeof(uint));
      ssbo_list_ranking_links_staging_buf_.push_update(); // glBufferSubData

      GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE);

      listranking_test_data_uploaded = true;
    }
    if (!listranking_test_data_validated && listranking_test_data_uploaded)
    { // list ranking test: validate uploaded nodes int the ssbo
      listranking_test_data_validated = true;

      ssbo_list_ranking_links_staging_buf_.read();
      uint* readback = ssbo_list_ranking_links_staging_buf_.data();
      for (size_t i = 0; i < listranking_test_nodes_prev_next.size(); ++i)
      {
        if (listranking_test_nodes_prev_next[i] != readback[i])
        {
          fprintf(stderr, "error: corrupted gpu buffer for testing list ranking.");
          listranking_test_data_validated = false;
        }
      }
    }


  }

  void GPUBufferPoolModule::sync_object(Object* ob)
  {

  }

  void GPUBufferPoolModule::end_sync()
  {
    // arr_buf_test_.resize(4096); // maybe needed for a few special buffers
  }


  void GPUBufferPoolModule::build_list_ranking_testing_data()
  { // TODO: make multiple lists instead of one long list
    listranking_test_nodes_prev_next.reinitialize(2 * NUM_ITEMS_BNPR_LIST_RANK_TEST);
    std::vector<int> nodes;
    nodes.resize((NUM_ITEMS_BNPR_LIST_RANK_TEST));
    for (int i = 0; i < NUM_ITEMS_BNPR_LIST_RANK_TEST; ++i)
      nodes[i] = i;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(nodes.begin(), nodes.end(), rng);

    Map<int, int> nodeToArr;
    Map<int, int> arrToNode;
    for (int i = 0; i < NUM_ITEMS_BNPR_LIST_RANK_TEST; ++i) {
      nodeToArr.add(i, nodes[i]);
      arrToNode.add(nodes[i], i);
    }

    for (int node = 0; node < NUM_ITEMS_BNPR_LIST_RANK_TEST; ++node) {
      int arr_id = nodeToArr.lookup(node);
      int arr_prev = arr_id > 0 ? (arr_id - 1) : arr_id; // head->prev == itself
      int arr_next = arr_id < NUM_ITEMS_BNPR_LIST_RANK_TEST - 1 ? arr_id + 1 : arr_id; // tail->next= itself
      listranking_test_nodes_prev_next[node*2] = arrToNode.lookup(arr_prev);
      listranking_test_nodes_prev_next[node*2+1] = arrToNode.lookup(arr_next);
    }
  }


}
