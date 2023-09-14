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
    int num_nodes = NUM_ITEMS_BNPR_LIST_RANK_TEST;

    listranking_test_nodes_prev_next.reinitialize(2 * num_nodes);

    listranking_test_nodes_rank.resize(num_nodes);
    listranking_test_nodes_head.resize(num_nodes);
    listranking_test_nodes_tail.resize(num_nodes);
    listranking_test_nodes_list_len.resize(num_nodes);

    std::vector<int> node_array;
    node_array.resize(num_nodes);

    std::vector<int> head, tail, seg_len, rank;
    head.resize(num_nodes);
    tail.resize(num_nodes);
    seg_len.resize(num_nodes);
    rank.resize(num_nodes);

    std::random_device rd;
    std::mt19937 rng(/*rd()*/(time(0)));
    int curr_head = 0;
    int curr_len =  std::max(1, (int)rng() % std::min(num_nodes, 2013));
    int curr_tail = curr_head + curr_len - 1;
    for (int arr_id = 0; arr_id < num_nodes; ++arr_id)
    {
        rank[arr_id] = arr_id - curr_head;
        seg_len[arr_id] = curr_len;
        head[arr_id] = curr_head;
        tail[arr_id] = curr_tail;
        if (curr_tail == arr_id)
        { // update to next sublist
          curr_head = arr_id + 1;
          curr_len = std::max(1, (int)rng() % std::min(std::max(1, num_nodes - curr_head), 2013));
          curr_tail = curr_head + curr_len - 1;
        }

        node_array[arr_id] = arr_id; // shuffle later
    }
    std::shuffle(node_array.begin(), node_array.end(), rng);

    Map<int, int> nodeToArr;
    Map<int, int> arrToNode;
    for (int node_id = 0; node_id < num_nodes; ++node_id) {
      uint arr_id = node_array[node_id];
      nodeToArr.add(node_id, arr_id);
      arrToNode.add(arr_id, node_id);
    }

    for (int node_id = 0; node_id < num_nodes; ++node_id) {
      uint arr_id = node_array[node_id];
      listranking_test_nodes_rank[node_id] = rank[arr_id];
      listranking_test_nodes_head[node_id] = arrToNode.lookup(head[arr_id]);
      listranking_test_nodes_tail[node_id] = arrToNode.lookup(tail[arr_id]);
      listranking_test_nodes_list_len[node_id] = seg_len[arr_id];
    }
    for (int node = 0; node < num_nodes; ++node) {
      int arr_id = nodeToArr.lookup(node);
      int arr_prev = (head[arr_id] == arr_id) ? arr_id : (arr_id - 1); // head->prev == itself
      int arr_next = (tail[arr_id] == arr_id) ? arr_id : (arr_id + 1); // tail->next= itself
      listranking_test_nodes_prev_next[node*2] = arrToNode.lookup(arr_prev);
      listranking_test_nodes_prev_next[node*2+1] = arrToNode.lookup(arr_next);
    }
  }
}
