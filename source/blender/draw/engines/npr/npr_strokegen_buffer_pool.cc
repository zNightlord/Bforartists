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

    UBO_ListRanking& ubo_list_ranking_tagging__ = ubo_list_ranking_tagging_;
    {
      ubo_list_ranking_tagging__.num_thread_groups = compute_num_groups(
        NUM_ITEMS_BNPR_LIST_RANK_TEST,
        GROUP_SIZE_BNPR_LIST_RANK_TEST
      );
      ubo_list_ranking_tagging__.num_nodes = NUM_ITEMS_BNPR_LIST_RANK_TEST;
      ubo_list_ranking_tagging__.subbuff_size = 4 * ((ubo_list_ranking_tagging__.num_nodes + 3) / 4);
      ubo_list_ranking_tagging__.num_tagging_iters = std::min(8u, ComputeTaggingIters(NUM_ITEMS_BNPR_LIST_RANK_TEST));
      // fprintf(stdout, "#tagging iters: %i ", num_iters); 
    }
    ubo_list_ranking_tagging_.push_update();

    UBO_BnprTreeScan &ubo_list_ranking_scan_infos = ubo_list_ranking_scan_infos_;
    {
      ubo_list_ranking_scan_infos.num_scan_items = NUM_ITEMS_BNPR_LIST_RANK_TEST;
      ubo_list_ranking_scan_infos.num_valid_scan_threads = compute_num_threads(
          ubo_list_ranking_scan_infos.num_scan_items, 2u);
      ubo_list_ranking_scan_infos.num_thread_groups = compute_num_groups(
          ubo_list_ranking_scan_infos.num_scan_items, GROUP_SIZE_BNPR_LIST_RANK_TEST, 2u);
      ubo_list_ranking_scan_infos.dummy = 0u;
    }
    ubo_list_ranking_scan_infos_.push_update();

    if (false == listranking_test_data_uploaded)
    { // list ranking test: build nodes on CPU & upload to ssbo
      build_list_ranking_testing_data(); // build listranking_test_nodes_prev_next

      ssbo_list_ranking_links_.resize(NUM_ITEMS_BNPR_LIST_RANK_TEST * 2);  
      ssbo_list_ranking_links_.clear_to_zero(); 
      GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE); 

      uint* data = ssbo_list_ranking_links_.data();
      ssbo_list_ranking_links_.read();
      memcpy(data, listranking_test_nodes_prev_next.data(), ssbo_list_ranking_links_.size() * sizeof(uint));
      ssbo_list_ranking_links_.push_update(); // glBufferSubData

      GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE); 

      listranking_test_data_uploaded = true; 
    }
    if (!listranking_test_data_validated && listranking_test_data_uploaded)
    { // list ranking test: validate uploaded nodes int the ssbo
      listranking_test_data_validated = true;

      ssbo_list_ranking_links_.read();
      uint* readback = ssbo_list_ranking_links_.data();
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
