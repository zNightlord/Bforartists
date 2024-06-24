/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#include "npr_strokegen_buffer_pool.hh"

#include "BLI_math_quaternion.hh"
#include "BLI_math_rotation.h"

#include <set>

#include "npr_strokegen_instance.hh"
#include "bnpr_defines.hh"
#include "DEG_depsgraph_query.hh"

#include <iostream>

namespace blender::npr::strokegen
{

  void GPUBufferPoolModule::on_begin_sync(
    const DRWView* drw_view, bool upload_list_ranking_test_data,
    bool update_view_matrices_for_contour_extraction, 
    bool update_view_matrices_for_dbg_view_rotation)
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
        GROUP_SIZE_BNPR_SCAN_SWEEP,
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

    UBO_ViewMatrices& ubo_view_matrices__ = ubo_view_matrices_;
    {
      ubo_view_matrices__.viewmat = drw_view->storage.viewmat;
      ubo_view_matrices__.viewinv = drw_view->storage.viewinv;
      ubo_view_matrices__.winmat = drw_view->storage.winmat;
      ubo_view_matrices__.wininv = drw_view->storage.wininv;

      ubo_view_matrices__.push_update();
    }

    UBO_ViewMatrices & ubo_view_matrices_cache__ = ubo_view_matrices_cache_;
    if (update_view_matrices_for_contour_extraction)
    {
      ubo_view_matrices_cache__.viewmat = drw_view->storage.viewmat;
      ubo_view_matrices_cache__.viewinv = drw_view->storage.viewinv;
      ubo_view_matrices_cache__.winmat = drw_view->storage.winmat;
      ubo_view_matrices_cache__.wininv = drw_view->storage.wininv;
      ubo_view_matrices_cache__.push_update();

      ubo_view_matrices_cache_2_.viewmat = drw_view->storage.viewmat;
      ubo_view_matrices_cache_2_.viewinv = drw_view->storage.viewinv;
      ubo_view_matrices_cache_2_.winmat = drw_view->storage.winmat;
      ubo_view_matrices_cache_2_.wininv = drw_view->storage.wininv;
      ubo_view_matrices_cache_2_.push_update();
    }else if (update_view_matrices_for_dbg_view_rotation){
      // Fixed view, or with custom rotations
      // fetch ui inputs
      float rot_ang;
      get_dbg_view_rot_angle(rot_ang);
      float4x4 v = ubo_view_matrices_cache_2_.viewmat;

      float r[4][4];
      float rot_axis[3] = { .0f, .0f, 1.0f };
      axis_angle_to_mat4(r, rot_axis, rot_ang);
      v = v * float4x4(r);

      float4x4 inv_v; // = math::invert(v); // invert_m4(v.ptr);
      invert_m4_m4(inv_v.ptr(), v.ptr());

      ubo_view_matrices_cache_.viewmat = v; // drw_view->storage.viewmat;
      ubo_view_matrices_cache_.viewinv = inv_v; // drw_view->storage.viewinv;
      ubo_view_matrices_cache_.winmat = ubo_view_matrices_cache_2_.winmat;
      ubo_view_matrices_cache_.wininv = ubo_view_matrices_cache_2_.wininv;
      ubo_view_matrices_cache_.push_update();
    }


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

    if (upload_list_ranking_test_data)
    {
      if (false == listranking_test_data_uploaded) {  // list ranking test: build nodes on CPU &
                                                      // upload to ssbo
        build_list_ranking_testing_data(true);        // build listranking_test_nodes_prev_next

        ssbo_list_ranking_links_staging_buf_.resize(NUM_ITEMS_BNPR_LIST_RANK_TEST * 2);
        ssbo_list_ranking_links_staging_buf_.clear_to_zero();
        GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE);

        uint *data = ssbo_list_ranking_links_staging_buf_.data();
        ssbo_list_ranking_links_staging_buf_.read();
        memcpy(data,
               listranking_test_nodes_prev_next.data(),
               ssbo_list_ranking_links_staging_buf_.size() * sizeof(uint));
        ssbo_list_ranking_links_staging_buf_.push_update();  // glBufferSubData

        GPU_memory_barrier(GPU_BARRIER_BUFFER_UPDATE);

        listranking_test_data_uploaded = true;
      }
      if (!listranking_test_data_validated &&
          listranking_test_data_uploaded) {  // list ranking test: validate uploaded nodes int the
                                             // ssbo
        listranking_test_data_validated = true;

        ssbo_list_ranking_links_staging_buf_.read();
        uint *readback = ssbo_list_ranking_links_staging_buf_.data();
        for (size_t i = 0; i < listranking_test_nodes_prev_next.size(); ++i) {
          if (listranking_test_nodes_prev_next[i] != readback[i]) {
            fprintf(stderr, "error: corrupted gpu buffer for testing list ranking.");
            listranking_test_data_validated = false;
          }
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


  void GPUBufferPoolModule::build_list_ranking_testing_data(bool test_loop_lists)
  { // TODO: make multiple lists instead of one long list
    int num_nodes = NUM_ITEMS_BNPR_LIST_RANK_TEST;

    listranking_test_nodes_prev_next.reinitialize(2 * num_nodes);

    listranking_test_nodes_rank.resize(num_nodes);
    listranking_test_nodes_head.resize(num_nodes);
    listranking_test_nodes_tail.resize(num_nodes);
    listranking_test_nodes_list_len.resize(num_nodes);
    listranking_test_head_nodes.clear();

    std::vector<int> node_array;
    node_array.resize(num_nodes);

    std::vector<int> head, tail, seg_len, rank;
    head.resize(num_nodes);
    tail.resize(num_nodes);
    seg_len.resize(num_nodes);
    rank.resize(num_nodes);

    std::random_device rd;
    // std::mt19937 rng(/*rd()*/ (time(0)));
    std::mt19937 rng((time(0)));
    std::srand(time(0));
    int curr_arr_head = 0;
    int curr_len = std::max(2, (int)std::rand() % num_nodes / 2);
    int curr_arr_tail = curr_arr_head + curr_len - 1;
    for (int arr_id = 0; arr_id < num_nodes; ++arr_id) {
      rank[arr_id] = arr_id - curr_arr_head;
      seg_len[arr_id] = curr_len;
      head[arr_id] = curr_arr_head;
      tail[arr_id] = curr_arr_tail;
      if (curr_arr_tail == arr_id) {  // update to next sublist
        curr_arr_head = arr_id + 1;
        curr_len = std::min(num_nodes - curr_arr_head,
                            std::max((int)std::rand() % (num_nodes / 2), 2));
          curr_arr_tail = curr_arr_head + curr_len - 1;
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

    std::vector<int> arr_looped_rank, arr_looped_head, arr_looped_tail;
    arr_looped_rank.resize(num_nodes);
    arr_looped_head.resize(num_nodes);
    arr_looped_tail.resize(num_nodes);
    if (test_loop_lists)
    {
      int arr_id = 0;
      while (arr_id < num_nodes)
      { // traverse sub-lists
        curr_arr_head = head[arr_id];
        curr_arr_tail = tail[arr_id];

        int arr_loop_shift = 0;
        int loop_head_node_id = arrToNode.lookup(curr_arr_head);
        int loop_len = seg_len[curr_arr_head];

        // scan list to find loop head
        for (int arr_id_scan = curr_arr_head; arr_id_scan <= curr_arr_tail; ++arr_id_scan)
        { // node with maximum address is the loop head
          int node_id = arrToNode.lookup(arr_id_scan);
          if (node_id > loop_head_node_id)
          {
            loop_head_node_id = node_id;
            arr_loop_shift = arr_id_scan - curr_arr_head;
          }
        }
        listranking_test_head_nodes.insert(loop_head_node_id);
        int loop_head_arr_id = nodeToArr.lookup(loop_head_node_id);

        // record mapping from sub-lists to shifted ones
        for (int arr_id_scan = curr_arr_head; arr_id_scan <= curr_arr_tail; ++arr_id_scan)
        {
          int local_rank = arr_id_scan - curr_arr_head;
          local_rank = (local_rank - arr_loop_shift + loop_len) % loop_len; // shift loop head node to the start
          arr_looped_rank[arr_id_scan] = local_rank;
          arr_looped_head[arr_id_scan] = loop_head_arr_id;
          arr_looped_tail[arr_id_scan] = curr_arr_head + ((arr_loop_shift - 1 + loop_len) % loop_len);
        }

        arr_id = curr_arr_tail + 1;
      }
    }


    for (int node_id = 0; node_id < num_nodes; ++node_id)
    {
      uint arr_id = nodeToArr.lookup(node_id);
      listranking_test_nodes_rank[node_id] = (!test_loop_lists) ? rank[arr_id] : arr_looped_rank[arr_id];
      listranking_test_nodes_head[node_id] = arrToNode.lookup(
        (!test_loop_lists) ? head[arr_id] : arr_looped_head[arr_id]
      );
      listranking_test_nodes_tail[node_id] = arrToNode.lookup(
        (!test_loop_lists) ? tail[arr_id] : arr_looped_tail[arr_id]
      );
      listranking_test_nodes_list_len[node_id] = seg_len[arr_id];
    }
    for (int node = 0; node < num_nodes; ++node) {
      int arr_id = nodeToArr.lookup(node);
      int arr_prev, arr_next;
      if (!test_loop_lists){
        arr_prev = (head[arr_id] == arr_id) ? arr_id : (arr_id - 1); // head->prev == itself
        arr_next = (tail[arr_id] == arr_id) ? arr_id : (arr_id + 1); // tail->next= itself
      }
      else{
        int curr_list_len = listranking_test_nodes_list_len[node];
        int curr_list_beg = head[arr_id]; // start addr of serialized loop list
        int list_rank = arr_id - curr_list_beg;
        arr_prev = curr_list_beg + (list_rank - 1 + curr_list_len) % curr_list_len; // head->prev == itself
        arr_next = curr_list_beg + (list_rank + 1 + curr_list_len) % curr_list_len; // tail->next == itself
      }
      listranking_test_nodes_prev_next[node*2]   = arrToNode.lookup(arr_prev);
      listranking_test_nodes_prev_next[node*2+1] = arrToNode.lookup(arr_next);
    }

    // Print linked lists
    if (num_nodes <= 64)
    {
      for (int head_node_id : listranking_test_head_nodes)
      {
        int list_len = listranking_test_nodes_list_len[head_node_id];
        int curr_node_id = head_node_id;
        for (int i = 0; i < list_len; ++i)
        {
          std::cout << curr_node_id;
          if ((i + 1) != list_len) std::cout << "->";
          curr_node_id = listranking_test_nodes_prev_next[curr_node_id * 2 + 1];
        }
        std::cout << std::endl;
      }
    }
  }

  void GPUBufferPoolModule::get_dbg_view_rot_angle(float& rot_ang) {
    const DRWContextState* draw_ctx = DRW_context_state_get();
    const Scene* scene_eval = DEG_get_evaluated_scene(draw_ctx->depsgraph);
    rot_ang = .01f * scene_eval->npr.npr_test_val_23;
  }

void GPUBufferPoolModule::should_dbg_rotate_view_matrices_cache(float& rot_ang, bool& dbg_rotate_view_matrix) {
    get_dbg_view_rot_angle(rot_ang);
    dbg_rotate_view_matrix = 1e-10f < math::abs(rot_angle_ubo_view_matrices_cache_ - rot_ang);
    rot_angle_ubo_view_matrices_cache_ = rot_ang;
  }
}
