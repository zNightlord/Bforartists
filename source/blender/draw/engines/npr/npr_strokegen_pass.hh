/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "draw_pass.hh"
#include "draw_manager.hh"

#include "npr_strokegen_shader.hh"
#include "bnpr_shader_shared.hh"
#include "gpu_batch_private.hh"
#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_texture_pool.hh"
#include "npr_strokegen_pass.hh"
#include "strokegen_mesh_raster_pass.hh"

#include <random>


namespace blender::npr::strokegen {
class Instance;

class StrokeGenPassModule // similar to "LineDrawingRenderPass"
{
private:
  /** Compute Passes */
  draw::PassSimple pass_comp_test = {"Strokegen Compute Test"};

  draw::PassSimple pass_extract_geom = {"StrokeGen Extract Geometry"};
  draw::PassSimple pass_process_contours = {"StrokeGen Process Contours"}; 
  draw::PassSimple pass_compress_contour_pixels = {"Generate Contour Pixel Mask"}; 

  draw::PassSimple pass_scan_test = {"Bnpr GPU Blelloch Scan Test"};
  draw::PassSimple pass_segscan_test = {"Bnpr GPU Blelloch SegScan Test"};
  draw::PassSimple pass_conv1d_test = {"Test GPU 1d conv on circular segments"};
  draw::PassSimple pass_listranking_test = {"List ranking test"};

  /** Draw Passes */
  StrokegenMeshRasterPass pass_draw_contour_edges = {"Draw Contour Edges"}; // Inherited from draw::PassMain 

  /** Dependent Modules */
  StrokeGenShaderModule &shaders_;
  GPUBufferPoolModule   &buffers_;
  GPUTexturePoolModule  &textures_;

  /** Geometry  */



public:
  StrokeGenPassModule(
      StrokeGenShaderModule &strokegen_shaders,
      GPUBufferPoolModule &strokegen_buffers,
      GPUTexturePoolModule &strokegen_textures
      )
    : shaders_(strokegen_shaders),
      buffers_(strokegen_buffers),
      textures_(strokegen_textures)
  {
  }

  ~StrokeGenPassModule() {}



  /* -------------------------------------------------------------------- */
  /** \name Passes
     * \{ */
  enum eType
  {
    SCAN_TEST = 0,
    SEGSCAN_TEST,
    SEGLOOPCONV_TEST,
    LIST_RANKING_TEST,

    GEOM_EXTRACTION,
    CONTOUR_PROCESS, 
    INDIRECT_DRAW_CONTOUR_EDGES,
    COMPRESS_CONTOUR_PIXELS,
  };

  PassSimple& get_compute_pass(eType passType);
  PassMain &get_contour_edge_draw_pass(eType passType);
  /** \} */



  /* -------------------------------------------------------------------- */
  /** \name Sync with Draw Module
     * \{ */
  void on_begin_sync();
  void on_end_sync(); 
  /** \} */



  /* -------------------------------------------------------------------- */
  /** \name Rebuild Render Passes
     * \{ */
  bool boostrap_before_extract_first_batch;
  int num_total_mesh_verts;
  int num_total_mesh_edges;
  int num_total_mesh_tris;
  void init_per_mesh_pass();
  void append_subpass_extract_contour_edges(GPUBatch *gpu_batch_line_adj,
                                            ResourceHandle &rsc_handle,
                                            gpu::Batch *edge_batch,
                                            int num_edges,
                                            gpu::GPUIndexBufType ib_type);
  void append_subpass_merge_vbo(GPUBatch *gpu_batch_surf, int batch_resource_index, int num_verts);
  void append_subpass_merge_line_adj_ibo(GPUBatch *gpu_batch_line_adj,
                                         int num_added_edges,
                                         gpu::GPUIndexBufType ib_type);
  void append_per_mesh_pass(Object* ob, GPUBatch* gpu_batch_line_adj, GPUBatch* gpu_batch_surf, ResourceHandle& rsc_handle, const DRWView* drw_view);
  void append_subpass_meshing_merge_verts(int num_verts_in, bool debug = false);
  void append_subpass_meshing_edge_adjacency(int num_edges_in, bool debug = false); 

  void rebuild_pass_process_contours();
  void append_subpass_fill_dispatch_args_contour_edges(bool all_contour_edges);
  void append_subpass_process_contour_edges();

  void rebuild_pass_contour_edge_drawcall();
  void rebuild_pass_compress_contour_pixels(bool debug = false); 



  bool test_scan; 
  void rebuild_pass_scan_test();
  void rebuild_pass_segscan_test();

  void rebuild_pass_conv_test();
  void rebuild_pass_list_ranking_pointer_jumping(
      PassSimple& pass_listranking, int num_splice_iters,
      const int num_jump_iters, int jumping_info_offset, bool loop_breaking_pass, bool loop_ranking_pass
      );

  bool test_list_ranking; 
  bool test_looped_pass_list_ranking;
  enum ListRankingPassType {
    Test = 0,
    ContourEdgeLinking = 1
  };
  void append_subpass_list_ranking(ListRankingPassType passType, PassSimple& pass_listranking, bool looped_pass_list_ranking);

  void rebuild_pass_list_ranking_fill_args(PassSimple& pass_listranking, bool per_anchor, bool per_spliced, int splicing_or_relinking_iter, int group_size_x, bool custom_pass);
  void print_list_ranking_nodes(int head_node_id, uint* computed_ranks, uint* computed_topo, uint* computed_links) const;
  /** \} */




  /* -------------------------------------------------------------------- */
  /** \name Scan Tester
   * \{ */
  template<typename T>
  void validate_pass_scan_test(bool (*equals)(const T &, const T &));
  template<typename T>
  bool validate_inter_block_exclusive_scan(
      const T *bufferInputVals,
      const T *bufferPrefixSum,
      bool (*equals)(const T &, const T &),
      uint num_scan_items,
      uint blk_size,
      uint numBlocks
      );
  template<typename T>
  bool validate_exclusive_scan(
      const T *bufferInputVals,
      const T *bufferPrefixSum,
      bool (*equals)(const T &, const T &),
      uint num_scan_items
      );
  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name Segment Scan Validation
   * \{ */
  template<typename T>
  void validate_segscan(
    bool (*func_equals)(const T&, const T&),
    uint (*func_get_hf)(const T&),
    T (*func_scan_op)(const T&, const T&),
    T zero_val,
    bool inclusive = false
  );
  template<typename T>
  bool validate_segscan_internal(
    const T *input,
    const T *output,
    const int blk_size, const int num_elems,
    uint (*func_get_hf)(const T&),
    bool (*func_equals)(const T &, const T &),
    T (*func_scan_op)(const T&, const T&),
    T zero_value,
    bool inclusive = false
  );
  /** \} */


  /* -------------------------------------------------------------------- */
  /** \name 1D Segmented Looped Convolution Validation
   * \{ */
  template<typename T>
  void validate_segloopconv1d(
    bool (*func_equals)(const T&, const T&),
    T (*func_conv_op)(const T&, const T&)
  );

  template<typename T>
  bool validate_segloopconv1d_internal(
    const T *input,
    const T *output,
    const int *seg_topo,
    int num_elems,
    int conv_radius,
    bool (*func_equals)(const T &, const T &),
    T (*func_conv_op)(const T&, const T&)
  );
  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name List Ranking Validation
   * \{ */
  bool validate_list_ranking() const;

  /** \} */


private:

};


template<typename T>
void StrokeGenPassModule::validate_pass_scan_test(bool (*equals)(const T &, const T &))
{
  SSBO_BnprScanData &buf_scan_inputs = buffers_.ssbo_in_scan_data_;
  buf_scan_inputs.read();
  T *data_scan_inputs = reinterpret_cast<T *>(buf_scan_inputs.data());

  SSBO_BnprScanData &buf_scan_output = buffers_.ssbo_out_scan_data_;
  buf_scan_output.read();
  T *data_scan_output = reinterpret_cast<T *>(buf_scan_output.data());

  bool valid_inter_block_scan = StrokeGenPassModule::validate_inter_block_exclusive_scan<T>(
      data_scan_inputs,
      data_scan_output,
      equals,

      buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items,
      GROUP_SIZE_BNPR_SCAN_TEST_SWEEP * 2u,
      buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups
      );
  if (!valid_inter_block_scan)
    fprintf(stderr, "bnpr: error: INTER-BLOCK scan test failed");

  bool valid_global_scan = StrokeGenPassModule::validate_exclusive_scan<T>(
      data_scan_inputs,
      data_scan_output,
      equals,
      buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items
      );
  if (!valid_global_scan)
    fprintf(stderr, "bnpr: error: GLOBAL scan test failed");
}


template<typename T>
bool StrokeGenPassModule::validate_inter_block_exclusive_scan(
    const T *const bufferInputVals,
    const T *const bufferPrefixSum,
    bool (*equals)(const T &, const T &),
    uint num_scan_items,
    uint blk_size,
    uint numBlocks
)
{
  Vector<int> failedElems(0);

  for (uint blk_id = 0u; blk_id < numBlocks; blk_id++) {
    for (uint blk_offset = 0u; blk_offset < blk_size - 1u; blk_offset++) {
      // index might go out of bound
      uint index = blk_id * blk_size + blk_offset;
      if (index >= num_scan_items - 2u)
        break;

      if (false == equals(
              bufferPrefixSum[index] + bufferInputVals[index],
              bufferPrefixSum[index + 1]
              )) {
        failedElems.append(index);
      }
    }
  }

  if (!failedElems.is_empty()) {
    return false;
  }

  return true;
}


template<typename T>
bool StrokeGenPassModule::validate_exclusive_scan(const T *bufferInputVals,
                                                  const T *bufferPrefixSum,
                                                  bool (* equals)(const T &, const T &),
                                                  uint num_scan_items)
{
  return validate_inter_block_exclusive_scan(
      bufferInputVals,
      bufferPrefixSum,
      equals,
      num_scan_items,
      num_scan_items,
      1
      );
}

template<typename T>
void StrokeGenPassModule::validate_segscan(
  bool (*func_equals)(const T&, const T&),
  uint (*func_get_hf)(const T&),
  T (*func_scan_op)(const T&, const T&),
  T zero_val,
  bool inclusive)
{
  SSBO_BnprScanData &buf_segscan_inputs = buffers_.ssbo_in_scan_data_;
  buf_segscan_inputs.read();
  T *data_segscan_inputs = reinterpret_cast<T *>(buf_segscan_inputs.data());

  SSBO_BnprScanData &buf_segscan_output = buffers_.ssbo_out_scan_data_;
  buf_segscan_output.read();
  T *data_segscan_output = reinterpret_cast<T *>(buf_segscan_output.data());

  bool succ = validate_segscan_internal(
    data_segscan_inputs,
    data_segscan_output,
    GROUP_SIZE_BNPR_SCAN_TEST_SWEEP * 2u,
    buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items,
    func_get_hf, func_equals, func_scan_op, zero_val,
    inclusive
  );

  if (!succ)
    fprintf(stderr, "strokegen error: segment scan test failed");
}

template<typename T>
bool StrokeGenPassModule::validate_segscan_internal(
  const T *input,
  const T *output,
  const int blk_size,
  const int num_elems,
  uint (*func_get_hf)(const T &),
  bool (*func_equals)(const T &, const T &),
  T (*func_scan_op)(const T&, const T&),
  T zero_value,
  const bool inclusive)
{
  T maxErrorFound = zero_value;
  T currScanValue = zero_value;

  for (int i = 0; i < num_elems; i++) {
    if (func_get_hf(input[i])) {
      currScanValue = zero_value;
    }

    if (inclusive) {
      // Inclusive scan
      currScanValue = func_scan_op(currScanValue, input[i]);
    }

    if (false == func_equals(currScanValue, output[i])) {
      return false;
    }

    if (!inclusive) {
      // Exclusive scan
      currScanValue = func_scan_op(currScanValue, input[i]);
    }
  }

  return true;
}

template <typename T>
void StrokeGenPassModule::validate_segloopconv1d(
  bool(* func_equals)(const T&, const T&),
  T(* func_conv_op)(const T&, const T&))
{
  SSBO_SegLoopConv1DData &buf_conv_inputs = buffers_.ssbo_in_segloopconv1d_data_;
  buf_conv_inputs.read();
  T *data_conv_inputs = reinterpret_cast<T *>(buf_conv_inputs.data());

  SSBO_SegLoopConv1DData &buf_conv_output = buffers_.ssbo_out_segloopconv1d_data_;
  buf_conv_output.read();
  T *data_conv_output = reinterpret_cast<T *>(buf_conv_output.data());

  SSBO_SegLoopConvDebugData &buf_dbg_data = buffers_.ssbo_debug_segloopconv1d_data_;
  buf_dbg_data.read();
  int* data_conv_debug = reinterpret_cast<int *>(buf_dbg_data.data());

  bool succ = validate_segloopconv1d_internal<T>(
    data_conv_inputs, data_conv_output, data_conv_debug,
    buffers_.ubo_segloopconv1d_.num_conv_items,
    NPR_SEGLOOPCONV1D_CONV_RADIUS,
    func_equals, func_conv_op
  );

  if (!succ)
    fprintf(stderr, "strokegen error: segloopconv1d test failed");
}

template <typename T>
bool StrokeGenPassModule::validate_segloopconv1d_internal(
  const T* input,
  const T* output,
  const int* seg_topo,
  int num_elems,
  int conv_radius,
  bool(* func_equals)(const T&, const T&),
  T(* func_conv_op)(const T&, const T&)
)
{
  T currOutput;

  for (int i = 0; i < num_elems; ++i)
  {
    const int seg_head = seg_topo[i * 4];
    const int seg_tail = seg_topo[i * 4 + 1];

    currOutput = input[i];
    { // convolution
      int neigh_id = i;
      for (int d = 1; d <= conv_radius; ++d)
      {
        neigh_id--;
        if (neigh_id < seg_head) neigh_id = seg_tail;

        currOutput = func_conv_op(currOutput, input[neigh_id]);
      }

      neigh_id = i;
      for (int d = 1; d <= conv_radius; ++d)
      {
        neigh_id++;
        if (neigh_id > seg_tail) neigh_id = seg_head;

        currOutput = func_conv_op(currOutput, input[neigh_id]);
      }
    }

    if (!func_equals(currOutput, output[i]))
      return false;
  }

  return true;
}
}
