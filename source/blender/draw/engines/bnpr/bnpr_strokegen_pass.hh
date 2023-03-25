/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "draw_pass.hh"
#include "draw_manager.hh"

#include "bnpr_shader.hh"
#include "bnpr_shader_shared.hh"
#include "bnpr_strokegen_buffer_pool.hh"
#include "bnpr_strokegen_texture_pool.hh"
#include "bnpr_strokegen_pass.hh"

namespace blender::bnpr {
class Instance;

class StrokeGenPassModule // similar to "LineDrawingRenderPass"
{
private:
  /** Compute Passes */
  draw::PassSimple pass_comp_test = {"Strokegen Compute Test"};
  draw::PassSimple pass_scan_test = {"Bnpr GPU Blelloch Scan Test"};
  draw::PassSimple pass_segscan_test = {"Bnpr GPU Blelloch SegScan Test"};
  draw::PassSimple pass_extract_geom = {"StrokeGen Extract Geometry"};
  
  /** Instance */
  ShaderModule &shaders_;
  GPUBufferPoolModule &buffers_;
  GPUTexturePoolModule &textures_;

  /** Geometry  */
  

public:
  StrokeGenPassModule(
      ShaderModule &strokegen_shaders,
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
    GEOM_EXTRACTION, 
  };

  PassSimple& get_compute_pass(eType passType);
  /** \} */


  
  /* -------------------------------------------------------------------- */
  /** \name Sync with Draw Module
     * \{ */
  void on_begin_sync();
  /** \} */


  
  /* -------------------------------------------------------------------- */
  /** \name Rebuild Render Passes
     * \{ */
  void reset_pass_extract_mesh_geom();
  void rebuild_sub_pass_extract_mesh_geom(Object* ob, GPUBatch* gpu_batch_line_adj);

  void rebuild_pass_scan_test();
  void rebuild_pass_segscan_test();
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
  /** \name Segment Scan Tester
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
    int blk_size,
    uint (*func_get_hf)(const T&),
    bool (*func_equals)(const T &, const T &),
    T (*func_scan_op)(const T&, const T&),
    T zero_value,
    bool inclusive = false
  );
  /** \} */

private:

};


template<typename T>
void StrokeGenPassModule::validate_pass_scan_test(bool (*equals)(const T &, const T &))
{
  SSBO_BnprScanData &buf_scan_inputs = buffers_.ssbo_bnpr_in_scan_data_;
  buf_scan_inputs.read();
  T *data_scan_inputs = reinterpret_cast<T *>(buf_scan_inputs.data());

  SSBO_BnprScanData &buf_scan_output = buffers_.ssbo_bnpr_out_scan_data_;
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
  SSBO_BnprScanData &buf_segscan_inputs = buffers_.ssbo_bnpr_in_scan_data_;
  buf_segscan_inputs.read();
  T *data_segscan_inputs = reinterpret_cast<T *>(buf_segscan_inputs.data());

  SSBO_BnprScanData &buf_segscan_output = buffers_.ssbo_bnpr_out_scan_data_;
  buf_segscan_output.read();
  T *data_segscan_output = reinterpret_cast<T *>(buf_segscan_output.data());

  bool succ = validate_segscan_internal(
    data_segscan_inputs,
    data_segscan_output,
    GROUP_SIZE_BNPR_SCAN_TEST_SWEEP * 2u,
    func_get_hf, func_equals, func_scan_op, zero_val,
    inclusive
  );

  if (!succ)
    fprintf(stderr, "bnpr: error: segment scan test failed");
}

template<typename T>
bool StrokeGenPassModule::validate_segscan_internal(
  const T *input,
  const T *output,
  const int blk_size,
  uint (*func_get_hf)(const T &),
  bool (*func_equals)(const T &, const T &),
  T (*func_scan_op)(const T&, const T&),
  T zero_value,
  const bool inclusive)
{
  T maxErrorFound = zero_value;
  T currScanValue = zero_value;

  for (int i = 0; i < blk_size; i++) {
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


}
