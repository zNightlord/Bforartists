/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "bnpr_shader_shared.hh"

namespace blender::npr::strokegen {
class Instance;

class GPUBufferPoolModule {
  friend class StrokeGenPassModule;

 private:
  /** Instance */
  Instance &instance;

  /** Compute Resources */
  SSBO_StrokeGenTest ssbo_bnpr_test_;

  SSBO_StrokeGenMeshPool     ssbo_bnpr_mesh_pool_;
  SSBO_StrokeGenMeshPoolArgs ssbo_bnpr_mesh_pool_args_;

  SSBO_BnprScanData       ssbo_in_scan_data_;
  SSBO_BnprScanData       ssbo_out_scan_data_;
  SSBO_BnprScanAggregates ssbo_scan_block_sum_;
  UBO_BnprTreeScan        ubo_bnpr_tree_scan_infos_;

  SSBO_SegLoopConv1DData     ssbo_in_segloopconv1d_data_;
  SSBO_SegLoopConv1DData     ssbo_out_segloopconv1d_data_;
  SSBO_SegLoopConvPatchTable ssbo_segloopconv1d_patch_table_;
  UBO_SegLoopConv1D          ubo_segloopconv1d_;

 public:
  GPUBufferPoolModule(Instance &inst) : instance(inst)
  {
  }
  ~GPUBufferPoolModule()
  {
  }

  void on_begin_sync();
  void sync_object(Object *ob);
  void end_sync();
};
}  // namespace blender::npr::strokegen
