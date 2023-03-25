/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#include "bnpr_strokegen_buffer_pool.hh"
#include "bnpr_instance.hh"
#include "bnpr_defines.hh"

namespace blender::bnpr
{
  void GPUBufferPoolModule::on_begin_sync()
  {
    UBO_BnprTreeScan& ubo_tree_scan = ubo_bnpr_tree_scan_infos_;
    {
      ubo_tree_scan.num_scan_items = NUM_ITEMS_BNPR_SCAN_TEST;
      ubo_tree_scan.num_valid_scan_threads = compute_num_threads(
        NUM_ITEMS_BNPR_SCAN_TEST,
        2u
      );
      ubo_tree_scan.num_thread_groups = compute_num_groups(
        NUM_ITEMS_BNPR_SCAN_TEST,
        GROUP_SIZE_BNPR_SCAN_TEST_SWEEP,
        2u
      );
      ubo_tree_scan.dummy = 0u;
    }
    ubo_bnpr_tree_scan_infos_.push_update();
  }

  void GPUBufferPoolModule::sync_object(Object* ob)
  {

  }

  void GPUBufferPoolModule::end_sync()
  {
    // arr_buf_test_.resize(4096); // maybe needed for a few special buffers
  }



}
