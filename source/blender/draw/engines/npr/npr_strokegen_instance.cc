/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * An instance contains all structures needed to do a complete render.
 */

#include <sstream>

#include "BKE_global.h"
#include "BKE_object.h"
#include "BLI_rect.h"
#include "DEG_depsgraph_query.h"
#include "DNA_ID.h"

#include "npr_strokegen_instance.hh"
#include "bnpr_defines.hh"

namespace blender::gpu {
class FrameBuffer;
}

namespace blender::npr::strokegen
{
  /* -------------------------------------------------------------------- */
  /** \name Initialization
   *
   * Initialization functions need to be called once at the start of a frame.
   * Active camera, render extent and enabled render passes are immutable until next init.
   * This takes care of resizing output buffers and view in case a parameter changed.
   * IMPORTANT: xxx.init() functions are NOT meant to acquire and allocate DRW resources.
   * Any attempt to do so will likely produce use after free situations.
   * \{ */
  void Instance::init(
    Depsgraph* depsgraph_,
    Manager* manager_,
    const View3D* v3d_,
    const RegionView3D* rv3d_,
    const DRWView* drw_view_,
    Object* camera_object_)
  {
    /* Init things static per render frame. (Not render graph related) */
    depsgraph = depsgraph_;
    manager = manager_;

    drw_view = drw_view_;
    v3d = v3d_;
    rv3d = rv3d_;

    camera_orig_object = camera_object_;

    info = "";
    frame_counter = 0;
  }

  void Instance::update_eval_members()
  {
    scene = DEG_get_evaluated_scene(depsgraph);
    camera_eval_object = (camera_orig_object) ?
                         DEG_get_evaluated_object(depsgraph, camera_orig_object) :
                         nullptr;
  }

  /** \} */






  /* -------------------------------------------------------------------- */
  /** \name Sync
   *
   * Sync will gather data from the scene that can change over a time step (i.e: motion steps).
   * IMPORTANT: xxx.sync() functions area responsible for creating DRW resources (i.e: DRWView) as
   * well as querying temp texture pool. All DRWPasses should be ready by the end on_finished_draw_viewport().
   * \{ */
  void Instance::begin_sync(Manager& manager, Texture& tex_prepass_depth)
  {
    /* Init draw passes and manager related stuff. (Begin render graph) */
    strokegen_passes.test_list_ranking = false;
    strokegen_passes.test_scan = false; 

    /* First setup resources */
    strokegen_buffers.on_begin_sync(
      drw_view, strokegen_passes.test_list_ranking, (frame_counter % 1024u) == 0u/*true*/
    );
    strokegen_textures.on_begin_sync(tex_prepass_depth); 

    /* Then setup render passes */
    strokegen_passes.on_begin_sync();
  }

  void Instance::end_sync(Manager&)
  {
    /* Post processing after all object. (Finish recording the commands) */
    strokegen_passes.on_end_sync();
  }

  void Instance::mesh_sync(
    Manager& manager, ObjectRef& object_ref, ResourceHandle& rsc_handle, GPUBatch** gpu_batch_surf)
  { /* Add object draw calls to passes. (Populate render graph) */
    Object *ob = object_ref.object;

    const bool is_renderable_type = ELEM(ob->type, OB_CURVES, OB_MESH, OB_LAMP);
    const int ob_visibility = DRW_object_visibility_in_active_context(ob);
    const bool partsys_is_visible = (ob_visibility & OB_VISIBLE_PARTICLES) != 0 &&
                                    (ob->type == OB_MESH);
    const bool object_is_visible = DRW_object_is_renderable(ob) &&
                                   (ob_visibility & OB_VISIBLE_SELF) != 0;

    if (!is_renderable_type || (!partsys_is_visible && !object_is_visible)) {
      return;
    }

    /* fclem: TODO cleanup. */
    ObjectRef ob_ref = DRW_object_ref_get(ob);
    BnprDrawData &bnpr_draw_data = sync.sync_object(ob);

    if (object_is_visible) {
      switch (ob->type) {
      case OB_MESH:
        sync.sync_mesh(ob, ob_ref, bnpr_draw_data, rsc_handle, drw_view);
        break;
      default:
        break;
      }
    }

    bnpr_draw_data.reset_recalc_flag();
  }


  /** \} */






  /* -------------------------------------------------------------------- */
  /** \name Rendering
   * \{ */

  void Instance::draw_viewport(Manager& manager, View& view, GPUTexture* pre_depth)
  {
    typedef StrokeGenPassModule::eType PType; 

    /* Submit passes here. (Execute render graph) */
    // GPU_storagebuf_copy_sub_from_vertbuf()
    GPU_texture_copy(strokegen_textures.tex_contour_raster_depth, pre_depth); 


    /* Geometry Extraction */
    manager.submit(strokegen_passes.get_compute_pass(PType::GEOM_EXTRACTION), view);
    manager.submit(strokegen_passes.get_compute_pass(PType::CONTOUR_PROCESS), view); 

    /* Draw Contour Edges */
    strokegen_textures.fb_contour_raster.bind();
    float fb_clear_col[4] = {0, 0, 0, 0}; 
    GPU_framebuffer_clear_color(strokegen_textures.fb_contour_raster, fb_clear_col);
    GPU_line_width(2.0f); // always snap to integer, see the opengl spec on line rasterization
    manager.submit(
      strokegen_passes.get_contour_edge_draw_pass(StrokeGenPassModule::INDIRECT_DRAW_CONTOUR_EDGES), view
    ); 
    GPU_line_width(1.0f);
    
    
    /* Pixel Extraction */
    manager.submit(strokegen_passes.get_compute_pass(PType::COMPRESS_CONTOUR_PIXELS), view); 










    /* GPU (Seg-)Scan Test ------------------------------------------------------------------------- */
    // manager.submit(strokegen_passes.get_compute_pass(ePassType::SCAN_TEST), view);
    manager.submit(strokegen_passes.get_compute_pass(PType::SEGSCAN_TEST), view);

    if (strokegen_passes.test_scan && frame_counter % 32 == 0)
    {
      // strokegen_passes.validate_pass_scan_test<BNPR_SCAN_TEST_DATA_TYPE>(
      //   [](const BNPR_SCAN_TEST_DATA_TYPE& a, const BNPR_SCAN_TEST_DATA_TYPE& b) {return a == b;}
      // );
      strokegen_passes.validate_segscan<SSBOData_SegScanTest>(
        [](const SSBOData_SegScanTest& a, const SSBOData_SegScanTest& b) { return a.val == b.val; },
        [](const SSBOData_SegScanTest& a) { return a.hf; },
        [](const SSBOData_SegScanTest& a, const SSBOData_SegScanTest& b)
        {
          return SSBOData_SegScanTest{a.val + b.val, 0};
        },
        SSBOData_SegScanTest{
          uint3(0u, 0u, 0u), 1u
        },
        false
      );

      strokegen_passes.validate_segloopconv1d<int>(
        [](const int& a, const int& b){ return a == b; },
        [](const int& a, const int& b){ return a + b; }
      );
    }
    frame_counter = (frame_counter + 1) % 100000000;



    /* GPU 1d looped segmented convolution Test ----------------------------------------- */
    manager.submit(strokegen_passes.get_compute_pass(PType::SEGLOOPCONV_TEST), view);



    /* GPU List Ranking Test ---------------------------------------------------------------- */
    if (strokegen_passes.test_list_ranking)
    {
      manager.submit(strokegen_passes.get_compute_pass(PType::LIST_RANKING_TEST), view);
      if (frame_counter % 64 == 0)
      {
        bool succ = strokegen_passes.validate_list_ranking();
      }
    }

  }

  void Instance::end_draw_viewport()
  {
    strokegen_textures.on_finished_draw_viewport(); 
  }

  /** \} */





}
