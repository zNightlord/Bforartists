/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * An instance contains all structures needed to do a complete render.
 */

#include <sstream>

#include "BKE_global.hh"
#include "BKE_object.hh"
#include "BLI_rect.h"
#include "DEG_depsgraph_query.hh"
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
    has_strokegen_enabled_mesh = false; // reset this switch each frame

    /* Init draw passes and manager related stuff. (Begin render graph) */
    strokegen_passes.test_list_ranking = false;
    strokegen_passes.test_scan = false;
    strokegen_passes.test_segloopconv = false; 

    /* First setup resources */
    int refresh_rate = std::max(1, strokegen_passes.meshing_params.seconds_sync_view_mat * 24000); 
    float rot_ang; bool dbg_rotate_view_matrix;
    strokegen_buffers.should_dbg_rotate_view_matrices_cache(rot_ang, dbg_rotate_view_matrix);
    strokegen_buffers.on_begin_sync(
      drw_view, strokegen_passes.test_list_ranking,
      (frame_counter % refresh_rate) == 0u,
      dbg_rotate_view_matrix/*true*/
    );
    strokegen_textures.on_begin_sync(tex_prepass_depth); 

    /* Then setup render passes */
    strokegen_passes.on_begin_sync((int)frame_counter);

    /* Setup gpu scene info */
    object_gpu_id_map.clear();
  }

  void Instance::end_sync(Manager&)
  {
    /* Post processing after all object. (Finish recording the commands) */
    strokegen_passes.on_end_sync();
    strokegen_buffers.on_end_sync(); 

    depsgraph_last_update_ = DEG_get_update_count(depsgraph);
  }

  void Instance::mesh_sync(
    Manager& manager, ObjectRef& object_ref, ResourceHandle& rsc_handle, gpu::Batch** gpu_batch_surf)
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

    const PerObjectStrokegenSettings &ui_input = ob->strokegen;
    if (ui_input.curve_type == 0)
      return; 

    /* fclem: TODO cleanup. */
    ObjectRef ob_ref = DRW_object_ref_get(ob);

    if (object_is_visible) {
      switch (ob->type) {
      case OB_MESH:
        sync.sync_mesh(ob, ob_ref, rsc_handle, drw_view);
        break;
      default:
        break;
      }
    }

  }


  /** \} */






  /* -------------------------------------------------------------------- */
  /** \name Rendering
   * \{ */

  void Instance::draw_viewport(Manager& manager, View& view, GPUTexture* pre_depth)
  {
    typedef StrokeGenPassModule::eType PType; 
    /* Submit passes here. (Execute render graph) */

    // We always need to clear the frame-buffers
    manager.submit(strokegen_passes.get_clear_framebuffer_pass());

    // Make sure buffer content is properly set
    if (false == has_strokegen_enabled_mesh)
      return; 

    GPU_texture_copy(strokegen_textures.tex_contour_raster_depth, pre_depth);
    
    for (int i = 0; i < strokegen_passes.get_num_passes_extract_geom(); ++i)
    {
      /* Geometry Extraction */
      manager.submit(strokegen_passes.get_compute_pass(PType::GEOM_EXTRACTION, i), view);

      /* Draw depth for remeshed surfaces */
      PassMain &render_pass_remesh_depth = strokegen_passes.get_render_pass(PType::INDIRECT_DRAW_REMESHED_DEPTH, i);
      render_pass_remesh_depth.framebuffer_set(&strokegen_textures.fb_remeshed_depth);
      strokegen_textures.fb_remeshed_depth.bind();
      manager.submit(render_pass_remesh_depth, view); 
    }


    /* Process Contour Edges */
    manager.submit(strokegen_passes.get_compute_pass(PType::CONTOUR_PROCESS), view); 



    /* Draw Contour Edges */
    strokegen_textures.fb_contour_raster.bind();
    GPU_line_width(2.0f); // always snap to integer, see the opengl spec on line rasterization
    manager.submit(strokegen_passes.get_render_pass(StrokeGenPassModule::INDIRECT_DRAW_DBG_VNOR), view);
    GPU_line_width(4.0f); 
    // manager.submit(strokegen_passes.get_render_pass(StrokeGenPassModule::INDIRECT_DRAW_CONTOUR_EDGES), view);
    // manager.submit(strokegen_passes.get_render_pass(StrokeGenPassModule::INDIRECT_DRAW_CONTOUR_2D_SAMPLES), view);
    GPU_line_width(1.0f);
    
    /* Pixel Extraction */
    manager.submit(strokegen_passes.get_compute_pass(PType::COMPRESS_CONTOUR_PIXELS), view); 







    /* GPU (Seg-)Scan Test ------------------------------------------------------------------------- */
    const bool dbg_segscan = false; 
    manager.submit(strokegen_passes.get_compute_pass(dbg_segscan ? PType::SEGSCAN_TEST : PType::SCAN_TEST), view);

    if (strokegen_passes.test_scan && frame_counter % 32 == 0)
    {
      // validate scan operators
      if (dbg_segscan) {
        strokegen_passes.validate_segscan<SSBOData_SegScanTest, SSBOData_SegScanTestEncoded>(
            SSBOData_SegScanTestDecodeFunc,
          [](const SSBOData_SegScanTest& a, const SSBOData_SegScanTest& b) { return a.val ==
          b.val; },
          [](const SSBOData_SegScanTest& a) { return a.hf; },
          [](const SSBOData_SegScanTest& a, const SSBOData_SegScanTest& b)
          {
            return SSBOData_SegScanTest{a.val + b.val, 0};
          },
          SSBOData_SegScanZeroValue,
          false
        );
      }else {
        strokegen_passes.validate_pass_scan_test<BNPR_SCAN_TEST_DATA_TYPE>(
            [](const BNPR_SCAN_TEST_DATA_TYPE &a, const BNPR_SCAN_TEST_DATA_TYPE &b) {
              return a == b;
            });  
      }
      
    }
    frame_counter = strokegen_frame_id_next(frame_counter);



    /* GPU 1d looped segmented convolution Test ----------------------------------------- */
    manager.submit(strokegen_passes.get_compute_pass(PType::SEGLOOPCONV_TEST), view);
    if (frame_counter % 32 == 0 && strokegen_passes.test_segloopconv) {
      // validate segloopconv1d
      strokegen_passes.validate_segloopconv1d<NPR_SEGLOOPCONV1D_TEST_DATA_TYPE>(
          [](const NPR_SEGLOOPCONV1D_TEST_DATA_TYPE &a,
             const NPR_SEGLOOPCONV1D_TEST_DATA_TYPE &b) { return a == b; },
          [](const NPR_SEGLOOPCONV1D_TEST_DATA_TYPE &a,
             const NPR_SEGLOOPCONV1D_TEST_DATA_TYPE &b) { return a + b; }
      );
    }


    /* GPU List Ranking Test ---------------------------------------------------------------- */
    if (strokegen_passes.test_list_ranking)
    {
      manager.submit(strokegen_passes.get_compute_pass(PType::LIST_RANKING_TEST), view);
      if (frame_counter % 8 == 0)
      {
        bool succ = strokegen_passes.validate_list_ranking();
      }
    }

  }

  void Instance::end_draw_viewport()
  {
    strokegen_textures.on_finished_draw_viewport(); 
  }

  int Instance::get_recalc_flags(const ObjectRef &ob_ref) const
  {
    auto get_flags = [&](const ObjectRuntimeHandle &runtime) {
      int flags = 0;
      SET_FLAG_FROM_TEST(
          flags, runtime.last_update_transform > depsgraph_last_update_, ID_RECALC_TRANSFORM);
      SET_FLAG_FROM_TEST(
          flags, runtime.last_update_geometry > depsgraph_last_update_, ID_RECALC_GEOMETRY);
      SET_FLAG_FROM_TEST(
          flags, runtime.last_update_shading > depsgraph_last_update_, ID_RECALC_SHADING);
      return flags;
    };

    int flags = get_flags(*ob_ref.object->runtime);
    if (ob_ref.dupli_parent) {
      flags |= get_flags(*ob_ref.dupli_parent->runtime);
    }

    return flags;
  }

  /** \} */


}
