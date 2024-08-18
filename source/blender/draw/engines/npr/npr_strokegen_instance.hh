/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * An renderer instance that contains all data to render a full frame.
*/

#pragma once

#include "DRW_engine.hh"
#include "DRW_render.hh"
#include "ED_view3d.hh"

#include "npr_strokegen_shader.hh"
#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_pass.hh"
#include "npr_strokegen_texture_pool.hh"
#include "npr_strokegen_sync.hh"
#include "draw_pass.hh"
#include "gpu_framebuffer_private.hh"

namespace blender::npr::strokegen
{
  using namespace draw;

  class Instance
  {
  public:
    /** Gfx Resources */
    StrokeGenShaderModule shaders; // singleton class for handling GPUShader(s)
    StrokegenSyncModule sync;
    GPUBufferPoolModule   strokegen_buffers;
    GPUTexturePoolModule  strokegen_textures;
    StrokeGenPassModule   strokegen_passes;
    bool has_strokegen_enabled_mesh; 

    /** Input data. */
    Depsgraph *depsgraph;
    Manager *manager;

    /** Evaluated IDs. */
    Scene *scene;
    Object *camera_eval_object;
    Object *camera_orig_object;
    /** Only available when rendering for viewport. */
    const DRWView *drw_view;
    const View3D *v3d;
    const RegionView3D *rv3d;

    /** Info string displayed at the top of the render / viewport. */
    std::string info = "";
    /** Debug mode from debug value. */
    // eDebugMode debug_mode = eDebugMode::DEBUG_NONE;
    uint frame_counter; // for debugging



   public:
    Instance()
        : shaders(*StrokeGenShaderModule::module_get()),
          sync(*this),
          strokegen_buffers(*this),
          strokegen_textures(*this),
          strokegen_passes(shaders, strokegen_buffers, strokegen_textures),
          has_strokegen_enabled_mesh(false)
    {
      strokegen_passes.test_looped_pass_list_ranking = true; // remember to also set flag at build_list_ranking_testing_data
    }


    /** Render Loop Events */
    void init(Depsgraph* depsgraph_, draw::Manager* manager_, const View3D* v3d_, const RegionView3D* rv3d_, const
              DRWView* drw_view_, Object* camera_object_);
    void update_eval_members();

    void begin_sync(Manager& manager, Texture& tex_prepass_depth);
    void end_sync(Manager& manager);

    void mesh_sync(Manager& manager, ObjectRef& object_ref, ResourceHandle& rsc_handle, gpu::Batch** gpu_batch_surf);

    void draw_viewport(Manager& manager, View& view, GPUTexture* pre_depth);
    void end_draw_viewport();

    /* Sync Object */
    uint64_t depsgraph_last_update_ = 0; 
    int get_recalc_flags(const ObjectRef &ob_ref) const;

  };
}


