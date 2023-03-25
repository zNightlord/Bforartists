/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * An renderer instance that contains all data to render a full frame.
*/

#pragma once

#include "BKE_image.h"
#include "DEG_depsgraph_query.h"
#include "DNA_shader_fx_types.h"
#include "DRW_engine.h"
#include "DRW_render.h"
#include "ED_view3d.h"
#include "GPU_capabilities.h"
#include "IMB_imbuf_types.h"

#include "draw_manager.hh"
#include "draw_pass.hh"
#include "bnpr_shader.hh"
#include "bnpr_sync.hh"
#include "bnpr_strokegen_buffer_pool.hh"
#include "bnpr_strokegen_texture_pool.hh"
#include "bnpr_strokegen_pass.hh"

namespace blender::bnpr
{
  using namespace draw;

  class Instance
  {
  private:

  public:
    /** Shading Modules */
    ShaderModule shaders; // singleton class for handling GPUShader(s)
    SyncModule sync;
    GPUBufferPoolModule   strokegen_buffers;
    GPUTexturePoolModule  strokegen_textures;
    StrokeGenPassModule   strokegen_passes;

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
    Instance() :
    shaders(*ShaderModule::module_get()),
    sync(*this),
    strokegen_buffers(*this),
    strokegen_textures(*this),
    strokegen_passes(shaders, strokegen_buffers, strokegen_textures)
    {
    };

    void init(Depsgraph* depsgraph_, draw::Manager* manager_, const View3D* v3d_, const RegionView3D* rv3d_, const
              DRWView* drw_view_, Object* camera_object_);
    void update_eval_members();

    void begin_sync(Manager& manager);
    void end_sync(Manager& manager);

    void object_sync(Manager& manager, ObjectRef& object_ref);


    void draw_viewport(Manager& manager, View& view);




  };
}


