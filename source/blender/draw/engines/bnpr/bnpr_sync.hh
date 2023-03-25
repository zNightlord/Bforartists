/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Module for synchronization between draw engine and outside
 */

#pragma once

#include "BKE_duplilist.h"
#include "BLI_ghash.h"
#include "BLI_map.hh"
#include "DNA_object_types.h"
#include "DRW_render.h"
#include "GPU_material.h"

#include "bnpr_shader_shared.hh"
#include "bnpr_sync_handles.hh"

namespace blender::bnpr
{

  class Instance;

  class SyncModule
  {
  private:
    Instance &inst_;

  public:
    SyncModule(Instance &inst) : inst_(inst) {};
    ~SyncModule(){};

    BnprDrawData &sync_object(Object *ob);
    WorldHandle &sync_world(::World *world) {};
    SceneHandle &sync_scene(::Scene *scene) {};

    void sync_mesh(Object *ob,
                   BnprDrawData &ob_draw_data,
                   draw::ResourceHandle res_handle,
                   const draw::ObjectRef &ob_ref);
  };



}


