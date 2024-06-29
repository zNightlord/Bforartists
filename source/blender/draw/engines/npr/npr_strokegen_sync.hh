/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Module for synchronization between draw engine and outside
 */

#pragma once

#include "BKE_duplilist.hh"
#include "BLI_ghash.h"
#include "BLI_map.hh"
#include "DNA_object_types.h"
#include "DRW_render.hh"
#include "GPU_material.hh"

#include "bnpr_shader_shared.hh"
#include "npr_sync_handles.hh"

namespace blender::npr::strokegen
{

  class Instance;

  class StrokegenSyncModule
  {
  private:
    Instance &inst_;


  public:
    StrokegenSyncModule(Instance &inst) : inst_(inst) {};
    ~StrokegenSyncModule(){};


    Map<ObjectKey, ObjectHandle> ob_handles = {};

    ObjectHandle sync_object(const ObjectRef& ob_ref);
    WorldHandle &sync_world(::World *world) {};
    SceneHandle &sync_scene(::Scene *scene) {};

    void sync_mesh(Object* ob,
                   const draw::ObjectRef& ob_ref,
                   draw::ResourceHandle& rsc_handle, const DRWView* drw_view);
  };



}


