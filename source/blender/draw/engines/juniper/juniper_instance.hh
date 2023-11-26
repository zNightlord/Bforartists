/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup juniper
 *
 */

#pragma once

#include "BKE_object.hh"
#include "DEG_depsgraph.hh"
#include "DNA_lightprobe_types.h"
#include "DRW_render.h"
#include "juniper_shaders.hh"
#include "juniper_sync.hh"
#include "juniper_pipeline.hh"
#include "juniper_buffers.hh"
#include "draw_manager.hh"
#include "draw_handle.hh"
#include "juniper_lights.hh"


namespace blender::juniper {

class JNPR {
public:
    /* Submodules */
    MaterialManager materials;
    Pipelines pipelines;
    Buffers buffers;
    LightManager lights;
    SceneSync sync_;

    /* Viewport info */
    std::stringstream info;

    /* Inputs */
    Depsgraph *depsgraph;
    draw::Manager *manager;


    draw::PassSimple composite_ps = {"JNPR.Composite"};
    JNPR();
    ~JNPR();

    void init();

    void object_sync(Object* ob);
    void composite(draw::View &view);

    void begin_sync();

    void render();
};

}  // namespace blender::juniper
