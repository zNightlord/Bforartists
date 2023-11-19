/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup juniper
 *
 */

#pragma once

#include "BKE_duplilist.h"
#include "BLI_ghash.h"
#include "BLI_map.hh"
#include "DNA_object_types.h"
#include "DRW_render.h"
#include "GPU_material.h"

namespace blender::juniper {

class JNPR;

class SceneSync {
private:
    JNPR &jnpr_;

public:
    SceneSync(JNPR &jnpr) : jnpr_(jnpr){};
    ~SceneSync() = default;

    void sync_object(Object *ob);

    void sync_mesh(Object *ob);

    void sync_light(Object *ob);
};

} /* namespace blender::juniper */