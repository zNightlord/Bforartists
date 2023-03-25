/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2016 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 */

#pragma once

#include "DRW_render.h"
#include "RE_engine.h"

#ifdef __cplusplus
extern "C" { // to satisfy .cc files
#endif

extern DrawEngineType draw_engine_bnpr_type;
extern RenderEngineType DRW_engine_viewport_bnpr_type;

// Also I spent some time reading the mesh extraction code,
// feel like I could add one or more new mesh buffers & extractors?
// In that way,
// 1) Never need to touch the DNA data;
// 2) DrawManager updates these buffers when populates caches for the DrawEngine;
// 3) DrawEngine uses a eevee-next-style sync module, in its function "sync_mesh":
// -- 3.1) duplicate a ref from the cache, and
// -- 3.2) issue ordered compute dispatches & indirect procedural draws;

#ifdef __cplusplus
}
#endif
