/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Converts the different renderable object types to drawcalls.
 */

#include "bnpr_engine.h"

#include "BKE_gpencil.h"
#include "BKE_object.h"
#include "DEG_depsgraph_query.h"
#include "DNA_curves_types.h"
#include "DNA_gpencil_types.h"
#include "DNA_modifier_types.h"
#include "DNA_particle_types.h"

#include "bnpr_instance.hh"
#include "bnpr_sync.hh"

#include "draw_cache_extract.hh"
#include "draw_cache_impl.h"

namespace blender::bnpr
{
  /* -------------------------------------------------------------------- */
  /** \name Draw Data
   *
   * \{ */

  static void draw_data_init_cb(struct DrawData* dd)
  {
    /* Object has just been created or was never evaluated by the engine. */
    dd->recalc = ID_RECALC_ALL; /* Tag given ID for an update in all the dependency graphs. */
  }


  BnprDrawData& SyncModule::sync_object(Object* ob)
  {
    DrawEngineType* owner = (DrawEngineType*)&DRW_engine_viewport_bnpr_type;
    struct DrawData* dd = DRW_drawdata_ensure(
      (ID*)ob, owner, sizeof(BnprDrawData), draw_data_init_cb, nullptr);

    BnprDrawData &dd_bnpr = *reinterpret_cast<BnprDrawData*>(dd); // draw-engine specific draw data.

    if (dd_bnpr.object_key.ob == nullptr)
    {
      dd_bnpr.object_key = ObjectKey(ob);
    }


    return dd_bnpr;
  }

  void SyncModule::sync_mesh(
    Object* ob, BnprDrawData& ob_draw_data,
    draw::ResourceHandle res_handle, const draw::ObjectRef& ob_ref
  )
  {
    bool mesh_is_manifold;
    GPUBatch *geobatch = DRW_cache_object_edge_detection_get(ob, &mesh_is_manifold);

    if (geobatch == nullptr) return;
    if (geobatch->elem == nullptr) return;

    // Old way to do this:
    // See "draw_subdiv_build_tris_buffer"
    // const char *defines = "#define XXX\n";
    // GPUShader *shader = get_strokegen_shader(...)
    //
    // eevvee_next way to do this:
    //  strokegen_passes.dispatch_extract_mesh_contour(ob);
    //  strokegen_passes.dispatch_XXX(...);
    //  ... ... ...
    inst_.strokegen_passes.rebuild_sub_pass_extract_mesh_geom(ob, geobatch);



  }
}
