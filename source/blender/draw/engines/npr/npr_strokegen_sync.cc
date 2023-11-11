/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Converts the different renderable object types to drawcalls.
 */

#include "bnpr_engine.h"

#include "DEG_depsgraph_query.h"
#include "DNA_curves_types.h"

#include "npr_strokegen_sync.hh"
#include "npr_strokegen_instance.hh"

#include "draw_cache_extract.hh"

namespace blender::npr::strokegen {
/* -------------------------------------------------------------------- */
/** \name Draw Data
 *
 * \{ */

static void draw_data_init_cb(struct DrawData *dd)
{
  /* Object has just been created or was never evaluated by the engine. */
  dd->recalc = ID_RECALC_ALL; /* Tag given ID for an update in all the dependency graphs. */
  }


  BnprDrawData& StrokegenSyncModule::sync_object(Object* ob)
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

  void StrokegenSyncModule::sync_mesh(
      Object* ob,
      const draw::ObjectRef& ob_ref,
      BnprDrawData& ob_draw_data,
      draw::ResourceHandle& rsc_handle,
      const DRWView* drw_view
      )
  {
    bool mesh_is_manifold;
    GPUBatch *gpu_batch_line_adj = DRW_cache_object_edge_detection_get(ob, &mesh_is_manifold);
    GPUBatch *gpu_batch_surf = DRW_cache_object_surface_get(ob_ref.object);

    if (gpu_batch_line_adj == nullptr) return;
    if (gpu_batch_line_adj->elem == nullptr) return;
    if (gpu_batch_surf == nullptr) return;
    if (gpu_batch_surf->elem == nullptr) return;

    // Old way to do this:
    // See "draw_subdiv_build_tris_buffer"
    // const char *defines = "#define XXX\n";
    // GPUShader *shader = get_strokegen_shader(...)
    //
    // eevvee_next way to do this:
    //  strokegen_passes.dispatch_extract_mesh_contour(ob);
    //  strokegen_passes.dispatch_XXX(...);
    //  ... ... ...


    if (inst_.strokegen_passes.boostrap_before_extract_first_batch)
    {
      inst_.strokegen_passes.append_per_mesh_pass(
          ob, gpu_batch_line_adj, gpu_batch_surf /**gpu_batch_surf*/, rsc_handle, drw_view);  // bootstrapping
      inst_.strokegen_passes.boostrap_before_extract_first_batch = false;  // switch off
    }
    inst_.strokegen_passes.append_per_mesh_pass(
        ob, gpu_batch_line_adj, gpu_batch_surf /**gpu_batch_surf*/, rsc_handle, drw_view);

  }



  }  // namespace blender::npr::strokegen
