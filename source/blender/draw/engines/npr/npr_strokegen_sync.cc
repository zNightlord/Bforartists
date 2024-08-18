/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Converts the different renderable object types to drawcalls.
 */


#include "DEG_depsgraph_query.hh"
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


  ObjectHandle StrokegenSyncModule::sync_object(const ObjectRef& ob_ref)
  {
    ObjectKey key(ob_ref.object);

    ObjectHandle &handle = ob_handles.lookup_or_add_cb(key, [&]() {
      ObjectHandle new_handle;
      new_handle.object_key = key;
      return new_handle;
    });

    handle.recalc = inst_.get_recalc_flags(ob_ref);

    return handle;
  }

  void StrokegenSyncModule::sync_mesh(
      Object* ob,
      const draw::ObjectRef& ob_ref,
      draw::ResourceHandle& rsc_handle,
      const DRWView* drw_view
      )
  {
    bool mesh_is_manifold;
    gpu::Batch *gpu_batch_line_adj = DRW_cache_object_edge_detection_get(ob, &mesh_is_manifold);
    gpu::Batch *gpu_batch_surf = DRW_cache_object_surface_get(ob_ref.object);

    if (gpu_batch_line_adj == nullptr)
      return;
    if (gpu_batch_line_adj->elem == nullptr)
      return;
    if (gpu_batch_surf == nullptr)
      return;
    if (gpu_batch_surf->elem == nullptr)
      return;

    // Old way to do this:
    // See "draw_subdiv_build_tris_buffer"
    // const char *defines = "#define XXX\n";
    // GPUShader *shader = get_strokegen_shader(...)
    //
    // eevvee_next way to do this:
    //  strokegen_passes.dispatch_extract_mesh_contour(ob);
    //  strokegen_passes.dispatch_XXX(...);
    //  ... ... ...

    const PerObjectStrokegenSettings &ui_input = ob->strokegen_settings;
    if (ui_input.curve_type == 0)
      return; 

    if (inst_.strokegen_passes.boostrap_before_extract_first_batch) {
      // bootstrapping
      inst_.strokegen_passes.append_per_mesh_pass(
          ob,
          gpu_batch_line_adj,
          gpu_batch_surf /**gpu_batch_surf*/,
          rsc_handle,
          drw_view
        );

      inst_.strokegen_passes.boostrap_before_extract_first_batch = false; // switch off
    }
    inst_.strokegen_passes.append_per_mesh_pass(
        ob,
        gpu_batch_line_adj,
        gpu_batch_surf /**gpu_batch_surf*/,
        rsc_handle,
        drw_view
        );
    inst_.strokegen_passes.append_pass_remeshed_surface_depth_drawcall();

    inst_.strokegen_passes.has_strokegen_enabled_mesh = true; 
  }



  }  // namespace blender::npr::strokegen
