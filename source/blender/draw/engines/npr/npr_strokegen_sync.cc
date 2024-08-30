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


    // Register object to gpu scene
    ObjectKey ob_key(ob);
    int gpu_obj_id = inst_.object_gpu_id_map.size();
    inst_.object_gpu_id_map.add(ob_key.hash(), gpu_obj_id);
    inst_.strokegen_buffers.sync_object(ob); // add object data to gpu buffers


    // Record per-object strokegen shaders
    if (inst_.strokegen_passes.boostrap_before_extract_first_batch) {
      // bootstrapping
      inst_.strokegen_passes.append_per_mesh_pass(
          ob, gpu_obj_id, 
          gpu_batch_line_adj,
          gpu_batch_surf /**gpu_batch_surf*/,
          rsc_handle, drw_view
      );

      inst_.strokegen_passes.boostrap_before_extract_first_batch = false; // switch off
    }
    inst_.strokegen_passes.append_per_mesh_pass(
        ob, gpu_obj_id, 
        gpu_batch_line_adj,
        gpu_batch_surf /**gpu_batch_surf*/,
        rsc_handle, drw_view
        );
    inst_.strokegen_passes.append_pass_remeshed_surface_depth_drawcall(
      STROKEGEN_SHADING_TYPE_TRANSPARENT == ob_ref.object->strokegen.surface_shading_type
    );

    inst_.has_strokegen_enabled_mesh = true;
  }


}  // namespace blender::npr::strokegen
