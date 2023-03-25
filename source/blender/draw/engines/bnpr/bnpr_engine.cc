/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2016 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 *
 * Simple engine for drawing color and/or depth.
 * When we only need simple flat shaders.
 */

#include "DRW_render.h"

#include "BKE_object.h"
#include "BKE_paint.h"

#include "BLI_alloca.h"

#include "GPU_shader.h"

#include "bnpr_engine.h"

#include "bnpr_instance.hh"
#include "bnpr_shader.hh"
#include "intern/mallocn_intern.h"


using namespace blender;



#define bnpr_ENGINE "BLENDER_bnpr"

/* *********** LISTS *********** */

/* GPUViewport.storage
 * Is freed every time the viewport engine changes. */
typedef struct bnpr_StorageList {
  struct bnpr_PrivateData *g_data;
} bnpr_StorageList;

typedef struct bnpr_PassList {
  struct DRWPass *depth_pass[2];
  struct DRWPass *depth_pass_pointcloud[2];
  struct DRWPass *depth_pass_cull[2];
} bnpr_PassList;

// Per-engine data
// Sent by draw manager to the draw-engine.
typedef struct bnpr_Data {
  void *engine_type;
  DRWViewportEmptyList *fbl;
  DRWViewportEmptyList *txl;
  bnpr_PassList *psl;
  bnpr_StorageList *stl;

  bnpr::Instance *instance;
  char info[GPU_INFO_SIZE];

} bnpr_Data;

/* *********** STATIC *********** */

typedef struct bnpr_PrivateData {
  DRWShadingGroup *depth_shgrp[2];
  DRWShadingGroup *depth_shgrp_cull[2];
  DRWShadingGroup *depth_hair_shgrp[2];
  DRWShadingGroup *depth_curves_shgrp[2];
  DRWShadingGroup *depth_pointcloud_shgrp[2];
  bool use_material_slot_selection;
} bnpr_PrivateData; /* Transient data */



static bool check_bnpr_support()
{
  return GPU_shader_storage_buffer_objects_support();
}


static void bnpr_engine_init(void *vedata)
{
  if (!check_bnpr_support()) { return; }

  bnpr_Data *ved = reinterpret_cast<bnpr_Data *>(vedata);
  if (ved->instance == nullptr) {
    ved->instance = new bnpr::Instance();
  }

  draw::Manager *drw_mgr = DRW_manager_get();

  const DRWContextState *ctx_state = DRW_context_state_get();
  View3D *v3d = ctx_state->v3d;
  RegionView3D *rv3d = ctx_state->rv3d;

  Object *camera = nullptr;
  if (v3d && rv3d && rv3d->persp == RV3D_CAMOB)
    camera = v3d->camera;

  const DRWView* default_drw_view = DRW_view_default_get();


  ved->instance->init(
    ctx_state->depsgraph,
    drw_mgr,
    ctx_state->v3d, rv3d,
    default_drw_view,
    camera
  );
}



static void bnpr_draw_scene_legacy(void *vedata)
{
  bnpr_PassList *psl = ((bnpr_Data *)vedata)->psl;

  DRW_draw_pass(psl->depth_pass[0]);
  DRW_draw_pass(psl->depth_pass_pointcloud[0]);
  DRW_draw_pass(psl->depth_pass_cull[0]);
  DRW_draw_pass(psl->depth_pass[1]);
  DRW_draw_pass(psl->depth_pass_pointcloud[1]);
  DRW_draw_pass(psl->depth_pass_cull[1]);
}

static void bnpr_draw_scene(void *vedata)
{
  bnpr_Data *ved = reinterpret_cast<bnpr_Data *>(vedata);
  if (!check_bnpr_support()) {
    STRNCPY(ved->info, "Error: No shader storage buffer support, required by bnpr.");
    return;
  }


  // TODO: not sure which is better, frame-buffer or texture list?
  DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
  // DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

  const DRWView *default_view = DRW_view_default_get();
  const DRWView *active_view = DRW_view_get_active();
  draw::Manager *manager = DRW_manager_get();
  draw::View view("DefaultView", default_view);

  // draw passes
  ved->instance->draw_viewport(*manager, view);
  // display error msg at the top of the render viewport
  STRNCPY(ved->info, ved->instance->info.c_str());


  /* Reset view for other following engines. */
  DRW_view_set_active(nullptr);
}




static void bnpr_cache_init_legacy(void *vedata)
{
  bnpr_PassList *psl = static_cast<bnpr_Data*>(vedata)->psl;
  bnpr_StorageList *stl = static_cast<bnpr_Data*>(vedata)->stl;
  DRWShadingGroup *grp;

  const DRWContextState *draw_ctx = DRW_context_state_get();

  if (!stl->g_data) {
    /* Alloc transient pointers */
    stl->g_data = (bnpr_PrivateData* )MEM_callocN(sizeof(*stl->g_data), __func__);
  }

  stl->g_data->use_material_slot_selection = DRW_state_is_material_select();

  /* Twice for normal and in front objects. */
  for (int i = 0; i < 2; i++) {
    DRWState clip_state = static_cast<DRWState>(
      (draw_ctx->sh_cfg == GPU_SHADER_CFG_CLIPPED) ? DRW_STATE_CLIP_PLANES : 0
    );
    DRWState infront_state = static_cast<DRWState>(
      (DRW_state_is_select() && (i == 1)) ? DRW_STATE_IN_FRONT_SELECT : 0
    );
    DRWState state = DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL;

    blender::bnpr::ShaderModule* shaderModule =
      blender::bnpr::ShaderModule::module_get();

    GPUShader *sh = DRW_state_is_select() ?
                        shaderModule->static_shader_get(blender::bnpr::DEPTH_CONSERVATIVE) :
                        shaderModule->static_shader_get(blender::bnpr::DEPTH);

    DRW_PASS_CREATE(psl->depth_pass[i], state | clip_state | infront_state);
    stl->g_data->depth_shgrp[i] = grp = DRW_shgroup_create(sh, psl->depth_pass[i]);
    DRW_shgroup_uniform_block(grp, "globalsBlock", G_draw.block_ubo);

    sh = DRW_state_is_select() ?
        shaderModule->static_shader_get(blender::bnpr::POINTCLOUD_DEPTH_CONSERVATIVE) :
        shaderModule->static_shader_get(blender::bnpr::POINTCLOUD_DEPTH);
    DRW_PASS_CREATE(psl->depth_pass_pointcloud[i], state | clip_state | infront_state);
    stl->g_data->depth_pointcloud_shgrp[i] = grp = DRW_shgroup_create(
        sh, psl->depth_pass_pointcloud[i]);
    DRW_shgroup_uniform_block(grp, "globalsBlock", G_draw.block_ubo);

    stl->g_data->depth_hair_shgrp[i] = grp = DRW_shgroup_create(
        shaderModule->static_shader_get(blender::bnpr::DEPTH), psl->depth_pass[i]);
    DRW_shgroup_uniform_block(grp, "globalsBlock", G_draw.block_ubo);

    stl->g_data->depth_curves_shgrp[i] = grp = DRW_shgroup_create(
        shaderModule->static_shader_get(blender::bnpr::CURVES_DEPTH), psl->depth_pass[i]);
    DRW_shgroup_uniform_block(grp, "globalsBlock", G_draw.block_ubo);

    sh = DRW_state_is_select() ? shaderModule->static_shader_get(blender::bnpr::DEPTH_CONSERVATIVE) :
                                 shaderModule->static_shader_get(blender::bnpr::DEPTH);
    state |= DRW_STATE_CULL_BACK;
    DRW_PASS_CREATE(psl->depth_pass_cull[i], state | clip_state | infront_state);
    stl->g_data->depth_shgrp_cull[i] = grp = DRW_shgroup_create(sh, psl->depth_pass_cull[i]);
    DRW_shgroup_uniform_block(grp, "globalsBlock", G_draw.block_ubo);
  }
}
static void bnpr_cache_init(void *vedata)
{
  if (!check_bnpr_support()) return;

  draw::Manager* drwmgr = DRW_manager_get();
  reinterpret_cast<bnpr_Data *>(vedata)->instance->begin_sync(*drwmgr);


  bnpr_cache_init_legacy(vedata);
}









/* TODO(fclem): DRW_cache_object_surface_material_get needs a refactor to allow passing NULL
 * instead of gpumat_array. Avoiding all this boilerplate code. */
static struct GPUBatch **bnpr_object_surface_material_get(Object *ob)
{
  const int materials_len = DRW_cache_object_material_count_get(ob);
  struct GPUMaterial **gpumat_array =
    static_cast<GPUMaterial**>(BLI_array_alloca(gpumat_array, materials_len));
  memset(gpumat_array, 0, sizeof(*gpumat_array) * materials_len);

  return DRW_cache_object_surface_material_get(ob, gpumat_array, materials_len);
}

static void bnpr_cache_populate_particles(void *vedata, Object *ob)
{
  // do nothing here.
}

static void bnpr_cache_populate_legacy(void *vedata, Object *ob)
{
  bnpr_StorageList *stl = ((bnpr_Data *)vedata)->stl;

  /* TODO(fclem): fix selection of smoke domains. */

  if (!DRW_object_is_renderable(ob) || (ob->dt < OB_SOLID)) {
    return;
  }

  const DRWContextState *draw_ctx = DRW_context_state_get();
  if (ob != draw_ctx->object_edit) {
    bnpr_cache_populate_particles(vedata, ob);
  }

  const bool do_in_front = (ob->dtx & OB_DRAW_IN_FRONT) != 0;
  if (ob->type == OB_CURVES) {
    DRW_shgroup_curves_create_sub(ob, stl->g_data->depth_curves_shgrp[do_in_front], NULL);
  }

  /* Make flat object selectable in ortho view if wireframe is enabled. */
  if ((draw_ctx->v3d->overlay.flag & V3D_OVERLAY_WIREFRAMES) ||
      (draw_ctx->v3d->shading.type == OB_WIRE) || (ob->dtx & OB_DRAWWIRE) || (ob->dt == OB_WIRE)) {
    int flat_axis = 0;
    bool is_flat_object_viewed_from_side = ((draw_ctx->rv3d->persp == RV3D_ORTHO) &&
                                            DRW_object_is_flat(ob, &flat_axis) &&
                                            DRW_object_axis_orthogonal_to_view(ob, flat_axis));

    if (is_flat_object_viewed_from_side) {
      /* Avoid losing flat objects when in ortho views (see T56549) */
      struct GPUBatch *geom = DRW_cache_object_all_edges_get(ob);
      if (geom) {
        DRW_shgroup_call(stl->g_data->depth_shgrp[do_in_front], geom, ob);
      }
      return;
    }
  }

  const bool use_sculpt_pbvh = BKE_sculptsession_use_pbvh_draw(ob, draw_ctx->rv3d) &&
                               !DRW_state_is_image_render();
  const bool do_cull = (draw_ctx->v3d &&
                        (draw_ctx->v3d->shading.flag & V3D_SHADING_BACKFACE_CULLING));

  DRWShadingGroup *shgrp = NULL;

  if (ob->type == OB_POINTCLOUD) {
    shgrp = stl->g_data->depth_pointcloud_shgrp[do_in_front];
  }
  else {
    shgrp = (do_cull) ? stl->g_data->depth_shgrp_cull[do_in_front] :
                        stl->g_data->depth_shgrp[do_in_front];
  }

  if (use_sculpt_pbvh) {
    DRW_shgroup_call_sculpt(shgrp, ob, false, false, false, false, false);
  }
  else {
    if (stl->g_data->use_material_slot_selection && BKE_object_supports_material_slots(ob)) {
      struct GPUBatch **geoms = bnpr_object_surface_material_get(ob);
      if (geoms) {
        const int materials_len = DRW_cache_object_material_count_get(ob);
        for (int i = 0; i < materials_len; i++) {
          if (geoms[i] == NULL) {
            continue;
          }
          const short material_slot_select_id = i + 1;
          DRW_select_load_id(ob->runtime.select_id | (material_slot_select_id << 16));
          DRW_shgroup_call(shgrp, geoms[i], ob);
        }
      }
    }
    else {
      struct GPUBatch *geom = DRW_cache_object_surface_get(ob);
      if (geom) {
        DRW_shgroup_call(shgrp, geom, ob);
      }
    }
  }
}
static void bnpr_cache_populate(void *vedata, Object *object)
{
  if (!check_bnpr_support()) return;

  bnpr_Data* ved = reinterpret_cast<bnpr_Data *>(vedata);

  draw::Manager* drw_mgr = DRW_manager_get();
  draw::ObjectRef ref {
    object,
    DRW_object_get_dupli(object),
    DRW_object_get_dupli_parent(object)
  };

  ved->instance->object_sync(*drw_mgr, ref);



  bnpr_cache_populate_legacy(vedata, object);
}


static void bnpr_cache_finish_legacy(void *vedata)
{
  bnpr_StorageList *stl = ((bnpr_Data *)vedata)->stl;

  UNUSED_VARS(stl);
}
static void bnpr_cache_finish(void *vedata)
{
  if (!check_bnpr_support()) return;

  bnpr_Data* ved = reinterpret_cast<bnpr_Data*>(vedata);
  draw::Manager* drw_mgr = DRW_manager_get();

  ved->instance->end_sync(*drw_mgr);



  bnpr_cache_finish_legacy(vedata);
}

static void bnpr_instance_free(void *instance) {
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  delete reinterpret_cast<bnpr::Instance *>(instance);
}

static void bnpr_engine_free(void)
{
  bnpr::ShaderModule::module_free();
}

static void bnpr_render_to_image(void *vedata, struct RenderEngine *engine,
                                    struct RenderLayer *layer,
                                    const struct rcti *rect) {
  UNUSED_VARS(vedata, engine, layer);
}






static const DrawEngineDataSize bnpr_data_size = DRW_VIEWPORT_DATA_SIZE(bnpr_Data);

DrawEngineType draw_engine_bnpr_type = {
    NULL,
    NULL,
    N_("bnpr"),
    &bnpr_data_size,
    bnpr_engine_init,
    &bnpr_engine_free,
    &bnpr_instance_free,
    &bnpr_cache_init,
    &bnpr_cache_populate,
    &bnpr_cache_finish,
    &bnpr_draw_scene,
    NULL,
    NULL,
    bnpr_render_to_image,
    NULL,
};

RenderEngineType DRW_engine_viewport_bnpr_type = {
  nullptr,
  nullptr,
  "bnpr_VIEW",
  N_("bnpr"),
  RE_INTERNAL | RE_USE_PREVIEW | RE_USE_STEREO_VIEWPORT | RE_USE_GPU_CONTEXT,
  nullptr,
  &DRW_render_to_image,
  nullptr,
  nullptr,
  nullptr,
  nullptr,
  nullptr,
  nullptr,
  nullptr, // TODO: impl this
  &draw_engine_bnpr_type,
  {nullptr, nullptr, nullptr},
};

#undef bnpr_ENGINE
