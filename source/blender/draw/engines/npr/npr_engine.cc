/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_object.h"
#include "BKE_report.h"
#include "DEG_depsgraph_query.h"
#include "DNA_camera_types.h"
#include "DRW_render.h"
#include "ED_view3d.h"
#include "GPU_capabilities.h"
#include "draw_manager.hh"
#include "draw_pass.hh"

#include "npr_engine.h" /* Own include. */
#include "npr_instance.hh"


/* -------------------------------------------------------------------- */
/** \name Interface with legacy C DRW manager
 * \{ */

using namespace blender;

struct NPR_Data {
  DrawEngineType *engine_type;
  DRWViewportEmptyList *fbl;
  DRWViewportEmptyList *txl;
  DRWViewportEmptyList *psl;
  DRWViewportEmptyList *stl;
  npr::Instance *instance;

  char info[GPU_INFO_SIZE];
};

static void npr_engine_init(void *vedata)
{
  /* TODO(fclem): Remove once it is minimum required. */
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }

  NPR_Data *ved = reinterpret_cast<NPR_Data *>(vedata);
  if (ved->instance == nullptr) {
    ved->instance = new npr::Instance();
  }

  ved->instance->init();
}

static void npr_cache_init(void *vedata)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  reinterpret_cast<NPR_Data *>(vedata)->instance->begin_sync();
}

static void npr_cache_populate(void *vedata, Object *object)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  draw::Manager *manager = DRW_manager_get();

  draw::ObjectRef ref;
  ref.object = object;
  ref.dupli_object = DRW_object_get_dupli(object);
  ref.dupli_parent = DRW_object_get_dupli_parent(object);

  reinterpret_cast<NPR_Data *>(vedata)->instance->object_sync(*manager, ref);
}

static void npr_cache_finish(void *vedata)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  reinterpret_cast<NPR_Data *>(vedata)->instance->end_sync();
}

static void npr_draw_scene(void *vedata)
{
  NPR_Data *ved = reinterpret_cast<NPR_Data *>(vedata);
  if (!GPU_shader_storage_buffer_objects_support()) {
    STRNCPY(ved->info, "Error: No shader storage buffer support");
    return;
  }
  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
  draw::Manager *manager = DRW_manager_get();

  const DRWView *default_view = DRW_view_default_get();
  draw::View view("DefaultView", default_view);

  ved->instance->draw_viewport(*manager, view, dtxl->depth, dtxl->color);
}

static void npr_instance_free(void *instance)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  delete reinterpret_cast<npr::Instance *>(instance);
}

static void npr_view_update(void *vedata)
{
  NPR_Data *ved = reinterpret_cast<NPR_Data *>(vedata);
}

static void npr_id_update(void *vedata, struct ID *id)
{
  UNUSED_VARS(vedata, id);
}

/* RENDER */

static bool npr_render_framebuffers_init(void)
{
  /* For image render, allocate own buffers because we don't have a viewport. */
  const float2 viewport_size = DRW_viewport_size_get();
  const int2 size = {int(viewport_size.x), int(viewport_size.y)};

  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

  /* When doing a multi view rendering the first view will allocate the buffers
   * the other views will reuse these buffers */
  if (dtxl->color == nullptr) {
    BLI_assert(dtxl->depth == nullptr);
    eGPUTextureUsage usage = GPU_TEXTURE_USAGE_GENERAL;
    dtxl->color = GPU_texture_create_2d(
        "txl.color", size.x, size.y, 1, GPU_RGBA16F, usage, nullptr);
    dtxl->depth = GPU_texture_create_2d(
        "txl.depth", size.x, size.y, 1, GPU_DEPTH24_STENCIL8, usage, nullptr);
  }

  if (!(dtxl->depth && dtxl->color)) {
    return false;
  }

  DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();

  GPU_framebuffer_ensure_config(
      &dfbl->default_fb,
      {GPU_ATTACHMENT_TEXTURE(dtxl->depth), GPU_ATTACHMENT_TEXTURE(dtxl->color)});

  GPU_framebuffer_ensure_config(&dfbl->depth_only_fb,
                                {GPU_ATTACHMENT_TEXTURE(dtxl->depth), GPU_ATTACHMENT_NONE});

  GPU_framebuffer_ensure_config(&dfbl->color_only_fb,
                                {GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(dtxl->color)});

  return GPU_framebuffer_check_valid(dfbl->default_fb, nullptr) &&
         GPU_framebuffer_check_valid(dfbl->color_only_fb, nullptr) &&
         GPU_framebuffer_check_valid(dfbl->depth_only_fb, nullptr);
}

#ifdef DEBUG
/* This is just to ease GPU debugging when the frame delimiter is set to Finish */
#  define GPU_FINISH_DELIMITER() GPU_finish()
#else
#  define GPU_FINISH_DELIMITER()
#endif

static void write_render_color_output(struct RenderLayer *layer,
                                      const char *viewname,
                                      GPUFrameBuffer *fb,
                                      const struct rcti *rect)
{
  RenderPass *rp = RE_pass_find_by_name(layer, RE_PASSNAME_COMBINED, viewname);
  if (rp) {
    GPU_framebuffer_bind(fb);
    GPU_framebuffer_read_color(fb,
                               rect->xmin,
                               rect->ymin,
                               BLI_rcti_size_x(rect),
                               BLI_rcti_size_y(rect),
                               4,
                               0,
                               GPU_DATA_FLOAT,
                               rp->rect);
  }
}

static void write_render_z_output(struct RenderLayer *layer,
                                  const char *viewname,
                                  GPUFrameBuffer *fb,
                                  const struct rcti *rect,
                                  float4x4 winmat)
{
  RenderPass *rp = RE_pass_find_by_name(layer, RE_PASSNAME_Z, viewname);
  if (rp) {
    GPU_framebuffer_bind(fb);
    GPU_framebuffer_read_depth(fb,
                               rect->xmin,
                               rect->ymin,
                               BLI_rcti_size_x(rect),
                               BLI_rcti_size_y(rect),
                               GPU_DATA_FLOAT,
                               rp->rect);

    int pix_num = BLI_rcti_size_x(rect) * BLI_rcti_size_y(rect);

    /* Convert GPU depth [0..1] to view Z [near..far] */
    if (DRW_view_is_persp_get(nullptr)) {
      for (float &z : MutableSpan(rp->rect, pix_num)) {
        if (z == 1.0f) {
          z = 1e10f; /* Background */
        }
        else {
          z = z * 2.0f - 1.0f;
          z = winmat[3][2] / (z + winmat[2][2]);
        }
      }
    }
    else {
      /* Keep in mind, near and far distance are negatives. */
      float near = DRW_view_near_distance_get(nullptr);
      float far = DRW_view_far_distance_get(nullptr);
      float range = fabsf(far - near);

      for (float &z : MutableSpan(rp->rect, pix_num)) {
        if (z == 1.0f) {
          z = 1e10f; /* Background */
        }
        else {
          z = z * range - near;
        }
      }
    }
  }
}

static void npr_render_to_image(void *vedata,
                                struct RenderEngine *engine,
                                struct RenderLayer *layer,
                                const struct rcti *rect)
{
  /* TODO(fclem): Remove once it is minimum required. */
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }

  if (!npr_render_framebuffers_init()) {
    RE_engine_report(engine, RPT_ERROR, "Failed to allocate GPU buffers");
    return;
  }

  GPU_FINISH_DELIMITER();

  /* Setup */

  DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();
  const DRWContextState *draw_ctx = DRW_context_state_get();
  Depsgraph *depsgraph = draw_ctx->depsgraph;

  NPR_Data *ved = reinterpret_cast<NPR_Data *>(vedata);
  if (ved->instance == nullptr) {
    ved->instance = new npr::Instance();
  }

  /* TODO(sergey): Shall render hold pointer to an evaluated camera instead? */
  Object *camera_ob = DEG_get_evaluated_object(depsgraph, RE_GetCamera(engine->re));

  /* Set the perspective, view and window matrix. */
  float4x4 winmat, viewmat, viewinv;
  RE_GetCameraWindow(engine->re, camera_ob, winmat.ptr());
  RE_GetCameraModelMatrix(engine->re, camera_ob, viewinv.ptr());
  viewmat = math::invert(viewinv);

  DRWView *view = DRW_view_create(viewmat.ptr(), winmat.ptr(), nullptr, nullptr, nullptr);
  DRW_view_default_set(view);
  DRW_view_set_active(view);

  /* Render */
  do {
    if (RE_engine_test_break(engine)) {
      break;
    }

    ved->instance->init(camera_ob);

    DRW_manager_get()->begin_sync();

    npr_cache_init(vedata);
    auto npr_render_cache = [](void *vedata,
                               struct Object *ob,
                               struct RenderEngine * /*engine*/,
                               struct Depsgraph * /*depsgraph*/) {
      npr_cache_populate(vedata, ob);
    };
    DRW_render_object_iter(vedata, engine, depsgraph, npr_render_cache);
    npr_cache_finish(vedata);

    DRW_manager_get()->end_sync();

    /* Also we weed to have a correct FBO bound for #DRW_curves_update */
    // GPU_framebuffer_bind(dfbl->default_fb);
    // DRW_curves_update(); /* TODO(@pragma37): Check this once curves are implemented */

    npr_draw_scene(vedata);

    /* Perform render step between samples to allow
     * flushing of freed GPUBackend resources. */
    GPU_render_step();
    GPU_FINISH_DELIMITER();
  } while (false); /* Sample count goes here */

  const char *viewname = RE_GetActiveRenderView(engine->re);
  write_render_color_output(layer, viewname, dfbl->default_fb, rect);
  write_render_z_output(layer, viewname, dfbl->default_fb, rect, winmat);
}

static void npr_render_update_passes(RenderEngine *engine, Scene *scene, ViewLayer *view_layer)
{
  if (view_layer->passflag & SCE_PASS_COMBINED) {
    RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_COMBINED, 4, "RGBA", SOCK_RGBA);
  }
  if (view_layer->passflag & SCE_PASS_Z) {
    RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_Z, 1, "Z", SOCK_FLOAT);
  }
}

extern "C" {

static const DrawEngineDataSize npr_data_size = DRW_VIEWPORT_DATA_SIZE(NPR_Data);

DrawEngineType draw_engine_npr_type = {
    nullptr,
    nullptr,
    N_("NPR"),
    &npr_data_size,
    &npr_engine_init,
    nullptr,
    &npr_instance_free,
    &npr_cache_init,
    &npr_cache_populate,
    &npr_cache_finish,
    &npr_draw_scene,
    &npr_view_update,
    &npr_id_update,
    &npr_render_to_image,
    nullptr,
};

RenderEngineType DRW_engine_viewport_npr_type = {
    nullptr,
    nullptr,
    "BLENDER_NPR",
    N_("NPR"),
    RE_INTERNAL | RE_USE_STEREO_VIEWPORT | RE_USE_GPU_CONTEXT,
    nullptr,
    &DRW_render_to_image,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    &npr_render_update_passes,
    &draw_engine_npr_type,
    {nullptr, nullptr, nullptr},
};
}

/** \} */
