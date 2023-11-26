/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_global.h"
#include "BLI_rect.h"

#include "BKE_editmesh.h"
#include "BKE_modifier.h"
#include "BKE_object.hh"
#include "BKE_paint.hh"
#include "BKE_particle.h"
#include "BKE_pbvh.h"
#include "BKE_report.h"
#include "DEG_depsgraph_query.hh"
#include "DNA_fluid_types.h"
#include "ED_view3d.hh"
#include "GPU_framebuffer.h"

#include "GPU_compute.h"

#include "DNA_object_types.h"
#include "DRW_render.h"
#include "draw_manager.hh"
#include "draw_pass.hh"

#include "RE_pipeline.h"

#include "juniper_engine.h" /* Own include. */

#include "juniper_instance.hh"
#include "draw_debug.hh"

#include "NOD_socket.hh"
#include "WM_types.hh"
#include "WM_api.hh"

#define JNPR_ENGINE "JNPR"

using namespace blender;

namespace blender::juniper {

using namespace draw;

/* Sync any object type */
void JNPR::object_sync(Object *ob)
{
    if (!DRW_object_is_renderable(ob)) {
        return;
    }

    if (ob->sculpt) {
        // TODO(Late): add pbvh support
        return;
    }

    if (ob->type == OB_MESH && ob->modifiers.first != nullptr) {
        LISTBASE_FOREACH(ModifierData *, md, &ob->modifiers) {
            if (md->type != eModifierType_ParticleSystem) {
                continue;
            }

            // TODO(Late): Add particle/hair support
        }
    }

    if (!(DRW_object_visibility_in_active_context(ob) & OB_VISIBLE_SELF)) {
        return;
    }

    if ((ob->dt < OB_SOLID) && !DRW_state_is_scene_render()) {
        return;
    }

    this->sync_.sync_object(ob);
}


void JNPR::composite(draw::View &view)
{
  composite_ps.init();
  composite_ps.shader_set(materials.static_shader_get(FILM_COMPOSITE));
  composite_ps.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_ALWAYS);
  composite_ps.bind_texture(SH_TEX_COLOR_SLOT, buffers.color_tx);
  composite_ps.bind_texture(SH_TEX_DEPTH_SLOT, buffers.depth_tx);
  composite_ps.bind_texture(SH_TEX_NORMAL_SLOT, buffers.normal_tx);

  lights.bind_resources(&composite_ps);

  composite_ps.draw_procedural(GPU_PRIM_TRIS, 1, 3);
  DRW_manager_get()->submit(composite_ps, view);
}

JNPR::JNPR() :
  materials(*this),
  pipelines(*this),
  buffers(*this),
  lights(*this),
  sync_(*this)
{

}

void JNPR::init() {
    manager = DRW_manager_get();
}

void JNPR::begin_sync() {
    DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
    int2 size = int2(GPU_texture_width(dtxl->color), GPU_texture_height(dtxl->color));

    // Render buffers must be acquired first before use in shaders/pipelines
    buffers.acquire(size);

    materials.begin_sync();
    pipelines.sync();
    lights.begin_sync();
}

JNPR::~JNPR() {

}

void JNPR::render() {
  const DefaultFramebufferList *dfbl = DRW_viewport_framebuffer_list_get();

  float clear_col[4] = {0.0f,
                        0.0f,
                        0.0f,
                        0.0f};

  draw::View drw_view("Main View", DRW_view_default_get());

// Draw prepass
  GPU_framebuffer_bind(buffers.prepass_fb_);
  GPU_framebuffer_clear_color(buffers.prepass_fb_, clear_col);
  GPU_framebuffer_clear_depth(buffers.prepass_fb_, 1.0f);
  pipelines.test.prepass(drw_view);

// Copy Depth buffer to new texture
  GPU_framebuffer_bind(buffers.zcopy_fb_);
  GPU_framebuffer_blit(buffers.prepass_fb_, 0,
                       buffers.zcopy_fb_, 0,
                       eGPUFrameBufferBits::GPU_DEPTH_BIT);

// Draw main shading pass
  GPU_framebuffer_bind(buffers.shading_fb_);
  GPU_framebuffer_clear_color(buffers.shading_fb_, clear_col);
  pipelines.test.render(drw_view);

// Draw composite pass (default_fb is displayed in viewport)
  GPU_framebuffer_bind(dfbl->default_fb);
  GPU_framebuffer_clear_color(dfbl->default_fb, clear_col);
  this->composite(drw_view);

  buffers.release();
  // Reset view for other engines (e.g. Overlay)
  DRW_view_set_active(nullptr);
}


} /* namespace blender::juniper */

struct JNPR_Data {
  DrawEngineType *engine_type;
  DRWViewportEmptyList *fbl;
  DRWViewportEmptyList *txl;
  DRWViewportEmptyList *psl;
  DRWViewportEmptyList *stl;
  juniper::JNPR *instance;

  char info[GPU_INFO_SIZE];
};

static void juniper_engine_init(void *vedata)
{
  JNPR_Data *ved = reinterpret_cast<JNPR_Data *>(vedata);
  if (ved->instance == nullptr) {
    ved->instance = new juniper::JNPR();
  }

  const DRWContextState *ctx_state = DRW_context_state_get();
  Depsgraph *depsgraph = ctx_state->depsgraph;
  Scene *scene = ctx_state->scene;
  View3D *v3d = ctx_state->v3d;
  const ARegion *region = ctx_state->region;
  RegionView3D *rv3d = ctx_state->rv3d;

  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
  int2 size = int2(GPU_texture_width(dtxl->color), GPU_texture_height(dtxl->color));

  const DRWView *default_view = DRW_view_default_get();

#if false
  /* Copied from eevee_engine_init to get view rect. */
  Object *camera = nullptr;
  /* Get render borders. */
  rcti rect;
  BLI_rcti_init(&rect, 0, size[0], 0, size[1]);
  if (v3d) {
    if (rv3d && (rv3d->persp == RV3D_CAMOB)) {
      camera = v3d->camera;
    }

    if (v3d->flag2 & V3D_RENDER_BORDER) {
      if (camera) {
        rctf viewborder;
        /* TODO(fclem) Might be better to get it from DRW. */
        ED_view3d_calc_camera_border(scene, depsgraph, region, v3d, rv3d, &viewborder, false);
        float viewborder_sizex = BLI_rctf_size_x(&viewborder);
        float viewborder_sizey = BLI_rctf_size_y(&viewborder);
        rect.xmin = floorf(viewborder.xmin + (scene->r.border.xmin * viewborder_sizex));
        rect.ymin = floorf(viewborder.ymin + (scene->r.border.ymin * viewborder_sizey));
        rect.xmax = floorf(viewborder.xmin + (scene->r.border.xmax * viewborder_sizex));
        rect.ymax = floorf(viewborder.ymin + (scene->r.border.ymax * viewborder_sizey));
      }
      else {
        rect.xmin = v3d->render_border.xmin * size[0];
        rect.ymin = v3d->render_border.ymin * size[1];
        rect.xmax = v3d->render_border.xmax * size[0];
        rect.ymax = v3d->render_border.ymax * size[1];
      }
    }
  }
#endif

  ved->instance->init();
  // TODO hack remove me
  ved->instance->depsgraph = depsgraph;
}

// This is passed JNPR_Data->instance directly
static void juniper_instance_free(void *instance)
{
    delete reinterpret_cast<juniper::JNPR *>(instance);
}

static void juniper_cache_init(void *vedata)
{
    JNPR_Data *ved = reinterpret_cast<JNPR_Data *>(vedata);
    ved->instance->begin_sync();
}

static void juniper_cache_populate(void *vedata, Object *ob)
{

    JNPR_Data *ved = reinterpret_cast<JNPR_Data *>(vedata);
    ved->instance->sync_.sync_object(ob);
}

static void juniper_cache_finish(void *vedata)
{
    JNPR_Data *ved = reinterpret_cast<JNPR_Data *>(vedata);
    ved->instance->lights.end_sync();
}

static void juniper_draw_scene(void *vedata)
{
  JNPR_Data *ved = reinterpret_cast<JNPR_Data *>(vedata);
  ved->instance->render();
}

static void juniper_render_to_image(void *vedata,
                                    struct RenderEngine *engine,
                                    struct RenderLayer *layer,
                                    const struct rcti *rect)
{
    /* Create new instance for rendering */
    juniper::JNPR* render_instance = new juniper::JNPR();
    // TODO: New render view from scene camera
//    draw::View render_view("Render View", DRW_view_default_get());
//    render_instance->composite(render_view);

    delete render_instance;
}

static void juniper_render_update_passes(RenderEngine *engine, Scene *scene, ViewLayer *view_layer)
{

  if (view_layer->passflag & SCE_PASS_COMBINED) {
    RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_COMBINED, 4, "RGBA", SOCK_RGBA);
  }
  if (view_layer->passflag & SCE_PASS_Z) {
    RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_Z, 1, "Z", SOCK_FLOAT);
  }
}

static void juniper_update_script_node(struct RenderEngine *engine,
                                     struct bNodeTree *ntree,
                                     struct bNode *node)
{
  if (node->id) {
    GPU_material_library_text_add((struct Text*) node->id, true);
  }
  update_node_declaration_and_sockets_p(ntree, node);
  Material* ma = (Material *) ntree->owner_id;

  if (ma) {
    ma->script_update_index++;
  }

  DEG_id_tag_update(&ma->id, ID_RECALC_SHADING | ID_RECALC_PARAMETERS);
  WM_main_add_notifier(NC_MATERIAL | ND_SHADING, ma);
}

extern "C" {

static const DrawEngineDataSize juniper_data_size = DRW_VIEWPORT_DATA_SIZE(JNPR_Data);

DrawEngineType draw_engine_juniper_type = {
    /*next*/ nullptr,
    /*prev*/ nullptr,
    /*idname*/ N_("Juniper"),
    /*vedata_size*/ &juniper_data_size,
    /*engine_init*/ &juniper_engine_init,
    /*engine_free*/ nullptr,
    /*instance_free*/ juniper_instance_free,
    /*cache_init*/ &juniper_cache_init,
    /*cache_populate*/ &juniper_cache_populate,
    /*cache_finish*/ &juniper_cache_finish,
    /*draw_scene*/ &juniper_draw_scene,
    /*view_update*/ nullptr,
    /*id_update*/ nullptr,
    /*render_to_image*/ &juniper_render_to_image,
    /*store_metadata*/ nullptr,
};

RenderEngineType DRW_engine_viewport_juniper_type = {
    /*next*/ nullptr,
    /*prev*/ nullptr,
    /*idname*/ JNPR_ENGINE,
    /*name*/ N_("Juniper"),
    /*flag*/ RE_INTERNAL | RE_USE_PREVIEW | RE_USE_GPU_CONTEXT,
    /*update*/ nullptr,
    /*render*/ &DRW_render_to_image,
    /*render_frame_finish*/ nullptr,
    /*draw*/ nullptr,
    /*bake*/ nullptr,
    /*view_update*/ nullptr,
    /*view_draw*/ nullptr,
    /*update_script_node*/ juniper_update_script_node,
    nullptr,  //&juniper_render_update_passes,
    &draw_engine_juniper_type,
    {nullptr, nullptr, nullptr},
};
} /* extern "C" */
