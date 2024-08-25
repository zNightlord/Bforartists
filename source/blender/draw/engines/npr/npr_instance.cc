#include "npr_instance.hh"

namespace blender::npr {

void Instance::init(Object *camera_ob)
{
  GPUTexture *viewport_tx = DRW_viewport_texture_list_get()->color;
  resolution_ = int2(GPU_texture_width(viewport_tx), GPU_texture_height(viewport_tx));

  if (prepass_sh_ == nullptr) {
    prepass_sh_ = GPU_shader_create_from_info_name("npr_prepass_mesh");
    deferred_pass_sh_ = GPU_shader_create_from_info_name("npr_deferred_pass");
  }

  if (strokegen_inst_ == nullptr) {
    strokegen_inst_ = new strokegen::Instance();

    Manager *drw_mgr = DRW_manager_get();

    const DRWContextState *ctx_state = DRW_context_state_get();
    View3D *v3d = ctx_state->v3d;
    RegionView3D *rv3d = ctx_state->rv3d;

    Object *camera = camera_ob;
    if (camera_ob == nullptr && v3d && rv3d && rv3d->persp == RV3D_CAMOB)
      camera = v3d->camera;  // use default cam if not specified

    const DRWView *default_drw_view = DRW_view_default_get();

    strokegen_inst_->init(
        ctx_state->depsgraph, drw_mgr, ctx_state->v3d, rv3d, default_drw_view, camera);
  }
}


void Instance::begin_sync()
{
  eGPUTextureUsage usage_fb_tx = GPU_TEXTURE_USAGE_SHADER_READ | GPU_TEXTURE_USAGE_ATTACHMENT;

  depth_tx_.ensure_2d(GPU_DEPTH24_STENCIL8, resolution_, usage_fb_tx);
  line_data_tx_.ensure_2d(GPU_R16F, resolution_, usage_fb_tx);
  line_color_tx_.ensure_2d(GPU_RGBA16F, resolution_, usage_fb_tx);

  prepass_ps_.init();
  prepass_ps_.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS);
  prepass_ps_.clear_color_depth_stencil(float4(0), 1.0f, 0);
  prepass_ps_.shader_set(prepass_sh_);

  deferred_pass_ps_.init();
  deferred_pass_ps_.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_BLEND_ALPHA);

  deferred_pass_ps_.shader_set(deferred_pass_sh_);
  deferred_pass_ps_.framebuffer_set(&deferred_pass_fb_);
  deferred_pass_ps_.bind_texture("depth_tx", &depth_tx_);
  deferred_pass_ps_.bind_texture("contour_tx", &(strokegen_inst_->strokegen_textures.tex_contour_raster));
  deferred_pass_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);

  strokegen_inst_->begin_sync(*DRW_manager_get(), depth_tx_);
}


void Instance::end_sync()
{
  strokegen_inst_->end_sync(*DRW_manager_get());
}


void Instance::object_sync(Manager &manager, ObjectRef &ob_ref)
{
  Object *ob = ob_ref.object;

  if (!DRW_object_is_renderable(ob)) {
    return;
  }
  if (!(DRW_object_visibility_in_active_context(ob) & OB_VISIBLE_SELF)) {
    return;
  }
  if ((ob->dt < OB_SOLID) && !DRW_state_is_scene_render()) {
    return;
  }

  if (ob->type == OB_MESH)
  { // Note: Manager::resource_handle() has side effect,
    // it will create a new resource every time gets called
    // so we need to cache the resource handle here
    ResourceHandle handle = manager.resource_handle(ob_ref);
    gpu::Batch *batch = DRW_cache_object_surface_get(ob_ref.object);

    strokegen_inst_->mesh_sync(manager, ob_ref, handle, &batch);

    prepass_ps_.draw(batch, handle);
  }
}

void Instance::draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
{
  view.sync(DRW_view_default_get());

  prepass_fb_.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx_)
                     // , GPU_ATTACHMENT_TEXTURE(id_tx_)
                     // , GPU_ATTACHMENT_TEXTURE(normal_tx_)
                     // , GPU_ATTACHMENT_TEXTURE(tangent_tx_)
  );
  prepass_fb_.bind();

  manager.submit(prepass_ps_, view);


  strokegen_inst_->draw_viewport(manager, view, depth_tx_);


  deferred_pass_fb_.ensure(GPU_ATTACHMENT_NONE,
                           GPU_ATTACHMENT_TEXTURE(color_tx),
                           GPU_ATTACHMENT_TEXTURE(line_color_tx_),
                           GPU_ATTACHMENT_TEXTURE(line_data_tx_));
  deferred_pass_fb_.bind();

  manager.submit(deferred_pass_ps_);

  GPU_texture_copy(depth_tx, depth_tx_);

  strokegen_inst_->end_draw_viewport(); 
}

void Instance::draw_viewport(Manager &manager,
                             View &view,
                             GPUTexture *depth_tx,
                             GPUTexture *color_tx)
{
  this->draw(manager, view, depth_tx, color_tx);
}


}  // namespace blender::npr
