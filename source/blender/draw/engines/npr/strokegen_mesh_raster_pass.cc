#include "strokegen_mesh_raster_pass.hh"

#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_texture_pool.hh"

void npr::strokegen::StrokegenMeshRasterPass::init_pass(
    npr::strokegen::StrokeGenShaderModule &shader_module,
    GPUTexturePoolModule &texture_module
){
  PassMain::init(); 

  // From "blender::workbench::SceneState::init"
  const DRWContextState *context = DRW_context_state_get();
  View3D *v3d = context->v3d;
  RegionView3D *rv3d = context->rv3d;

  Vector<float4> clip_planes; 
  clip_planes.clear();
  DRWState new_clip_state = RV3D_CLIPPING_ENABLED(v3d, rv3d) ? DRW_STATE_CLIP_PLANES :
                                                               DRW_STATE_NO_DRAW;
  if (new_clip_state & DRW_STATE_CLIP_PLANES) {
    int plane_len = (RV3D_LOCK_FLAGS(rv3d) & RV3D_BOXCLIP) ? 4 : 6;
    for (auto i : IndexRange(plane_len)) {
      clip_planes.append(rv3d->clip[i]);
    }
  }
  
  shader_set(shader_module.static_shader_get(eShaderType::INDIRECT_DRAW_CONTOUR_EDGES));

  DRWState drw_state = DRW_STATE_NO_DRAW;
  drw_state |= (DRW_STATE_WRITE_COLOR); 
  drw_state |= (DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_WRITE_DEPTH); // no z-write, lequal
  drw_state |= (DRW_STATE_STENCIL_ALWAYS); 
  state_set(drw_state, clip_planes.size());

  framebuffer_set(&texture_module.fb_contour_raster); 
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_subpass(GPUBufferPoolModule& buffers, GPUTexturePoolModule& textures)
{
  draw::PassMain::Sub* subpass = &sub("strokegen raster pass");
  subpass->bind_ssbo(0, buffers.ssbo_bnpr_mesh_pool_);
  subpass->bind_ssbo(1, buffers.ssbo_bnpr_mesh_pool_counters_);
  float2 fb_res = textures.get_contour_raster_screen_res();
  float2 fb_res_inv = float2(1.0f / fb_res.x, 1.0f / fb_res.y); 
  subpass->push_constant("pcs_screen_size_inv_", fb_res_inv); 
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, buffers.ssbo_bnpr_mesh_pool_draw_args_); 
}
