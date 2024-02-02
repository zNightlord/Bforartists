#include "strokegen_mesh_raster_pass.hh"

#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_pass.hh"
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
  drw_state |= (
    (draw_settings.draw_hidden_lines ? DRW_STATE_DEPTH_ALWAYS : DRW_STATE_DEPTH_LESS_EQUAL)
    | DRW_STATE_WRITE_DEPTH
  ); // z-write, lequal
  drw_state |= (DRW_STATE_STENCIL_ALWAYS); 
  state_set(drw_state, clip_planes.size());

  framebuffer_set(&texture_module.fb_contour_raster); 
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_contour_subpass(
  StrokeGenShaderModule& shaders, GPUBufferPoolModule& buffers, GPUTexturePoolModule& textures)
{
  draw::PassMain::Sub* subpass = &sub("strokegen raster pass");
  subpass->shader_set(shaders.static_shader_get(eShaderType::INDIRECT_DRAW_CONTOUR_EDGES));

  subpass->bind_ssbo(0, buffers.reused_ssbo_bnpr_mesh_pool_());
  subpass->bind_ssbo(1, buffers.ssbo_bnpr_mesh_pool_counters_);
  subpass->bind_ssbo(2, buffers.ssbo_contour_edge_rank_);
  subpass->bind_ssbo(3, buffers.ssbo_contour_edge_list_len_);
  subpass->bind_ssbo(4, buffers.ssbo_contour_edge_list_head_);
  float2 fb_res = textures.get_contour_raster_screen_res();
  float2 fb_res_inv = float2(1.0f / fb_res.x, 1.0f / fb_res.y); 
  subpass->push_constant("pcs_screen_size_inv_", fb_res_inv);

  // debug info
  subpass->bind_ssbo(5, buffers.ssbo_contour_to_contour_); 


  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE); 
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, buffers.ssbo_bnpr_mesh_pool_draw_args_); 
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_dbg_lines_subpass(
    npr::strokegen::StrokeGenShaderModule& shaders,
    npr::strokegen::GPUBufferPoolModule& buffers)
{
  draw::PassMain::Sub *subpass = &sub("draw debug lines");
  subpass->shader_set(shaders.static_shader_get(eShaderType::INDIRECT_DRAW_DBG_LINES));

  subpass->bind_ssbo(0, buffers.ssbo_dbg_lines_);
  subpass->bind_ubo(0, buffers.ubo_view_matrices_);

  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, buffers.ssbo_bnpr_vert_debug_draw_args_);
}
