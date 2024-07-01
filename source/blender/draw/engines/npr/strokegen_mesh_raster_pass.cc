#include "strokegen_mesh_raster_pass.hh"

#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_pass.hh"
#include "npr_strokegen_texture_pool.hh"

void npr::strokegen::StrokegenMeshRasterPass::init_pass(
    npr::strokegen::StrokeGenShaderModule& shader_module,
    npr::strokegen::GPUTexturePoolModule& texture_module,
    npr::strokegen::StrokegenMeshRasterPass::Usage usage
    ){
  PassMain::init(); 

  this->usage = usage; 

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

  DRWState drw_state = DRW_STATE_NO_DRAW;
  if (usage == Usage::DRAW_CONTOUR_EDGES || usage == Usage::DBG_LINES) {
    drw_state |= (DRW_STATE_WRITE_COLOR);
    drw_state |= ((draw_settings.draw_hidden_lines ? DRW_STATE_DEPTH_ALWAYS :
                                                     DRW_STATE_DEPTH_LESS_EQUAL)
                  /*| DRW_STATE_WRITE_DEPTH*/);  // lequal
    drw_state |= (DRW_STATE_STENCIL_ALWAYS);
    state_set(drw_state, clip_planes.size());  
  }
  if (usage == Usage::DRAW_CONTOUR_2D_SAMPLES) {
    drw_state |= (DRW_STATE_WRITE_COLOR);
    drw_state |= ((draw_settings.draw_hidden_lines ? DRW_STATE_DEPTH_ALWAYS :
                                                     DRW_STATE_DEPTH_LESS_EQUAL) |
                  DRW_STATE_WRITE_DEPTH);  // z-write, lequal
    drw_state |= (DRW_STATE_STENCIL_ALWAYS);
    state_set(drw_state, clip_planes.size());
  }
  if (usage == Usage::REMESHED_SURFACE_DEPTH) {
    drw_state |= (DRW_STATE_WRITE_COLOR);
    drw_state |= (DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_WRITE_DEPTH);  // z-write, lequal
    drw_state |= (DRW_STATE_STENCIL_ALWAYS);
    state_set(drw_state, clip_planes.size());  
  }
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_contour_subpass(
  StrokeGenShaderModule& shaders, GPUBufferPoolModule& buffers, GPUTexturePoolModule& textures)
{
  draw::PassMain::Sub* subpass = &sub("bnpr_geom_draw_contour_edges");
  subpass->shader_set(shaders.static_shader_get(eShaderType::INDIRECT_DRAW_CONTOUR_EDGES));

  subpass->bind_ssbo(0, buffers.reused_ssbo_bnpr_mesh_pool_());
  subpass->bind_ssbo(1, buffers.ssbo_bnpr_mesh_pool_counters_);
  subpass->bind_ssbo(2, buffers.ssbo_contour_snake_rank_);
  subpass->bind_ssbo(3, buffers.ssbo_contour_snake_list_len_);
  subpass->bind_ssbo(4, buffers.ssbo_contour_snake_list_head_);
  subpass->bind_ssbo(5, buffers.ssbo_contour_to_contour_);
  subpass->bind_ssbo(6, buffers.ssbo_contour_snake_seg_rank_);
  subpass->bind_ssbo(7, buffers.ssbo_contour_snake_seg_len_);
  subpass->bind_ssbo(8, buffers.ssbo_contour_snake_flags_);
  float2 fb_res = textures.get_contour_raster_screen_res();
  float2 fb_res_inv = float2(1.0f / fb_res.x, 1.0f / fb_res.y); 
  subpass->push_constant("pcs_screen_size_inv_", fb_res_inv);


  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE); 
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, buffers.ssbo_bnpr_contour_mesh_draw_args_); 
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_contour_2d_subpass(
    StrokeGenShaderModule &shaders, GPUBufferPoolModule &buffers, GPUTexturePoolModule &textures
)
{
  draw::PassMain::Sub *subpass = &sub("bnpr_geom_draw_contour_edges_2d");
  subpass->shader_set(shaders.static_shader_get(INDIRECT_DRAW_CONTOUR_2D_SAMPLES));

  subpass->bind_ssbo(0, buffers.ssbo_bnpr_mesh_pool_counters_);
  subpass->bind_ssbo(1, buffers.reused_ssbo_stroke_mesh_pool_());
  subpass->push_constant("pcs_screen_size_", textures.get_contour_raster_screen_res()); 

  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_TRIS, buffers.ssbo_bnpr_2d_sample_draw_args_);
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_dbg_lines_subpass(
    StrokeGenShaderModule & shaders,
    GPUBufferPoolModule & buffers,
    int line_type)
{
  draw::PassMain::Sub *subpass = &sub("bnpr_geom_draw_debug_lines");
  subpass->shader_set(shaders.static_shader_get(INDIRECT_DRAW_DBG_LINES));

  subpass->bind_ssbo(0, buffers.ssbo_dbg_lines_);
  subpass->bind_ssbo(1, buffers.ssbo_bnpr_mesh_pool_counters_);
  subpass->bind_ubo(0, buffers.ubo_view_matrices_);
  subpass->push_constant("pcs_line_type_", line_type); 

  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, buffers.ssbo_bnpr_vert_debug_draw_args_);
}

void npr::strokegen::StrokegenMeshRasterPass::append_draw_remeshed_surface_depth_subpass(
    npr::strokegen::StrokeGenShaderModule& shaders,
    npr::strokegen::GPUBufferPoolModule& buffers,
    npr::strokegen::GPUTexturePoolModule& textures)
{
  draw::PassMain::Sub *subpass = &sub("draw contour per obj z");

  subpass->shader_set(shaders.static_shader_get(INDIRECT_DRAW_CONTOUR_MESH_DEPTH));

  subpass->bind_ssbo(0, buffers.reused_ssbo_face_to_vert_draw_depth_());
  subpass->bind_ssbo(1, buffers.ssbo_vbo_full_);
  subpass->bind_ubo(0, buffers.ubo_view_matrices_cache_);

  subpass->barrier(GPU_BARRIER_COMMAND | GPU_BARRIER_SHADER_STORAGE);
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_TRIS, buffers.ssbo_bnpr_contour_mesh_draw_args_); 
}
