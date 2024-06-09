#pragma once

#include "draw_manager.hh"
#include "draw_pass.hh"
#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_shader.hh"
#include "npr_strokegen_texture_pool.hh"

using namespace blender;

namespace blender::npr::strokegen {

struct SurfaceDebugContext; 


class StrokegenMeshRasterPass : public draw::PassMain {
 public:
  enum Usage {
    DRAW_CONTOUR_EDGES,
    DRAW_CONTOUR_2D_SAMPLES, 
    REMESHED_SURFACE_DEPTH,
    DBG_LINES,
  } usage;

  StrokegenMeshRasterPass(const char *name = "unknown strokegen render pass") : draw::PassMain(name)
  {
    usage = DRAW_CONTOUR_EDGES; // setup in init() func
  }

  StrokegenMeshRasterPass(const StrokegenMeshRasterPass& pass) : draw::PassMain(pass.debug_name)
  {
    // don't copy anything
  }

  

  struct DrawSettings {
    bool draw_hidden_lines; 
  } draw_settings;
  void init_pass(StrokeGenShaderModule& shader_module, GPUTexturePoolModule& texture_module, Usage usage);
  void append_draw_contour_subpass(StrokeGenShaderModule &shaders,
                                   GPUBufferPoolModule &buffers,
                                   GPUTexturePoolModule &textures);

  void append_draw_contour_2d_subpass(StrokeGenShaderModule &shaders,
                                      GPUBufferPoolModule &buffers,
                                      GPUTexturePoolModule &textures); 

  void append_draw_dbg_lines_subpass(StrokeGenShaderModule & shaders,
                                     GPUBufferPoolModule & buffers, int line_type);

  void append_draw_remeshed_surface_depth_subpass(
      npr::strokegen::StrokeGenShaderModule& shaders,
      npr::strokegen::GPUBufferPoolModule& buffers, npr::strokegen::GPUTexturePoolModule& textures);

  
};

}

