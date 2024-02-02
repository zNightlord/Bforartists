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
  StrokegenMeshRasterPass(const char *name) : draw::PassMain(name)
  {
  }

  struct DrawSettings {
    bool draw_hidden_lines; 
  } draw_settings;
  void init_pass(StrokeGenShaderModule& shader_module, GPUTexturePoolModule& texture_module);
  void append_draw_contour_subpass(StrokeGenShaderModule &shaders,
                                   GPUBufferPoolModule &buffers,
                                   GPUTexturePoolModule &textures);

  void append_draw_dbg_lines_subpass(npr::strokegen::StrokeGenShaderModule& shaders,
                                     npr::strokegen::GPUBufferPoolModule& buffers);
};

}

