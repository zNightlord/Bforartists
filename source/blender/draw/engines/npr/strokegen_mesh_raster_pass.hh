#pragma once

#include "DNA_camera_types.h"
#include "DRW_render.h"
#include "draw_manager.hh"
#include "draw_pass.hh"
#include "bnpr_shader_shared.hh"
#include "gpu_shader_private.hh"
#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_shader.hh"
#include "npr_strokegen_texture_pool.hh"

using namespace blender;

namespace blender::npr::strokegen {

class StrokegenMeshRasterPass : public draw::PassMain {
 public:
  StrokegenMeshRasterPass(const char *name) : draw::PassMain(name)
  {
  }
  
  void init_pass(StrokeGenShaderModule& shader_module, GPUTexturePoolModule& texture_module);
  void append_draw_subpass(GPUBufferPoolModule& buffers);
};

}

