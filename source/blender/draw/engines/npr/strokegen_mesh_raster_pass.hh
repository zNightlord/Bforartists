#pragma once

#include "DNA_camera_types.h"
#include "DRW_render.h"
#include "draw_manager.hh"
#include "draw_pass.hh"
#include "bnpr_shader_shared.hh"

using namespace blender;

class StrokegenMeshRasterPass : draw::PassMain
{

public:
  struct IndirectDrawGfxResource
  {
    npr::strokegen::SSBO_StrokeGenMesh* ssbo_mesh;
    draw::StorageBuffer<DrawCommand, true> *ssbo_draw_args; 
  };

  void append_draw_subpass(IndirectDrawGfxResource &gfx_rsc); 

};
