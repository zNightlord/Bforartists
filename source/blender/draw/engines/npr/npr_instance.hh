#pragma once

#include "BKE_object.hh"
#include "BKE_report.hh"
#include "DEG_depsgraph_query.hh"
#include "DNA_camera_types.h"
#include "DRW_render.hh"
#include "draw_manager.hh"
#include "draw_pass.hh"

#include "npr_strokegen_instance.hh"

namespace blender::npr {

using namespace draw;

class Instance {
 public:
  View view = {"DefaultView"};
  int2 resolution_;

  PassMain prepass_ps_ = {"PrePass"};
  GPUShader *prepass_sh_ = nullptr;

  Texture id_tx_ = {"PrePass.ID"};
  // Texture normal_tx_ = {"PrePass.Normal"};
  // Texture tangent_tx_ = {"PrePass.Tangent"};
  Texture depth_tx_ = {"PrePass.Depth"};

  /* Line Width stored in red channel */
  Texture line_data_tx_ = {"PrePass.Line Data"};
  Texture line_color_tx_ = {"PrePass.Line Color"};

  Framebuffer prepass_fb_ = {"PrePass"};

  PassSimple deferred_pass_ps_ = {"Deferred Pass"};
  GPUShader *deferred_pass_sh_ = nullptr;

  Framebuffer deferred_pass_fb_ = {"Deferred Pass"};

  strokegen::Instance* strokegen_inst_; /* Strokegen Instance */

  ~Instance()
  {
    DRW_SHADER_FREE_SAFE(prepass_sh_);
    DRW_SHADER_FREE_SAFE(deferred_pass_sh_);

    delete strokegen_inst_;
  }

  
  void init(Object *camera_ob = nullptr);

  void begin_sync();

  void end_sync();

  void object_sync(Manager &manager, ObjectRef &ob_ref);

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx);

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx);
};

}  // namespace blender::npr
