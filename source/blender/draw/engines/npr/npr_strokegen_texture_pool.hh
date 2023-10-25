/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "bnpr_shader_shared.hh"

namespace blender::npr::strokegen
{
  class Instance;

  class GPUTexturePoolModule
  {
  private:
    /** Instance */
    Instance &instance_;

  public:
    /** Pooled Textures */
    TextureFromPool tex_contour_raster;
    TextureFromPool tex_contour_raster_depth; 

    /** Framebuffers */
    Framebuffer fb_contour_raster;
    float2 get_contour_raster_screen_res() const;

  public:
    GPUTexturePoolModule(Instance &inst) :
      instance_(inst),
      tex_contour_raster("contour edge gbuffer"){};
    ~GPUTexturePoolModule() {};

    void on_begin_sync(Texture& tex_prepass_depth); 
    void sync(Object* object);;
    void on_finished_draw_viewport();;

  };
}
