/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "bnpr_shader_shared.hh"

namespace blender::strokegen
{
  class StrokegenInstance;

  class GPUTexturePoolModule
  {
  private:
    /** Instance */
    StrokegenInstance &instance_;

    /** Compute Resources */
    TextureFromPool strokegen_tex_test_;


  public:
    GPUTexturePoolModule(StrokegenInstance &inst) :
      instance_(inst),
      strokegen_tex_test_("StrokegenTexture_Test"){};
    ~GPUTexturePoolModule() {};

    void sync(Object* object) {};
    void end_sync() {};

  };
}
