/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#include "npr_strokegen_texture_pool.hh"
#include "npr_strokegen_instance.hh"

namespace blender::npr::strokegen
{

void GPUTexturePoolModule::on_begin_sync(Texture& tex_prepass_depth)
{
  // Query screen res. From blender::workbench::SceneState::init
  GPUTexture *viewport_tx = DRW_viewport_texture_list_get()->color;
  int2 screen_res = int2(GPU_texture_width(viewport_tx), GPU_texture_height(viewport_tx));

  // Acquire pooled texture. From blender::workbench::OpaquePass::draw
  eGPUTextureUsage usage_fb_tx = GPU_TEXTURE_USAGE_SHADER_READ | GPU_TEXTURE_USAGE_ATTACHMENT;
  tex_contour_raster.acquire(screen_res, eGPUTextureFormat::GPU_RGBA32F, usage_fb_tx);

  int2 prez_res = int2(GPU_texture_width(tex_prepass_depth), GPU_texture_height(tex_prepass_depth)); 
  tex_contour_raster_depth.acquire(prez_res, eGPUTextureFormat::GPU_DEPTH24_STENCIL8, usage_fb_tx);

  GPUAttachment depth_att = GPU_ATTACHMENT_TEXTURE(tex_contour_raster_depth); 
  GPUAttachment col_att = GPU_ATTACHMENT_TEXTURE(tex_contour_raster);

  fb_contour_raster.ensure(depth_att, col_att); 
}

void GPUTexturePoolModule::sync(Object *object)
{

}

void GPUTexturePoolModule::on_finished_draw_viewport()
{
  tex_contour_raster.release();
  tex_contour_raster_depth.release(); 
}


}
