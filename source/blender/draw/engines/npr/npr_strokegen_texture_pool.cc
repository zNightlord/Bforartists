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

float2 GPUTexturePoolModule::get_contour_raster_screen_res() const
{
  float2 fb_res = float2((float)GPU_texture_width(tex_contour_raster),
                         (float)GPU_texture_height(tex_contour_raster));
  return fb_res; 
}

void GPUTexturePoolModule::on_begin_sync(Texture& tex_prepass_depth)
{
  // Query screen res. From blender::workbench::SceneState::init
  GPUTexture *viewport_tx = DRW_viewport_texture_list_get()->color;
  int2 screen_res = int2(GPU_texture_width(viewport_tx), GPU_texture_height(viewport_tx));

  eGPUTextureUsage usage_fb_tx = GPU_TEXTURE_USAGE_SHADER_READ | GPU_TEXTURE_USAGE_ATTACHMENT;

  // Acquire pooled texture. From blender::workbench::OpaquePass::draw
  {
    tex_contour_raster.acquire(screen_res, eGPUTextureFormat::GPU_RGBA32F, usage_fb_tx);

    int2 prez_res = int2(GPU_texture_width(tex_prepass_depth), GPU_texture_height(tex_prepass_depth)); 
    tex_contour_raster_depth.acquire(prez_res, eGPUTextureFormat::GPU_DEPTH24_STENCIL8, usage_fb_tx);

    GPUAttachment depth_att = GPU_ATTACHMENT_TEXTURE(tex_contour_raster_depth); 
    GPUAttachment col_att = GPU_ATTACHMENT_TEXTURE(tex_contour_raster);
    fb_contour_raster.ensure(depth_att, col_att);
  }

  {
    int2 perobjz_res = int2(GPU_texture_width(tex_prepass_depth), GPU_texture_height(tex_prepass_depth));
    tex_contour_per_obj_z_col.acquire(
        perobjz_res, eGPUTextureFormat::GPU_R32F, usage_fb_tx);
    tex_remeshed_surf_depth_.acquire(perobjz_res, eGPUTextureFormat::GPU_DEPTH24_STENCIL8, usage_fb_tx);

    GPUAttachment depth_att = GPU_ATTACHMENT_TEXTURE(/*tex_remeshed_surf_depth_*/ tex_contour_raster_depth);
    GPUAttachment col_att = GPU_ATTACHMENT_TEXTURE(tex_contour_per_obj_z_col);
    fb_remeshed_depth.ensure(depth_att, col_att);
  }

  {
    eGPUTextureUsage tex_pix_mark_usage = GPU_TEXTURE_USAGE_SHADER_READ | GPU_TEXTURE_USAGE_SHADER_WRITE; 
    int2 compressed_rect_size = CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE; 
    int2 pix_marks_res = (screen_res + compressed_rect_size - int2(1, 1)) / compressed_rect_size; 
    tex2d_contour_pix_marks_.acquire(pix_marks_res, eGPUTextureFormat::GPU_R32UI, tex_pix_mark_usage);
    tex2d_contour_dbg_.acquire(screen_res, eGPUTextureFormat::GPU_RGBA32F, tex_pix_mark_usage);

    GPUAttachment col_att = GPU_ATTACHMENT_TEXTURE(tex2d_contour_dbg_);
    fb_contour_dbg.ensure(GPU_ATTACHMENT_NONE, col_att);
  }

}

void GPUTexturePoolModule::sync(Object *object)
{

}

void GPUTexturePoolModule::on_finished_draw_viewport()
{
  tex_contour_raster.release();
  tex_contour_raster_depth.release();
  tex_contour_per_obj_z_col.release();
  tex_remeshed_surf_depth_.release(); 
  tex2d_contour_pix_marks_.release();
  tex2d_contour_dbg_.release(); 
}


}
