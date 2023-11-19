
#include "juniper_buffers.hh"
#include "juniper_instance.hh"

namespace blender::juniper {

void Buffers::acquire(int2 extent) {
  eGPUTextureFormat color_format = GPU_RGBA16F;
  eGPUTextureFormat normal_format = GPU_RGBA16F;

  // Render to double resolution
  // extent *= 2;
  depth_tx.acquire(extent, GPU_DEPTH24_STENCIL8);
  depth_copy_tx.acquire(extent, GPU_DEPTH24_STENCIL8);
  color_tx.acquire(extent, color_format);
  normal_tx.acquire(extent, normal_format);

  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();

  // For now we use default depth tex.
  prepass_fb_.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx), GPU_ATTACHMENT_TEXTURE(normal_tx));
  shading_fb_.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx), GPU_ATTACHMENT_TEXTURE(color_tx));
  zcopy_fb_.ensure(GPU_ATTACHMENT_TEXTURE(depth_copy_tx));
}

void Buffers::release() {
  depth_tx.release();
  depth_copy_tx.release();
  color_tx.release();
  normal_tx.release();
}

}

