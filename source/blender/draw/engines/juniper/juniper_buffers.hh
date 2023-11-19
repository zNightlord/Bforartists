
#pragma once

#include "DRW_render.h"
#include "GPU_common_types.h"
#include "GPU_texture.h"
#include "juniper_defines.hh"
#include "DRW_gpu_wrapper.hh"

namespace blender::juniper {

class JNPR;

class Buffers {
public:
    draw::TextureFromPool depth_tx;
    draw::TextureFromPool depth_copy_tx;
    draw::TextureFromPool color_tx;
    draw::TextureFromPool normal_tx;

    draw::Framebuffer prepass_fb_ = {"Prepass"};
    draw::Framebuffer zcopy_fb_ = {"Z-Copy"};
    draw::Framebuffer shading_fb_;
private:
    JNPR &jnpr_;

public:
    Buffers(JNPR &jnpr) : jnpr_(jnpr){};

    void acquire(int2 extent);
    void release();
};

}