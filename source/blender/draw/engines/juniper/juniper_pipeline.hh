#pragma once

#include "DRW_render.h"
#include "draw_shader_shared.h"
#include "draw_pass.hh"
#include "juniper_shaders.hh"

namespace blender::juniper {
class JNPR;

class NPRPipeline {
private:
    JNPR &jnpr_;

    draw::PassMain npr_opaque_ps_ = {"Test draw pass"};
    draw::PassMain::Sub *npr_opaque_sub_ = nullptr;

    draw::PassMain npr_prepass_ps_ = {"Juniper Prepass"};
    draw::PassMain::Sub *npr_prepass_sub_ = nullptr;
public:
    NPRPipeline(JNPR &jnpr) : jnpr_(jnpr){};

    void sync();

    draw::PassMain::Sub *material_add(Material *blender_mat, GPUMaterial *gpumat, eMaterialPassType type) {
      switch (type) {
        case MAT_PASS_PREPASS:
          return &npr_prepass_sub_->sub(GPU_material_get_name(gpumat));
        case MAT_PASS_SHADING:
          return &npr_opaque_sub_->sub(GPU_material_get_name(gpumat));
        default:
          fprintf(stderr, "Invalid eMaterialPassType %d\n", type);
          return nullptr;
      }
    }
    void render(draw::View &view);
    void prepass(draw::View &view);

};

class Pipelines {
public:
    NPRPipeline test;

    void sync()
    {
        test.sync();
    }

    draw::PassMain::Sub *
    material_add(Object *ob, Material *blender_mat, GPUMaterial *gpumat, eMaterialPassType pass_type)
    {
        (void)ob; // Warning
        return test.material_add(blender_mat, gpumat, pass_type);
    }

    Pipelines(JNPR &jnpr) : test(jnpr){}
};

} // namespace blender::juniper
