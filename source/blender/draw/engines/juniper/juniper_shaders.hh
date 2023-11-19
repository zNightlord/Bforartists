/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "DRW_render.h"

#include "BLI_map.hh"
#include "BLI_vector.hh"
#include "GPU_material.h"
#include <array>
#include <string>

#include "BLI_string_ref.hh"
#include "DRW_render.h"
#include "GPU_shader.h"
#include "draw_manager.hh"
#include "draw_pass.hh"

namespace blender::juniper {

/* -------------------------------- */
//          STATIC SHADERS

enum eStaticShaderType {
    /* Builtin pipeline shaders */
    FILM_COMPOSITE,
    MAX_SHADER_TYPE,
};

/* -------------------------------- */
//          MATERIAL SHADERS

enum eMaterialPassType {
    MAT_PASS_PLACEHOLDER, // Tag for failed/compiling shaders
    MAT_PASS_PREPASS,
    MAT_PASS_SHADING,
};

struct MaterialPass {
    GPUMaterial *gpumat; // GPU shader
    draw::PassMain::Sub *sub_pass; // draw pass this material renders from
};

/** Shading pipeline material generated from shader node-tree */
struct NPRMaterial {
    bool is_alpha_blend; // Unused for now TODO alpha blend material pipeline
    MaterialPass pass_pre, pass_main, pass_depth; // Depth-only pass for Shadow
};

/** Unique shader variations should output a unique value (geometry type, pass type, etc.) */
// TODO add non-mesh variations
uint64_t options_from_shader_type(eMaterialPassType pass_type);

// Keep in sync with options_from_shader_type
eMaterialPassType shader_type_from_options(uint64_t options);

/** Unique key identifying a material */
struct MaterialKey {
    Material *mat; // Blender Material ID
    uint64_t options; // Unused, for specialisations (e.g. Mesh/Hair or Opaque/Blend )

    MaterialKey(Material *mat_) : mat(mat_)
    {
        options = 0;
    }

    uint64_t hash() const
    {
        // Hack - assume NPRMaterial size is less than total number of variations.
        BLI_assert(options < sizeof(*mat));
        return uint64_t(mat) + options;
    }

    /* Comparison operators for Map */
    bool operator<(const MaterialKey &k) const {
        return (mat < k.mat) || (options < k.options);
    }

    bool operator==(const MaterialKey &k) const {
        return (mat == k.mat) && (options == k.options);
    }
};

struct ShaderKey {
    GPUShader *sh;
    uint64_t options;

    ShaderKey(GPUMaterial *gpumat, eMaterialPassType pass_type)
    {
        sh = GPU_material_get_shader(gpumat);
        options = options_from_shader_type(pass_type);
    }

    uint64_t hash() const
    {
        return uint64_t(sh) + options;
    }

    bool operator<(const ShaderKey &k) const
    {
        return (sh == k.sh) ? (options < k.options) : (sh < k.sh);
    }

    bool operator==(const ShaderKey &k) const
    {
        return (sh == k.sh) && (options == k.options);
    }
};

class JNPR;

/** Manager class storing both static (pipeline) shaders and Material (node-tree) shaders */
class MaterialManager {
 public:
    int64_t compile_shader_queue = 0;
 private:
    JNPR &jnpr_;
    std::array<GPUShader *, MAX_SHADER_TYPE> static_shaders_;
    Material* default_surface_material_;

    // Cache for re-using Materials (e.g. with same Nodetree layout)
    Map<MaterialKey, NPRMaterial> material_map_;
    // Cache for re-using individual sub-pass/shaders (e.g. same pre-pass layout)
    Map<ShaderKey, draw::PassMain::Sub *> shader_map_;
    // Valid until next object_materials_get call.
    Vector<NPRMaterial *> object_materials;
    bNodeTree *default_surface_ntree_;
 public:
    MaterialManager(JNPR &jnpr);
    ~MaterialManager();

    void begin_sync();

    GPUShader *static_shader_get(eStaticShaderType type);
    GPUShader *static_shader_override_create(const char* shader_name);

    NPRMaterial &create_npr_material(Object *ob, Material *blender_mat);
    Vector<NPRMaterial *> &object_materials_get(Object* ob);

    MaterialPass material_pass_get(Object *ob, Material *blend_mat, eMaterialPassType type);

    void material_create_info_amend(GPUMaterial *pMaterial, GPUCodegenOutput *pOutput);

    GPUMaterial *material_shader_get(Material *blend_mat, bNodeTree *ntree, eMaterialPassType type, bool deferred);

};

} /* namespace blender::juniper */
