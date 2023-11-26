/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "gpu_shader_create_info.hh"

#include "juniper_instance.hh"
#include "juniper_shaders.hh" // Own include
#include "BKE_node.h"

// For finding shader source overrides
#include "BKE_main.h"
#include "DNA_text_types.h"
#include "BKE_text.h"
#include "gpu_shader_create_info.hh"

#include "NOD_shader.h"
#include "intern/depsgraph.hh"

// For overriding builtin shader sources within blender
#include "gpu_shader_create_info_private.hh"

namespace blender::juniper {

using namespace blender;

MaterialManager::MaterialManager(JNPR &jnpr) : jnpr_(jnpr) {
    // Initialize static shaders to null.
    for (GPUShader *&shader : static_shaders_) {
        shader = nullptr;
    }

    // Initialise default nodetree to error value (default output color)
    // TODO: Material property based default surface
    default_surface_ntree_ = ntreeAddTree(nullptr, "Shader Nodetree", ntreeType_Shader->idname);
    bNode *output = nodeAddStaticNode(nullptr, default_surface_ntree_, SH_NODE_OUTPUT_MATERIAL);
    bNodeSocket *col_sock = static_cast<bNodeSocket *>(BLI_listbase_string_or_index_find(
            &output->inputs, nullptr, 0,4));
    float default_col[4] = {1.0, 0.0, 0.0, 1.0};
    copy_v4_v4((float*)col_sock->default_value, default_col);
    nodeSetActive(default_surface_ntree_, output);
}

MaterialManager::~MaterialManager() {
    // Free static shaders.
    for (GPUShader *&shader : static_shaders_) {
        DRW_SHADER_FREE_SAFE(shader);
    }

  ntreeFreeEmbeddedTree(default_surface_ntree_);
  MEM_SAFE_FREE(default_surface_ntree_);
}

void MaterialManager::begin_sync() {
    compile_shader_queue = 0;

    material_map_.clear();
    shader_map_.clear();
}

/* ------------------------------- */
//         STATIC SHADERS



// Checks if blendfile has a text block matching a source GLSL filename, and overrides the source.
GPUShader *MaterialManager::static_shader_override_create(const char* shader_name) {
  using namespace blender::gpu::shader;
  ShaderCreateInfo* cinfo_orig = (ShaderCreateInfo *) gpu_shader_create_info_get(shader_name);
  ShaderCreateInfo cinfo_new("static_shader_override");
  cinfo_new.additional_info(shader_name);

  // Overrides match original source file name, scan through texts to find those
  const char *frag_filename = cinfo_orig->fragment_source_.c_str();
  const char *vert_filename = cinfo_orig->vertex_source_.c_str();
  bool has_custom = false;

  if (jnpr_.depsgraph != nullptr) {
    deg::Depsgraph* dg = (deg::Depsgraph *)(jnpr_.depsgraph);
    Text* frag_source_txt = nullptr;
    Text* vert_source_txt = nullptr;

    LISTBASE_FOREACH(Text*, text, &dg->bmain->texts) {
      if (strcmp(text->id.name + 2, frag_filename) == 0) {
        frag_source_txt = text;
        has_custom = true;
      }
      if (strcmp(text->id.name + 2, vert_filename) == 0) {
        vert_source_txt = text;
        has_custom = true;
      }
    }

    // If we have overridden, keep the original source, but override main function
    // (BLENDER_REQUIRE needs original file intact for dependency linking)
    if (frag_source_txt) {
      cinfo_new.define("SRC_OVERRIDE_FRAG");
      size_t dummy;
      char* frag_src = txt_to_buf(frag_source_txt, &dummy);
      cinfo_new.fragment_source_generated = std::string(frag_src);
      MEM_SAFE_FREE(frag_src);
    }
    if (vert_source_txt) {
      cinfo_new.define("SRC_OVERRIDE_VERT");
      size_t dummy;
      char* vert_src = txt_to_buf(vert_source_txt, &dummy);
      cinfo_new.vertex_source_generated = std::string(vert_src);
      MEM_SAFE_FREE(vert_src);
    }
  }

  if (has_custom) {
    GPUShader* custom_shader = GPU_shader_create_from_info(reinterpret_cast<const GPUShaderCreateInfo *>(&cinfo_new));
    if (custom_shader) {
      return custom_shader;
    }
  }

  // Custom shader compilation failed, use builtin compile-time shader
  return GPU_shader_create_from_info_name(shader_name);
}

GPUShader *MaterialManager::static_shader_get(eStaticShaderType type) {
    if (static_shaders_[type] == nullptr) {
        const char *shader_name;
        switch (type) {
            case FILM_COMPOSITE:
                shader_name = "juniper_film_composite";
                break;
            case MAX_SHADER_TYPE:
                /* This should be an error */
                fprintf(stderr, "JNPR: error: Could not compile static shader.\"%s\"\n", shader_name);
                break;
        }

        // static_shaders_[type] = GPU_shader_create_from_info_name(shader_name);
      static_shaders_[type] = this->static_shader_override_create(shader_name);
    }


    return static_shaders_[type];
}


/* ----------------------------- */
//          MATERIALS

/** Unique shader variations should output a unique value (geometry type, pass type, etc.) */
// TODO add non-mesh variations
uint64_t options_from_shader_type(eMaterialPassType pass_type)
{
  return pass_type;
}

// Keep in sync with options_from_shader_type
eMaterialPassType shader_type_from_options(uint64_t options)
{
  // Lower bits of shader uuid (up to size of GPUShader struct)
  return static_cast<eMaterialPassType>(options & ((1u << 3u) - 1u));
}

void MaterialManager::material_create_info_amend(GPUMaterial *gpumat, GPUCodegenOutput *codegen_) {
    using namespace blender::gpu::shader;
    GPUCodegenOutput &codegen = *codegen_;

    ShaderCreateInfo &info = *reinterpret_cast<ShaderCreateInfo *>(codegen.create_info);

    eMaterialPassType pass_type = shader_type_from_options(GPU_material_uuid_get(gpumat));

    /* WORKAROUND: Replace by new ob info. */
    int64_t ob_info_index = info.additional_infos_.first_index_of_try("draw_object_infos");
    if (ob_info_index != -1) {
        info.additional_infos_[ob_info_index] = "draw_object_infos_new";
    }

    /* WORKAROUND: Add new ob attr buffer. */
    if (GPU_material_uniform_attributes(gpumat) != nullptr) {
        info.additional_info("draw_object_attribute_new");
    }

    std::stringstream attr_load;
    attr_load << "void attrib_load()\n";
    attr_load << "{\n";
    attr_load << ((codegen.attr_load) ? codegen.attr_load : "");
    attr_load << "}\n\n";

    {
        std::stringstream vert_gen, frag_gen;

        {
            vert_gen << "vec3 nodetree_displacement()\n";
            vert_gen << "{\n";
            vert_gen << ((codegen.displacement) ? codegen.displacement : "return vec3(0);\n");
            vert_gen << "}\n\n";
        }
        
        if (codegen.displacement) {
            /* Bump displacement. Needed to recompute normals after displacement. */
            info.define("MAT_DISPLACEMENT_BUMP");

            frag_gen << "vec3 nodetree_displacement()\n";
            frag_gen << "{\n";
            frag_gen << codegen.displacement;
            frag_gen << "}\n\n";
        }

        frag_gen << "vec4 nodetree_jnpr_color()\n";
        frag_gen << "{\n";
        frag_gen << ((codegen.jnpr_color) ? codegen.jnpr_color : "return vec4(0.0, 1.0, 0.0, 1.0);\n");
        frag_gen << "}\n\n";

        vert_gen << attr_load.str();
        info.vertex_source_generated = vert_gen.str();
        info.fragment_source_generated = frag_gen.str();
    }

    switch (pass_type) {
      case MAT_PASS_PREPASS:
        info.additional_info("juniper_surface_prepass");
        break;
      case MAT_PASS_SHADING:
        info.additional_info("juniper_surface_test");
        break;
    }

}

static void codegen_cb(void* thunk, GPUMaterial *mat, GPUCodegenOutput* codegen)
{
    reinterpret_cast<MaterialManager *>(thunk)->material_create_info_amend(mat, codegen);
}

GPUMaterial* MaterialManager::material_shader_get(Material *blend_mat, bNodeTree* ntree, eMaterialPassType type, bool deferred) {
    return DRW_shader_from_material(
            blend_mat,
            ntree,
            options_from_shader_type(type),
            false,
            deferred,
            codegen_cb,
            this);
}

MaterialPass MaterialManager::material_pass_get(Object *ob, Material *blend_mat, eMaterialPassType type)
{
    bNodeTree *ntree = blend_mat->nodetree ? blend_mat->nodetree : default_surface_ntree_;

    MaterialPass matpass = MaterialPass();

    // TODO shader variations
    matpass.gpumat = material_shader_get(blend_mat, ntree, type, true);

    switch (GPU_material_status(matpass.gpumat)) {
      case GPU_MAT_SUCCESS:
        break;
      case GPU_MAT_FAILED: // TODO add separate error shader
      case GPU_MAT_QUEUED:
      default:
        compile_shader_queue++;
        blend_mat = BKE_material_default_surface();
        matpass.gpumat = material_shader_get(blend_mat, default_surface_ntree_, type, false);
        break;
      // case GPU_MAT_FAILED:
      // default:
      //   break;
    }

    jnpr_.manager->register_layer_attributes(matpass.gpumat);

    ShaderKey shader_key(matpass.gpumat, type);

    draw::PassMain::Sub *shader_sub = shader_map_.lookup_or_add_cb(shader_key, [&]() {
        return jnpr_.pipelines.material_add(ob, blend_mat, matpass.gpumat, type);
    });

    if (shader_sub != nullptr) {
        matpass.sub_pass = &shader_sub->sub(GPU_material_get_name(matpass.gpumat));
        matpass.sub_pass->material_set(*jnpr_.manager, matpass.gpumat);
        // Bind light groups from material
        int4 light_groups = blend_mat->light_group_bits;
        matpass.sub_pass->push_constant("light_groups", light_groups);
    } else{
        matpass.sub_pass = nullptr;
    }

    return matpass;
}

NPRMaterial &MaterialManager::create_npr_material(Object *ob, Material *blender_mat)
{
    MaterialKey key(blender_mat);

    NPRMaterial &mat = material_map_.lookup_or_add_cb(key, [&]() {
        NPRMaterial mat = NPRMaterial();
        mat.pass_main = material_pass_get(ob, blender_mat, MAT_PASS_SHADING);
        mat.pass_pre = material_pass_get(ob, blender_mat, MAT_PASS_PREPASS);
        return mat;
    });

    return mat;
}

Vector<NPRMaterial *> &MaterialManager::object_materials_get(Object *ob) {
    object_materials.clear();

    const int num_materials = DRW_cache_object_material_count_get(ob);
    for (auto i : IndexRange(num_materials)) {
        // Material slot indices start at 1
        Material *blender_mat = BKE_object_material_get(ob, i + 1);
        if (blender_mat == nullptr) {
          blender_mat = BKE_material_default_surface();
        }

        NPRMaterial &mat = create_npr_material(ob, blender_mat);
        object_materials.append(&mat);
    }

    return object_materials;
}

}
