/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Shader module that manage shader libraries, deferred compilation,
 * and static shader usage.
 */

#include "GPU_capabilities.h"

#include "gpu_shader_create_info.hh"

#include "bnpr_shader.hh"

namespace blender::bnpr {

/* -------------------------------------------------------------------- */
/** \name Module
 *
 * \{ */

ShaderModule *ShaderModule::g_shader_module = nullptr;

ShaderModule *ShaderModule::module_get()
{
  if (g_shader_module == nullptr) {
    /* TODO(@fclem) thread-safety. */
    g_shader_module = new ShaderModule();
  }
  return g_shader_module;
}

void ShaderModule::module_free()
{
  if (g_shader_module != nullptr) {
    /* TODO(@fclem) thread-safety. */
    delete g_shader_module;
    g_shader_module = nullptr;
  }
}

ShaderModule::ShaderModule()
{
  for (GPUShader *&shader : shaders_) {
    shader = nullptr;
  }

#ifdef DEBUG
  /* Ensure all shader are described. */
  for (auto i : IndexRange(MAX_SHADER_TYPE)) {
    const char *name = static_shader_create_info_name_get(eShaderType(i));
    if (name == nullptr) {
      std::cerr << "bnpr: Missing case for eShaderType(" << i
                << ") in static_shader_create_info_name_get().";
      BLI_assert(0);
    }
    const GPUShaderCreateInfo *create_info = GPU_shader_create_info_get(name);
    BLI_assert_msg(create_info != nullptr, "bnpr: Missing create info for static shader.");
  }
#endif
}

ShaderModule::~ShaderModule()
{
  for (GPUShader *&shader : shaders_) {
    DRW_SHADER_FREE_SAFE(shader);
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Static shaders
 *
 * \{ */

const char *ShaderModule::static_shader_create_info_name_get(eShaderType shader_type)
{ // returns the name of ShaderCreateInfo
  switch (shader_type) {
    case DEPTH:
      return "basic_depth_mesh";
    case POINTCLOUD_DEPTH:
      return "basic_depth_pointcloud";
    case CURVES_DEPTH:
      return "basic_depth_curves";
    case DEPTH_CONSERVATIVE:
      return "basic_depth_mesh_conservative";
    case POINTCLOUD_DEPTH_CONSERVATIVE:
      return "basic_depth_pointcloud_conservative";
    case COMPUTE_TEST:
      return "bnpr_strokegen_test_xxx";

    case COMPUTE_GEOM_EXTRACT:
      return "bnpr_geom_extract";
    
    case SCAN_TEST_UPSWEEP:
      return "bnpr_scan_test_upsweep";
    case SCAN_TEST_AGGREGATE:
      return "bnpr_scan_test_aggregate";
    case SCAN_TEST_DWSWEEP:
      return "bnpr_scan_test_dwsweep";
    case SEGSCAN_TEST_UPSWEEP:
      return "bnpr_segscan_test_upsweep";
    case SEGSCAN_TEST_AGGREGATE:
      return "bnpr_segscan_test_aggregate";
    case SEGSCAN_TEST_DWSWEEP:
      return "bnpr_segscan_test_dwsweep";

    /* To avoid compiler warning about missing case. */
    case MAX_SHADER_TYPE:
      return "";
  }
  return "";
}

GPUShader *ShaderModule::static_shader_get(eShaderType shader_type)
{
  if (shaders_[shader_type] == nullptr) {
    const char *shader_name = static_shader_create_info_name_get(shader_type);

    shaders_[shader_type] = GPU_shader_create_from_info_name(shader_name);

    if (shaders_[shader_type] == nullptr) {
      fprintf(stderr, "bnpr: error: Could not compile static shader \"%s\"\n", shader_name);
    }
    BLI_assert(shaders_[shader_type] != nullptr);
  }
  return shaders_[shader_type];
}


/** \} */

}  // namespace blender::bnpr
