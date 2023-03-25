/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * Shader module that manage shader libraries, deferred compilation,
 * and static shader usage.
 */

#pragma once

#include <array>
#include <string>

#include "BLI_string_ref.hh"
#include "DRW_render.h"
#include "GPU_material.h"
#include "GPU_shader.h"


namespace blender::bnpr {

/* Keep alphabetical order and clean prefix. */
enum eShaderType {
  COMPUTE_TEST = 0,
  DEPTH,
  POINTCLOUD_DEPTH,
  CURVES_DEPTH,
  DEPTH_CONSERVATIVE,
  POINTCLOUD_DEPTH_CONSERVATIVE,

  COMPUTE_GEOM_EXTRACT, 
 
  SCAN_TEST_UPSWEEP,
  SCAN_TEST_AGGREGATE,
  SCAN_TEST_DWSWEEP,

  SEGSCAN_TEST_UPSWEEP,
  SEGSCAN_TEST_AGGREGATE,
  SEGSCAN_TEST_DWSWEEP,

  MAX_SHADER_TYPE,
};

/**
 * Shader module. shared between instances.
 */
class ShaderModule {
 private:
  std::array<GPUShader *, MAX_SHADER_TYPE> shaders_;

  /** Shared shader module across all engine instances. */
  static ShaderModule *g_shader_module;

 public:
  ShaderModule();
  ~ShaderModule();

  GPUShader *static_shader_get(eShaderType shader_type);
  // TODO: GPUMaterial ? (see impl in eevee)


  /** Only to be used by Instance constructor. */
  static ShaderModule *module_get();
  static void module_free();

 private:
  const char *static_shader_create_info_name_get(eShaderType shader_type);
};

}  // namespace blender::bnpr
