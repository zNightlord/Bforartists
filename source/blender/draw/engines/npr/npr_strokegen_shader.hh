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


namespace blender::npr::strokegen {

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

  CONV1D_TEST_BUILD_PATCH,
  CONV1D_TEST_CONVOLUTION,

  LISTRANKING_INIT_ANCHORS,
  LISTRANKING_COMPACT_ANCHORS,
  LISTRANKING_SPLICE_OUT_NODES, 
  LISTRANKING_FILL_DISPATCH_ARGS,
  LISTRANKING_SUBLIST_POINTER_JUMPING, 
  LISTRANKING_UPLOAD_CPU_DATA, 

  MAX_SHADER_TYPE,
};

/**
 * Shader module. shared between instances.
 */
class StrokeGenShaderModule {
 private:
  std::array<GPUShader *, MAX_SHADER_TYPE> shaders_;

  /** Shared shader module across all engine instances. */
  static StrokeGenShaderModule *g_shader_module;

 public:
  StrokeGenShaderModule();
  ~StrokeGenShaderModule();

  GPUShader *static_shader_get(eShaderType shader_type);
  // TODO: GPUMaterial ? (see impl in eevee)


  /** Only to be used by Instance constructor. */
  static StrokeGenShaderModule *module_get();
  static void module_free();

 private:
  const char *static_shader_create_info_name_get(eShaderType shader_type);
};

}  // namespace blender::bnpr
