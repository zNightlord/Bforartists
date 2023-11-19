/* SPDX-FileCopyrightText: 2005 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup gpu
 *
 * Parsing of and code generation using GLSL shaders in gpu/shaders/material. */

#pragma once

#include "GPU_material.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_FUNCTION_NAME 64
#define MAX_PARAMETER 36

struct GSet;

typedef enum {
  FUNCTION_QUAL_IN,
  FUNCTION_QUAL_OUT,
  FUNCTION_QUAL_INOUT,
} GPUFunctionQual;

// Additional data for GPUFunctions registered at runtime
typedef struct GPUFunctionRuntimeMeta {
    char paramname[MAX_PARAMETER][MAX_FUNCTION_NAME];
} GPUFunctionRuntimeMeta;

typedef struct GPUFunction {
  char name[MAX_FUNCTION_NAME];
  eGPUType paramtype[MAX_PARAMETER];
  GPUFunctionQual paramqual[MAX_PARAMETER];
  // Optional - used for user-defined functions
  GPUFunctionRuntimeMeta* runtimemeta;
  int totparam;
  /* TODO(@fclem): Clean that void pointer. */
  void *source; /* GPUSource */
} GPUFunction;

GPUFunction *gpu_material_library_use_function(struct GSet *used_libraries, const char *name);
GPUFunction *gpu_material_library_lookup_function(const char *name);

// Add and parse new code to material library
void gpu_material_library_runtime_add(const char* name, char *source);
// Remove runtime code from material library
void gpu_material_library_runtime_remove(const char* name);

#ifdef __cplusplus
}
#endif
