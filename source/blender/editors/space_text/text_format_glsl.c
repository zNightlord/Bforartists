/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup sptext
 */

#include <string.h>

#include "BLI_blenlib.h"

#include "DNA_space_types.h"
#include "DNA_text_types.h"

#include "BKE_text.h"

#include "text_format.h"

/* *** Local Functions (for format_line) *** */

static int txtfmt_glsl_find_builtin_func(const char *string)
{
  int i, len;

  /* Keep aligned args for readability. */
  /* clang-format off */

  /* list is from
   * https://registry.khronos.org/OpenGL/specs/gl/GLSLangSpec.4.60.pdf
   */
  if        (STR_LITERAL_STARTSWITH(string, "acosh",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "acos",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "asinh",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "asin",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atanh",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atan",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "cosh",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "cos",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "degrees",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "radians",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sin",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sinh",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "tan",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "tanh",         len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "abs",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "ceil",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "clamp",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dFdx",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dFdy",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "exp2",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "exp",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "floor",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "fma",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "fract",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "fwidth",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "inversesqrt",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isinf",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isnan",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "log2",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "log",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "max",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "min",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mix",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mod",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "modf",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "noise",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "pow",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "roundEven",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "round",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sign",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "smoothstep",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sqrt",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "step",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "trunc",        len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "floatBitsToInt",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "frexp",                 len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "intBitsToFloat",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "ldexp",                 len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "packDouble2x32",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "packHalf2x16",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "packUnorm",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "unpackDouble2x32",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "unpackHalf2x16",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "unpackUnorm",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "cross",                 len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "distance",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dot",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "equal",                 len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "faceforward",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "length",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "normalize",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "notEqual",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "reflect",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "refract",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "all",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "any",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "greaterThanEqual",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "greaterThan",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "lessThanEqual",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "lessThan",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "not",                   len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "EmitStreamVertex",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "EmitVertex",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "EndPrimitive",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "EndStreamPrimitive",    len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "interpolateAtCentroid", len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "interpolateAtOffset",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "interpolateAtSample",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "texelFetchOffset",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "texelFetch",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureGatherOffsets",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureGatherOffset",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureGather",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureGradOffset",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureGrad",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureLodOffset",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureLod",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureOffset",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProjGradOffset", len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProjGrad",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProjLodOffset",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProjOffset",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProjLod",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureProj",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureQueryLevels",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureQueryLod",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureSamples",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "textureSize",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "texture",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "determinant",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "groupMemoryBarrier",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "inverse",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "matrixCompMult",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "outerProduct",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "transpose",            len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "bitCount",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bitfieldExtract",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bitfieldInsert",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bitfieldReverse",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "findLSB",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "findMSB",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uaddCarry",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "umulExtended",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usubBorrow",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicAdd",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicAnd",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicCompSwap",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicExchange",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicMax",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicMin",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicOr",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageAtomicXor",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageLoad",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageSamples",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageSize",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageStore",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "atomicAdd",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicAnd",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicCompSwap",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicCounter",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicCounterDecrement",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicCounterIncrement",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicExchange",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicMax",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicMin",                   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicOr",                    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomicXor",                   len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "barrier",                     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "groupMemoryBarrier",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "memoryBarrierAtomicCounter",  len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "memoryBarrierBuffer",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "memoryBarrierImage",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "memoryBarrierShared",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "memoryBarrier",               len)) { i = len;

  } else                                                                         { i = 0;
  }

  /* clang-format on */

  /* If next source char is an identifier (eg. 'i' in "definite") no match */
  if (i == 0 || text_check_identifier(string[i])) {
    return -1;
  }

  /* If next source char is *not* a ( then no match */
  /* TODO: Check for potential whitespace between function name and parenthesis */
  if (string[i] != '(') {
    return -1;
  }
  return i;
}

static int txtfmt_glsl_find_reserved(const char *string)
{
  int i, len;

  /* Keep aligned args for readability. */
  /* clang-format off */

  /* list is from...
   * https://registry.khronos.org/OpenGL/specs/gl/GLSLangSpec.4.60.pdf
   */
  if        (STR_LITERAL_STARTSWITH(string, "const",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uniform",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "buffer",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "shared",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "attribute",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "varying",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "volatile",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "restrict",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "readonly",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "writeonly",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "atomic_uint",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "layout",        len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "centroid",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "flat",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "smooth",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "noperspective", len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "patch",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sample",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "invariant",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "precise",       len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "break",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "continue",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "do",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "for",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "while",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "switch",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "case",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "default",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "if",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "else",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "subroutine",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "discard",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "return",        len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "in",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "out",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "inout",         len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "lowp",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mediump",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "highp",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "precision",     len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "gl_CullDistance",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_FragCoord",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_FragDepth",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_FrontFacing",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_GlobalInvocationID",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_HelperInvocation",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_InstanceID",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_InvocationID",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_Layer",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_LocalInvocationID",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_LocalInvocationIndex", len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_NumSamples",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_NumWorkGroups",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_PatchVerticesIn",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_PointCoord",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_PointSize",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_Position",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_PrimitiveID",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_PrimitiveIDIn",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_SampleID",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_SampleMask",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_SampleMaskIn",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_SamplePosition",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_TessCoord",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_TessLevelInner",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_TessLevelOuter",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_VertexID",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_ViewportIndex",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_WorkGroupID",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "gl_WorkGroupSize",        len)) { i = len;
  } else                                                                     { i = 0;
  }

  /* clang-format on */

  /* If next source char is an identifier (eg. 'i' in "definite") no match */
  if (i == 0 || text_check_identifier(string[i])) {
    return -1;
  }
  return i;
}

/* Checks the specified source string for a GLSL special name. This name must
 * start at the beginning of the source string and must be followed by a non-
 * identifier (see text_check_identifier(char)) or null character.
 *
 * If a special name is found, the length of the matching name is returned.
 * Otherwise, -1 is returned. */

static int txtfmt_glsl_find_typename(const char *string)
{
  int i, len;

  /* Keep aligned args for readability. */
  /* clang-format off */

  /* OSL shader types */
  if        (STR_LITERAL_STARTSWITH(string, "int",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "void",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bool",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "true",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "false",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "float",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "double",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "struct",     len)) { i = len;

    /* Longer identifiers with postfixes must be checked first (to avoid false early matches) */
  } else if (STR_LITERAL_STARTSWITH(string, "mat2x2",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat3x2",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat4x2",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat2x3",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat3x3",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat4x3",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat2x4",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat3x4",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat4x4",     len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "dmat2x2",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat3x2",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat4x2",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat2x3",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat3x3",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat4x3",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat2x4",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat3x4",    len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dmat4x4",    len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "vec2",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "vec3",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "vec4",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "ivec2",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "ivec3",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "ivec4",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bvec2",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bvec3",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "bvec4",      len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "uint",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uvec2",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uvec3",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uvec4",      len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "dvec2",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dvec3",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "dvec4",      len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "mat2",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat3",       len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "mat4",       len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "sampler1DArrayShadow",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler1DShadow",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler1DArray",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler1D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler1DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler1D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler1DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler1D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DArrayShadow",     len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DShadow",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DArray",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler2D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler2DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler2DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler2D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler2D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DRectShadow",      len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DRect",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler2DRect",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler2DRect",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DMSArray",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "sampler2DMS",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler2DMSArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler2DMS",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler2DMSArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler2DMS",             len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "sampler3D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isampler3D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usampler3D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "samplerCubeArrayShadow",   len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "samplerCubeShadow",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "samplerCubeArray",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "samplerCube",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isamplerCubeArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isamplerCube",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usamplerCubeArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usamplerCube",             len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "samplerBuffer",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "isamplerBuffer",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "usamplerBuffer",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "image1DArray",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "image1D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage1DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage1D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage1DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage1D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "image2DArray",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "image2D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage2DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage2D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage2DArray",          len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage2D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "image2DRect",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage2DRect",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage2DRect",           len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "image2DMSArray",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "image2DMS",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage2DMSArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage2DMS",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage2DMSArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage2DMS",             len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "image3D",                len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimage3D",               len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimage3D",               len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "imageCubeArray",         len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "imageCube",              len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimageCubeArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimageCube",             len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimageCubeArray",        len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimageCube",             len)) { i = len;

  } else if (STR_LITERAL_STARTSWITH(string, "imageBuffer",            len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "iimageBuffer",           len)) { i = len;
  } else if (STR_LITERAL_STARTSWITH(string, "uimageBuffer",           len)) { i = len;

  } else                                                          { i = 0;
  }

  /* clang-format on */

  /* If next source char is an identifier (eg. 'i' in "definite") no match */
  if (i == 0 || text_check_identifier(string[i])) {
    return -1;
  }
  return i;
}

/* matches py 'txtfmt_osl_find_decorator' */
static int txtfmt_glsl_find_preprocessor(const char *string)
{
  if (string[0] == '#') {
    int i = 1;
    /* White-space is ok '#  foo'. */
    while (text_check_whitespace(string[i])) {
      i++;
    }
    while (text_check_identifier(string[i])) {
      i++;
    }
    return i;
  }
  return -1;
}

static char txtfmt_glsl_format_identifier(const char *str)
{
  char fmt;

  /* Keep aligned args for readability. */
  /* clang-format off */

  if        (txtfmt_glsl_find_typename(str) != -1)     { fmt = FMT_TYPE_SPECIAL;
  } else if (txtfmt_glsl_find_builtin_func(str) != -1) { fmt = FMT_TYPE_DIRECTIVE;
  } else if (txtfmt_glsl_find_reserved(str)     != -1) { fmt = FMT_TYPE_KEYWORD;
  } else if (txtfmt_glsl_find_preprocessor(str) != -1) { fmt = FMT_TYPE_RESERVED;
  } else                                              { fmt = FMT_TYPE_DEFAULT;
  }

  /* clang-format on */

  return fmt;
}

static void txtfmt_glsl_format_line(SpaceText *st, TextLine *line, const bool do_next)
{
  FlattenString fs;
  const char *str;
  char *fmt;
  char cont_orig, cont, find, prev = ' ';
  int len, i;

  /* Get continuation from previous line */
  if (line->prev && line->prev->format != NULL) {
    fmt = line->prev->format;
    cont = fmt[strlen(fmt) + 1]; /* Just after the null-terminator */
    BLI_assert((FMT_CONT_ALL & cont) == cont);
  }
  else {
    cont = FMT_CONT_NOP;
  }

  /* Get original continuation from this line */
  if (line->format != NULL) {
    fmt = line->format;
    cont_orig = fmt[strlen(fmt) + 1]; /* Just after the null-terminator */
    BLI_assert((FMT_CONT_ALL & cont_orig) == cont_orig);
  }
  else {
    cont_orig = 0xFF;
  }

  len = flatten_string(st, &fs, line->line);
  str = fs.buf;
  if (!text_check_format_len(line, len)) {
    flatten_string_free(&fs);
    return;
  }
  fmt = line->format;

  while (*str) {
    /* Handle continuations */
    if (cont) {
      /* C-Style comments */
      if (cont & FMT_CONT_COMMENT_C) {
        if (*str == '*' && *(str + 1) == '/') {
          *fmt = FMT_TYPE_COMMENT;
          fmt++;
          str++;
          *fmt = FMT_TYPE_COMMENT;
          cont = FMT_CONT_NOP;
        }
        else {
          *fmt = FMT_TYPE_COMMENT;
        }
      }
      str += BLI_str_utf8_size_safe(str) - 1;
    }
    /* Not in a comment already... */
    else {
      /* Deal with comments first */
      if (*str == '/' && *(str + 1) == '/') {
        /* fill the remaining line */
        text_format_fill(&str, &fmt, FMT_TYPE_COMMENT, len - (int)(fmt - line->format));
      }
        /* C-Style (multi-line) comments */
      else if (*str == '/' && *(str + 1) == '*') {
        cont = FMT_CONT_COMMENT_C;
        *fmt = FMT_TYPE_COMMENT;
        fmt++;
        str++;
        *fmt = FMT_TYPE_COMMENT;
      }
        /* White-space (all white-space has been converted to spaces). */
      else if (*str == ' ') {
        *fmt = FMT_TYPE_WHITESPACE;
      }
        /* Numbers (digits not part of an identifier and periods followed by digits) */
      else if ((prev != FMT_TYPE_DEFAULT && text_check_digit(*str)) ||
               (*str == '.' && text_check_digit(*(str + 1)))) {
        *fmt = FMT_TYPE_NUMERAL;
      }
        /* Punctuation */
      else if ((*str != '#') && text_check_delim(*str)) {
        *fmt = FMT_TYPE_SYMBOL;
      }
        /* Identifiers and other text (no previous white-space or delimiters. so text continues). */
      else if (prev == FMT_TYPE_DEFAULT) {
        str += BLI_str_utf8_size_safe(str) - 1;
        *fmt = FMT_TYPE_DEFAULT;
      }
      /* Not white-space, a digit, punctuation, or continuing text.
       * Must be new, check for special words. */
      else {
        /* Keep aligned arguments for readability. */
        /* clang-format off */

        /* Special vars(v) or built-in keywords(b) */
        /* keep in sync with `txtfmt_glsl_format_identifier()`. */
        if        ((i = txtfmt_glsl_find_typename(str))     != -1) { prev = FMT_TYPE_SPECIAL;
        } else if ((i = txtfmt_glsl_find_builtin_func(str)) != -1) { prev = FMT_TYPE_DIRECTIVE;
        } else if ((i = txtfmt_glsl_find_reserved(str))     != -1) { prev = FMT_TYPE_KEYWORD;
        } else if ((i = txtfmt_glsl_find_preprocessor(str)) != -1) { prev = FMT_TYPE_RESERVED;
        }

        /* clang-format on */
        if (i > 0) {
          if (prev == FMT_TYPE_DIRECTIVE) { /* can contain utf8 */
            text_format_fill(&str, &fmt, prev, i);
          }
          else {
            text_format_fill_ascii(&str, &fmt, prev, i);
          }
        }
        else {
          str += BLI_str_utf8_size_safe(str) - 1;
          *fmt = FMT_TYPE_DEFAULT;
        }
      }
    }
    prev = *fmt;
    fmt++;
    str++;
  }

  /* Terminate and add continuation char */
  *fmt = '\0';
  fmt++;
  *fmt = cont;

  /* If continuation has changed and we're allowed, process the next line */
  if (cont != cont_orig && do_next && line->next) {
    txtfmt_glsl_format_line(st, line->next, do_next);
  }

  flatten_string_free(&fs);
}

void ED_text_format_register_glsl(void)
{
  static TextFormatType tft = {NULL};
  static const char *ext[] = {"glsl", NULL};

  tft.format_identifier = txtfmt_glsl_format_identifier;
  tft.format_line = txtfmt_glsl_format_line;
  tft.ext = ext;

  ED_text_format_register(&tft);
}
