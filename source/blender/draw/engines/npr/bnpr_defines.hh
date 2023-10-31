/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * List of defines that are shared with the GPUShaderCreateInfos. We do this to avoid
 * dragging larger headers into the createInfo pipeline which would cause problems.
 */

#ifndef GPU_SHADER
#  pragma once
#endif

/* -------------------------------------------------------------------- */
/** Empty shader test
 * \{ */
#define GROUP_SIZE_STROKEGEN_TEST 512u
/** \} */



/* -------------------------------------------------------------------- */
/** Geometry Extraction from GPUBatch(es)
 * \{ */
#define GROUP_SIZE_STROKEGEN_GEOM_EXTRACT 256u
#define GROUP_SIZE_FILL_ARGS 32u
#define GROUP_SIZE_X_CONTOUR_PIXEL_COMPRESS 8u
#define GROUP_SIZE_Y_CONTOUR_PIXEL_COMPRESS 8u

#define MAX_GPU_HASH_TABLE_SIZE ((2048 * 2048 * 16)) /* we just cannot afford more than this */
 /** \} */



/* -------------------------------------------------------------------- */
/** \Scan Test
 * \{ */
#define GROUP_SIZE_BNPR_SCAN_TEST_SWEEP 1024u
#define GROUP_SIZE_BNPR_SCAN_TEST_AGGRG 1024u

#define NUM_ITEMS_BNPR_SCAN_TEST 1973581u

#ifndef GPU_SHADER
  #define BNPR_SCAN_TEST_DATA_TYPE uint3
  #define BNPR_SCAN_TEST_DATA_TYPE_STR "uvec3"

  // remember to update SSBOData_SegScanTest if this changes
  #define BNPR_SEG_SCAN_TEST_STRUCT_TYPE SSBOData_SegScanTest
  #define BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR "SSBOData_SegScanTest"
#endif
/** \} */



/* -------------------------------------------------------------------- */
/** \Segmented Convolution Test
 * \{ */
#define NUM_ITEMS_SEGLOOPCONV1D_TEST 1973581u
#define GROUP_SIZE_SEGLOOPCONV1D_TEST 512u

#define NPR_SEGLOOPCONV1D_CONV_RADIUS 32u
#ifndef GPU_SHADER
#  define NPR_SEGLOOPCONV1D_CONV_RADIUS_STR "32u"
#endif

#define NPR_SEGLOOPCONV1D_TEST_DATA_TYPE uint
#ifndef GPU_SHADER
#  define NPR_SEGLOOPCONV1D_TEST_DATA_TYPE_STR "uint"
#endif
/** \} */



/* -------------------------------------------------------------------- */
/** \List Ranking Test
 * \{ */
#define GROUP_SIZE_BNPR_LIST_RANK_TEST 1024
#define GROUP_SIZE_BNPR_LIST_RANK_TEST_STR "1024u"
#define NUM_ITEMS_BNPR_LIST_RANK_TEST 4194304u
/*4194304u*/ /*1973588u*/
#define MAX_NUM_JUMPS_BNPR_LIST_RANK_TEST 32u

#define NUM_ITERS_BNPR_LIST_RANK_SPLICE 3u
#define BNPR_LIST_RANK_CACHED_ANCHOR_COUNTERS ((1u + NUM_ITERS_BNPR_LIST_RANK_SPLICE))
/* must align with 4 bytes */
#define BNPR_LIST_RANK_ANCHOR_COUNTER_BUFFER_SIZE ((((BNPR_LIST_RANK_CACHED_ANCHOR_COUNTERS + 3u) >> 2u) << 2u))


#define NUM_ITERS_BNPR_LIST_RANK_RELINK NUM_ITERS_BNPR_LIST_RANK_SPLICE

/** \} */
