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
#define GROUP_SIZE_STROKEGEN_GEOM_EXTRACT 512u
#define GROUP_SIZE_FILL_ARGS 32u
#define GROUP_SIZE_X_CONTOUR_PIXEL_COMPRESS 8u
#define GROUP_SIZE_Y_CONTOUR_PIXEL_COMPRESS 8u

#define MAX_GPU_HASH_TABLE_SIZE ((2048 * 2048 * 16)) /* we just cannot afford more than this */

#define MAX_NUM_STROKEGEN_OBJECTS ((65536u))

#define MAX_NUM_EDGES_PER_BATCH ((2048 * 2048 * 4))

/* https://math.stackexchange.com/questions/1541125/total-number-of-edges-in-a-triangle-mesh-with-n-vertices */
// have to make bigger since verts in pos_nor are not merged & vbo is ridiculously huge
// #define MAX_NUM_VERTS_PER_BATCH (((MAX_NUM_EDGES_PER_BATCH / 3) / 4) * 4)
#define MAX_NUM_VERTS_PER_BATCH (((MAX_NUM_EDGES_PER_BATCH / 2) / 4) * 4)

#define MAX_NUM_FACES_PER_BATCH (((((MAX_NUM_EDGES_PER_BATCH / 3) * 2) / 4) * 4)

#define MAX_NUM_CONTOUR_EDGES_PER_BATCH ((8 * 1024 * 1024)) /* about MAX_NUM_EDGES_PER_BATCH ^ 0.8f */

#define MAX_REMESHING_ITERS 4u /* max iters for the main remeshing loop */

/* max consecutive operations in one remesh iter */
#define MAX_CONSEQ_EDGE_SPLITS 8u 
#define MAX_CONSEQ_EDGE_COLLAPSES 8u 
#define MAX_CONSEQ_EDGE_FLIPS 8u
#define MAX_CONSEQ_FACE_SPLITS 8u
#define MAX_REMESH_OPS_PER_ITER ((MAX_CONSEQ_EDGE_SPLITS + MAX_CONSEQ_EDGE_COLLAPSES + MAX_CONSEQ_EDGE_FLIPS))

#define MAX_REMESHING_OPS (((MAX_REMESHING_ITERS) * (MAX_REMESH_OPS_PER_ITER)))
/** \} */

#define FRAME_COUNTER_CLAMP 100000000 // must be even number
#define MAX_TEMPORAL_FRAMES 2u // Match against define in shader
#define MAX_TEMPORAL_TRACKED_OBJECTS 64u // Match against define in shader

/* -------------------------------------------------------------------- */
/** \Scan Test
 * \{ */
/* Match against bnpr_shader_shared.hh --- */
#define GROUP_SIZE_BNPR_SCAN_SWEEP 1024u
#define GROUP_SIZE_BNPR_SCAN_AGGRG 1024u
// for recursive-scan should be as large as possible, or make each thread scan multiple items
/* --------------------------------------- */
#define NUM_ITEMS_BNPR_SCAN_TEST 1973581u

#ifndef GPU_SHADER
  #define BNPR_SCAN_TEST_DATA_TYPE float
#endif
/** \} */



/* -------------------------------------------------------------------- */
/** \Segmented Convolution Test
 * \{ */
#define NUM_ITEMS_SEGLOOPCONV1D_TEST 1973581u
#define GROUP_SIZE_SEGLOOPCONV1D_TEST 512u

#define NPR_TEST_SEGLOOPCONV1D_CONV_RADIUS 32u
#ifndef GPU_SHADER
#  define NPR_TEST_SEGLOOPCONV1D_CONV_RADIUS_STR "32u"
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
