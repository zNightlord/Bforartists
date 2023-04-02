/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 *
 * List of defines that are shared with the GPUShaderCreateInfos. We do this to avoid
 * dragging larger headers into the createInfo pipeline which would cause problems.
 */

#pragma once

/* -------------------------------------------------------------------- */
/** Geometry Extraction from GPUBatch(es)
 * \{ */
#define GROUP_SIZE_STROKEGEN_GEOM_EXTRACT 256u
/** \} */



/* -------------------------------------------------------------------- */
/** \Scan Test
 * \{ */
#define GROUP_SIZE_STROKEGEN_TEST 512u
#define GROUP_SIZE_BNPR_SCAN_TEST_SWEEP 1024u
#define GROUP_SIZE_BNPR_SCAN_TEST_AGGRG 1024u

#define NUM_ITEMS_BNPR_SCAN_TEST 1973581u

#ifndef GPU_SHADER
# define BNPR_SCAN_TEST_DATA_TYPE uint3
# define BNPR_SCAN_TEST_DATA_TYPE_STR "uvec3"

// remember to update SSBOData_SegScanTest if this changes
# define BNPR_SEG_SCAN_TEST_STRUCT_TYPE SSBOData_SegScanTest
# define BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR "SSBOData_SegScanTest"
#endif
/** \} */






