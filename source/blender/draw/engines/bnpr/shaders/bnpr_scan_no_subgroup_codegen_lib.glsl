#ifndef BNPR_SCAN_NO_SUBGROUP_CODEGEN_LIB
#define BNPR_SCAN_NO_SUBGROUP_CODEGEN_LIB

#pragma BLENDER_REQUIRE(bnpr_scan_no_subgroup_lib.glsl)

/* -------------------------------------------------------------------- */
/** \name Tree Scan LDS Cache
 * \{ */
DECLARE_TREE_SCAN_CACHE

#ifdef IS_TREE_SEG_SCAN

	DECLARE_TREE_SCAN_CACHE_HF

#endif
/** \} */



/* -------------------------------------------------------------------- */
/** \name Utility Functions
 * \{ */
DECLARE_TREE_SCAN_INDEXING_FUNCTION

DECLARE_TREE_SCAN_FUNC_CLEAN_SCAN_DATA

DECLARE_TREE_SCAN_FUNC_CLEAN_SEG_SCAN_DATA

/** \} */




#ifndef IS_TREE_SEG_SCAN
/* -------------------------------------------------------------------- */
/** \name Scan Functions
 * \{ */
DECLARE_TREE_SCAN_FUNC_BLOCK

DECLARE_TREE_SCAN_FUNC_AGGREGATE

/** \} */
#endif





#ifdef IS_TREE_SEG_SCAN
/* -------------------------------------------------------------------- */
/** \name Segmented Tree Scan Functions
 * \{ */

DECLARE_TREE_SEGSCAN_FUNC_DWSWEEP_FILL_CACHE

DECLARE_TREE_SEGSCAN_FUNC_AGGREGATE_FILL_CACHE

DECLARE_TREE_SEGSCAN_FUNC_UPSWEEP

DECLARE_TREE_SEGSCAN_FUNC_DWSWEEP

/** \} */
#endif





#endif
