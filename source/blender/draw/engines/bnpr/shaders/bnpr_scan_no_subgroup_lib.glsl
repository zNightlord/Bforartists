#pragma BLENDER_REQUIRE(bnpr_hlsl_support_lib.glsl)

#ifndef BNPR_SCAN_NO_SUBGROUP_INCLUDED
#define BNPR_SCAN_NO_SUBGROUP_INCLUDED

/* Scan Operators */
uint u32_add(uint a, uint b)
{
	return a + b;
}
uvec2 uvec2_add(uvec2 a, uvec2 b)
{
	return a + b;
}
uvec3 uvec3_add(uvec3 a, uvec3 b)
{
	return a + b;
}
uvec4 uvec4_add(uvec4 a, uvec4 b)
{
	return a + b;
}
float f32_add(float a, float b)
{
	return a + b;
}





/* Inputs ---------------------                                                       */
/* Example:                                                                           */
/* (- basic -)                                                                        */
/* #define SCAN_DATA_TYPE uvec2			    // Input data type for scan operation        */
/* #define SCAN_OP(a, b) (a + b)				  // Scan operation                            */
/* #define SCAN_ZERO_VAL                // Zero value for scan operator               */
/* #define SCAN_BLOCK_SIZE 256			      // Typically thread_group_size               */
/* #define SCAN_FUNCTION_TAG ScanTest   // Alias name for this set of scan ops        */
/* ---------------------------------------                                            */


/* Macro expansion, for details, see
/* ---------------------------------------
/* https://stackoverflow.com/questions/1489932/how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg */
#define CAT(x, y) CAT_(x, y)
#define CAT_(x, y) x ## y


/* Type & Type conversion & Scan OP
/* --------------------------------------- */
#define T SCAN_DATA_TYPE
#define tag SCAN_FUNCTION_TAG
#define OP SCAN_OP


/* thread group size provided in glsl
/* --------------------------------------- */
#ifdef SCAN_BLOCK_SIZE
# undef SCAN_BLOCK_SIZE
#endif
#define SCAN_BLOCK_SIZE ((gl_WorkGroupSize.x))

#define DATA_SIZE       ((2u * SCAN_BLOCK_SIZE))


/* Padding Macros for Eliminating Bank Conflicts
/* --------------------------------------------------------- */
#define NUM_BANKS       32
#define LOG_NUM_BANKS   5
#define OFFSET_BANK_CONFLICT_FREE(x) (((x) >> LOG_NUM_BANKS))


// Tree Scan LDS Caches
// ------------------------------------------------------------
#define TREE_SCAN_CACHE CAT(TreeScanCache, tag)
#define TREE_SCAN_CACHE_SIZE (DATA_SIZE + DATA_SIZE / NUM_BANKS)
#define TREE_SCAN_CACHE_HF CAT(TreeScanCacheHF, tag)


#define DECLARE_TREE_SCAN_CACHE \
	shared T TREE_SCAN_CACHE[TREE_SCAN_CACHE_SIZE]; \

#define DECLARE_TREE_SCAN_CACHE_HF \
	shared uint TREE_SCAN_CACHE_HF[TREE_SCAN_CACHE_SIZE]; \



// Tree Scan
// ------------------------------------------------------------
struct TreeScanIndices
{
  uvec2 global_x2; /* index of element in global compute buffer */
  uvec2 lds_x2;    /* index of element 0/1 in shared memory */
};

/**
 * \brief Returns (Global_scanAddr0, Global_scanAddr1, LDS_scanAddr0, LDS_scanAddr1)
 *			LDS_scanAddr0/1 : index of element 0/1 in shared memory
 *			Global_scanAddr0/1 : index of element in global compute buffer
 */
#define DECLARE_TREE_SCAN_INDEXING_FUNCTION                                                             \
TreeScanIndices GetTreeScanIndices(					                                                        \
	uint groupIdx, uint gIdx					                                                              \
){  																				                                         \
	const uint groupOffset = (DATA_SIZE) * gIdx.x;									                             \
                                                                                                        \
	uint ai = groupIdx;                /*   0   1   2   3 ... 255  => ai	*/	                          \
	/* ------ + 1 * 512 ------- (Suppose gIdx.x == 1)						*/		                          \
	uint scanAddrA = groupOffset + ai; /* 512 513 514 515 ... 767  => scanAddrA	*/	                    \
                                                                                                        \
	uint bi = ai + DATA_SIZE / 2; /* 256 257 258 259 ... 511   => bi			*/	                          \
	uint scanAddrB = groupOffset + bi; /* 768 641 642 643 ... 1151  => scanAddrB*/	                    \
                                                                                                        \
	return TreeScanIndices(uvec2(scanAddrA, scanAddrB), uvec2(ai, bi));						                 \
}                                                                                                       \


#define _FUNC_CLEAN_SCAN_DATA CAT(TreeScanCleanData, tag)
/**
 * \brief Clear item that not mapped to actual scanned data
 */
#define DECLARE_TREE_SCAN_FUNC_CLEAN_SCAN_DATA \
void _FUNC_CLEAN_SCAN_DATA( \
    TreeScanIndices scan_ids,       \
    uint num_scanned_items,         \
    inout T scan_data_A,            \
    inout T scan_data_B             \
){ \
  if (scan_ids.global_x2.x >= num_scanned_items)        \
    scan_data_A = SCAN_ZERO_VAL;                        \
  if (scan_ids.global_x2.y >= num_scanned_items)        \
    scan_data_B = SCAN_ZERO_VAL;                        \
} \


#define _FUNC_CLEAN_SEG_SCAN_DATA CAT(TreeScanCleanData, tag)
/**
 * \brief Clear item that not mapped to actual scanned data
 */
#define DECLARE_TREE_SCAN_FUNC_CLEAN_SEG_SCAN_DATA \
void _FUNC_CLEAN_SEG_SCAN_DATA( \
	TreeScanIndices scan_ids,       \
	uint num_scanned_items,         \
	inout uint hf_A,                \
	inout T scan_data_A,            \
	inout uint hf_B,                \
	inout T scan_data_B             \
){ \
	if (scan_ids.global_x2.x >= num_scanned_items)        \
	{ \
		scan_data_A = SCAN_ZERO_VAL;                       \
		hf_A = 1u;                                         \
	} \
	if (scan_ids.global_x2.y >= num_scanned_items)        \
	{ \
		scan_data_B = SCAN_ZERO_VAL;                       \
		hf_B = 1u;                                         \
	} \
} \


#define _FUNC_TREE_SCAN_BLOCK CAT(TreeScanBlockExc, tag)
/**
 * \brief Declares a block-wise exclusive tree scan function.
 */
#define DECLARE_TREE_SCAN_FUNC_BLOCK \
void _FUNC_TREE_SCAN_BLOCK( \
    uint groupIdx,                    \
    uint gIdx,                        \
    T initialData_A,                  \
    T initialData_B,                  \
    out T scanRes_A,                  \
    out T scanRes_B                   \
){ \
  TreeScanIndices scanAddrs = GetTreeScanIndices(groupIdx, gIdx);						           \
  uint ai = scanAddrs.lds_x2.x;                                                        \
  uint bi = scanAddrs.lds_x2.y;                                                        \
  /* Bank Offset == index >> bits_banks(5 in Nvidia card) */                         \
  uint aiOffset = OFFSET_BANK_CONFLICT_FREE(ai);                                       \
  uint biOffset = OFFSET_BANK_CONFLICT_FREE(bi);                                       \
                                                                                       \
  /*  Store data into LDS with memory bank offset                                      \
  ---------------------------------------------------------------------                \
  about 'tailvalue':                                                                   \
  in prefix sum, last elem is going to be erased                                       \
  but we will need it later, so cache it here */                                       \
  TREE_SCAN_CACHE[ai + aiOffset] = initialData_A;                                      \
  TREE_SCAN_CACHE[bi + biOffset] = initialData_B;                                      \
  /* about LDS memory layout:                                                          \
  Interleaved storage,                                                                 \
  that is, ith(i % 32 == 0) is not used;                                               \
  e.g:                                                                                 \
  [0, 31]  X [32, 63] X  [64, 95]  X [96, 127]  -- Input CBuffer                       \
      + 0________+1___________+2___________+3 ... -- + OFFSET_BANK...(x)               \
  [0, 31] 32 [33, 64] 65 [66, 97] 98 [99, 130]  -- TREE_SCAN_CACHE                   */\
  \
  \
  \
  /* //////////////////////////////////////////////////////////////////////// */       \
  /* Scan --- Phase II        Up-Sweeping                                     */       \
  /* Work Indices:                                                            */       \
  /* offset = 2^k                                                             */       \
  /* a(i, k) = (2^k) * (2i + 1) - 1 = (2*gidx)*offset + offset - 1            */       \
  /* b(i, k) = a(i, k) + 2^k = a(i, k) + offset                               */       \
  /* i ~ groupIdx, k ~ iteration, all start from 0.                           */       \
  uint offset = 1;     /* Step Length == 2^k */                                        \
  uint d = DATA_SIZE / 2u; /* [0, ... , d]th threads are dispatched */           \
  for (; d > 0; d >>= 1){                                                       \
      barrier();                                                                \
      if (groupIdx < d){                                                        \
          ai = offset * (2 * groupIdx + 1) - 1;                                 \
          bi = offset * (2 * groupIdx + 2) - 1;                                 \
          ai += OFFSET_BANK_CONFLICT_FREE(ai);                                  \
          bi += OFFSET_BANK_CONFLICT_FREE(bi);                                  \
                                                                                \
          TREE_SCAN_CACHE[bi] = OP(TREE_SCAN_CACHE[ai], TREE_SCAN_CACHE[bi]);		\
      }                                                                         \
      offset *= 2;                                                              \
  }                                                                             \
  \
  \
  \
  /* ////////////////////////////////////////////////////////////////////////*/ \
  /* Phase III */                                                               \
  if (groupIdx == 0)                                                            \
  {                                                                             \
      /* Zero out last elem, prepare for up-sweeping */                         \
      uint lastIndex = DATA_SIZE - 1 + OFFSET_BANK_CONFLICT_FREE(DATA_SIZE - 1);\
      TREE_SCAN_CACHE[lastIndex] = SCAN_ZERO_VAL;                               \
  }                                                                             \
  \
  \
  \
  /* ///////////////////////////////////////////////////////////////////////// */ \
  /* Phase IV                 Down-Sweeping                                    */ \
  /* Util this point,                                                          */ \
  /* d == 0,                                                                   */ \
  /* offset == GROUP_SIZE * 2 == DATA_SIZE                                     */ \
  /* This is actually "rolling back + mirror" version of Phase I,              */ \
  /* So this execution code is a mirrored loop                                 */ \
  for (d = 1; d < DATA_SIZE; d *= 2){                                              \
      offset >>= 1;                                                                \
      barrier();                                           \
      if (groupIdx < d){                                                           \
          /* So the indexing function is the same, (rolling back)                  \
          just the roles of ai & bi are switched                              */   \
          ai = offset * (2 * groupIdx + 1) - 1;                                    \
          bi = offset * (2 * groupIdx + 2) - 1;                                    \
          ai += OFFSET_BANK_CONFLICT_FREE(ai);                                     \
          bi += OFFSET_BANK_CONFLICT_FREE(bi);                                     \
          /* swap */                                                               \
          T aiValOld = TREE_SCAN_CACHE[ai];                                        \
          TREE_SCAN_CACHE[ai] = TREE_SCAN_CACHE[bi];                               \
          TREE_SCAN_CACHE[bi] = OP(aiValOld, TREE_SCAN_CACHE[bi]);                 \
      }                                                                            \
  }                                                                                \
  barrier();                                               \
  \
  \
  \
  T pSumAtAi = TREE_SCAN_CACHE[groupIdx + aiOffset];                               \
  T pSumAtBi = TREE_SCAN_CACHE[groupIdx + SCAN_BLOCK_SIZE + biOffset];             \
  \
  \
  scanRes_A = pSumAtAi;                                     \
  scanRes_B = pSumAtBi;                                     \
} \


#define _FUNC_TREE_SCAN_AGGREGATE CAT(TreeScanBlockAggregate, tag)
/**
 \brief Second step for tree scan.
 *  Apply exclusive scan on
 *  inclusive sums from each scanned data-block.
*/
#define DECLARE_TREE_SCAN_FUNC_AGGREGATE \
void _FUNC_TREE_SCAN_AGGREGATE( \
  uint groupIdx,        \
  uint gIdx,            \
  T aggregateA,         \
  T aggregateB,         \
  out T aggSumA,        \
  out T aggSumB         \
) \
{				                    \
  _FUNC_TREE_SCAN_BLOCK	        \
  (			                    \
    groupIdx,                   \
    gIdx,                       \
    aggregateA,                 \
    aggregateB,                 \
    aggSumA, /*out*/            \
    aggSumB  /*out*/            \
  );			                    \
} \






/* -------------------------------------------------------------------- */
/** \name Segmented Tree Scan
 *
 * \{ */
#define _FUNC_TREE_SEG_SCAN_DWSWEEP_FILL_CACHE CAT(TreeSegScanExc_DwSweep_FillLDS_, tag)

#define DECLARE_TREE_SEGSCAN_FUNC_DWSWEEP_FILL_CACHE	\
void _FUNC_TREE_SEG_SCAN_DWSWEEP_FILL_CACHE(	\
	uint groupId, 					\
	TreeScanIndices scanAddrs, 	  \
	/* --- LDS inputs --- */ 	     \
	uint encodedHFs_A,	           \
	T partialSum_A,	              \
	uint encodedHFs_B,	           \
	T partialSum_B,	              \
	T aggregateScanRes 				  \
) {																																\
	uint ai = scanAddrs.lds_x2.x;																			            \
	uint bi = scanAddrs.lds_x2.y;																			            \
	\
	/* Bank Offset == index >> bits_banks(5 in Nvidia card) */										      \
	uint aiOffset = OFFSET_BANK_CONFLICT_FREE(ai);													            \
	uint biOffset = OFFSET_BANK_CONFLICT_FREE(bi);													            \
	\
	/*  Store data into LDS with memory bank offset								 	\
	--------------------------------------------------------------------- 	\
	about 'tailvalue':																		   \
	in prefix sum, last elem is going to be erased								   \
	but we will need it later, so cache it here                          */	\
	uint cacheAddrAi = ai + aiOffset;																              \
	uint cacheAddrBi = bi + biOffset;																              \
	TREE_SCAN_CACHE[cacheAddrAi] = partialSum_A;													              \
	TREE_SCAN_CACHE_HF[cacheAddrAi] = encodedHFs_A;													           \
	/* Different from normal down-sweep that zeros out last elem, */								     \
	/* We use output from prev inter-block scan kernel instead */									     \
	TREE_SCAN_CACHE[cacheAddrBi] =                                                                 \
		(groupId == SCAN_BLOCK_SIZE - 1)                                                            \
			? aggregateScanRes : partialSum_B;	                                                     \
	TREE_SCAN_CACHE_HF[cacheAddrBi] = encodedHFs_B;													           \
} \


#define _FUNC_TREE_SEG_SCAN_AGGREGATE_FILL_CACHE CAT(TreeSegScanExc_Aggregate_FillLDS_, tag)

#define DECLARE_TREE_SEGSCAN_FUNC_AGGREGATE_FILL_CACHE \
void _FUNC_TREE_SEG_SCAN_AGGREGATE_FILL_CACHE( \
	uint groupId, 					\
	TreeScanIndices scanAddrs, \
	/* --- LDS inputs --- */ 	\
	uint firstInitialHFAi,	   \
	uint firstInitialHFBi	   \
){ \
	uint ai = scanAddrs.lds_x2.x;																			            \
	uint bi = scanAddrs.lds_x2.y;																			            \
	\
	/* Bank Offset == index >> bits_banks(5 in Nvidia card) */										      \
	uint aiOffset = OFFSET_BANK_CONFLICT_FREE(ai);													            \
	uint biOffset = OFFSET_BANK_CONFLICT_FREE(bi);													            \
	uint cacheAddrAi = ai + aiOffset;																               \
	uint cacheAddrBi = bi + biOffset;																               \
	if (groupId == 0u)                                                                              \
	{                                                                                               \
		/* Zero out last elem, prepare for up-sweeping */                                          \
		uint lastIndex = DATA_SIZE - 1 + OFFSET_BANK_CONFLICT_FREE(DATA_SIZE - 1);                   \
		TREE_SCAN_CACHE[lastIndex] = SCAN_ZERO_VAL;                                                  \
	}                                                                                               \
	barrier();                                                                                      \
	                                                                                                \
	/* Compared to normal seg-scan,						  */                                            \
	/* need to encode original hfs differently here */                                            \
	TREE_SCAN_CACHE_HF[cacheAddrAi] = tree_seg_scan_encode_upsweep_hfs(                             \
		TREE_SCAN_CACHE_HF[cacheAddrAi], firstInitialHFAi                                            \
	);                                                                                              \
	TREE_SCAN_CACHE_HF[cacheAddrBi] = tree_seg_scan_encode_upsweep_hfs(                             \
		TREE_SCAN_CACHE_HF[cacheAddrBi], firstInitialHFBi                                            \
	);                                                                                              \
}


#define _FUNC_TREE_SEG_SCAN_UPSWEEP CAT(TreeSegScanExc_UpSweep_, tag)

#define DECLARE_TREE_SEGSCAN_FUNC_UPSWEEP	\
void _FUNC_TREE_SEG_SCAN_UPSWEEP(	\
	uint groupIdx,				        \
	uint gIdx,                      \
	/* --- scan inputs --- */ 	  \
	uint headFlagAi,	              \
	T initialDataAi,	              \
	uint headFlagBi,	              \
	T initialDataBi,	              \
	/* --- block partial sums --- */ \
 	out uint headFlagPartialSum_A,  \
	out T partialSum_A,		        \
	out uint headFlagPartialSum_B,  \
	out T partialSum_B		        \
) \
{ \
	/* -------------------------------------------------------	*/									      \
	/* nAddr:													*/									                     \
	/* .x: Global_scanAddr0, .y: Global_scanAddr1, 				*/									         \
	/* .z: LDS_scanAddr0, .w: LDS_scanAddr1						*/									            \
	TreeScanIndices scanAddrs = GetTreeScanIndices(groupIdx, gIdx);									      \
	uint ai = scanAddrs.lds_x2.x;																			            \
	uint bi = scanAddrs.lds_x2.y;																			            \
																																	\
	/* Bank Offset == index >> bits_banks(5 in Nvidia card) */										      \
	uint aiOffset = OFFSET_BANK_CONFLICT_FREE(ai);													            \
	uint biOffset = OFFSET_BANK_CONFLICT_FREE(bi);													            \
																																	\
	/*  Store data into LDS with memory bank offset								 \
	--------------------------------------------------------------------- \
	about 'tailvalue':																		\
	in prefix sum, last elem is going to be erased								\
	but we will need it later, so cache it here                          */	\
	uint cacheAddrAi = ai + aiOffset;																               \
	uint cacheAddrBi = bi + biOffset;																               \
	TREE_SCAN_CACHE[cacheAddrAi] = initialDataAi;													            \
	TREE_SCAN_CACHE_HF[cacheAddrAi] = headFlagAi;													            \
	TREE_SCAN_CACHE[cacheAddrBi] = initialDataBi;													            \
	TREE_SCAN_CACHE_HF[cacheAddrBi] = headFlagBi;													            \
	/* about LDS memory layout:																		               \
	Interleaved storage,																			                  \
	that is, ith(i % 32 == 0) is not used;															            \
	e.g:																							                        \
	[0, 31]  X [32, 63] X  [64, 95]  X [96, 127]  -- Input CBuffer									      \
		+ 0________+1___________+2___________+3 ... -- + OFFSET_BANK...(x)							   \
	[0, 31] 32 [33, 64] 65 [66, 97] 98 [99, 130]  -- TREE_SCAN_CACHE			*/					      \
                                                                                                	\
                                                                                                	\
	/* //////////////////////////////////////////////////////////////////////// */					\
	/* Scan --- Phase II        Up-Sweeping                                     */					\
	/* Work Indices:                                                            */					\
	/* offset = 2^k                                                             */					\
	/* a(i, k) = (2^k) * (2i + 1) - 1 = (2*gidx)*offset + offset - 1            */					\
	/* b(i, k) = a(i, k) + 2^k = a(i, k) + offset                               */					\
	/* i ~ groupIdx, k ~ iteration, all start from 0.                           */					\
	uint offset = 1u; /* Step Length == 2^k */														            \
	uint d = DATA_SIZE / 2u; /* [0, ... , d]th threads are dispatched */								   \
                                                                                                   \
	bool activeThread;																				                  \
	for (; d > 0; d >>= 1)																			                  \
	{																								                        \
		activeThread = groupIdx < d;																                  \
																																	\
		ai = offset * (2u * groupIdx + 1u) - 1u;														            \
		bi = offset * (2u * groupIdx + 2u) - 1u;														            \
		ai += OFFSET_BANK_CONFLICT_FREE(ai);														               \
		bi += OFFSET_BANK_CONFLICT_FREE(bi);														               \
																																	\
		barrier();															                                       \
		bool isSegHeadAtBi = (0u != (1u & TREE_SCAN_CACHE_HF[bi]));												\
		if (activeThread && (!isSegHeadAtBi))														               \
		{																							                        \
			TREE_SCAN_CACHE[bi] = OP(TREE_SCAN_CACHE[ai], TREE_SCAN_CACHE[bi]);						      \
		}																							                        \
		barrier();															                                       \
	                                                                                                \
		TREE_SCAN_CACHE_HF[bi] = activeThread																		   \
			? uint(isSegHeadAtBi || bool(TREE_SCAN_CACHE_HF[ai]))												   \
			: uint(isSegHeadAtBi);																		               \
                                                                                                   \
		offset *= 2u;																				                     \
	}																								                        \
                                                                                                   \
	barrier();																                                       \
                                                                                                   \
																																	\
	partialSum_A = TREE_SCAN_CACHE[cacheAddrAi];													               \
	partialSum_B = TREE_SCAN_CACHE[cacheAddrBi];													               \
	headFlagPartialSum_A = TREE_SCAN_CACHE_HF[cacheAddrAi];											         \
	headFlagPartialSum_B = TREE_SCAN_CACHE_HF[cacheAddrBi];											         \
} \


#define _FUNC_TREE_SEG_SCAN_DWSWEEP CAT(TreeSegScanExc_DwSweep_, tag)

#define DECLARE_TREE_SEGSCAN_FUNC_DWSWEEP	\
void _FUNC_TREE_SEG_SCAN_DWSWEEP( \
	uint groupIdx,	                  \
	uint gIdx,		                  \
	/* --- scan results --- */ \
	out T scanResult_A,	            \
	out T scanResult_B	            \
) \
{																									            \
	/* Addressing & Data Loading											*/					         \
	/* scanAddrs:																*/					         \
	/* -- .x: Global_scanAddr0, .y: Global_scanAddr1, 			*/				            \
	/* -- .z: LDS_scanAddr0,    .w: LDS_scanAddr1					*/					         \
	TreeScanIndices scanAddrs = GetTreeScanIndices(groupIdx, gIdx);                     \
	\
	/* Bank Offset == index >> bits_banks(5 in Nvidia card) */                        \
	uint ai = scanAddrs.lds_x2.x;																			\
	uint bi = scanAddrs.lds_x2.y;																			\
	uint aiOffset = OFFSET_BANK_CONFLICT_FREE(ai);													\
	uint biOffset = OFFSET_BANK_CONFLICT_FREE(bi);													\
	\
	/*  Store data into LDS with memory bank offset												\
	/* instead should call "_FUNC_TREE_SEG_SCAN_FILL_CACHE"                            \
	* ---------------------------------------------------------------------				\
	* about 'tailvalue':																				   \
	* in prefix sum, last elem is going to be erased												\
	* but we will need it later, so cache it here */					                     \
	uint cacheAddrAi = ai + aiOffset;																   \
	uint cacheAddrBi = bi + biOffset;																   \
	/**TREE_SCAN_CACHE[cacheAddrAi] = partialSum_A;		*/										\
	/**TREE_SCAN_CACHE_HF[cacheAddrAi] = encodedHFs_A;	*/										\
	/**TREE_SCAN_CACHE[cacheAddrBi] = partialSum_B;		*/										\
	/**TREE_SCAN_CACHE_HF[cacheAddrBi] = encodedHFs_B;	*/										\
	\
	\
	/* ///////////////////////////////////////////////////////////////////////// */	\
	/* Phase IV                 Down-Sweeping                                    */	\
	/* Util this point,                                                          */	\
	/* d == 0,                                                                   */	\
	/* offset == GROUP_SIZE * 2 == DATA_SIZE                                     */	\
	/* This is actually "rolling back + mirror" version of Phase I,              */	\
	/* So this execution code is a mirrored loop                                 */	\
	uint offset = DATA_SIZE;																		      \
	uint d = 0u;																						      \
	bool activeThread;																				      \
	for (d = 1u; d < DATA_SIZE; d *= 2u)																\
	{																								            \
		offset >>= 1;																				         \
		/* So the indexing function is the same, (rolling back)	                    \
		 * just the roles of ai & bi are switched */                                   \
		ai = offset * (2u * groupIdx + 1u) - 1u;														\
		bi = offset * (2u * groupIdx + 2u) - 1u;														\
		uint aiNext = ai + 1u + OFFSET_BANK_CONFLICT_FREE(ai + 1u);						      \
		ai += OFFSET_BANK_CONFLICT_FREE(ai);														   \
		bi += OFFSET_BANK_CONFLICT_FREE(bi);														   \
	                                                                                    \
		activeThread = groupIdx < d;																      \
		                                                                                 \
		barrier();															                           \
		T valAi = TREE_SCAN_CACHE[ai];																   \
		T valBi = TREE_SCAN_CACHE[bi];																   \
		                                                                                 \
		barrier();															                           \
		if (activeThread) /* swap */																      \
			TREE_SCAN_CACHE[ai] = valBi;															      \
		                                                                                 \
		barrier();															                           \
		uint origHFAiNext = /*DECODE_ORIG_HF(TREE_SCAN_CACHE_HF[aiNext]);*/            \
			tree_seg_scan_decode_upsweep_hfs_get_origHF(TREE_SCAN_CACHE_HF[aiNext]);	   \
		uint currHFAi = /**DECODE_CURR_HF(TREE_SCAN_CACHE_HF[ai]);*/							\
			tree_seg_scan_decode_upsweep_hfs_get_sumHF(TREE_SCAN_CACHE_HF[ai]);           \
		                                                                                 \
		if (activeThread)																			         \
		{																							            \
			TREE_SCAN_CACHE[bi] = (origHFAiNext == 1u)                                    \
 				? SCAN_ZERO_VAL : ((currHFAi == 1u) ? valAi : OP(valAi, valBi));           \
		}																							            \
		                                                                                 \
		barrier();															                           \
		/* Clear partial sum hf, keep original flag */										   \
		TREE_SCAN_CACHE_HF[ai] &= 0x00000002u;														   \
	}																								            \
	                                                                                    \
	barrier();																                           \
	                                                                                    \
	scanResult_A = TREE_SCAN_CACHE[cacheAddrAi];                                        \
	scanResult_B = TREE_SCAN_CACHE[cacheAddrBi];                                        \
} \



/** \} */




#endif

