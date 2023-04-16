#pragma BLENDER_REQUIRE(npr_strokegen_hlsl_support_lib.glsl)

#ifndef BNPR_SEGLOOPCONV1D_LIB_INCLUDED
#define BNPR_SEGLOOPCONV1D_LIB_INCLUDED

void func_device_store_loopconv1d_patch_id(uint stAddr, uint val)
{
	ssbo_segloopconv1d_patch_table_[stAddr] = val;
}

/* Macro expansion */
#ifndef CAT
	#define CAT_(x, y) x ## y
	#define CAT(x, y) CAT_(x, y)
#endif





/* Inputs /////////////////////////////////////////////////////////////////////
 *-Common------------------------------*                                                
 * LOOPCONV1D_TAG              Test    *
 * LOOPCONV1D_MAX_RADIUS       32u     *
 * LOOPCONV1D_GROUP_SIZE       256u    *
*/
#define tag                         LOOPCONV1D_TAG
#define MAX_CONV_RADIUS             LOOPCONV1D_MAX_RADIUS
#define GROUP_SIZE_CONV             ((gl_WorkGroupSize.x))
/*-Convolution------------------------------------*
 * #if _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION
 *    DATA_TYPE_LOOPCONV1D                   float
 *    FUNC_DEVICE_LOAD_LOOPCONV1D_PATCH_ID
 *    FUNC_DEVICE_LOAD_LOOPCONV1D_DATA
 * #endif 
*/ 
#if _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION
	#define T_CONV                      DATA_TYPE_LOOPCONV1D
	#define DEVICE_LOAD_CONV_DATA       FUNC_DEVICE_LOAD_LOOPCONV1D_DATA
	#define DEVICE_LOAD_CONV_PATCH_ID   FUNC_DEVICE_LOAD_LOOPCONV1D_PATCH_ID
#endif
/*-Build Conv Acc Cache----------------------------*
 * #if _KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE
 *    FUNC_DEVICE_STORE_LOOPCONV1D_PATCH_ID
 * #endif
*/
#if _KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE
	#define DEVICE_STORE_CONV_PATCH_ID  FUNC_DEVICE_STORE_LOOPCONV1D_PATCH_ID
#endif






/* Common Defines ///////////////////////////////////////////////// */
#define NUM_PATCHES_PER_GROUP ((2u * MAX_CONV_RADIUS))

#define _FUNC_GET_PATCH_TABLE_BUFFER_INDEX CAT(GetPatchTableBufferIndex_, tag)
/* Calculate global buffer address for item in per-thread-group patch table */
uint _FUNC_GET_PATCH_TABLE_BUFFER_INDEX(uint blockIdx, uint patchId)
{
	return (blockIdx * (NUM_PATCHES_PER_GROUP)) + patchId; 
}



/* 1D Convolution on Circular Segments  //////////////////////////////////////////////////////////// */
#ifdef _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION

	/* Avaliable Interfaces /////////////////////////////////////////////// */
	#define _FUNC_SETUP_SEGLOOP1DCONV      CAT(SetupSegmentedConvolution_, tag) /* Setup LDS caches needed for convolution */
	#define _FUNC_LOAD_CONV_DATA_LDS_LEFT  CAT(LoadLDSConvData_AtLeft_, tag) /* load neighbor item with ngtv circular offset */
	#define _FUNC_LOAD_CONV_DATA_LDS_RIGHT CAT(LoadLDSConvData_AtRight_, tag) /* load neighbor item with pstv circular offset */


	/* 1) Main Cache ////////////////////////////////////////////
	 * Maintains per-elem data to be convoluted, 
	 * with padding on both sides */
	#define CONV_LDS_LEN ((GROUP_SIZE_CONV + 2 * MAX_CONV_RADIUS))
	shared T_CONV CAT(LDS_ConvData_, tag)[CONV_LDS_LEN];
	/* Main Cache Layout:   
	 * Suppose convolution radius == 4, group size == 6: 
	 * --- Padding is applied to LDS_ConvData_tag[]: ---
	 *	0   1   2   3   4	5   6   7   8   9  10  11  12  13   Main_Cache[] ("LDS_ConvData_tag[]")
	 *	|<-Padding->|				                |<--Padding->|
	 *  radius==4		   groupsize==6           radius==4 
	 */
	
	
	/* 2) Patch Cache ////////////////////////////////////////////
	 * Patching for loop segments */
	/*	Padding is not enough for convolution when
	 *	circular segment(loop) exists.
	 *	FIRST-TAIL and LAST-HEAD can cause a long "jump"
	 *	when convolution goes across them, and the jump
	 * may lead to elements far away from LDS_ConvData_tag.
	 *
	 * FIRST-TAIL t and LAST-HEAD h have special patches,
	 * stored in LDS_ConvPatch_tag.
	 * For details, see "https://zhuanlan.zhihu.com/p/263566817" 
	 */
	shared T_CONV CAT(LDS_ConvPatch_, tag)[NUM_PATCHES_PER_GROUP];
	/* Patch Cache Layout:  (let MAX_CONV_RADIUS==3)
	/* ElemIdGl |1stSegHead, +1, +2, lastSegTail-2, -1, -0, <- LDS_ConvPatch_tag[0~5]*/
	/* patchId  |         0,  1,  2,             3,  4,  5   */
	
	
	
	/* Data Indexing ////////////////////////////////////////////
	 * ElemIdGl: index used in Global Memory, usually id of the convolved element
	 * ElemIdLc: Main Cache index ("LDS_ConvData_tag[]")
	 * ElemIdPatch: Patch Cache index ("LDS_ConvPatch_tag[]")
	 */
	/* gl_LocalInvocationID.x -> ElemIdLc */
	int CAT(GroupIdx_To_ElemIdLc_, tag)(uint blockId, uint groupIdx)
	{
		return (blockId == 0u)
			? int(groupIdx)
			: int(groupIdx + MAX_CONV_RADIUS);
	}
	
	/* gl_LocalInvocationID.x -> ElemIdGl */
	int CAT(GroupIdx_To_ElemIdGl_, tag)(uint blockId, uint groupIdx)
	{ /* essentially global invoc id.x */
		return int(groupIdx + blockId * GROUP_SIZE_CONV);
	}
	
	
	/* ElemIdLc -> ~Gl*/
	int CAT(ElemIdLc_To_Gl_, tag)(int elemIdLc, uint blockId)
	{
		int elemIdGlStart =
			(blockId == 0u) ? 0
			: (int(blockId * GROUP_SIZE_CONV) - int(MAX_CONV_RADIUS));
		
		return elemIdLc + elemIdGlStart;
	}
	/* ElemIdGl -> ~Lc*/
	int CAT(ElemIdGl_To_Lc_, tag)(int elemIdGl, uint blockId)
	{
		int elemIdGlStart =
			(blockId == 0u) ? 0
			: int(blockId * GROUP_SIZE_CONV) - int(MAX_CONV_RADIUS);
		
		return elemIdGl - elemIdGlStart;
	}
	
	
	/* If elem is indexed using right half of the patch cache */
	bool CAT(IsElemIdLc_RightPatch_, tag)(int elemIdLc, uint blockId)
	{
		return (elemIdLc >= (
			blockId != 0u
				? int(CONV_LDS_LEN)
				: int(CONV_LDS_LEN - MAX_CONV_RADIUS)
		));
	}
	/* If elem is indexed using left half of the patch cache */
	bool CAT(IsElemIdLc_LeftPatch_, tag)(int elemIdLc, uint blockId)
	{
		return elemIdLc < 0;
	}
	
	
	/* elemIdGl of curr elem & segment tail  -> elemIdPatch */
	uint CAT(RightPatchElemId_LastHead_, tag)(
		uint segTailIdGl, uint elemIdGl
	) {
		/* elemIdGl | segTail-2, -1, -0, (let MAX_CONV_RADIUS==3) */
		/* patchId  |		    3,  4,  5                            */
		uint dist = (segTailIdGl - elemIdGl + 1u);
		return 2u * MAX_CONV_RADIUS - dist;
	}
	/* elemIdGl of curr elem & segment head -> elemIdPatch */
	uint CAT(LeftPatchElemId_FirstTail_, tag)(
		uint segHeadIdGl, uint elemIdGl
	) {
		return elemIdGl - segHeadIdGl;
	}
	
	
	/* Move element along its loop segment.
	 * Make sure offset<=MAX_CONV_RADIUS */
	void CAT(MoveConvElemId_, tag)(
		bool moveLeft, uint offset,
		uint blockId, uint groupIdx,
		uint segLen, uint segHeadId,
		out uint elemIdGl, out int elemIdLc
	) {
		elemIdGl = CAT(GroupIdx_To_ElemIdGl_, tag)(blockId, groupIdx);
		
		offset = offset % segLen;
		
		elemIdGl -= segHeadId;
		elemIdGl += (moveLeft ? (segLen - offset) : offset);
		elemIdGl = elemIdGl % segLen;
		elemIdGl += segHeadId;
		
		elemIdLc = CAT(ElemIdGl_To_Lc_, tag)(elemIdGl, blockId);
	}
	
	/* Setup Patch LDS Cache ------------------*/
	/* Load elemIdGl at patch_cache[patchIdLc] */
	void CAT(PatchData_LoadDevice_StoreLDS_, tag)(
		uint blockId, uint patchIdLc, uint elemCount
	) {
		if (patchIdLc < NUM_PATCHES_PER_GROUP)
		{ 
			uint patchBufferAddress = _FUNC_GET_PATCH_TABLE_BUFFER_INDEX(blockId, patchIdLc);
			uint patchIdGl = DEVICE_LOAD_CONV_PATCH_ID(patchBufferAddress);
			CAT(LDS_ConvPatch_, tag)[patchIdLc] = DEVICE_LOAD_CONV_DATA(patchIdGl);
		}
	}
	
	/* Setup Main LDS Cache ----------------------------------------------*/
	/* Load conv data to main cache, for thread@groupIdx in group@blockId */
	void CAT(ConvData_LoadDevice_StoreLDS_, tag)(
		uint blockId, uint groupIdx, out T_CONV convData
	) {
		int elemIdGl = CAT(GroupIdx_To_ElemIdGl_, tag)(blockId, groupIdx);
		int elemIdLc = CAT(GroupIdx_To_ElemIdLc_, tag)(blockId, groupIdx);
		
		convData = DEVICE_LOAD_CONV_DATA(elemIdGl);
		
		CAT(LDS_ConvData_, tag)[elemIdLc] = convData;
	}
	/* Load PADDING conv data to main cache, for thread@groupIdx in group@blockId */
	void CAT(Padding_LoadDevice_StoreLDS_, tag)(uint blockId, uint groupIdx)
	{
		/* Determine active threads for padding */
		bool leftPadding = (0 < blockId) && (groupIdx < MAX_CONV_RADIUS);    /* workers for left padding */
		bool rightPadding = (GROUP_SIZE_CONV - MAX_CONV_RADIUS) <= groupIdx; /* workers for right padding */
		
		int paddingIdLc =
			leftPadding ? (groupIdx) : (
				(blockId == 0) ? (groupIdx + MAX_CONV_RADIUS) : (groupIdx + 2 * MAX_CONV_RADIUS) /* only padding at right side */
			);
		
		uint paddingIdGl = CAT(ElemIdLc_To_Gl_, tag)(paddingIdLc, blockId);
		
		/*[branch]*/if (leftPadding || rightPadding)
		{
			CAT(LDS_ConvData_, tag)[paddingIdLc] = DEVICE_LOAD_CONV_DATA(paddingIdGl);
		}
	}

	/**
	 * \brief Setup everything for a 1d segmented convolution
	 */
	void _FUNC_SETUP_SEGLOOP1DCONV(
		uint3 gIdx, uint groupIdx, uint elemCount,
		out T_CONV convData
	){
		CAT(PatchData_LoadDevice_StoreLDS_, tag)(
			gIdx.x, groupIdx, elemCount
		);
	
		CAT(ConvData_LoadDevice_StoreLDS_, tag)(
			gIdx.x, groupIdx, /*out*/convData
		);
	
		CAT(Padding_LoadDevice_StoreLDS_, tag)(
			gIdx.x, groupIdx
		);
		GroupMemoryBarrierWithGroupSync();
	}
	
	/**
	 * \brief Load convolution data with left offset
	 * \param offset assert(offset <= MAX_CONV_RADIUS)
	 */
	T_CONV _FUNC_LOAD_CONV_DATA_LDS_LEFT(
		uint offset,
		uint blockId, uint groupIdx : SV_GroupIndex,
		uint segLen, uint segHeadId
	)
	{
		T_CONV convData;
		
		uint elemIdGl = 0;
		int elemIdLc = 0;
		
		CAT(MoveConvElemId_, tag)(
			true, offset,
			blockId, groupIdx,
			segLen, segHeadId,
			// out -----------
			elemIdGl, elemIdLc
		);
		
		bool patch = CAT(IsElemIdLc_RightPatch_, tag)(elemIdLc, blockId);
		[branch] if (patch)
		{
			uint patchId = CAT(RightPatchElemId_LastHead_, tag)(
				segHeadId + segLen - 1,
				elemIdGl
			);
			convData = CAT(LDS_ConvPatch_, tag)[patchId];
		}
		else
		{
			convData = CAT(LDS_ConvData_, tag)[elemIdLc];
		}
		
		return convData;
	}
	
	/**
	 * \brief Load convolution data with left offset
	 * \param offset assert(offset <= MAX_CONV_RADIUS)
	 */
	T_CONV _FUNC_LOAD_CONV_DATA_LDS_RIGHT(
	  uint offset,
	  uint blockId, uint groupIdx,
	  uint segLen, uint segHeadId
	)
	{
		T_CONV convData;
		
		uint elemIdGl = 0;
		int elemIdLc = 0;
		
		CAT(MoveConvElemId_, tag) (
			false, offset,
			blockId, groupIdx,
			segLen, segHeadId,
			// out -----------
			elemIdGl, elemIdLc
		);
		
		bool patch = CAT(IsElemIdLc_LeftPatch_, tag) (elemIdLc, blockId);
		[branch] if (patch)
		{
			uint patchId = CAT(LeftPatchElemId_FirstTail_, tag) (
		      segHeadId, elemIdGl
			);
			convData = CAT(LDS_ConvPatch_, tag) [patchId];
		}
		else
		{
			convData = CAT(LDS_ConvData_, tag) [elemIdLc];
		}
		
		return convData;
	}
#endif /* _KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION */





/* Build Patch Table for Given set of segmened loops ///////////////////////////////// */
#ifdef _KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE
	
	/* Avaliable Interfaces /////////////////////////////////////////////// */
	#define _FUNC_COMPUTE_PATCH_TABLE   CAT(ComptueBlockPatchElemIds_, tag)
	#define _LDS_PATCH_ELEMIDGL         CAT(LDS_PatchElemIds_, tag)


	#ifdef NULL_LEFT_TAIL
		#undef NULL_LEFT_TAIL
	#endif
	#define NULL_LEFT_TAIL (0xffffffff)
	
	#ifdef NULL_RIGHT_HEAD
		#undef NULL_RIGHT_HEAD
	#endif
	#define NULL_RIGHT_HEAD 0
	
	#ifdef INVALID_PATCH_EDGE_ID
		#undef INVALID_PATCH_EDGE_ID
	#endif
	#define INVALID_PATCH_EDGE_ID 0xffffffff


	shared uint CAT(LDS_LeftTailGroupId_, tag)  = NULL_LEFT_TAIL;
	shared uint CAT(LDS_RightHeadGroupId_, tag) = NULL_RIGHT_HEAD;
	shared uint _LDS_PATCH_ELEMIDGL[MAX_CONV_RADIUS * 2];
	
	uint CAT(MoveElemIdAlongLoop_, tag)(uint elemId, int offset, uint loopHeadElemId, uint loopLen)
	{
		bool moveLeft = offset < 0;
	
		uint d = abs(offset);
		d = d % loopLen;
	
		elemId -= loopHeadElemId;
		elemId += (moveLeft ? (loopLen - d) : d);
		elemId = (elemId % loopLen);
		elemId += loopHeadElemId;
	
		return elemId;
	}

	bool CAT(IsActiveThreadForStorePatchTable_, tag)(uint groupIdx)
	{
		return groupIdx < NUM_PATCHES_PER_GROUP;
	}
	

	/**
	 * \brief _FUNC_COMPUTE_PATCH_TABLE
	 * Computes patching info for current Thread Group,
	 * and stores at group shared memory _LDS_PATCH_ELEMIDGL
	 * 
	 * For details, 
	 * see comments on Patch Cache above
	 * and "https://zhuanlan.zhihu.com/p/263566817"
	 * 
	 * \param segHead head element id
	 * \param segTail tail element id
	 * \param segLen segment length
	 */
	void _FUNC_COMPUTE_PATCH_TABLE(
		uint blockId, uint groupIdx, 
		bool isSegHead, bool isSegTail, 
		uint segHead, uint segTail, uint segLen
	){
		barrier();
	
		/* ------------Loading Extra Neighboring Data---------------- */
		/* Step 1. Vote for right - most head && left - most tail */
		if (isSegHead)
		{
			atomicMax(CAT(LDS_RightHeadGroupId_, tag), groupIdx + 1u);
		}
		if (isSegTail)
		{
			atomicMin(CAT(LDS_LeftTailGroupId_, tag), groupIdx);
		}
		barrier();
	
		/* Step 2. Compute patching data address */
		bool foundLeftTail = 
			(CAT(LDS_LeftTailGroupId_, tag) != NULL_LEFT_TAIL);
		bool foundRightHead = 
			(CAT(LDS_RightHeadGroupId_, tag) != NULL_RIGHT_HEAD);
	
		if ((foundLeftTail && (groupIdx == CAT(LDS_LeftTailGroupId_, tag))) ||
			((!foundLeftTail) && (groupIdx == GROUP_SIZE_CONV - 1u)))
		{
			/* Patch at LEFT */
			/*[unroll]*/
			for (uint i = 0u; i < MAX_CONV_RADIUS; ++i)
			{
				_LDS_PATCH_ELEMIDGL[i] =
					CAT(MoveElemIdAlongLoop_, tag)(
						segHead,
						int(i), /*offset*/
						segHead, segLen
					);
			}
		}
	
		if ((foundRightHead && ((groupIdx + 1u) == CAT(LDS_RightHeadGroupId_, tag))) 
			|| ((!foundRightHead) && (groupIdx == 0u)))
		{
			uint patchStart = 
				CAT(MoveElemIdAlongLoop_, tag)(
					segTail,
					-int(MAX_CONV_RADIUS - 1), /* offset */
					segHead, segLen
				);
			/* Patch at RIGHT */
			/*[unroll]*/
			for (uint i = 0u; i < MAX_CONV_RADIUS; ++i)
			{
				_LDS_PATCH_ELEMIDGL[i + MAX_CONV_RADIUS] =
					CAT(MoveElemIdAlongLoop_, tag)(
						patchStart,
						int(i), /* offset */
						segHead, segLen
					);
			}
		}
		barrier();
		
		uint patchIdLc = groupIdx;
		if (CAT(IsActiveThreadForStorePatchTable_, tag)(patchIdLc))
		{
			uint patchBufferAddress = _FUNC_GET_PATCH_TABLE_BUFFER_INDEX(blockId, patchIdLc);
			DEVICE_STORE_CONV_PATCH_ID(patchBufferAddress, _LDS_PATCH_ELEMIDGL[patchIdLc]);
		}
	}

#endif 
/* _KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE */



#endif 
/* BNPR_SEGLOOPCONV1D_LIB_INCLUDED */