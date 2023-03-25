#pragma BLENDER_REQUIRE(bnpr_scan_no_subgroup_codegen_lib.glsl)


/* input buffers:
 * -----------------------------------------------
 * BnprScanDataBuf        bnpr_in_scan_data_buf_
 * BnprScanDataBuf        bnpr_out_scan_data_buf_
 * BnprScanBlockSumBuf    bnpr_scan_block_sum_buf_
*/

#define T_To_Uint(x) x
#define Uint_To_T(x) x

#if defined(_KERNEL_MULTI_COMPILE__TREE_SCAN_UPSWEEP)
void main()
{
	const uint groupId = gl_LocalInvocationID.x;

	TreeScanIndices scan_ids = GetTreeScanIndices(groupId, gl_WorkGroupID.x);

	T scanval_A, scanval_B;
	{ /* init & store random scan input vals */
		scanval_A = T(
		wang_hash(scan_ids.global_x2.y + scan_ids.lds_x2.x) % 12u
		);
		scanval_B = T(
		wang_hash(scan_ids.global_x2.x + scan_ids.lds_x2.y) % 12u
		);
		/* ------------------------------------------------------------------------ */
		/* avoid invalid loads */
		_FUNC_CLEAN_SCAN_DATA(
		scan_ids, ubo_bnpr_tree_scan_infos_.num_scan_items,
		scanval_A, scanval_B /* <- inout */
		);
		/* ------------------------------------------------------------------------ */

		bnpr_in_scan_data_buf_[scan_ids.global_x2.x] = /**floatBitsToUint*/(scanval_A);
		bnpr_in_scan_data_buf_[scan_ids.global_x2.y] = /**floatBitsToUint*/(scanval_B);
	}



	/* execute block-wise exlusive scan */
	T scanRes_ai, scanRes_bi;
	_FUNC_TREE_SCAN_BLOCK(
	groupId,
	gl_WorkGroupID.x,
	scanval_A,
	scanval_B,
	/* -out- */
	scanRes_ai,
	scanRes_bi
	);

	/* store scan results */
	bnpr_out_scan_data_buf_[scan_ids.global_x2.x] = /**floatBitsToUint*/(scanRes_ai);
	bnpr_out_scan_data_buf_[scan_ids.global_x2.y] = /**floatBitsToUint*/(scanRes_bi);

	/* store block aggregate */
	if (groupId == gl_WorkGroupSize.x - 1)
	{
		bnpr_scan_block_sum_buf_[gl_WorkGroupID.x] = /**floatBitsToUint*/(SCAN_OP(scanRes_bi, scanval_B));
	}
}
#endif






#if defined(_KERNEL_MULTI_COMPILE__TREE_SCAN_AGGREGATE)
void main()
{
	const uint groupId =  gl_LocalInvocationID.x;
	const uint gIdx =     gl_WorkGroupID.x;

	TreeScanIndices scanAddrs = GetTreeScanIndices(groupId, 0);
	T aggregate_A = /**uintBitsToFloat*/(bnpr_scan_block_sum_buf_[scanAddrs.global_x2.x]);
	T aggregate_B = /**uintBitsToFloat*/(bnpr_scan_block_sum_buf_[scanAddrs.global_x2.y]);

	T aggregateSum_A, aggregateSum_B;
	_FUNC_TREE_SCAN_AGGREGATE(
	groupId,
	gIdx,
	aggregate_A,
	aggregate_B,
	aggregateSum_A,
	aggregateSum_B
	);

	bnpr_scan_block_sum_buf_[scanAddrs.global_x2.x] = /**floatBitsToUint*/(aggregateSum_A);
	bnpr_scan_block_sum_buf_[scanAddrs.global_x2.y] = /**floatBitsToUint*/(aggregateSum_B);


}
#endif






#if defined(_KERNEL_MULTI_COMPILE__TREE_SCAN_DWSWEEP)
void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint gIdx = gl_WorkGroupID.x;

	TreeScanIndices scan_ids = GetTreeScanIndices(groupId, gl_WorkGroupID.x);

	T block_scan_res_A, block_scan_res_B;
	block_scan_res_A = /**uintBitsToFloat*/(bnpr_out_scan_data_buf_[scan_ids.global_x2[0]]);
	block_scan_res_B = /**uintBitsToFloat*/(bnpr_out_scan_data_buf_[scan_ids.global_x2[1]]);

	T aggregate_sum = /**uintBitsToFloat*/(bnpr_scan_block_sum_buf_[gIdx]);

	bnpr_out_scan_data_buf_[scan_ids.global_x2[0]] = /**floatBitsToUint*/(SCAN_OP(aggregate_sum, block_scan_res_A));
	bnpr_out_scan_data_buf_[scan_ids.global_x2[1]] = /**floatBitsToUint*/(SCAN_OP(aggregate_sum, block_scan_res_B));
}
#endif



#if defined(_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_UPSWEEP)
void main()
{
	const uint groupId =  gl_LocalInvocationID.x;
	const uint idx = gl_GlobalInvocationID.x;

	TreeScanIndices scan_ids = GetTreeScanIndices(groupId, gl_WorkGroupID.x);

	T scanval_A, scanval_B;
	uint hf_A, hf_B;
	{ /* init & store random scan input vals */
		hf_A = 1u & uint(wang_hash(scan_ids.global_x2.x + scan_ids.lds_x2.x) % 128u == 0u);
		if (idx == 0) hf_A = 0u;
			scanval_A = T(
			wang_hash(scan_ids.global_x2.y + scan_ids.lds_x2.x) % 12u
		);
		hf_B = 1u & uint(wang_hash(scan_ids.global_x2.y + scan_ids.lds_x2.y) % 128u == 0u);
			scanval_B = T(
			wang_hash(scan_ids.global_x2.x + scan_ids.lds_x2.y) % 12u
		);

		/* avoid invalid loads */
		/* ----------------------------------------------------------- */
		_FUNC_CLEAN_SEG_SCAN_DATA(
			scan_ids, ubo_bnpr_tree_scan_infos_.num_scan_items,
			hf_A, scanval_A, hf_B, scanval_B /* <- inout */
		);
		/* ----------------------------------------------------------- */

		bnpr_in_scan_data_buf_[scan_ids.global_x2.x] = SEGSCAN_STRUCT_TYPE(scanval_A, hf_A);
		bnpr_in_scan_data_buf_[scan_ids.global_x2.y] = SEGSCAN_STRUCT_TYPE(scanval_B, hf_B);
	}



	/* execute block-wise exlusive scan */
	/* ----------------------------------------------------------- */
	uint headFlagPartialSum_A, headFlagPartialSum_B;
	T scanRes_ai, scanRes_bi;
	_FUNC_TREE_SEG_SCAN_UPSWEEP(
		groupId,
		gl_WorkGroupID.x,
		hf_A, 						scanval_A,
		hf_B, 						scanval_B,
		/* -out- */
		headFlagPartialSum_A, 	scanRes_ai,
		headFlagPartialSum_B, 	scanRes_bi
	);
	/* ----------------------------------------------------------- */

	
	/* store scan results */
	/* ----------------------------------------------------------- */
	bnpr_out_scan_data_buf_[scan_ids.global_x2.x] = SEGSCAN_STRUCT_TYPE(
	scanRes_ai,
		tree_seg_scan_encode_upsweep_hfs(headFlagPartialSum_A, hf_A)
	);
	bnpr_out_scan_data_buf_[scan_ids.global_x2.y] = SEGSCAN_STRUCT_TYPE(
		scanRes_bi,
		tree_seg_scan_encode_upsweep_hfs(headFlagPartialSum_B, hf_B)
	);
	/* ----------------------------------------------------------- */


	/* store block aggregate */
	/* ----------------------------------------------------------- */
	if (groupId == gl_WorkGroupSize.x - 1)
	{
		uint debug_hf = uint((wang_hash(gl_WorkGroupID.x * 17u) % 12u) == 0u); /* debug only */

		bnpr_scan_block_sum_buf_[gl_WorkGroupID.x] = SEGSCAN_STRUCT_TYPE(
			/* different from ordinary scan, we store exclusive sum here */
			scanval_B,
			tree_seg_scan_encode_upsweep_hfs(
				headFlagPartialSum_B, /* or sum of block hfs */
				TREE_SCAN_CACHE_HF[0]  /* original hf of block */
			)
		);
	}
	/* ----------------------------------------------------------- */

}
#endif


#if defined(_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_AGGREGATE)
void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint gIdx = gl_WorkGroupID.x;

	TreeScanIndices scan_ids = GetTreeScanIndices(groupId, 0);

	SEGSCAN_STRUCT_TYPE aggregate_A = bnpr_scan_block_sum_buf_[scan_ids.global_x2.x];
	T partialSumTreeAi = aggregate_A.val;
	uint partialOrTreeAi = tree_seg_scan_decode_upsweep_hfs_get_sumHF(aggregate_A.hf);
	uint firstInitialHFAi = tree_seg_scan_decode_upsweep_hfs_get_origHF(aggregate_A.hf);

	SEGSCAN_STRUCT_TYPE aggregate_B = bnpr_scan_block_sum_buf_[scan_ids.global_x2.y];
	T partialSumTreeBi = aggregate_B.val;
	uint partialOrTreeBi = tree_seg_scan_decode_upsweep_hfs_get_sumHF(aggregate_B.hf);
	uint firstInitialHFBi = tree_seg_scan_decode_upsweep_hfs_get_origHF(aggregate_B.hf);



	T upsweep_res_sum_A, upsweep_res_sum_B;
	uint upsweep_res_hf_A, upsweep_res_hf_B;
	_FUNC_TREE_SEG_SCAN_UPSWEEP(
		groupId,
		gl_WorkGroupID.x,
		partialOrTreeAi, partialSumTreeAi,
		partialOrTreeBi, partialSumTreeBi,
		/* -out- */
		upsweep_res_hf_A, upsweep_res_sum_A,
		upsweep_res_hf_B, upsweep_res_sum_B
	);



	_FUNC_TREE_SEG_SCAN_AGGREGATE_FILL_CACHE(
	groupId, scan_ids,
		/* --- LDS inputs --- */
		firstInitialHFAi, firstInitialHFBi
	);



	T scan_res_A, scan_res_B;
	_FUNC_TREE_SEG_SCAN_DWSWEEP
	(
		groupId, gIdx,
		/* --- out --- */
		scan_res_A, scan_res_B
	);


	/* store scan results */
	bnpr_scan_block_sum_buf_[scan_ids.global_x2.x] = SEGSCAN_STRUCT_TYPE(scan_res_A, 0); /* no hf needed */
	bnpr_scan_block_sum_buf_[scan_ids.global_x2.y] = SEGSCAN_STRUCT_TYPE(scan_res_B, 0); /* no hf needed */
}
#endif


#if defined(_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_DWSWEEP)
void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint gIdx = gl_WorkGroupID.x;

	TreeScanIndices scan_ids = GetTreeScanIndices(groupId, gl_WorkGroupID.x);

	SEGSCAN_STRUCT_TYPE block_scan_res_A, block_scan_res_B, aggregate_scan_res;
	block_scan_res_A = /**uintBitsToFloat*/(bnpr_out_scan_data_buf_[scan_ids.global_x2.x]);
	block_scan_res_B = /**uintBitsToFloat*/(bnpr_out_scan_data_buf_[scan_ids.global_x2.y]);
	aggregate_scan_res = /**uintBitsToFloat*/(bnpr_scan_block_sum_buf_[gIdx]);

	_FUNC_TREE_SEG_SCAN_DWSWEEP_FILL_CACHE(
	    groupId, scan_ids,
	    /* --- block partial sums --- */
	    block_scan_res_A.hf, block_scan_res_A.val,
	    block_scan_res_B.hf, block_scan_res_B.val,
	    /* --- scanned block aggregate --- */
	    aggregate_scan_res.val
	);

	T scan_res_A, scan_res_B;
	_FUNC_TREE_SEG_SCAN_DWSWEEP
	(
		groupId, gIdx,
		/* --- out --- */
		scan_res_A, scan_res_B
	);

	bnpr_out_scan_data_buf_[scan_ids.global_x2.x] = SEGSCAN_STRUCT_TYPE(scan_res_A, block_scan_res_A.hf);
	bnpr_out_scan_data_buf_[scan_ids.global_x2.y] = SEGSCAN_STRUCT_TYPE(scan_res_B, block_scan_res_B.hf);
}
#endif
