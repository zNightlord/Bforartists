
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)


#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION)
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
	const uint contour_id = idx; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = contour_id < num_contours; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION__PASS_0)
    uint head_contour_id = ssbo_contour_edge_list_head_in_[contour_id];
    uint list_len        = ssbo_contour_edge_list_len_in_[contour_id]; 
    uint rank            = ssbo_contour_edge_rank_in_[contour_id]; 

	bool is_head_contour = rank == 0; 
	// TODO: can be bad if ssbo_contour_to_contour has invalid links
	// should rectify the links in the future
	uint prev_contour_id = ssbo_contour_to_contour_[2u*contour_id]; 
	bool prev_equ_self = prev_contour_id == contour_id; 
    
    uint output_addr = head_contour_id + rank; 
    if (valid_thread)
    {
		// transfer topology data from list ranking
        ssbo_contour_edge_rank_out_[output_addr]      = rank; 
        ssbo_contour_edge_list_len_out_[output_addr]  = list_len;
        ssbo_contour_edge_list_head_out_[output_addr] = head_contour_id;

		// transfer packed data from the intermediate buffer
		uvec4 enc_data[2];
		Load4(ssbo_contour_edge_transfer_data_, contour_id*2u,    enc_data[0]); 
		Load4(ssbo_contour_edge_transfer_data_, contour_id*2u+1u, enc_data[1]);

        uvec3 vpos_0_enc, vpos_1_enc; 
		vpos_0_enc = enc_data[0].xyz;
		vpos_1_enc = uvec3(enc_data[0].w, enc_data[1].xy);
		
		ContourFlags cf = decode_contour_flags(enc_data[1].z); 
		if (is_head_contour && prev_equ_self)
			cf.looped_curve = false; // will broadcast this next pass
		else 
			cf.looped_curve = true; 

		Store3(ssbo_contour_edge_vpos_out_, output_addr*2u,    vpos_0_enc);
        Store3(ssbo_contour_edge_vpos_out_, output_addr*2u+1u, vpos_1_enc);
		store_contour_flags(output_addr, cf); 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION__PASS_1)
    uint head_contour_id = ssbo_contour_edge_list_head_out_[contour_id];

	ContourFlags cf = load_contour_flags(contour_id);
	ContourFlags cf_head = load_contour_flags(head_contour_id); 
	if (valid_thread && !cf_head.looped_curve)
	{ /* broadcast loop flag to all contour edges in the curve */
		cf.looped_curve = false;
		store_contour_flags(contour_id, cf); /* race condition does not hurt */

		/* copy to segloopconv1d input buffer */
		ssbo_in_segloopconv1d_data_[contour_id] = encode_contour_flags(cf); 
	}

	if (idx == 0u)
		ssbo_segloopconv1d_info_.num_conv_items = num_contours; 
#endif

}
#endif



#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION)
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
	const uint contour_id = idx; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = contour_id < num_contours; 

	ContourFlags cf = load_contour_flags(contour_id);
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf);
	bool is_curve_tail = cct.tail_contour_id == contour_id;
	bool is_curve_head = cct.head_contour_id == contour_id; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__SETUP)
	uint next_contour_id = is_curve_tail ? cct.head_contour_id : contour_id + 1u;
	ContourFlags cf_next = load_contour_flags(next_contour_id);
	
	bool is_scan_seg_head = cf.seg_head 	 || is_curve_head;  
	bool is_scan_seg_tail = cf_next.seg_head || is_curve_tail;  

	if (valid_thread)
	{
		if (cf_next.seg_head)
		{
			cf.seg_tail = true; 
			store_contour_flags(contour_id, cf); 
		}

		uint scan_item_id = contour_id; 
		uint scan_val = 1u; 

		uint hf = is_scan_seg_head ? 1u : 0u; 
		scan_data_buf_0_[scan_item_id] = segscan_uint_hf_encode(SSBOData_SegScanType_uint(scan_val, hf)); 

		scan_item_id = num_contours - 1u - contour_id; 
		hf = is_scan_seg_tail ? 1u : 0u; 
		scan_data_buf_1_[scan_item_id] = segscan_uint_hf_encode(SSBOData_SegScanType_uint(scan_val, hf)); 
	}

	if (idx.x == 0) /* I dont care, this is cheap, just keep this up to date here */
	{
		ssbo_tree_scan_infos_contour_segmentation_.num_scan_items = num_contours; 
		ssbo_tree_scan_infos_contour_segmentation_.num_valid_scan_threads = compute_num_threads(
			num_contours, 2u
		);
		ssbo_tree_scan_infos_contour_segmentation_.num_thread_groups = compute_num_groups(
			num_contours, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
		);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__FINISH)
#define REVERSE_ID(id) ((num_contours - 1u - id))
	SSBOData_SegScanType_uint scan_res_step_0 = segscan_uint_hf_decode(scan_output_buf_0_[contour_id]);
	SSBOData_SegScanType_uint scan_res_step_1 = segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(contour_id)]);

	uint dist_to_seg_head = scan_res_step_0.val; // == contour_id - head_id
	uint dist_to_seg_tail = scan_res_step_1.val; // == tail_id - contour_id

	uint seg_len = dist_to_seg_head + dist_to_seg_tail + 1u; 
	uint seg_head = contour_id - dist_to_seg_head;	
	uint seg_tail = contour_id + dist_to_seg_tail; 
	uint seg_rank = dist_to_seg_head; 

	// Fix seg rank & seg len for looped curve
	if (cf.looped_curve)
	{
		ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf); 
		bool first_seg_in_loop = seg_head == cct.head_contour_id;
		bool last_seg_in_loop  = seg_tail == cct.tail_contour_id;
		bool is_self_loop = seg_len == cct.len; 

		if (first_seg_in_loop && !is_self_loop)
		{
			ContourFlags cf_curv_head = load_contour_flags(cct.head_contour_id);
			bool downflow = !(cf_curv_head.seg_head);

			if (downflow)
			{ // prev half segment, at the end of loop
				SSBOData_SegScanType_uint scan_res_step_0_tailseg = 
					segscan_uint_hf_decode(scan_output_buf_0_[cct.tail_contour_id]);
				// SSBOData_SegScanType_uint scan_res_step_1_tailseg = 
				// 	segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(cct.tail_contour_id)]); // should equal to 0
				
				uint tail_seg_len = scan_res_step_0_tailseg.val + 1u; 
				seg_rank += tail_seg_len; // fix seg rank
				seg_len += tail_seg_len; // fix seg len
			}
		}
		if (last_seg_in_loop && !is_self_loop)
		{
			ContourFlags cf_curv_tail = load_contour_flags(cct.tail_contour_id);
			bool overflow = !(cf_curv_tail.seg_tail);

			if (overflow)
			{ // next half segment, at the head of loop
				// SSBOData_SegScanType_uint scan_res_step_0_headseg = 
				// 	segscan_uint_hf_decode(scan_output_buf_0_[cct.head_contour_id]); // should equal to 0
				SSBOData_SegScanType_uint scan_res_step_1_headseg = 
					segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(cct.head_contour_id)]);

				uint head_seg_len = scan_res_step_1_headseg.val + 1u;
				seg_len += head_seg_len; // fix seg len
			}
		}
	}

	if (valid_thread)
	{ // Store to buffer
		// Update seg_rank & seg_len
		ssbo_contour_edge_seg_rank_[contour_id] = seg_rank;
		ssbo_contour_edge_seg_len_[contour_id] = seg_len;

		// TODO: Update seg_head/tail for curve head/tail in next pass
		// cannot do it here due to data racing
	}

#endif


}
#endif









#if defined(_KERNEL_MULTICOMPILE__CALC_CONTOUR_EDGE_SCREEN_POSITIONS)
/*
 * ubo_view_matrices_
 * pcs_screen_size_
*/

vec2 ndc_to_screen_uv(vec4 pos_ndc)
{
	pos_ndc.xyz /= pos_ndc.w;
	pos_ndc.xy = pos_ndc.xy * 0.5f + 0.5f; 
	return pos_ndc.xy * pcs_screen_size_.xy; /* TODO: flip y or not? */
}

void main()
{
	const uint idx = gl_GlobalInvocationID.x; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours; 

	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 world_to_view = ubo_view_matrices_.viewmat;
	mat4 mat_camera_proj = ubo_view_matrices_.winmat; 

	vec4 vpos_ws[2]; /* v0, v1 on edge */
	vec3 vnor_ws[2]; 
	vec4 vpos_ndc[2]; 
	vec2 vpos_uv[2]; 

	if (valid_thread)
	{ /* read vertex pos transformed to world space */
	  /* Note: wpos_and_edgeid will be overwrite, for saving space */
        uvec3 vpos_0_enc, vpos_1_enc;
        Load3(ssbo_contour_edge_vpos_, contour_id*2u,    vpos_0_enc);
        Load3(ssbo_contour_edge_vpos_, contour_id*2u+1u, vpos_1_enc);
        vec3 vpos_0, vpos_1;
        vpos_0 = uintBitsToFloat(vpos_0_enc);
        vpos_ws[0] = vec4(vpos_0.xyz, 1.0f);
		vpos_ndc[0] = mat_camera_proj * vec4((world_to_view * vpos_ws[0]).xyz, 1.0f); 
        vpos_1 = uintBitsToFloat(vpos_1_enc);
        vpos_ws[1] = vec4(vpos_1.xyz, 1.0f);
		vpos_ndc[1] = mat_camera_proj * vec4((world_to_view * vpos_ws[1]).xyz, 1.0f); 

		/* write to mesh pool TODO: these should be removed */
		uint addr_st = mesh_pool_addr__zwhclip(contour_id); 
		buf_strokegen_mesh_pool[addr_st+0] = floatBitsToUint(vpos_ndc[0].z); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_ndc[0].w); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_ndc[1].z); 
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_ndc[1].w); 
		
        addr_st = mesh_pool_addr__edgedir(contour_id); 
		vpos_uv[0].xy = ndc_to_screen_uv(vpos_ndc[0]);
		vpos_uv[1].xy = ndc_to_screen_uv(vpos_ndc[1]);
		vec2 edge_dir = (vpos_uv[1].xy - vpos_uv[0].xy); 
		float edge_dir_len = length(edge_dir);
		vec2 edge_dir_norm = edge_dir / edge_dir_len; 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(edge_dir_norm.x);
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(edge_dir_norm.y);

		addr_st = mesh_pool_addr__edgeuv(contour_id); 
		vpos_uv[0].xy /= pcs_screen_size_.xy; 
		vpos_uv[1].xy /= pcs_screen_size_.xy; 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(vpos_uv[0].x); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_uv[0].y); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_uv[1].x);
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_uv[1].y);
	}



}
#endif