
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)



#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION)
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
	const uint contour_id = idx; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = contour_id < num_contours; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION__PASS_0)
	if (idx == 0u)
	{
		uint num_snakes = ssbo_list_ranking_addressing_counters_[0]; 
		ssbo_bnpr_mesh_pool_counters_.num_contour_verts = num_snakes; 
		ssbo_segloopconv1d_info_.num_conv_items 		= num_snakes; 
	}

    uint head_vtx_addr; bool looped_curve; 
	uint encoded_info = ssbo_list_ranking_list_head_info_[contour_id*2u + 1u];
	head_vtx_addr = (encoded_info >> 1u); 
	looped_curve  = ((encoded_info & 1u) == 1u); 

	uint num_edges_in_curve = ssbo_contour_edge_list_len_in_[contour_id]; 
	uint num_verts_in_curve = looped_curve ? num_edges_in_curve : num_edges_in_curve + 1u; 

	uint edge_rank_in_curve = ssbo_contour_edge_rank_in_[contour_id]; 
	uint vtx_rank_in_curve = edge_rank_in_curve; 

	bool is_head_edge = edge_rank_in_curve == 0; 
	bool is_tail_edge = edge_rank_in_curve == num_edges_in_curve - 1; 
    
    uint vtx_addr = head_vtx_addr + edge_rank_in_curve; 
	bool additional_output_tail_vtx = is_tail_edge && !looped_curve; 
	additional_output_tail_vtx = additional_output_tail_vtx;
    if (valid_thread)
    {
		// transfer topology data from list ranking
        ssbo_contour_snake_rank_[vtx_addr]      = vtx_rank_in_curve; 
        ssbo_contour_snake_list_len_[vtx_addr]  = num_verts_in_curve;
        ssbo_contour_snake_list_head_[vtx_addr] = head_vtx_addr;
		if (additional_output_tail_vtx)
		{ // last edge outputs its end vertex
			ssbo_contour_snake_rank_[vtx_addr + 1u]      = vtx_rank_in_curve + 1; 
			ssbo_contour_snake_list_len_[vtx_addr + 1u]  = num_verts_in_curve;
			ssbo_contour_snake_list_head_[vtx_addr + 1u] = head_vtx_addr;
		}

		// transfer packed edge data from the intermediate buffer
		ContourEdgeTransferData cetd = load_contour_edge_transfer_data(contour_id); 

		ContourFlags cf = cetd.cf; 
		init_contour_looped_curve(looped_curve, cf);
		init_contour_curve_head(is_head_edge, cf);
		init_contour_curve_tail(is_tail_edge, cf);

		vec2 cusp_func_v = cetd.cusp_funcs; 
		set_contour_cusp_flags(.0f < cusp_func_v[0], cf); 

		Store3(ssbo_contour_snake_vpos_, vtx_addr, floatBitsToUint(cetd.vpos_ws[0]));
		store_contour_flags(vtx_addr, cf); 
		
		ssbo_in_segloopconv1d_data_[vtx_addr] = encode_contour_flags(cf); // copy to segloopconv1d input buffer
		
		
		if (additional_output_tail_vtx)
		{
        	Store3(ssbo_contour_snake_vpos_, vtx_addr+1u, floatBitsToUint(cetd.vpos_ws[1]));
			cf.seg_head = false; 
			init_contour_curve_head(false, cf);
			init_contour_curve_tail(true, cf);
			set_contour_cusp_flags(.0f < cusp_func_v[1], cf); 
			store_contour_flags(vtx_addr+1u, cf);

			ssbo_in_segloopconv1d_data_[vtx_addr+1u] = encode_contour_flags(cf); // copy to segloopconv1d input buffer
		}
    }
#endif

}
#endif



#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION)
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
	const uint contour_id = idx; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	bool valid_thread = contour_id < num_contours; 

	ContourFlags cf = load_contour_flags(contour_id);
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf);
	bool is_curve_head = cct.head_contour_id == contour_id; 
	bool is_curve_tail = cct.tail_contour_id == contour_id; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__SETUP)

	uint prev_contour_id = is_curve_head ? cct.tail_contour_id : contour_id - 1u;
	ContourFlags cf_prev = load_contour_flags(prev_contour_id);

	uint next_contour_id = is_curve_tail ? cct.head_contour_id : contour_id + 1u;
	ContourFlags cf_next = load_contour_flags(next_contour_id);
	
	// cf.seg_head:= segment head due to cusp variation
	bool is_scan_seg_head = cf.seg_head 	 || is_curve_head;  
	bool is_scan_seg_tail = cf_next.seg_head || is_curve_tail;  

	if (valid_thread)
	{
		if (cf_next.seg_head)
		{
			set_contour_seg_tail(true, cf); 
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

	if (idx.x == 0) 
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
		ssbo_contour_snake_seg_rank_[contour_id] = seg_rank;
		ssbo_contour_snake_seg_len_[contour_id] = seg_len;

		// TODO: Update seg_head/tail for curve head/tail in next pass
		// cannot do it here due to data racing
	}

#endif


}
#endif







#if defined(_KERNEL_MULTICOMPILE__CALC_CONTOUR_EDGES_RENDER_DATA)
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

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours; 

	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 world_to_view = ubo_view_matrices_.viewmat;
	mat4 mat_camera_proj = ubo_view_matrices_.winmat; 

	ContourFlags cf = load_contour_flags(contour_id); 
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf); // TODO: encode tail flag
	uint next_contour_id = move_contour_id_along_loop(cct, contour_id, 1.0f); 
	bool curve_tail = cct.tail_contour_id == contour_id; 
	
	vec3 vpos_0, vpos_1;
	{ // read vertex pos
		uvec3 vpos_0_enc, vpos_1_enc;
		Load3(ssbo_contour_snake_vpos_, contour_id,    	vpos_0_enc);
		Load3(ssbo_contour_snake_vpos_, next_contour_id, vpos_1_enc);
		vpos_0 = uintBitsToFloat(vpos_0_enc);
		vpos_1 = uintBitsToFloat(vpos_1_enc);
	}

#if defined(_KERNEL_MULTICOMPILE__CALC_CONTOUR_EDGES_DRAW_DATA)
	if (valid_thread)
	{ 
		if (curve_tail && !cf.looped_curve) vpos_1 = vpos_0; /* omit at the end of a curve */
		
		vec4 vpos_ws[2]; /* v0, v1 on edge */
		vec3 vnor_ws[2]; 
		vec4 vpos_ndc[2]; 
		vec2 vpos_uv[2]; 

        vpos_ws[0] = vec4(vpos_0.xyz, 1.0f);
		vpos_ndc[0] = mat_camera_proj * vec4((world_to_view * vpos_ws[0]).xyz, 1.0f); 
        vpos_ws[1] = vec4(vpos_1.xyz, 1.0f);
		vpos_ndc[1] = mat_camera_proj * vec4((world_to_view * vpos_ws[1]).xyz, 1.0f); 

		/* write to draw data */
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
#endif
}
#endif






#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE)
/*
	ssbo_bnpr_mesh_pool_counters_
	ssbo_tree_scan_infos_2d_resampler_
	
	// segscan
 	uint ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_[]
	uint ssbo_tree_scan_output_2d_resampler_accumulate_curvlen_[] ( == ssbo_contour_arc_len_param_[])
	// scan
	uint ssbo_tree_scan_input_2d_resampler_alloc_samples_[] (<- exclusive to above 2 ssbos)
	uint ssbo_tree_scan_output_2d_resampler_alloc_samples_[] ( == ssbo_contour_to_start_sample_[])
	// segscan
	uint ssbo_tree_scan_input_2d_resample_contour_idmapping_[] (<- exclusive to above 2 ssbos)
	uint ssbo_tree_scan_output_2d_resample_contour_idmapping_[] (== ssbo_2d_sample_to_contour_[])

	uint ssbo_contour_2d_resample_raster_data_[] 
	uint ssbo_contour_2d_sample_geometry_[]
	vec2 pcs_screen_size_
	float pcs_sample_rate_
*/

vec2 get_raster_resolution()
{
	return pcs_screen_size_.xy / pcs_sample_rate_; 
}

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY)
struct Contour2DSampleSegmentationInfo
{ // 2d samples are segmented if
	// (at different contour curves)
	uint contour_curve_head_id; 
	// || (at different contour segments)
	uint contour_seg_head; 
	// || (arclen is not consecutive) <= due to partially clipped contour curve(s)
	float sample_arc_len_param; 
}
Contour2DSampleSegmentationInfo load_2d_sample_seg_key(uint sample_id, uint contour_id)
{
	Contour2DSampleSegmentationInfo info; 
	
	uint contour_id = ssbo_2d_sample_to_contour_[sample_id]; 
	ContourFlags cf = load_contour_flags(contour_id); 
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf);

	info.contour_curve_head_id = ssbo_contour_snake_list_head_[contour_id]; 

	uint contour_seg_rank = ssbo_contour_snake_seg_rank_[contour_id];
	info.contour_seg_head = move_contour_id_along_loop(cct, contour_id, -float(contour_seg_rank)); 

	info.sample_arc_len_param  = load_ssbo_contour_2d_sample_geometry__curv_arclen_param(sample_id, num_samples);
}
#endif

void main()
{
	const uint idx = gl_GlobalInvocationID.x; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__PREP_ARCLEN_PARAM)
	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours; 

	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 world_to_view = ubo_view_matrices_.viewmat;
	mat4 mat_camera_proj = ubo_view_matrices_.winmat; 

	ContourFlags cf = load_contour_flags(contour_id); 
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf); 
	uint next_contour_id = move_contour_id_along_loop(cct, contour_id, 1.0f); 
	bool non_looped_curve_tail = (cct.tail_contour_id == contour_id) && (!cf.looped_curve); 
	
	vec3 vpos_0, vpos_1;
	{ // read vertex pos
		uvec3 vpos_0_enc, vpos_1_enc;
		Load3(ssbo_contour_snake_vpos_, contour_id,    	vpos_0_enc);
		Load3(ssbo_contour_snake_vpos_, next_contour_id, vpos_1_enc);
		vpos_0 = uintBitsToFloat(vpos_0_enc);
		vpos_1 = uintBitsToFloat(vpos_1_enc);
	}

	const vec2 raster_resolution = get_raster_resolution(); 
	LineRasterResult raster_result = raster_line_segment(
		vec4(vpos_0.xyz, 1.0f), vec4(vpos_1.xyz, 1.0f), 
		raster_resolution, world_to_view, mat_camera_proj
	);

	Contour2DResampleRasterData c2rd; 
	c2rd.begend_uvs 	= raster_result.beg_from_p0 ? raster_result.begend_uvs : raster_result.begend_uvs.zwxy; 
	c2rd.has_samples 	= !(raster_result.is_line_clip_rejected || non_looped_curve_tail) && valid_thread; 
	if (!c2rd.has_samples)
		c2rd.begend_uvs.zw = c2rd.begend_uvs.xy;

	if (valid_thread)
	{
		// accumulate 2d-curve-len for each curve
		uint hf = cct.head_contour_id == contour_id ? 1 : 0; 
		
		float scan_val = length(raster_resolution * (c2rd.begend_uvs.zw - c2rd.begend_uvs.xy)); 
		if (!c2rd.has_samples) scan_val = 0.0f; // play it safe
		
		uvec2 segscan_input_enc = segscan_float_hf_encode(SSBOData_SegScanType_float(scan_val, hf)); 
		Store2(ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_, contour_id, segscan_input_enc);

		// Cache raster result of each contour edge
		uvec3 resample_data = encode_contour_2d_resample_data(c2rd);
		Store3(ssbo_contour_2d_resample_raster_data_, contour_id, resample_data);
	}

	if (idx.x == 0) 
	{
		ssbo_tree_scan_infos_2d_resampler_.num_scan_items = num_contours; 
		ssbo_tree_scan_infos_2d_resampler_.num_valid_scan_threads = compute_num_threads(
			num_contours, 2u
		);
		ssbo_tree_scan_infos_2d_resampler_.num_thread_groups = compute_num_groups(
			num_contours, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
		);

		ssbo_bnpr_mesh_pool_counters_.num_2d_samples = 0u; 
	}

	if (valid_thread)
	{
		vec4 dbg_col = vec4(.0f); 
		ivec2 dbg_pos = ivec2(c2rd.begend_uvs.xy * pcs_screen_size_.xy);
		// imageStore(tex2d_contour_dbg_, dbg_pos, dbg_col);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__ALLOC_SAMPLES)
	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours;
	
	ContourFlags cf = load_contour_flags(contour_id); 
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf); 

	// Arc-length parameterized 3d contours
	float len_2d, arc_len_param; {
		SSBOData_SegScanType_float scan_input = segscan_float_hf_decode(uvec2(
			ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_[contour_id * 2u], 
			ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_[contour_id * 2u + 1u]
		));
		len_2d = scan_input.val; 

		SSBOData_SegScanType_float scan_output = segscan_float_hf_decode(uvec2(
			/* == ssbo_tree_scan_output_2d_resampler_accumulate_curvlen_ */
			ssbo_contour_arc_len_param_[2u * contour_id], 
			ssbo_contour_arc_len_param_[2u * contour_id + 1u]
		));
		arc_len_param = scan_output.val; 
	}

	// Generate 2d sample at integer arc-len 
	uint num_samples = 0; {
		vec2 arc_len_range = vec2(arc_len_param, arc_len_param + len_2d); 
		
		const float range_low  = ceil(arc_len_range.x); 
		const float range_high = floor(arc_len_range.y); 
		float inside_count = max(.0f, range_high - range_low + 1.0f);
		// don't count the last sample if it is at the range boundary
		bool range_high_touch_sample = arc_len_range.y == range_high;
		bool last_edge_in_curve = (!cct.looped_curve) && ((contour_id + 1u) == cct.tail_contour_id);  
		if (range_high_touch_sample && (!last_edge_in_curve)) inside_count -= 1.0f; 
		
		num_samples = uint(inside_count + 1e-10f);
		if (len_2d == 0.0f) num_samples = 0u; 
	}

	if (valid_thread)
	{ // allocate samples for each curve, scan to preserve orders
		uint scan_val = num_samples; 
		ssbo_tree_scan_input_2d_resampler_alloc_samples_[contour_id] = scan_val;

		ssbo_contour_arc_len_param_[2u * contour_id] = floatBitsToUint(arc_len_param);
		ssbo_contour_arc_len_param_[2u * contour_id + 1u] = floatBitsToUint(len_2d);
	}

	if (valid_thread)
	{
		Contour2DResampleRasterData c2rd; 
		{
			uvec3 resample_data; 
			Load3(ssbo_contour_2d_resample_raster_data_, contour_id, resample_data);
			c2rd = decode_contour_2d_resample_data(resample_data);
		}

		vec4 dbg_col = vec4(float(num_samples), .0f, .0f, 1.0f); 
		ivec2 dbg_pos = ivec2(c2rd.begend_uvs.xy * pcs_screen_size_.xy);
		// imageStore(tex2d_contour_dbg_, dbg_pos, dbg_col);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__ALLOC_SAMPLES_FINISH)
	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours;

	uint num_samples = ssbo_tree_scan_input_2d_resampler_alloc_samples_[contour_id]; 

	if (contour_id == num_contours - 1u)
	{
		uint alloc_sample_offset = /* == ssbo_tree_scan_output_2d_resampler_alloc_samples_ */
			ssbo_contour_to_start_sample_[contour_id]; 

		ssbo_bnpr_mesh_pool_counters_.num_2d_samples = alloc_sample_offset + num_samples;
		/* fill dispatch args after this */
	}

	if (valid_thread)
	{
		Contour2DResampleRasterData c2rd; 
		{
			uvec3 resample_data; 
			Load3(ssbo_contour_2d_resample_raster_data_, contour_id, resample_data);
			c2rd = decode_contour_2d_resample_data(resample_data);
		}

		uint alloc_sample_offset = /* == ssbo_tree_scan_output_2d_resampler_alloc_samples_ */
			ssbo_contour_to_start_sample_[contour_id]; 
		vec4 dbg_col = vec4(float(alloc_sample_offset), .0f, .0f, 1.0f); 
		ivec2 dbg_pos = ivec2(c2rd.begend_uvs.xy * pcs_screen_size_.xy);
		// imageStore(tex2d_contour_dbg_, dbg_pos, dbg_col);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__CLEAR_BUFFER)
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	const uint sample_id = idx; 
	bool valid_thread = sample_id < num_samples;

	if (valid_thread)
		ssbo_tree_scan_input_2d_resample_contour_idmapping_[sample_id] = 
			segscan_uint_hf_encode(SSBOData_SegScanType_uint(0x3fffffffu, 0u));
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__SETUP_SEGSCAN)
	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
	const uint contour_id = idx; 
	bool valid_thread = contour_id < num_contours;

	uint num_2d_samples = ssbo_tree_scan_input_2d_resampler_alloc_samples_[contour_id];
	uint alloc_sample_offset = ssbo_contour_to_start_sample_[contour_id]; // == ssbo_tree_scan_output_2d_resampler_alloc_samples_

	if (valid_thread && 0 < num_2d_samples)
	{ // Seed at segment heads
		ssbo_tree_scan_input_2d_resample_contour_idmapping_[alloc_sample_offset] = 
			segscan_uint_hf_encode(SSBOData_SegScanType_uint(contour_id, 1u));
	}

	if (idx == 0u)
	{
		uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples;
		ssbo_tree_scan_infos_2d_resampler_.num_scan_items = num_samples; 
		ssbo_tree_scan_infos_2d_resampler_.num_valid_scan_threads = compute_num_threads(
			num_samples, 2u
		);
		ssbo_tree_scan_infos_2d_resampler_.num_thread_groups = compute_num_groups(
			num_samples, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
		);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__FINISH)
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	const uint sample_id = idx; 
	bool valid_thread = sample_id < num_samples;

	if (valid_thread)
	{
		SSBOData_SegScanType_uint segscan_input = 
			segscan_uint_hf_decode(ssbo_tree_scan_input_2d_resample_contour_idmapping_[sample_id]);
		SSBOData_SegScanType_uint segscan_output = 
			// == ssbo_tree_scan_output_2d_resample_contour_idmapping_[]
			segscan_uint_hf_decode(ssbo_2d_sample_to_contour_[sample_id]); 
		
		segscan_output.val = min(segscan_output.val, segscan_input.val); // exclusive scan into inclusive
		ssbo_2d_sample_to_contour_[sample_id] = (segscan_output.val);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_POSITION)
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	const uint sample_id = idx; 
	bool valid_thread = sample_id < num_samples;
	vec2 raster_resolution = get_raster_resolution(); 

	uint contour_id 			= ssbo_2d_sample_to_contour_[sample_id]; 
	uint beg_sample_id 			= ssbo_contour_to_start_sample_[contour_id]; 
	float contour_arc_len_param = uintBitsToFloat(ssbo_contour_arc_len_param_[contour_id * 2u]); 
	Contour2DResampleRasterData c2rd; 
	{
		uvec3 resample_data; 
		Load3(ssbo_contour_2d_resample_raster_data_, contour_id, resample_data);
		c2rd = decode_contour_2d_resample_data(resample_data);
	}


	// Evaluate Position ---
	// contour edge is parameterized as P = C + (s + k) * V, 
	// C: contour snake 2d pos, 
	vec2 C = c2rd.begend_uvs.xy * raster_resolution;
	// V: normalized screen dir, 
	vec2 V = normalize(c2rd.begend_uvs.zw * raster_resolution - c2rd.begend_uvs.xy * raster_resolution);
	// s: offset parameter because we sample at integer arc-lens
	float s = ceil(contour_arc_len_param) - contour_arc_len_param; 
	// k: this is the kth sample from the contour edge
	float k = float(sample_id - beg_sample_id);

	vec2 P = C + (s + k) * V;
	float sample_arc_len_param = ceil(contour_arc_len_param) + k;

	if (valid_thread)
	{
		store_ssbo_contour_2d_sample_geometry__position(sample_id, P); 
		store_ssbo_contour_2d_sample_geometry__curv_arclen_param(sample_id, sample_arc_len_param, num_samples);
	}	

	if (valid_thread)
	{
		vec4 dbg_col = vec4(1.0f); 
		imageStore(tex2d_contour_dbg_, ivec2(P * pcs_sample_rate_), dbg_col);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY)

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__STEP_0)
	// Prep for evaluating topo
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	const uint sample_id = idx; 
	bool valid_thread = sample_id < num_samples;
	vec2 raster_resolution = get_raster_resolution(); 

	Contour2DSampleSegmentationInfo sample_si = load_2d_sample_seg_key(sample_id, contour_id); 
	Contour2DSampleSegmentationInfo next_sample_si = load_2d_sample_seg_key(sample_id + 1u, contour_id);
	#endif


#endif

}
#endif