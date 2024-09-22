
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_brush_toolbox_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_debug_view_lib.glsl)


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

		ssbo_contour_snake_to_temporal_record_[vtx_addr] = cetd.temporal_rec_id; 
		ssbo_contour_snake_to_object_id_[vtx_addr] = cetd.obj_id;
		
		CuspSegmentDenoiseData seg_denoise_data; // setup segloopconv1d input buffer
		seg_denoise_data.cf = cf; 
		seg_denoise_data.vpos_ws = cetd.vpos_ws[0]; 
		uvec4 seg_denoise_data_enc = encode_cusp_segment_denoise_data(seg_denoise_data); 
		Store4(ssbo_in_segloopconv1d_data_, vtx_addr, seg_denoise_data_enc);
		
		
		if (additional_output_tail_vtx)
		{
        	Store3(ssbo_contour_snake_vpos_, vtx_addr+1u, floatBitsToUint(cetd.vpos_ws[1]));
			
			cf.seg_head = false; 
			init_contour_curve_head(false, cf);
			init_contour_curve_tail(true, cf);
			set_contour_cusp_flags(.0f < cusp_func_v[1], cf); 
			store_contour_flags(vtx_addr+1u, cf);
		
			// TODO: fix this for non-looped curves, for now we're just copying...
			// but this requires 2 rec_ids for both verts in the ContourEdgeTransferData
			ssbo_contour_snake_to_temporal_record_[vtx_addr+1u] = cetd.temporal_rec_id; 
			ssbo_contour_snake_to_object_id_[vtx_addr+1u] = cetd.obj_id;

			// setup segloopconv1d input buffer
			seg_denoise_data.cf = cf; 
			seg_denoise_data.vpos_ws = cetd.vpos_ws[1]; 
			seg_denoise_data_enc = encode_cusp_segment_denoise_data(seg_denoise_data);
			Store4(ssbo_in_segloopconv1d_data_, (vtx_addr+1u), seg_denoise_data_enc);
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
	bool is_curve_head = cct.head_id == contour_id; 
	bool is_curve_tail = cct.tail_id == contour_id; 

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__SETUP)

	uint prev_contour_id = is_curve_head ? cct.tail_id : contour_id - 1u;
	ContourFlags cf_prev = load_contour_flags(prev_contour_id);

	uint next_contour_id = is_curve_tail ? cct.head_id : contour_id + 1u;
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
	ContourFlags cf_curv_head = load_contour_flags(cct.head_id);
	SSBOData_SegScanType_uint scan_res_step_1_headseg = segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(cct.head_id)]);
	ContourFlags cf_curv_tail = load_contour_flags(cct.tail_id);
	SSBOData_SegScanType_uint scan_res_step_0_tailseg = segscan_uint_hf_decode(scan_output_buf_0_[cct.tail_id]); 
	FixLoopedJumps(
		/*inout*/seg_rank, seg_len, 
		cct, seg_head, seg_tail, cf_curv_head, cf_curv_tail, 
		scan_res_step_0_tailseg, scan_res_step_1_headseg
	); 

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
	bool curve_tail = cct.tail_id == contour_id; 
	// TODO: encapsulate vpos load
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

vec2 get_raster_resolution()
{
	return pcs_screen_size_.xy / pcs_sample_rate_; 
}

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__STEP_0)
struct Contour2DSampleSegmentationInfo
{ // 2d samples are segmented if
	// (at different contour curves)
	uint contour_curve_head_id; 
	// || (at different contour segments)
	uint contour_seg_head; 
	// || (contour rank is not consecutive) <= due to partially clipped contour curve(s)
	float contour_seg_rank; 

	vec2 uv; 
	ContourFlags cf; 
}; 
Contour2DSampleSegmentationInfo load_2d_sample_seg_key(uint sample_id, uint num_samples)
{
	Contour2DSampleSegmentationInfo info; 
	
	uint contour_id = ssbo_2d_sample_to_contour_[sample_id]; 
	ContourFlags cf = load_contour_flags(contour_id); 
	ContourCurveTopo cct = load_contour_curve_topo(contour_id, cf);

	info.contour_curve_head_id = ssbo_contour_snake_list_head_[contour_id]; 

	info.contour_seg_rank = ssbo_contour_snake_seg_rank_[contour_id];
	info.contour_seg_head = move_contour_id_along_loop(cct, contour_id, -float(info.contour_seg_rank)); 

	info.uv = load_ssbo_contour_2d_sample_geometry__position(sample_id);
	info.cf = cf; 

	return info; 
}

bool is_2d_sample_curve_head(Contour2DSampleSegmentationInfo si, Contour2DSampleSegmentationInfo si_prev)
{
	return si.contour_curve_head_id != si_prev.contour_curve_head_id; 
}
void is_2d_sample_seg_head(
	Contour2DSampleSegmentationInfo si, Contour2DSampleSegmentationInfo si_prev, 
	out bool seg_head_contour, out bool seg_head_clipped)
{
	bool break_contour_seg = false; 
	if (si.contour_seg_head != si_prev.contour_seg_head) break_contour_seg = true; // different contour segments
	// non-consecutive contour ranks - a segment jump from curve end to start  
	else if (si.contour_seg_rank < si_prev.contour_seg_rank) break_contour_seg = true; 

	float dist = length(get_raster_resolution() * (si.uv - si_prev.uv)); 
	bool break_arc_len = 1.1f < dist && (si.contour_curve_head_id == si_prev.contour_curve_head_id); 

	seg_head_contour = break_contour_seg; 
	seg_head_clipped = break_arc_len; 
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
	bool non_looped_curve_tail = (cct.tail_id == contour_id) && (!cf.looped_curve); 
	
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
	
	bool contour_edge_clipped = raster_result.is_line_clip_rejected || raster_result.is_line_clipped; 
	if (valid_thread && contour_edge_clipped)
	{ // Mark at curve head. race cond' wont hurt here, note in the future you need to clear this flag before reuse it 
		ContourFlags cf_curve_head = load_contour_flags(cct.head_id); 
		set_contour_flags_curve_clipped(true, cf_curve_head); 
		store_contour_flags(cct.head_id, cf_curve_head); 
	}

	if (valid_thread)
	{
		// accumulate 2d-curve-len for each curve
		uint hf = cct.head_id == contour_id ? 1 : 0; 
		
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
		bool last_edge_in_curve = (!cct.looped_curve) && ((contour_id + 1u) == cct.tail_id);  
		if (range_high_touch_sample && (!last_edge_in_curve)) inside_count -= 1.0f; 
		
		num_samples = uint(inside_count + 1e-10f);
		if (len_2d == 0.0f) num_samples = 0u; 
	}

	if (valid_thread)
	{ // allocate samples for each curve, use scan to preserve orders
		uint scan_val = num_samples; 
		ssbo_tree_scan_input_2d_resampler_alloc_samples_[contour_id] = scan_val;

		ssbo_contour_arc_len_param_[2u * contour_id] = floatBitsToUint(arc_len_param);
		ssbo_contour_arc_len_param_[2u * contour_id + 1u] = floatBitsToUint(len_2d);
	}

	if (valid_thread)
	{ // do a little extra here, broadcast curve-clipped flag from the header contour snake
		ContourFlags cf_curve_head = load_contour_flags(cct.head_id); 
		set_contour_flags_curve_clipped(cf_curve_head.curve_clipped, cf); 
		store_contour_flags(contour_id, cf); 
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

	// record total number of allocated 2d samples
	if (contour_id == num_contours - 1u) { 
		uint alloc_sample_offset = /* == ssbo_tree_scan_output_2d_resampler_alloc_samples_ */
			ssbo_contour_to_start_sample_[contour_id]; 
		ssbo_bnpr_mesh_pool_counters_.num_2d_samples = alloc_sample_offset + num_samples;
		/* fill dispatch args after this */
	}

	// debug draw
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
	Contour2DResampleRasterData c2rd; {
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
		store_ssbo_contour_2d_sample_geometry__position(sample_id, P.xy / raster_resolution.xy); 
		store_ssbo_contour_2d_sample_geometry__curv_arclen_param(sample_id, sample_arc_len_param, num_samples);

		ContourFlags cf = load_contour_flags(contour_id); // init sample flags
		cf.curve_head = false;
		cf.curve_tail = false;
		cf.seg_head = false;
		cf.seg_tail = false;
		// other fields keep the same as the parent contour snake
		store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
	}	
#endif

#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY)
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	const uint sample_id = idx; 
	bool valid_thread = sample_id < num_samples;
	vec2 raster_resolution = get_raster_resolution(); 

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__STEP_0) || defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SETMENTATION__PREP_SEGTAILS) || defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__SETUP_SEGSCAN) || defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__FINISH_SEGSCAN)
		// Detect curve topo first, then build sub-segments
		bool segment_by_seg = (0 < pcs_segment_by_seg_); 
	#endif

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__STEP_0)
		// Detect curve head and seg head
		// Let's say we have a contour curve,
		// with 2 segments A and B, and B is partially clipped.(from the view frustum)
		// B is split into B1, BX and B2. 
		//                                                
		//             ___b5b6______b4________b3              
		//            b7         BX            \              
		//           /                          b2            
		//          b8                           |            
		//          |                            | 
		// +======= b9 ======================== b1 =======+           b9                          b1                                          
		// |         \                        b0/         |            \                        b0/                                           
		// |      B1  b10        A         a5/  B2        |             b10        A         a5/                                              
		// |           \a0_a1_a2||a3___a4_/               |  ==>         \a0_a1_a2||a3___a4_/                                                 
		// |                  a3:=curve head              |                     a3:=curve head                                                
		// |                                              |          
		// +================== Viewport ==================+    Mem Layout:a3 a4 a5 b0 b1 b9 b10 a0 a1 a2                                         
		//                                                                        |     |      |        
		// After 2d resampling, 
		// contour edges in BX do not have 2d samples.                                                   
		// Hence the 2d samples are laied out as follows:
		// [samples from edge a3-5, b0-b1, b9- b10, a0-a2] (in data buffers, the whole curve starts at a3) 
		// curve head sample in a3,   detected by comparing the contour-curve-head-id with prev sample in memory
		// seg head samples in a0,b0, detected by comparing the contour-seg-head-id with prev sample in memory
		// seg head sample  in b9,    detected by comparing the arc-len-param with prev sample(b1) in memory
		// tails in a2, a5, b10, b1 can be inferred from heads in a3, b0, a0, b9
		//
		// The clipping complicated the topology of resampled curves:  
		// - Breaks can happen within adjacent segments (gap from b1 to b9)
		uint contour_id = ssbo_2d_sample_to_contour_[sample_id]; 
		Contour2DSampleSegmentationInfo sample_si = load_2d_sample_seg_key(sample_id, num_samples); 
		ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id); 

		if (!segment_by_seg)
		{
			uint prev_sample_id = sample_id == 0u ? 0u : sample_id - 1u; 
			Contour2DSampleSegmentationInfo prev_sample_si = load_2d_sample_seg_key(prev_sample_id, num_samples);

			bool is_curve_head = sample_id == 0 || is_2d_sample_curve_head(sample_si, prev_sample_si); 
			if (valid_thread)
			{
				cf.curve_head = is_curve_head; 
				store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
			}
		} else {
			// curve segments are built so that we can divide them to build the sub-segments
			ContourCurveTopo cct = load_contour_2d_sample_curve_topo(sample_id, cf, num_samples);

			uint prev_sample_id = move_contour_id_along_loop(cct, sample_id, -1.0f);
			Contour2DSampleSegmentationInfo prev_sample_si = load_2d_sample_seg_key(prev_sample_id, num_samples);

			bool seg_head_contour = false; bool seg_head_clipped = false; 
			is_2d_sample_seg_head(
				sample_si, prev_sample_si, 
				/*out*/seg_head_contour, seg_head_clipped
			);
			
			bool is_seg_entire_curve = false; 
			uint contour_seg_len = ssbo_contour_snake_seg_len_[contour_id];
			uint contour_curve_len = ssbo_contour_snake_list_len_[contour_id]; 
			is_seg_entire_curve = cf.looped_curve && contour_seg_len == contour_curve_len; 

			if (valid_thread)
			{
				cf.seg_head = seg_head_contour || seg_head_clipped; 
				init_seg_head_contour(seg_head_contour/* is_seg_head */, cf); 
				init_seg_head_clipped(seg_head_clipped, cf); 
				init_no_segmentation_on_contour_curve(is_seg_entire_curve, cf);
				store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
			}

			// Debug draw 
			if (valid_thread)
			{
				vec2 dbg_pix = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
				vec4 dbg_col = vec4(.0f, .0f, 1.0f, 1.0f); 
				if (seg_head_contour) {
					dbg_col.r = sample_si.contour_seg_head; 
					dbg_col.g = prev_sample_si.contour_seg_head; 
					dbg_col.b = sample_si.contour_seg_rank; 
				}
				// imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
			}
		}

		if (valid_thread && !segment_by_seg)
		{ // in the initial segmentation pass we also store the sample's pointer to the scene object
			uint object_id = ssbo_contour_snake_to_object_id_[contour_id];
			store_ssbo_contour_2d_sample_topology__object_id(sample_id, object_id, num_samples); 
		}
		
		if (sample_id == 0) // Prepare for segmented convolution
			ssbo_segloopconv1d_info_.num_conv_items = num_samples; 
	#endif

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SETMENTATION__PREP_SEGTAILS)
		// Given cf.seg_head, detect curve tail and seg tail
		ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id);
		
		if (!segment_by_seg)
		{
			uint next_sample_id = sample_id == num_samples - 1u ? num_samples - 1u : sample_id + 1u;
			ContourFlags cf_next = load_ssbo_contour_2d_sample_topology__flags(next_sample_id); 

			bool is_curve_tail = cf_next.curve_head || (sample_id == num_samples - 1u);
			if (valid_thread)
			{
				cf.curve_tail = is_curve_tail; 
				store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
			}
		}
		else
		{ // Default path unless we are initializing the 2d curve topo
			ContourCurveTopo cct = load_contour_2d_sample_curve_topo(sample_id, cf, num_samples); 
			uint next_sample_id = move_contour_id_along_loop(cct, sample_id, +1.0f); 
			ContourFlags cf_next = load_ssbo_contour_2d_sample_topology__flags(next_sample_id); 

			bool is_seg_tail   = cf_next.seg_head; 
			if (valid_thread)
			{
				cf.seg_tail = is_seg_tail; 
				set_seg_tail_clipped(cf_next.seg_head_clipped, cf);
				store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
			}
		}
		
		if (valid_thread)
		{
			vec2 dbg_pix = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
			vec4 dbg_col = vec4(1.0f); 
			// imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
		}
	#endif

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__SETUP_SEGSCAN)
		// Setup segscan inputs
		if (valid_thread)
		{
			ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id);

			bool segment_by_seg = (0 < pcs_segment_by_seg_); 
			uint scan_item_id = sample_id; 
			uint scan_val = 1u; 

			uint hf = (segment_by_seg ? (cf.seg_head || cf.curve_head) : cf.curve_head) ? 1u : 0u; 
			ssbo_tree_scan_input_2d_sample_segmentation_0_[scan_item_id] = 
				segscan_uint_hf_encode(SSBOData_SegScanType_uint(scan_val, hf)); 

			scan_item_id = num_samples - 1u - sample_id; 
			hf = (segment_by_seg ? (cf.seg_tail || cf.curve_tail) : cf.curve_tail) ? 1u : 0u; 
			ssbo_tree_scan_input_2d_sample_segmentation_1_[scan_item_id] = 
				segscan_uint_hf_encode(SSBOData_SegScanType_uint(scan_val, hf)); 
		}

		if (idx.x == 0) 
		{
			ssbo_tree_scan_infos_2d_resampler_.num_scan_items = num_samples; 
			ssbo_tree_scan_infos_2d_resampler_.num_valid_scan_threads = compute_num_threads(
				num_samples, 2u
			);
			ssbo_tree_scan_infos_2d_resampler_.num_thread_groups = compute_num_groups(
				num_samples, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
			);
		}

		if (valid_thread)
		{
			vec2 dbg_pix = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
			vec4 dbg_col = vec4(1.0f); 
			// imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
		}
	#endif

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__FINISH_SEGSCAN)
		#define REVERSE_ID(id) ((num_samples - 1u - id))
		bool segment_by_seg_ = (0 < pcs_segment_by_seg_); 

		SSBOData_SegScanType_uint scan_res_step_0 = 
			segscan_uint_hf_decode(ssbo_tree_scan_output_2d_sample_segmentation_0_[sample_id]);
		SSBOData_SegScanType_uint scan_res_step_1 = 
			segscan_uint_hf_decode(ssbo_tree_scan_output_2d_sample_segmentation_1_[REVERSE_ID(sample_id)]);

		uint dist_to_head = scan_res_step_0.val; 
		uint dist_to_tail = scan_res_step_1.val; 

		uint scanseg_head = sample_id - dist_to_head; 
		uint scanseg_tail = sample_id + dist_to_tail;
		uint scanseg_rank = dist_to_head; 
		uint scanseg_len = dist_to_head + dist_to_tail + 1u; 


		if (!segment_by_seg_)
		{ // segment by curves
			if (valid_thread)
			{
				// Update seg_rank & seg_len
				store_ssbo_contour_2d_sample_topology__curve_rank(sample_id, scanseg_rank, num_samples);
				store_ssbo_contour_2d_sample_topology__curve_len(sample_id,  scanseg_len,  num_samples);
			}
		}else{
			// segment by sub-segs within each curve
			// special care for looped curves
			ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id);
			ContourCurveTopo cct = load_contour_2d_sample_curve_topo(sample_id, cf, num_samples);

			ContourFlags cf_curv_head = load_ssbo_contour_2d_sample_topology__flags(cct.head_id);
			SSBOData_SegScanType_uint scan_res_step_1_headseg = segscan_uint_hf_decode(
				ssbo_tree_scan_output_2d_sample_segmentation_1_[REVERSE_ID(cct.head_id)]
			);
			ContourFlags cf_curv_tail = load_ssbo_contour_2d_sample_topology__flags(cct.tail_id);
			SSBOData_SegScanType_uint scan_res_step_0_tailseg = segscan_uint_hf_decode(
				ssbo_tree_scan_output_2d_sample_segmentation_0_[cct.tail_id]
			); 
			FixLoopedJumps(
				/*inout*/scanseg_rank, scanseg_len, 
				cct, scanseg_head, scanseg_tail, cf_curv_head, cf_curv_tail, 
				scan_res_step_0_tailseg, scan_res_step_1_headseg
			);

			if (valid_thread)
			{
				// Update seg_rank & seg_len
				store_ssbo_contour_2d_sample_topology__seg_rank(sample_id, scanseg_rank, num_samples);
				store_ssbo_contour_2d_sample_topology__seg_len(sample_id,  scanseg_len,  num_samples);
			}

			if (valid_thread) // debug draw
			{
				ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id);
				uint curve_len = load_ssbo_contour_2d_sample_topology__curve_len(sample_id, num_samples); 
				bool is_looped_samples = is_2d_sample_curve_looped(
					cf.looped_curve && cf.no_segmentation_on_contour_curve, 
					cf.curve_clipped, scanseg_len == curve_len
				); 
				
				vec4 dbg_col = cf.curve_clipped ? vec4(1, 0, 0, 1) : vec4(1.0f); 
				dbg_col.rgb = rand_col_rgb(scanseg_len, scanseg_len + 17); 
				if (!is_looped_samples) dbg_col = vec4(.0f); 

				vec2 dbg_pix = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
				// imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
			}
		}
	#endif

	#if defined(_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__REMOVE_FAKE_CORNERS)
		#define PI 3.1415926535614f
		ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id); 

		if (cf.is_corner) 
		{ // Detect fake corners
			ContourCurveTopo cct = load_contour_2d_sample_curve_topo(sample_id, cf, num_samples); 
			uint corner_seg_rank = load_ssbo_contour_2d_sample_topology__seg_rank(sample_id, num_samples); 
			uint corner_seg_len  = load_ssbo_contour_2d_sample_topology__seg_len(sample_id, num_samples); 

			uint seg_tail_id = move_contour_id_along_loop(cct, sample_id, +float(corner_seg_len - 1u - corner_seg_rank)); 
			uint seg_head_id = move_contour_id_along_loop(cct, sample_id, -float(corner_seg_rank));
			// Note: need to consider clipped gaps between adjacent segments!!! 
			vec2 p = load_ssbo_contour_2d_sample_geometry__position(sample_id); 
			vec2 pp = load_ssbo_contour_2d_sample_geometry__position(seg_head_id); 
			vec2 pn = load_ssbo_contour_2d_sample_geometry__position(seg_tail_id); 

			vec2 vp = pcs_screen_size_.xy * (pp - p); 
			vec2 vpdir = normalize(vp); 
			vec2 vn = pcs_screen_size_.xy * (pn - p); 
			vec2 vndir = normalize(vn); 

			float angle = acos(dot(vpdir, vndir));
			bool fake_corner = (angle > PI * .4f); 

			if (valid_thread)
			{
				float angle_degree = angle * 180.0f / PI; 
				vec2 dbg_pix = pcs_screen_size_.xy * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
				vec4 dbg_col = vec4(.0f); 
				dbg_col.r = fake_corner ? 1.0f : .0f;
				dbg_col.g =  .0f;   
				imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
			}
			if (fake_corner)
				cf.is_corner = false; 
		}

		cf.seg_head = /* cf.is_corner ||  */cf.seg_head_contour || cf.seg_head_clipped; 
		if (valid_thread)
			store_ssbo_contour_2d_sample_topology__flags(sample_id, cf); 
	#endif
#endif

}
#endif





#if defined(_KERNEL_MULTICOMPILE__CALC_CONTOUR_2D_STROKE_RENDER_DATA)
/*
vec2 pcs_screen_size_
float pcs_stroke_width_
ssbo_contour_2d_sample_geometry
ssbo_contour_2d_sample_topology
*/
void main()
{
	const uint sample_id = gl_GlobalInvocationID.x; 
	const uint num_samples = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 
	bool valid_thread = sample_id < num_samples; 
	
	// Load object info
	const uint object_id = load_ssbo_contour_2d_sample_topology__object_id(sample_id, num_samples); 
	StrokegenObjectInfo object_info = load_strokegen_object_info(object_id); 

	// Load topology
	ContourFlags cf = load_ssbo_contour_2d_sample_topology__flags(sample_id); 
	ContourCurveTopo cct = load_contour_2d_sample_curve_topo(sample_id, cf, num_samples); 
	uint seg_rank = load_ssbo_contour_2d_sample_topology__seg_rank(sample_id, num_samples);  
	uint seg_len = load_ssbo_contour_2d_sample_topology__seg_len(sample_id, num_samples); 

	bool is_looped_samples = is_2d_sample_curve_looped(
		cf.looped_curve && cf.no_segmentation_on_contour_curve, 
		cf.curve_clipped, seg_len == cct.len
	); 
	float arc_len_param = float(seg_rank) / float(seg_len); 

	// Load Geometry
	vec2 pos = vec2(pcs_screen_size_.xy) * load_ssbo_contour_2d_sample_geometry__position(sample_id); 
	vec2 tangent = load_ssbo_contour_2d_sample_geometry__tangent(sample_id, num_samples); 

	// Calc stroke width
	float stylized_width_arc_len = 1.0f - abs(arc_len_param - .5f) * 2.0f; 
	if (is_looped_samples) stylized_width_arc_len = 0.7f; 
	
	float stylized_width_stk_len = 1.0f; 
	if (seg_len < 24u) stylized_width_stk_len = float(seg_len) / 24.0f; 
	
	float object_stk_width = object_info.stroke_width; 

	float stk_width = pcs_stroke_width_ * stylized_width_arc_len * stylized_width_stk_len * object_stk_width; 

	// Calc stroke mesh
	mat3x2 verts = compute_wing_quad_verts(
		pos, pcs_screen_size_.xy, 
		vec2(tangent.y, -tangent.x), stk_width
	);
	for (uint i = 0; i < 3; i++)
		verts[i] /= vec2(pcs_screen_size_.xy); // back to 01-uv space


	vec4 col = vec4(.0f, .0f, .0f, 1.0f); 
	// col.rgb = rand_col_rgb(seg_len, seg_len); 
	if (cf.occluded/* occluded_filtered */) col.a = .0f; 
	// if (is_looped_samples) col = vec4(1.0f, 0.0f, 0.0f, 1.0f); 

	if (valid_thread) 
	{
		store_ssbo_stroke_mesh_pool__skeletal_VB(sample_id, verts);
		store_ssbo_stroke_mesh_pool__skeletal_color(sample_id, col, num_samples);
		SpineTopoInfo sti; 
		sti.looped = is_looped_samples && !cf.seg_tail_clipped; 
		sti.is_tail_spine = cf.curve_tail || cf.seg_tail_clipped; 
		sti.head_sample_id = cct.head_id; 
		store_ssbo_stroke_mesh_pool__skeletal_topo_info(sample_id, sti, num_samples); 
	} 

	if (valid_thread)
	{
		vec2 dbg_pix = pcs_screen_size_.xy * verts[0]; 
		vec4 dbg_col = vec4(1.0f, 1.0f, 1.0f, 1.0f);  
		imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
		
		dbg_pix = pcs_screen_size_.xy * verts[1]; 
		dbg_col = vec4(.0f, 1.0f, .0f, 1.0f);  
		imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);

		dbg_pix = pcs_screen_size_.xy * verts[2]; 
		dbg_col = vec4(1.0f, 0.0f, 1.0f, 1.0f);  
		imageStore(tex2d_contour_dbg_, ivec2(dbg_pix), dbg_col);
	}
}
#endif


