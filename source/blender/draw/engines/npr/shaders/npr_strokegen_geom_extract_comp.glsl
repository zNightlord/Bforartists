
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_allocation_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_brush_toolbox_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)

 
/* all counters are cleared in _KERNEL_MULTICOMPILE__GEOM_EXTRACT kernel ------------------- */
#if defined(_KERNEL_MULTICOMPILE_BOOSTRAP_GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
	#define GLOBAL_COMPACTION_COUNTER__MESH_CONTOUR_EDGES ssbo_bnpr_mesh_pool_counters_.num_contour_edges
	#define GLOBAL_COMPACTION_COUNTER__MESH_VERTS ssbo_bnpr_mesh_pool_counters_.num_verts
	#define GLOBAL_COMPACTION_COUNTER__MESH_EDGES ssbo_bnpr_mesh_pool_counters_.num_edges
	#define GLOBAL_COMPACTION_COUNTER__MESH_FACES ssbo_bnpr_mesh_pool_counters_.num_faces
#endif



#if defined(_KERNEL_MULTICOMPILE_BOOSTRAP_GEOM_EXTRACT)
void main()
{
	/* bootstrapping pass */
	if (gl_GlobalInvocationID.x == 0u)
	{
		ssbo_bnpr_mesh_pool_counters_.num_contour_edges 		 = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_verts         		 = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_edges         		 = 0;
		ssbo_bnpr_mesh_pool_counters_.num_faces         		 = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_contour_verts = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_filtered_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_filtered_verts = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_dbg_vnor_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_dbg_edge_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_draw_faces = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_frags = 0;
		
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_verts         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_edges         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_faces         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_verts = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_verts = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_vnor_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_general_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_edge_lines = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_draw_faces = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_frags = 0;
	}
}
#endif


#if defined(_KERNEL_MULTICOMPILE__COPY_VBO)
 void main()
 {
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 

	const uint VertID = idx; 
	const uint NumVerts = pcs_meshbatch_num_verts_;
	bool valid_thread = VertID < NumVerts;
	if (!valid_thread) return; 

	const uint ResourceID = pcs_rsc_handle_; 

	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 model_to_world = drw_matrix_buf[ResourceID].model;  

	vec3 vpos_ls, vpos_ws; 

	uint ld_base_addr = VertID * 3u; 
	vpos_ls.x = ssbo_meshbatch_vbo_[ld_base_addr + 0];
	vpos_ls.y = ssbo_meshbatch_vbo_[ld_base_addr + 1]; 
	vpos_ls.z = ssbo_meshbatch_vbo_[ld_base_addr + 2];

	vpos_ws = (model_to_world * vec4(vpos_ls, 1.0f)).xyz; 

	uint st_base_addr = (VertID) * 3u;  
	ssbo_vbo_full_[st_base_addr+0] = vpos_ws.x; 
	ssbo_vbo_full_[st_base_addr+1] = vpos_ws.y;
	ssbo_vbo_full_[st_base_addr+2] = vpos_ws.z;

	if (idx == 0)
	{ /* cache counters until last extracted mesh */
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges   = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_verts           = ssbo_bnpr_mesh_pool_counters_.num_verts; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_edges           = ssbo_bnpr_mesh_pool_counters_.num_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_faces           = ssbo_bnpr_mesh_pool_counters_.num_faces; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_verts   = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_edges  = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_verts  = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_vnor_lines  = ssbo_bnpr_mesh_pool_counters_.num_dbg_vnor_lines;
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_general_lines = ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines;
		ssbo_bnpr_mesh_pool_counters_prev_.num_dbg_edge_lines  = ssbo_bnpr_mesh_pool_counters_.num_dbg_edge_lines;
		ssbo_bnpr_mesh_pool_counters_prev_.num_draw_faces 	   = ssbo_bnpr_mesh_pool_counters_.num_draw_faces;
		ssbo_bnpr_mesh_pool_counters_prev_.num_frags 		   = ssbo_bnpr_mesh_pool_counters_.num_frags;

		/* reset for each mesh */
		ssbo_bnpr_mesh_pool_counters_.num_draw_faces = 0u; 
	}
 }
#endif


#if defined(_KERNEL_MULTICOMPILE__COPY_EDGE_ADJ_IBO)
 void main()
 {
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    bool valid_thread = (idx < pcs_edge_count_); 
	if (!valid_thread) return; 

    const uint EdgeID = idx.x; 
	uvec4 vid; 
    load_and_decode_ibo__edge_adj(EdgeID, 4u, /*out*/vid); 

	uint st_base_addr = (EdgeID) * 4u;  
	ssbo_edge_to_vert_[st_base_addr+0] = vid[0];
	ssbo_edge_to_vert_[st_base_addr+1] = vid[1]; 
	ssbo_edge_to_vert_[st_base_addr+2] = vid[2];
	ssbo_edge_to_vert_[st_base_addr+3] = vid[3];
 }
#endif






#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT)

vec3 ld_vbo(uint vert)
{
	uint base_addr = vert * 3; /* vbo is aligned to vec4 per vert pos */
	return vec3(buf_vbo[base_addr], buf_vbo[base_addr+1], buf_vbo[base_addr+2]); 
}

uint get_num_edges()
{
	return pcs_num_edges + ssbo_dyn_mesh_counters_.num_edges; 
}

vec2 ndc_to_screen_uv(vec4 pos_ndc)
{
	pos_ndc.xyz /= pos_ndc.w;
	pos_ndc.xy = pos_ndc.xy * 0.5f + 0.5f; 
	return pos_ndc.xy * pcs_screen_size_.xy; 
}


void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint idx 	   = gl_GlobalInvocationID.x;
	
	const uint NumEdges = get_num_edges();
	bool valid_thread = idx < NumEdges; 
	const uint wedge_id = idx;
	const uint resource_id = pcs_rsc_handle; 

	EdgeFlags ef = load_edge_flags(wedge_id); 
	AdjWedgeInfo w[4] = {
		decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 0u]),
		decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 1u]),
		decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 2u]),
		decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 3u])
	}; 
	
	
	/* Extract Contour Edges -------------------------------------------------------------- */
	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 model_to_world = drw_matrix_buf[resource_id].model; 
	mat4 world_to_model = drw_matrix_buf[resource_id].model_inverse; 
	mat4 view_to_world = ubo_view_matrices_.viewinv;
	mat4 world_to_view = ubo_view_matrices_.viewmat;
	mat4 mat_camera_proj = ubo_view_matrices_.winmat; 
	
	bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

	uvec4 vids; 
	for (uint i = 0; i < 4; ++i)
		vids[i] = buf_ibo[4 * wedge_id + i]; 
	
	vec3 v[4]; 
	for (uint i = 0; i < 4; ++i)
		v[i] = ld_vbo(vids[i]);  

	VertFlags vf[4]; 
	for (uint i = 0; i < 4; ++i)
		vf[i] = load_vert_flags(vids[i]);
	
	bool is_contour = false; 
	float face_orient_123, face_orient_013; 
	if (0 == pcs_extract_interpo_contour_)
	{
		is_contour = is_contour_edge(
			v[0], v[1], v[2], v[3], cam_pos_ws, ef, 
			/*out*/face_orient_123, face_orient_013
		);
	}else
		is_contour = is_interp_contour_edge__after_tessellation(
			vf[0], vf[1], vf[2], vf[3], ef, 
			/*out*/ face_orient_123, face_orient_013
		);

	if (false == valid_thread) is_contour = false; 

	uint compacted_idx = compact_contour_edge(is_contour, groupId); 
	uint local_contour_edge_id = calc_local_contour_edge_id(compacted_idx);
	
	barrier();

	uint ifrontface = is_back_face(face_orient_123) ? 1 : 0; 
	if (is_contour)
	{ /* write world pos to output buffer */
		/* Note: wpos_and_edgeid will be overwrite in the next pass, for saving space */
		uint base_addr = local_contour_edge_id * 2u; 
		ssbo_contour_temp_data_[base_addr+0] = wedge_id; 
		
		PerContourWedgeInfo pcwi; 
		pcwi.is_border = ef.border; 
		pcwi.wedge_id = wedge_id;
		pcwi.ifrontface = ifrontface; 
		ssbo_contour_temp_data_[base_addr+1] = encode_per_contour_wedge_info(pcwi); 
	}

	if (valid_thread)
	{
		PerWedgeContourInfo peci; 
		peci.is_border = ef.border; 
		peci.is_contour = is_contour;
		peci.contour_id = is_contour ? local_contour_edge_id : NULL_EDGE; 
		peci.ifrontface = ifrontface; 
		ssbo_edge_to_contour_[wedge_id] = encode_per_wedge_contour_info(peci); 
	}



	/* Extract Triangle Index List -------------------------------------------------------------- */
    bvec2 gen_face = bvec2(valid_thread, valid_thread);
	for (uint iface = 0; iface < 2; ++iface)
	{
		uvec2 iwedges_fi = (iface == 0u) ? uvec2(1, 2) : uvec2(3, 0); 
		for (uint iiw = 0u; iiw < 2u; ++iiw)
		{ // Only wedge with the highest id generates the face 
			uint iwedge = iwedges_fi[iiw]; 
			if (wedge_id < w[iwedge].wedge_id)
				gen_face[iface] = false;
		}
	}

	uint num_gen_faces = uint(gen_face.x) + uint(gen_face.y); 
	if (!valid_thread) num_gen_faces = 0; 
	
	uint faces_offset = alloc_draw_face(groupId, num_gen_faces); 

	barrier(); 

	uint st_face_addr = 4u * faces_offset;
	for (uint iface = 0; iface < 2; ++iface)
		if (gen_face[iface])
		{
			uvec3 iverts_face = mark__face_to_winded_verts(iface); 

			ssbo_face_to_vert_draw_depth_[st_face_addr + 0u] = vids[iverts_face[0]];
			ssbo_face_to_vert_draw_depth_[st_face_addr + 1u] = vids[iverts_face[1]];
			ssbo_face_to_vert_draw_depth_[st_face_addr + 2u] = vids[iverts_face[2]];
			ssbo_face_to_vert_draw_depth_[st_face_addr + 3u] = wedge_id;

			st_face_addr += 4u; 
		}



	/* debug view */
	if (pcs_edge_visualize_mode_ > 0)
	{
		// Check topo consistency
		bool valid_ee, valid_ev, valid_ve;
		validate_wedge_topo(wedge_id, /*out*/ valid_ee, valid_ev, valid_ve);  
		
		bool dbg_line = valid_thread && (!ef.del_by_collapse) && (!ef.dupli) && (!ef.del_by_split); 
		vec3 dbg_col = vec3(1.0f); 

		/* visualize edges with invalid topology */
		if (pcs_edge_visualize_mode_ == 1) 
		{
			uint w0 = ((w[0].wedge_id)); 
			uint w1 = ((w[1].wedge_id)); 
			uint w2 = ((w[2].wedge_id)); 
			uint w3 = ((w[3].wedge_id)); 
		
			if (/* w0 == w2 || w1 == w3 ||  */w0 == w3 || w1 == w2) 
				valid_ee = false; 
			else valid_ee = true; 

			dbg_line = dbg_line && (!valid_ee); 
		}
		if (pcs_edge_visualize_mode_ == 2) 
			dbg_line = dbg_line && (!valid_ev); 
		if (pcs_edge_visualize_mode_ == 3) 
			dbg_line = dbg_line && (!valid_ve); 
		if (pcs_edge_visualize_mode_ == 4) 
			dbg_line = dbg_line && (ef.selected); 

		/* visualize selected edges */ 
		if (pcs_edge_visualize_mode_ == 5) 
			dbg_line = valid_thread && (!ef.dupli) && (ef.selected) && ef.temp_dbg_draw_edge; 

		if (pcs_edge_visualize_mode_ == 6) 
		{
			dbg_line = dbg_line && (vf[1].contour && vf[3].contour); 
		}

		if (pcs_edge_visualize_mode_ == 7)
		{
			bool sel_border = ef.sel_border;
			dbg_line = dbg_line && sel_border; 
		}

		if (pcs_edge_visualize_mode_ == 8)
		{
			bool border_from_vtx_flags = vf[1].border_eval && vf[3].border_eval; 
			dbg_line = dbg_line && border_from_vtx_flags; 
		}

		if (pcs_edge_visualize_mode_ == 9)
		{
			dbg_line = dbg_line && (0 < ef.crease_level); 
		}

		if (10 <= pcs_edge_visualize_mode_ && pcs_edge_visualize_mode_ < 14)
		{
			LoopSubdEdgeTreeUpNode node = decode_loop_subd_tree_node(ssbo_subd_edge_tree_node_up_[wedge_id]); 

			uint match_code = pcs_edge_visualize_mode_ - 10; // 0, 1, 2, 3
			dbg_line = valid_thread && (!ef.dupli) && (!ef.del_by_split) && (ef.selected); 
			// dbg_line = dbg_line && (match_code == node.code); 
			// dbg_line = dbg_line && ((match_code / 2u) == (node.code / 2u)); 

			dbg_col = node.code == 0 ? vec3(1.0f, 0.0f, 0.0f) : 
				node.code == 1 ? vec3(0.0f, 1.0f, 0.0f) : 
				node.code == 2 ? vec3(0.0f, .5f, 1.0f) : 
				node.code == 3 ? vec3(1.0f, 1.0f, 0.0f) : vec3(0.0f);
			dbg_col *= .5f; 
		}


		uint dbg_line_idx = compact_dbg_edge(dbg_line, groupId); 
		dbg_line_idx += get_debug_line_offset(DBG_LINE_TYPE__EDGES); 
		if (dbg_line)
		{
			vec3 dbg_vpos_0 = v[1]; 
			vec3 dbg_vpos_1 = v[3]; 
		 
			DebugVertData dvd_0 = DebugVertData(v[1], dbg_col); 
			DebugVertData dvd_1 = DebugVertData(v[3], dbg_col);
			store_debug_line_data(dbg_line_idx, dvd_0, dvd_1); 
		}
	}

}
#endif




#if defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
	
	#if defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS__CONTOUR_EDGES)
	if (idx == 0u)
	{
		const uint num_draws = ssbo_bnpr_mesh_pool_counters_.num_contour_verts; // last edge of non-loop curve has zero length
		ssbo_bnpr_contour_mesh_draw_args_.vertex_len 		= 2u * num_draws;   		/*#verts*/
		ssbo_bnpr_contour_mesh_draw_args_.instance_len 	= 1;  						/*#instances*/
		ssbo_bnpr_contour_mesh_draw_args_.vertex_first 	= 0;  						/*ibo offset*/
		ssbo_bnpr_contour_mesh_draw_args_.base_index 		= 0;  						/*vbo offset*/
		ssbo_bnpr_contour_mesh_draw_args_.instance_first_indexed = 0; 					/*instance offset*/
	}
	#endif

	#if defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS__CONTOUR_SAMPLES)
	if (idx == 0u)
	{
		const uint num_wings = ssbo_bnpr_mesh_pool_counters_.num_2d_samples; 		
		ssbo_bnpr_2d_sample_draw_args_.vertex_len 		= POINTS_PER_WING_QUAD * num_wings;  /*#verts*/
		ssbo_bnpr_2d_sample_draw_args_.instance_len 	= 1;  						/*#instances*/
		ssbo_bnpr_2d_sample_draw_args_.vertex_first 	= 0;  						/*ibo offset*/
		ssbo_bnpr_2d_sample_draw_args_.base_index 		= 0;  						/*vbo offset*/
		ssbo_bnpr_2d_sample_draw_args_.instance_first_indexed = 0; 					/*instance offset*/
	}
	#endif

	#if defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS_DEPTH)
	if (idx == 0u)
	{
		const uint num_draws = ssbo_bnpr_mesh_pool_counters_.num_draw_faces; // last edge of non-loop curve has zero length
		ssbo_bnpr_contour_mesh_draw_args_.vertex_len 		= 3u * num_draws;   		/*#verts*/
		ssbo_bnpr_contour_mesh_draw_args_.instance_len 	= 1;  						/*#instances*/
		ssbo_bnpr_contour_mesh_draw_args_.vertex_first 	= 0;  						/*ibo offset*/
		ssbo_bnpr_contour_mesh_draw_args_.base_index 		= 0;  						/*vbo offset*/
		ssbo_bnpr_contour_mesh_draw_args_.instance_first_indexed = 0; 					/*instance offset*/
	}
	#endif
}
#endif





#if defined(_KERNEL_MULTICOMPILE__EXTRACT_MESH_CONTOUR_DATA)
vec2 ndc_to_screen_uv(vec4 pos_ndc)
{
	pos_ndc.xyz /= pos_ndc.w;
	pos_ndc.xy = pos_ndc.xy * 0.5f + 0.5f; 
	return pos_ndc.xy * pcs_screen_size_.xy; /* TODO: flip y or not? */
}

struct VECircContext_ContourLinking
{
	PerWedgeContourInfo pwci; /* pwci of current wedge ("wi") */
}; 
bool func_ve_circulator(CirculatorIterData iter, inout VECircContext_ContourLinking ctx)
{
	ctx.pwci = decode_per_wedge_contour_info(ssbo_edge_to_contour_[iter.awi_next.wedge_id]);

	if (ctx.pwci.is_contour) return false; /* found a contour, quit circulation */
	return true; 
}

void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 

	const uint NumContourEdgesCurr = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 

	const uint contour_id_local = idx; 
	const uint contour_id_global = calc_global_contour_edge_id(contour_id_local); 
	bool valid_thread = contour_id_global < NumContourEdgesCurr; 

	PerContourWedgeInfo pcwi; 
	uint v0, v1; 
	vec4 vpos_ws[2]; 
	{
		uint base_addr = contour_id_local * 2u; 
		
		uint wedge_id = ssbo_contour_temp_data_[base_addr+0];
		pcwi = decode_per_contour_wedge_info(ssbo_contour_temp_data_[base_addr+1]);  
		
		uvec2 iverts_frontface = mark__cwedge_to_verts(pcwi.ifrontface); 
		v0 = ssbo_edge_to_vert_[wedge_id*4 + iverts_frontface[0]];
		vpos_ws[0].x = (ssbo_vbo_full_[v0*3+0]); 
		vpos_ws[0].y = (ssbo_vbo_full_[v0*3+1]); 
		vpos_ws[0].z = (ssbo_vbo_full_[v0*3+2]); 
		vpos_ws[0].w = 1.0f; 

		v1 = ssbo_edge_to_vert_[wedge_id*4 + iverts_frontface[1]];
		vpos_ws[1].x = (ssbo_vbo_full_[v1*3+0]);
		vpos_ws[1].y = (ssbo_vbo_full_[v1*3+1]);
		vpos_ws[1].z = (ssbo_vbo_full_[v1*3+2]); 
		vpos_ws[1].w = 1.0f; 
	}

	/* Generate Contour Data ------------------------------------------------------------- */
	/* use local id to load temp data from prev pass, use global id to store */
	if (valid_thread)
	{ /* Note: wpos_and_edgeid will be overwrite, for saving space */
		{ /* Generate Geometry */ 
			vec2 maxcurv; vec2 cusp_func; 
			ld_vcurv_max_with_cusp(v0, /*out*/maxcurv[0], cusp_func[0]);
			ld_vcurv_max_with_cusp(v1, /*out*/maxcurv[1], cusp_func[1]);
			bool seg_head = sign(cusp_func[0]) != sign(cusp_func[1]); 
			ContourFlags cf = init_contour_flags(seg_head);

			ContourEdgeTransferData cetd; 
			cetd.vpos_ws[0] = vpos_ws[0].xyz;
			cetd.vpos_ws[1] = vpos_ws[1].xyz;
			cetd.cf = cf;
			cetd.cusp_funcs = cusp_func;

			{ /* Link contour edge to temporal record */
				uint old_edge_0 = ssbo_contour_vert_to_old_edge_[v0]; 
				uint rec_id_0 = load_ssbo_edge_to_new_temporal_record_(old_edge_0); 

				uint old_edge_1 = ssbo_contour_vert_to_old_edge_[v1];
				uint rec_id_1 = load_ssbo_edge_to_new_temporal_record_(old_edge_1);

				cetd.temporal_rec_id = rec_id_0; 
			}

			/* write to intermediate buffer, will be shuffled after list ranking */
			store_contour_edge_transfer_data_(contour_id_global, cetd); 
		}

		{ /* Generate Topology */
			/* build contour edge adjacency */
			bool is_border = pcwi.is_border; 
			bool backface_border = is_border && !is_border_edge_front_facing(pcwi.ifrontface); 

			VertWedgeListHeader vwlh_beg_vtx, vwlh_end_vtx; 
			vwlh_beg_vtx.wedge_id = pcwi.wedge_id; 
			vwlh_beg_vtx.ivert = mark__cwedge_to_beg_vert(pcwi.ifrontface); 
			vwlh_end_vtx.wedge_id = pcwi.wedge_id;
			vwlh_end_vtx.ivert = mark__cwedge_to_end_vert(pcwi.ifrontface);

			VECircContext_ContourLinking ctx; 

			/* Rotate backwards around beg vert of this edge */
			ctx.pwci.is_contour = true;
			ctx.pwci.contour_id = contour_id_local; 
			bool rot_fwd = false; 

			VE_CIRCULATOR(vwlh_beg_vtx, func_ve_circulator, ctx, rot_fwd)

			if (ctx.pwci.contour_id == 0x1fffffffu) // hitting bad border edges (e0==e3,e1==e2) this is bad and should not happen
				ctx.pwci.contour_id = contour_id_local; 

			bool back_face_T_junction = (backface_border && !ctx.pwci.is_border); 
			ssbo_contour_to_contour_[contour_id_global*2] = calc_global_contour_edge_id(
				(!back_face_T_junction) ? ctx.pwci.contour_id : /*break backface border chain at T-junction*/contour_id_local
			); 

			/* Rotate fowards around end vert of this edge */
			ctx.pwci.is_contour = true; 
			ctx.pwci.contour_id = contour_id_local; 
			rot_fwd = true; 

			VE_CIRCULATOR(vwlh_end_vtx, func_ve_circulator, ctx, rot_fwd)

			if (ctx.pwci.contour_id == 0x1fffffffu) // hitting bad border edges (e0==e3,e1==e2)
				ctx.pwci.contour_id = contour_id_local; 

			back_face_T_junction = (backface_border && !ctx.pwci.is_border); 
			ssbo_contour_to_contour_[contour_id_global*2+1] = calc_global_contour_edge_id(
				(!back_face_T_junction) ? ctx.pwci.contour_id : /*break backface border chain at T-junction*/contour_id_local
			); 
		}
	}



	/* Prepare Soft Raster Data ------------------------------------------------------------- */
	LineRasterResult line_raster_data; 
	uint num_raster_frags = 0u;
	if (valid_thread)
	{
		/* transform matrices, see "common_view_lib.glsl" */ 
		mat4 world_to_view = ubo_view_matrices_.viewmat;
		mat4 mat_camera_proj = ubo_view_matrices_.winmat; 

		line_raster_data = raster_line_segment(
			vpos_ws[0], vpos_ws[1], pcs_screen_size_.xy, 
			world_to_view, mat_camera_proj
		);
		num_raster_frags = line_raster_data.num_frags; 
	}
	
	uint frag_offset = alloc_raster_frags(groupId, num_raster_frags);
	barrier(); 

	if (valid_thread)
	{ /* store header fragment */
		store_contour_edge_raster_data(contour_id_global, line_raster_data, frag_offset); 
	}
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FILL_LIST_RANKING_INPUTS_FOR_CONTOUR_EDGES)
void main()
{
	uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges;
	ssbo_list_ranking_inputs_.num_nodes = num_contours; 
}
#endif





// Visibility Test. This happens after the extraction loop (collecting contour edges from each object)
#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS)

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__VISIBILITY_TEST)
float load_depth_at(ivec2 coord, mat4 mat_camera_proj_inv)
{
	vec2 uvcoords = vec2(coord + .5f) / pcs_screen_size_.xy; 

	float z_ndc = texture(tex_remeshed_surf_depth_, uvcoords).r;
	vec3 ndc = vec3(uvcoords, z_ndc) * 2.0 - 1.0;
	
	vec4 pos_vs = mat_camera_proj_inv * vec4(ndc, 1.0);
	pos_vs.xyz = pos_vs.xyz / pos_vs.w;

	return pos_vs.z; 
}

mat3 load_depth_3x3(ivec2 coord, mat4 mat_camera_proj_inv)
{
	mat3 z_samples_3x3 = mat3(1.0f);
	z_samples_3x3[0][0] = load_depth_at(coord + ivec2(-1, 1),  mat_camera_proj_inv);
	z_samples_3x3[1][0] = load_depth_at(coord + ivec2(0, 1),   mat_camera_proj_inv);
	z_samples_3x3[2][0] = load_depth_at(coord + ivec2(1, 1),   mat_camera_proj_inv);

	z_samples_3x3[0][1] = load_depth_at(coord + ivec2(-1, 0),  mat_camera_proj_inv);
	z_samples_3x3[1][1] = load_depth_at(coord,                 mat_camera_proj_inv);
	z_samples_3x3[2][1] = load_depth_at(coord + ivec2(1, 0),   mat_camera_proj_inv);

	z_samples_3x3[0][2] = load_depth_at(coord + ivec2(-1, -1), mat_camera_proj_inv);
	z_samples_3x3[1][2] = load_depth_at(coord + ivec2(0, -1),  mat_camera_proj_inv);
	z_samples_3x3[2][2] = load_depth_at(coord + ivec2(1, -1),  mat_camera_proj_inv);

	return z_samples_3x3; 
}

float frag_depth_test(mat3 z_vs_3x3, float z_frag, float z_tolerance)
{	
	float z = z_frag + z_tolerance; /* move a bit closer */

	float occ_counter = .0f; 	
	for (uint i = 0; i < 3; ++i)
		for (uint j = 0; j < 3; ++j)
			if (z_vs_3x3[i][j] < z)
				++occ_counter; 

	return occ_counter;
}
#endif



#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2) || defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_3)

/* interpolate edge data for contour #contour_id */
void interpo_contour_edge_transfer_data(
	uint contour_id, uint parent_contour_id, 
	float beg_ratio, float end_ratio, 
	bool only_update_visibility = false, bool seg_visible = true)
{
	ContourEdgeTransferData cetd_parent = load_contour_edge_transfer_data(parent_contour_id);
	ContourEdgeTransferData cetd; 
	cetd.vpos_ws[0]    = mix(cetd_parent.vpos_ws[0], 	cetd_parent.vpos_ws[1],    beg_ratio);
	cetd.vpos_ws[1]    = mix(cetd_parent.vpos_ws[0], 	cetd_parent.vpos_ws[1],    end_ratio);
	cetd.cf = cetd_parent.cf; // TODO: interpolate flags?
	if (!seg_visible) set_contour_flags_occluded(cetd.cf);
	cetd.cusp_funcs[0] = mix(cetd_parent.cusp_funcs[0], cetd_parent.cusp_funcs[1], beg_ratio);
	cetd.cusp_funcs[1] = mix(cetd_parent.cusp_funcs[0], cetd_parent.cusp_funcs[1], end_ratio); 
	cetd.temporal_rec_id = cetd_parent.temporal_rec_id;
	
	store_contour_edge_transfer_data_(contour_id, cetd); 
}

void set_contour_edge_transfer_data_occluded(uint contour_id)
{
	ContourEdgeTransferData cetd = load_contour_edge_transfer_data(contour_id);
	set_contour_flags_occluded(cetd.cf);
	store_contour_edge_transfer_data_(contour_id, cetd); 
}

ContourVisibilitySplitInfo load_contour_visibility_split_info(uint contour_id)
{
	ContourVisibilitySplitInfo cvsi; 
	uvec4 cvsi_enc; 
	Load4(ssbo_contour_visibility_split_info_, contour_id, cvsi_enc);
	cvsi = decode_contour_visibility_split_info(cvsi_enc);

	return cvsi; 
}
#endif


#define DEBUG_FRAGMENTS 1
#if defined(DEBUG_FRAGMENTS)
void store_frag_coord(uint frag_id, vec2 frag_coord, uint num_frags)
{
	uint offset = (frag_id + num_frags) * 2u;
	ssbo_frag_raster_data_[offset + 0u] = floatBitsToUint(frag_coord.x); 
	ssbo_frag_raster_data_[offset + 1u] = floatBitsToUint(frag_coord.y); 
}
vec2 load_frag_coord(uint frag_id, uint num_frags)
{
	uint offset = (frag_id + num_frags) * 2u;
	return vec2(
		uintBitsToFloat(ssbo_frag_raster_data_[offset + 0u]), 
		uintBitsToFloat(ssbo_frag_raster_data_[offset + 1u])
	); 
}
vec3 hash32(vec2 p) 
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}
vec3 hsl2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
vec3 rand_col_rgb(uint seed0, uint seed1)
{
    vec3 rnd = hash32(vec2(float(seed0 * 17), float(seed1 * 33))); 

    float hue = rnd.x;
    float saturation = 0.6 + rnd.z*0.4;
    float luminosity  = 0.6 + rnd.y*0.4;

    return hsl2rgb(vec3(hue, saturation, luminosity)); 
}
#endif

void main()
{
/* Use min-segscan to broadcast contour id ---------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__CLEAR_BUFFER)
	uint frag_id_global = gl_GlobalInvocationID.x;
	uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;
	bool valid_thread = frag_id_global < num_frags; 

	if (valid_thread)
		ssbo_tree_scan_input_contour_fragment_idmapping_[frag_id_global] = 
			segscan_uint_hf_encode(SSBOData_SegScanType_uint(0x3fffffffu, 0u));
	
	// debug only
	if (valid_thread)
		ssbo_frag_to_contour_[frag_id_global] = frag_id_global; 

#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__SETUP_SEGSCAN)
	uint contour_id = gl_GlobalInvocationID.x;
	uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges;
	bool valid_thread = contour_id < num_contours;

	/* Setup segscan data */
	if (valid_thread)
	{
		uint frag_offset = ssbo_contour_raster_data_[contour_id * 6u + 5u];
		bool contour_has_frags = frag_offset != 0xffffffffu; 
		
		if (contour_has_frags)
			ssbo_tree_scan_input_contour_fragment_idmapping_[frag_offset] = 
				segscan_uint_hf_encode(SSBOData_SegScanType_uint(contour_id, 1u)); 
	}

	if (contour_id == 0) 
	{
		uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;

		ssbo_tree_scan_infos_contour_segmentation_.num_scan_items = num_frags; 
		ssbo_tree_scan_infos_contour_segmentation_.num_valid_scan_threads = compute_num_threads(
			num_frags, 2u
		);
		ssbo_tree_scan_infos_contour_segmentation_.num_thread_groups = compute_num_groups(
			num_frags, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
		);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__FINISH_SEGSCAN)
	uint frag_id_global = gl_GlobalInvocationID.x;
	uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;
	bool valid_thread = frag_id_global < num_frags; 

	/* Finish segscan results */
	if (valid_thread)
	{
		// ssbo_tree_scan_output_contour_fragment_idmapping_ := ssbo_frag_to_contour_
		SSBOData_SegScanType_uint segscan_output = 
			segscan_uint_hf_decode(ssbo_frag_to_contour_[frag_id_global]); 
		SSBOData_SegScanType_uint segscan_input = 
			segscan_uint_hf_decode(ssbo_tree_scan_input_contour_fragment_idmapping_[frag_id_global]);
		
		// transform exclusive scan into inclusive
		segscan_output.val = min(segscan_output.val, segscan_input.val); 

		ssbo_frag_to_contour_[frag_id_global] = (segscan_output.val);
	}
#endif

/* Visibility Test ---------------------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__VISIBILITY_TEST)
	uint frag_id = gl_GlobalInvocationID.x;
	uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;
	bool valid_thread = frag_id < num_frags;

	mat4 mat_camera_proj_inv = ubo_view_matrices_.wininv; 

	uint contour_edge_id = ssbo_frag_to_contour_[frag_id]; 

	uint head_frag_id; 
	LineRasterResult line_raster_data = load_contour_edge_raster_data(contour_edge_id, /*out*/head_frag_id); 

	vec4 begend_frags = line_raster_data.begend_uvs.xyzw * pcs_screen_size_.xyxy;
	float linear_interp = 0;
    float linear_step = 0; // how much "factor" costs to go to neighbor frag on edge
	vec2 sampleTexel = calc_frag_screen_pos(
        begend_frags,
        head_frag_id,
        frag_id,
        line_raster_data.is_x_major_line,
        /* out */linear_interp, linear_step 
    );

	/* visibility test */
	mat3 z_vs_3x3 = load_depth_3x3(ivec2(sampleTexel.xy + 1e-10f), mat_camera_proj_inv); 
	float z_frag = - interpolate_frag_depth(
		line_raster_data.begend_wclips.x, 
		line_raster_data.begend_wclips.y, 
		linear_interp
	);
	bool visible = 1.0f < frag_depth_test(z_vs_3x3, z_frag, 0.01f * pcs_visibility_thresh_);
	

	
	if (valid_thread && visible)
		imageStore(tex2d_contour_dbg_, ivec2(sampleTexel), vec4(-z_frag, -z_vs_3x3[1][1], visible, .0f)); 

	if (valid_thread)
	{
		FragVisibilityTestResult fvtr;
		fvtr.visible = visible;
		fvtr.head_frag_id = head_frag_id; 
		ssbo_frag_raster_data_[frag_id] = encode_frag_visibility_test_result(fvtr); 

#if defined(DEBUG_FRAGMENTS)
		store_frag_coord(frag_id, sampleTexel, num_frags);	
#endif
	}

#endif

/* Segment fragments based on visibility --------------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__SETUP_SEGSCAN)
	uint frag_id = gl_GlobalInvocationID.x;
	uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;
	bool valid_thread = frag_id < num_frags;

	/* Setup segscan data */
	if (valid_thread)
	{
		FragVisibilityTestResult fvtr = decode_frag_visibility_test_result(ssbo_frag_raster_data_[frag_id]);
		
		/* Mark segment heads */
		uint prev_frag_id = max(1u, frag_id) - 1u; 
		FragVisibilityTestResult fvtr_prev = decode_frag_visibility_test_result(ssbo_frag_raster_data_[prev_frag_id]); 
		uint hf = 0u; 
		{
			bool head_frag_at_contour = (frag_id == 0u 				|| fvtr.head_frag_id == frag_id); // start of contour edge
			bool split_frag = (fvtr_prev.visible != fvtr.visible); // split frag if visibility changes
			hf = (head_frag_at_contour || split_frag) ? 1u : 0u; 
		}

		/* Mark segment tails */
		uint next_frag_id = min(num_frags - 1u, frag_id + 1u);
		FragVisibilityTestResult fvtr_next = decode_frag_visibility_test_result(ssbo_frag_raster_data_[next_frag_id]); 
		uint tf = 0u; 
		{
			bool tail_frag_at_contour = (frag_id == (num_frags - 1u) || fvtr_next.head_frag_id == next_frag_id); // end of contour edge
			bool split_frag = (fvtr.visible != fvtr_next.visible); // split frag if visibility changes
			tf = (tail_frag_at_contour || split_frag) ? 1u : 0u; 
		}

		/* 2 segscans to build segment topo */
		uint scan_item_id = frag_id; 		
		ssbo_tree_scan_input_contour_visibility_split_0_[scan_item_id] = segscan_uint_hf_encode(SSBOData_SegScanType_uint(1u, hf)); 
		scan_item_id = num_frags - 1u - frag_id; 
		ssbo_tree_scan_input_contour_visibility_split_1_[scan_item_id] = segscan_uint_hf_encode(SSBOData_SegScanType_uint(1u, tf)); 
	}

	if (frag_id == 0) 
	{
		uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;

		ssbo_tree_scan_infos_contour_segmentation_.num_scan_items = num_frags; 
		ssbo_tree_scan_infos_contour_segmentation_.num_valid_scan_threads = compute_num_threads(
			num_frags, 2u
		);
		ssbo_tree_scan_infos_contour_segmentation_.num_thread_groups = compute_num_groups(
			num_frags, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
		);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_0)
	uint contour_id = gl_GlobalInvocationID.x;
	uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges;
	bool valid_thread = contour_id < num_contours;

	/* Initialize per-contour vis-split info */
	if (valid_thread)
	{
		ContourVisibilitySplitInfo cvsi;
		init_contour_visibility_split_info(contour_id, /*inout*/cvsi);

		uvec4 cvsi_enc = encode_contour_visibility_split_info(cvsi); 
		Store4(ssbo_contour_visibility_split_info_, contour_id, cvsi_enc);
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_1)
	uint frag_id = gl_GlobalInvocationID.x;
	uint num_frags = ssbo_bnpr_mesh_pool_counters_.num_frags;
	bool valid_thread = frag_id < num_frags;
	const uint groupId = gl_LocalInvocationID.x; 

	FragVisibilityTestResult fvtr = decode_frag_visibility_test_result(ssbo_frag_raster_data_[frag_id]);

	uint contour_edge_id = ssbo_frag_to_contour_[frag_id]; 
	LineRasterResult line_raster_data = load_contour_edge_raster_data(contour_edge_id); 

	uint prev_contour_id = ssbo_contour_to_contour_[2u * contour_edge_id]; 
	uint next_contour_id = ssbo_contour_to_contour_[2u * contour_edge_id + 1u]; 


	/* Analyse segmented scan outputs */
	uint seg_len, seg_head, seg_tail, seg_rank;
	bool is_seg_head; 
	{
		SSBOData_SegScanType_uint segscan_output_0 = 
			segscan_uint_hf_decode(ssbo_tree_scan_output_contour_visibility_split_0_[frag_id]);
		SSBOData_SegScanType_uint segscan_output_1 =
			segscan_uint_hf_decode(ssbo_tree_scan_output_contour_visibility_split_1_[num_frags - 1u - frag_id]);

		uint dist_to_seg_head = segscan_output_0.val; // frag_id - seg_head_frag_id
		uint dist_to_seg_tail = segscan_output_1.val; // seg_tail_frag_id - frag_id
		
		seg_len = dist_to_seg_head + dist_to_seg_tail + 1u; 
		seg_head = frag_id - dist_to_seg_head;
		seg_tail = frag_id + dist_to_seg_tail;
		seg_rank = dist_to_seg_head; 
		is_seg_head = valid_thread && (seg_rank == 0u); 
	}

	/* Store mappings 
	 * frag_seg_head->prev/next_frag_seg_head at each RASTERIZED contour,  
	 * frag_seg_head->contour_id			  at each fragment 
	 * 
	 * so that mapping contour_id->prev/next_contour_id can be inferred using contour_id->frag_seg_head, 
	 * as long as the contour is rasterized. */
	
	uint prev_seg_head, next_seg_head; // frag_seg_head->prev/next_frag_seg_head (link between segment heads)
	/* However, 
	 * prev/next contours can be clipped away and have NO fragment mapped to it */
	bool prev_contour_not_rastered = false; 
	bool next_contour_not_rastered = false; 
	{
		/* Find previous segment head := (prev_seg_tail - prev_seg_tail_rank) */
		uint prev_seg_tail, prev_seg_tail_rank;

		bool has_prev_seg = valid_thread && (fvtr.head_frag_id < frag_id);
		if (has_prev_seg)
		{ /* find previous segment within the SAME contour */
			prev_seg_tail = max(1u, seg_head) - 1u;

			SSBOData_SegScanType_uint segscan_output_0_prev_seg_tail = 
				segscan_uint_hf_decode(ssbo_tree_scan_output_contour_visibility_split_0_[prev_seg_tail]);
			prev_seg_tail_rank = segscan_output_0_prev_seg_tail.val;
		}else
		{ /* Link to previous contour */
			if (prev_contour_id == contour_edge_id)
			{ /* point to self if there is nothing ahead */
				prev_seg_tail = seg_head; 
				prev_seg_tail_rank = 0u;
			}
			else
			{ /* find last seg of previous contour */
				uint prev_contour_head_frag_id; 
				LineRasterResult line_raster_data_prev = 
					load_contour_edge_raster_data(prev_contour_id, /*out*/prev_contour_head_frag_id);
				uint prev_contour_last_frag = prev_contour_head_frag_id + line_raster_data_prev.num_frags - 1u;
				prev_seg_tail = prev_contour_last_frag;

				SSBOData_SegScanType_uint segscan_output_0_prev_seg_tail = 
					segscan_uint_hf_decode(ssbo_tree_scan_output_contour_visibility_split_0_[prev_contour_last_frag]);
				prev_seg_tail_rank = segscan_output_0_prev_seg_tail.val;

				prev_contour_not_rastered = line_raster_data_prev.num_frags == 0u; /* prev_seg_head is invalid in this case */
			}
		}
		prev_seg_head = prev_seg_tail - prev_seg_tail_rank; 
		

		/* Find next segment */
		FragVisibilityTestResult fvtr_next_seg = decode_frag_visibility_test_result(ssbo_frag_raster_data_[seg_tail + 1u]);
		bool is_last_seg_at_contour = valid_thread && 
			((seg_tail == num_frags - 1u) || (fvtr_next_seg.head_frag_id == (seg_tail + 1u)));

		bool has_next_seg = valid_thread && !is_last_seg_at_contour; 
		if (has_next_seg) 
		  /* Find next segment within the SAME contour */
			next_seg_head = seg_tail + 1u;
		else
		{ /* Link to next contour */
			if (next_contour_id == contour_edge_id)
			{ /* point to self if there is nothing after */
				next_seg_head = seg_head; 
			}
			else
			{ /* find first seg of next contour */
				uint next_contour_head_frag_id; 
				LineRasterResult line_raster_data_next = 
					load_contour_edge_raster_data(next_contour_id, /*out*/next_contour_head_frag_id);
				next_seg_head = next_contour_head_frag_id;

				next_contour_not_rastered = line_raster_data_next.num_frags == 0u; 
			}
		}
	}
	
	
	/* Analyse contour edge */	
	bool zero_rastered_frags = 0u == line_raster_data.num_frags /* totally clipped */; 
	bool is_contour_segmented = seg_len != line_raster_data.num_frags; 
	
	bool is_frag_contour_header = fvtr.head_frag_id == frag_id; 
	bool add_new_contour = valid_thread 
		&& (is_contour_segmented) 
		&& (!is_frag_contour_header) /* 1st split just reuses the original contour edge */
		&& (!zero_rastered_frags)
		&& is_seg_head /* each seg head is responsible for this */
	;
	
	uint new_contour_id = compact_visibility_contour_split(add_new_contour, groupId); 
	barrier(); 

	if (is_seg_head && valid_thread)
	{ /* This should be ran for EVERY rasterized contour edge EXACTLY once */
		ContourVisibilitySplitInfo cvsi;
		cvsi.parent_contour_id = contour_edge_id;
		cvsi.is_visible = fvtr.visible;
		cvsi.is_new_contour = add_new_contour; 
		if (add_new_contour || (is_contour_segmented && is_frag_contour_header))
		{
			float beg_factor = float(seg_head - fvtr.head_frag_id) / float(line_raster_data.num_frags);
			float end_factor = float(seg_tail - fvtr.head_frag_id + 1.0f) / float(line_raster_data.num_frags);

			cvsi.begend_ratios = vec2(beg_factor, end_factor);
		}else
		{
			cvsi.begend_ratios = vec2(0.0f, 1.0f);
			cvsi.is_new_contour = false; 
		}
		cvsi.no_rastered_prev_contour = prev_contour_not_rastered;
		cvsi.prev_frag_seg_head_id = prev_seg_head; 
		cvsi.no_rastered_next_contour = next_contour_not_rastered;
		cvsi.next_frag_seg_head_id = next_seg_head; 
		cvsi.no_rastered_frags 	   = zero_rastered_frags; 
		cvsi.split_by_occlusion    = is_contour_segmented; 
		
		uint update_contour_id = add_new_contour ? new_contour_id : contour_edge_id;
		uvec4 cvsi_enc = encode_contour_visibility_split_info(cvsi); 
		Store4(ssbo_contour_visibility_split_info_, update_contour_id, cvsi_enc);

		if (add_new_contour)
		{ /* cache original contour links, mainly for tail segs when they reach to null edges(that is, not rastered) */
			ssbo_contour_to_contour_[2 * new_contour_id] = prev_contour_id; 
			ssbo_contour_to_contour_[2 * new_contour_id + 1] = next_contour_id;
		}

		/* Cache frag to contour mapping */
		ssbo_frag_seg_head_to_visibility_split_contour_[frag_id] = update_contour_id; 
	}

	#if defined(DEBUG_FRAGMENTS)
	if (valid_thread)
	{
		vec2 frag_coord = load_frag_coord(frag_id, num_frags);	
		imageStore(tex2d_contour_dbg_, ivec2(frag_coord), 
			fvtr.visible ? vec4(.0, 1.0, .0, .0) : vec4(1.0, .0, 1.0, .0));
	}
	#endif

#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2_3)
	uint contour_id = gl_GlobalInvocationID.x;
	uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges;

	bool valid_thread = contour_id < num_contours;
	if (!valid_thread) return;


	ContourVisibilitySplitInfo cvsi = load_contour_visibility_split_info(contour_id);
	uint split_type = calc_contour_visibility_states(cvsi);

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2)
	uint split_type_mask = BIT_TYPE_CONTOUR_VIS_SPLIT__NEW_EDGE; 
	
	if (0u != (split_type & split_type_mask))
	{ /* Generate geometry for new contour edges */
		interpo_contour_edge_transfer_data(
			contour_id, 
			cvsi.parent_contour_id, cvsi.begend_ratios[0], cvsi.begend_ratios[1], 
			false, cvsi.is_visible
		); 
	}

	
	split_type_mask = (BIT_TYPE_CONTOUR_VIS_SPLIT__NEW_EDGE | BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_SPLIT);
	
	if (0u != (split_type & split_type_mask))
	{ /* build new edge adjacency & update old edge adjacency */
		uint prev_contour_id; 
		if (false == cvsi.no_rastered_prev_contour)
			prev_contour_id = ssbo_frag_seg_head_to_visibility_split_contour_[cvsi.prev_frag_seg_head_id]; 
		else
			prev_contour_id = ssbo_contour_to_contour_[2u * contour_id + 0u];
		ssbo_contour_to_contour_[2u * contour_id + 0u] = prev_contour_id;

		ContourVisibilitySplitInfo cvsi_prev = load_contour_visibility_split_info(prev_contour_id);
		uint split_type_prev = calc_contour_visibility_states(cvsi_prev); 
		if (0u == (split_type_prev & split_type_mask))
			ssbo_contour_to_contour_[2u * prev_contour_id + 1u] = contour_id;


		uint next_contour_id; 
		if (false == cvsi.no_rastered_next_contour)
			next_contour_id = ssbo_frag_seg_head_to_visibility_split_contour_[cvsi.next_frag_seg_head_id]; 
		else /* ssbo_contour_to_contour_ has been cached  
			  * for new contour generated fromt the last segment in its old contour */
			next_contour_id = ssbo_contour_to_contour_[2u * contour_id + 1u];
		ssbo_contour_to_contour_[2u * contour_id + 1u] = next_contour_id;

		ContourVisibilitySplitInfo cvsi_next = load_contour_visibility_split_info(next_contour_id);
		uint split_type_next = calc_contour_visibility_states(cvsi_next);
		if (0u == (split_type_next & split_type_mask))
			ssbo_contour_to_contour_[2u * next_contour_id + 0u] = contour_id;
	}
#endif

#if defined(_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_3)
	/* Overwrite geometry for old edges by its 1st split new edge */
	uint split_type_mask = BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_SPLIT; 
	if (0u != (split_type & split_type_mask))
	{ /* Generate geometry for new contour edges */
		interpo_contour_edge_transfer_data(
			contour_id, 
			contour_id, cvsi.begend_ratios[0], cvsi.begend_ratios[1], 
			false, cvsi.is_visible
		); 
	}else if (!cvsi.is_visible){
		set_contour_edge_transfer_data_occluded(contour_id); 
	}
#endif
	

#endif


}

#endif