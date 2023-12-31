
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)


 
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
		ssbo_bnpr_mesh_pool_counters_.num_contour_edges_curr = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_filtered_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_.num_filtered_verts = 0; 

		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_verts         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_edges         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_faces         = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges_curr = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_edges = 0; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_verts = 0; 
	}
}
#endif


#if defined(_KERNEL_MULTICOMPILE__COMPACT_VBO)
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
	float4x4 model_to_world = drw_matrix_buf[ResourceID].model;  
	float4x4 world_to_model = drw_matrix_buf[ResourceID].model_inverse; 
	float4x4 view_to_world = ubo_view_matrices_.viewinv;
	float4x4 world_to_view = ubo_view_matrices_.viewmat;

	vec3 vpos_ls, vpos_ws; 

	uint ld_base_addr = VertID; /* TODO: posnor buf, verify this */ 
	SSBO_Data_PosNor posnor_packed = ssbo_meshbatch_vbo_[ld_base_addr]; 
	vpos_ls.x = posnor_packed.x;
	vpos_ls.y = posnor_packed.y; 
	vpos_ls.z = posnor_packed.z;

	vpos_ws = (model_to_world * vec4(vpos_ls, 1.0f)).xyz; 

	uint st_base_addr = (VertID) * 3u;  
	ssbo_vbo_full_[st_base_addr+0] = vpos_ws.x; 
	ssbo_vbo_full_[st_base_addr+1] = vpos_ws.y;
	ssbo_vbo_full_[st_base_addr+2] = vpos_ws.z;

	if (idx == 0)
	{ /* cache counters from last mesh extraction pass */
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_verts         = ssbo_bnpr_mesh_pool_counters_.num_verts; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_edges         = ssbo_bnpr_mesh_pool_counters_.num_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_faces         = ssbo_bnpr_mesh_pool_counters_.num_faces; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges_curr = ssbo_bnpr_mesh_pool_counters_.num_contour_edges_curr; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_edges = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
		ssbo_bnpr_mesh_pool_counters_prev_.num_filtered_verts = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
	}
 }
#endif


#if defined(_KERNEL_MULTICOMPILE__COMPACT_EDGE_ADJ_IBO)
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




bool is_back_face(float ndv)
{
	return ndv <= .0f; 
}

#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT)
/*    v0
 *   /  \
 *  /    \                
 * v1----v2             
 *  \    /                
 *   \  /                     
 *    v3    winding 012, 321                
*/
bool is_contour_edge(
	vec3 v0, vec3 v1, vec3 v2, vec3 v3, vec3 cam_pos, 
	out float face_orient_012, out float face_orient_321
)
{ /* impl based on overlay_outline_prepass_vert_no_geom.glsl */
	vec3 v10 = v0 - v1;
   	vec3 v12 = v2 - v1;
	vec3 v13 = v3 - v1;

	vec3 n0 = cross(v12, v10);
	vec3 n3 = cross(v13, v12);

	vec3 view_dir = cam_pos - v1; 

	face_orient_012 = dot(view_dir, n0);
  	face_orient_321 = dot(view_dir, n3);
	bool is_contour = 
		/* !is_back_face(face_orient_012) && !is_back_face(face_orient_321);  */
		is_back_face(face_orient_012) != is_back_face(face_orient_321); 
		/* (sign(face_orient_012) != sign(face_orient_321)); */

	/* convexity test */
	vec3 p0 = v0; 
	vec3 p1 = v1;
	vec3 p10 = p0 - p1;
	float p10_d_n3 = dot(normalize(p10), normalize(n3)); 
	bool concave_edge = p10_d_n3 > .01f; 
	/* Note: don't do this when linking contour edges */
	/* if (concave_edge) is_contour = false; */ /* must be hidden */
	
	return is_contour; 
}

vec3 ld_vbo(uint vert)
{
	uint base_addr = vert * 3; /* vbo is aligned to vec4 per vert pos */
	return vec3(buf_vbo[base_addr], buf_vbo[base_addr+1], buf_vbo[base_addr+2]); 
}

uint get_num_edges()
{
	return pcs_num_edges + ssbo_dyn_mesh_counters_.num_edges; 
}

void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint idx 	   = gl_GlobalInvocationID.x;
	
	const uint NumEdges = get_num_edges();
	bool valid_thread = idx < NumEdges; 
	const uint EdgeId = idx;

	const uint resource_id = pcs_rsc_handle; 
	
	/* transform matrices, see "common_view_lib.glsl" */ 
	mat4 model_to_world = drw_matrix_buf[resource_id].model; 
	mat4 world_to_model = drw_matrix_buf[resource_id].model_inverse; 
	mat4 view_to_world = ubo_view_matrices_.viewinv;
	mat4 world_to_view = ubo_view_matrices_.viewmat;
	
	bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
	vec3 cam_pos_loc = (world_to_model * vec4(cam_pos_ws, 1.0f)).xyz; 
	
	uvec4 vids; 
	for (uint i = 0; i < 4; ++i)
		vids[i] = buf_ibo[4 * EdgeId + i]; 
	vids = wing_verts_to_line_adj(vids); 
	
	vec3 v[4]; 
	for (uint i = 0; i < 4; ++i)
		v[i] = ld_vbo(vids[i]);  


	float face_orient_012, face_orient_321; 
	bool is_contour = is_contour_edge(
		v[0], v[1], v[2], v[3], cam_pos_ws
		, face_orient_012, face_orient_321/*out*/
	); 

	EdgeFlags ef = load_edge_flags(EdgeId); 
	is_contour = is_contour && (!ef.del_by_split) && (!ef.del_by_collapse) && (!ef.dupli); 
	bool is_border = ef.border; 

	/* debug view */
	if (pcs_dbg_wedge_flooding_ > 0)
	{
		// Check wfptr.is_seed; 
		WedgeFloodingPointer wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_out_[EdgeId]);
		
		// Check topo consistency
		bool valid_ee, valid_ev, valid_ve;
		validate_wedge_topo(EdgeId, /*out*/ valid_ee, valid_ev, valid_ve);  
		
		is_contour = (!ef.del_by_collapse) && (!ef.dupli) && (!ef.del_by_split) 
			&& ( 
				(!valid_ev)
				|| (!valid_ev)
				|| (!valid_ve)
			); 
	}

	bool rev_edge_dir = is_back_face(face_orient_012); 
	if (false == valid_thread) is_contour = false; 

	uint compacted_idx = compact_contour_edge(is_contour, groupId); 

	if (is_contour)
	{
		/* write world pos to output buffer */
		uint base_addr = mesh_pool_addr__wpos_and_edgeid(compacted_idx); 
		uint addr_p0 = rev_edge_dir ? base_addr + 3 : base_addr;  
		uint addr_p1 = rev_edge_dir ? base_addr : base_addr + 3; 
		buf_strokegen_mesh_pool[addr_p0+0] = floatBitsToUint(v[1].x); 
		buf_strokegen_mesh_pool[addr_p0+1] = floatBitsToUint(v[1].y); 
		buf_strokegen_mesh_pool[addr_p0+2] = floatBitsToUint(v[1].z); 
		buf_strokegen_mesh_pool[addr_p1+0] = floatBitsToUint(v[2].x); 
		buf_strokegen_mesh_pool[addr_p1+1] = floatBitsToUint(v[2].y); 
		buf_strokegen_mesh_pool[addr_p1+2] = floatBitsToUint(v[2].z); 
		
		PerContourWedgeInfo pcwi; 
		pcwi.is_border = is_border; 
		pcwi.wedge_id = EdgeId;
		pcwi.ifrontface = (!is_back_face(face_orient_012)) ? 1 : 0; /* see "line_adj_to_wing_verts" */ 
		buf_strokegen_mesh_pool[base_addr+6]  = encode_per_contour_wedge_info(pcwi); 
	}

	if (valid_thread)
	{
		PerWedgeContourInfo peci; 
		peci.is_border = is_border; 
		peci.is_contour = is_contour;
		peci.contour_id = is_contour ? compacted_idx : NULL_EDGE; 
		peci.ifrontface = (!is_back_face(face_orient_012)) ? 1 : 0; /* see "line_adj_to_wing_verts" */
		ssbo_edge_to_contour_[EdgeId] = encode_per_wedge_contour_info(peci); 
	}
}
#endif




#if defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
	if (idx == 0u)
	{
		const uint num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
		ssbo_bnpr_mesh_pool_draw_args_.vertex_len 		= 2u * num_contour_edges;  /*#verts*/
		ssbo_bnpr_mesh_pool_draw_args_.instance_len 	= 1;  						/*#instances*/
		ssbo_bnpr_mesh_pool_draw_args_.vertex_first 	= 0;  						/*ibo offset*/
		ssbo_bnpr_mesh_pool_draw_args_.base_index 		= 0;  						/*vbo offset*/
		ssbo_bnpr_mesh_pool_draw_args_.instance_first_indexed = 0; 					/*instance offset*/
	}
}
#endif





#if defined(_KERNEL_MULTICOMPILE_CALC_CONTOUR_EDGE_RASTER_DATA)
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

	const uint NumContourEdges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges_curr; 
	const uint NumContourEdgesTotal = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = idx < NumContourEdges; 
	const uint ContourEdgeIdx = idx + ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges; 

	/* transform matrices, see "common_view_lib.glsl" */ 
	float4x4 world_to_view = ubo_view_matrices_.viewmat;
	float4x4 mat_camera_proj = ubo_view_matrices_.winmat; 

	vec4 vpos_ws[2]; /* v0, v1 on edge */
	vec4 vpos_ndc[2]; 
	vec2 vpos_uv[2]; 
	if (valid_thread)
	{ /* read vertex pos transformed to world space */
		uint base_addr = mesh_pool_addr__wpos_and_edgeid(ContourEdgeIdx); 
		vpos_ws[0].x = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+0]); 
		vpos_ws[0].y = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+1]); 
		vpos_ws[0].z = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+2]); 
		vpos_ws[0].w = 1.0f; 
		vpos_ws[1].x = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+3]); 
		vpos_ws[1].y = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+4]); 
		vpos_ws[1].z = uintBitsToFloat(buf_strokegen_mesh_pool[base_addr+5]); 
		vpos_ws[1].w = 1.0f; 

		vpos_ndc[0] = mat_camera_proj * vec4((world_to_view * vpos_ws[0]).xyz, 1.0f); 
		vpos_ndc[1] = mat_camera_proj * vec4((world_to_view * vpos_ws[1]).xyz, 1.0f); 
		
		uint addr_st = mesh_pool_addr__zwhclip(ContourEdgeIdx); 
		buf_strokegen_mesh_pool[addr_st+0] = floatBitsToUint(vpos_ndc[0].z); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_ndc[0].w); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_ndc[1].z); 
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_ndc[1].w); 

		

		vpos_uv[0].xy = ndc_to_screen_uv(vpos_ndc[0]);
		vpos_uv[1].xy = ndc_to_screen_uv(vpos_ndc[1]);
		vec2 edge_dir = (vpos_uv[1].xy - vpos_uv[0].xy); 
		float edge_dir_len = length(edge_dir);
		vec2 edge_dir_norm = edge_dir / edge_dir_len; 

		addr_st = mesh_pool_addr__edgedir(ContourEdgeIdx); 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(edge_dir_norm.x);
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(edge_dir_norm.y);



		/* Elongate edges that are too short on screen */
		const float min_edge_len = 2.0f; 
		if (edge_dir_len < min_edge_len)
		{
			float elongation = (min_edge_len - edge_dir_len) * .5f; 
			vpos_uv[0].xy -= edge_dir_norm * elongation; 
			vpos_uv[1].xy += edge_dir_norm * elongation; 
		}

		addr_st = mesh_pool_addr__edgeuv(ContourEdgeIdx); 
		vpos_uv[0].xy /= pcs_screen_size_.xy; 
		vpos_uv[1].xy /= pcs_screen_size_.xy; 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(vpos_uv[0].x); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_uv[0].y); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_uv[1].x);
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_uv[1].y);

		
		/* build contour edge adjacency */
		PerContourWedgeInfo pcwi = decode_per_contour_wedge_info(buf_strokegen_mesh_pool[base_addr+6]);  
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
		ctx.pwci.contour_id = ContourEdgeIdx; 
		bool rot_fwd = false; 

		VE_CIRCULATOR(vwlh_beg_vtx, func_ve_circulator, ctx, rot_fwd)

		bool back_face_T_junction = (backface_border && !ctx.pwci.is_border); 
		ssbo_contour_to_contour_[ContourEdgeIdx*2] = 
			(!back_face_T_junction) ? ctx.pwci.contour_id : /*break backface border chain at T-junction*/ContourEdgeIdx; 


		/* Rotate fowards around end vert of this edge */
		ctx.pwci.is_contour = true;
		ctx.pwci.contour_id = ContourEdgeIdx; 
		rot_fwd = true; 

		VE_CIRCULATOR(vwlh_end_vtx, func_ve_circulator, ctx, rot_fwd)

		back_face_T_junction = (backface_border && !ctx.pwci.is_border); 
		ssbo_contour_to_contour_[ContourEdgeIdx*2+1] = 
			(!back_face_T_junction) ? ctx.pwci.contour_id : /*break backface border chain at T-junction*/ContourEdgeIdx; 
	}

	if (idx.x == 0) /* I dont care, this is cheap, just keep this up to date here */
		ssbo_list_ranking_inputs_.num_nodes = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
}
#endif




