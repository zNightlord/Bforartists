
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_decode_ibo_lib.glsl)


 
/* all counters are cleared in _KERNEL_MULTICOMPILE__GEOM_EXTRACT kernel ------------------- */
#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
	#define GLOBAL_COMPACTION_COUNTER__MESH_CONTOUR_EDGES ssbo_bnpr_mesh_pool_counters_.num_contour_edges
	#define GLOBAL_COMPACTION_COUNTER__MESH_VERTS ssbo_bnpr_mesh_pool_counters_.num_verts
	#define GLOBAL_COMPACTION_COUNTER__MESH_EDGES ssbo_bnpr_mesh_pool_counters_.num_edges
	#define GLOBAL_COMPACTION_COUNTER__MESH_FACES ssbo_bnpr_mesh_pool_counters_.num_faces
#endif


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
	out float face_orient_012, out float face_orient_321, 
	out vec3 edge_offset
)
{ /* impl based on overlay_outline_prepass_vert_no_geom.glsl */
	vec3 v10 = v0 - v1;
   	vec3 v12 = v2 - v1;
	vec3 v13 = v3 - v1;

	vec3 n0 = cross(v12, v10);
	vec3 n3 = cross(v13, v12);
	edge_offset = .5f * (normalize(n0), normalize(n3)); 

	vec3 view_dir = cam_pos - v1; 

	face_orient_012 = dot(view_dir, n0);
  	face_orient_321 = dot(view_dir, n3);
	bool is_contour = (sign(face_orient_012) != sign(face_orient_321)); 

	/* convexity test */
	vec3 p0 = v0; 
	vec3 p1 = v1;
	vec3 p10 = p0 - p1;
	float p10_d_n3 = dot(normalize(p10), normalize(n3)); 
	bool concave_edge = p10_d_n3 > .01f; 
	if (concave_edge) is_contour = false; /* must be hidden */
	
	return is_contour; 
}

vec3 ld_vbo(uint vert)
{
	uint base_addr = vert * 4; /* vbo is aligned to vec4 per vert pos */
	return vec3(buf_vbo[base_addr], buf_vbo[base_addr+1], buf_vbo[base_addr+2]); 
}

void main()
{
	const uint groupId = gl_LocalInvocationID.x;
	const uint idx 	   = gl_GlobalInvocationID.x;

	if (pcs_clear_compaction_counter_ != 0u)
	{ /* bootstrapping pass */
		GLOBAL_COMPACTION_COUNTER__MESH_CONTOUR_EDGES = 0; 
		GLOBAL_COMPACTION_COUNTER__MESH_VERTS = 0; 
		GLOBAL_COMPACTION_COUNTER__MESH_EDGES = 0;
		GLOBAL_COMPACTION_COUNTER__MESH_FACES = 0; 
		return; 
	}
	
	const uint EdgeId = idx;
	const uint NumEdges = pcs_num_verts / 4u;
	bool valid_thread = EdgeId < NumEdges; 

	const uint resource_id = pcs_rsc_handle; 
	
	/* transform matrices, see "common_view_lib.glsl" */ 
	float4x4 model_to_world = drw_matrix_buf[resource_id].model; 
	float4x4 world_to_model = drw_matrix_buf[resource_id].model_inverse; 
	float4x4 view_to_world = ubo_view_matrices_.viewinv;
	float4x4 world_to_view = ubo_view_matrices_.viewmat;
	
	bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
	vec3 cam_pos_loc = (world_to_model * vec4(cam_pos_ws, 1.0f)).xyz; 
	
	vec3 v[4]; 
#if defined(_KERNEL_MULTICOMPILE__INDEX_BUFFER_16BIT)
	for (uint i = 0; i < 2; ++i)
	{
		uint ibo_addr = 2 * EdgeId + i;
		uint ibo_data = buf_ibo[ibo_addr];
		/* decode 16 bit index */
		uint ibo_data_16h = (ibo_data >> 16u);
    	uint ibo_data_16l = (ibo_data & 0xFFFFu);
		/* fetch vertex pos */
		v[i*2]   = ld_vbo(ibo_data_16l);
		v[i*2+1] = ld_vbo(ibo_data_16h); 
	}
#else
	for (uint i = 0; i < 4; ++i)
	{
		uint ibo_addr = 4 * EdgeId + i;
		uint ibo_data = buf_ibo[ibo_addr];
		/* fetch vertex pos */
		v[i] = ld_vbo(ibo_data);  
	}
#endif
	
	float face_orient_012, face_orient_321; 
	vec3 edge_offset_dir; 
	bool is_contour = is_contour_edge(
		v[0], v[1], v[2], v[3], cam_pos_loc
		, face_orient_012, face_orient_321, edge_offset_dir /*out*/
	); 
	bool rev_edge_dir = face_orient_012 < .0f; 
	if (false == valid_thread) is_contour = false; 

	uint compacted_idx = compact_contour_edge(is_contour, groupId); 

	if (is_contour)
	{
		/* transform to world space */
		for (uint i = 0; i < 4; ++i)
			v[i] = (model_to_world * vec4(v[i], 1.0f)).xyz; 
		
		/* write world pos to output buffer */
		uint base_addr = mesh_pool_addr__wpos_and_edgeid(compacted_idx); 
		uint addr_p0 = rev_edge_dir ? base_addr + 3 : base_addr;  
		uint addr_p1 = rev_edge_dir ? base_addr : base_addr + 3; 
		buf_strokegen_mesh_pool[addr_p0+0]  = floatBitsToUint(v[1].x); 
		buf_strokegen_mesh_pool[addr_p0+1]  = floatBitsToUint(v[1].y); 
		buf_strokegen_mesh_pool[addr_p0+2]  = floatBitsToUint(v[1].z); 
		buf_strokegen_mesh_pool[addr_p1+0]  = floatBitsToUint(v[2].x); 
		buf_strokegen_mesh_pool[addr_p1+1]  = floatBitsToUint(v[2].y); 
		buf_strokegen_mesh_pool[addr_p1+2]  = floatBitsToUint(v[2].z); 
		buf_strokegen_mesh_pool[base_addr+6]  = EdgeId; 
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
		ssbo_bnpr_mesh_pool_draw_args_.vertex_len 		= 2u * num_contour_edges;  	/*#verts*/
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

void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 

	const uint ContourEdgeIdx = idx; 
	const uint NumContourEdges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = ContourEdgeIdx < NumContourEdges; 

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
		
		uint addr_st = mesh_pool_addr__zwhclip(ContourEdgeIdx, NumContourEdges); 
		buf_strokegen_mesh_pool[addr_st+0] = floatBitsToUint(vpos_ndc[0].z); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_ndc[0].w); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_ndc[1].z); 
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_ndc[1].w); 

		

		vpos_uv[0].xy = ndc_to_screen_uv(vpos_ndc[0]);
		vpos_uv[1].xy = ndc_to_screen_uv(vpos_ndc[1]);
		vec2 edge_dir = (vpos_uv[1].xy - vpos_uv[0].xy); 
		float edge_dir_len = length(edge_dir);
		vec2 edge_dir_norm = edge_dir / edge_dir_len; 

		addr_st = mesh_pool_addr__edgedir(ContourEdgeIdx, NumContourEdges); 
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

		addr_st = mesh_pool_addr__edgeuv(ContourEdgeIdx, NumContourEdges); 
		vpos_uv[0].xy /= pcs_screen_size_.xy; 
		vpos_uv[1].xy /= pcs_screen_size_.xy; 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(vpos_uv[0].x); 
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(vpos_uv[0].y); 
		buf_strokegen_mesh_pool[addr_st+2] = floatBitsToUint(vpos_uv[1].x);
		buf_strokegen_mesh_pool[addr_st+3] = floatBitsToUint(vpos_uv[1].y);
	}
}
#endif


#if defined(_KERNEL_MULTICOMPILE__COMPACT_VBO)
/*
 * float ssbo_meshbatch_vbo_[]
 * int pcs_meshbatch_num_verts_
 * 
 * float ssbo_vbo_full_[]
 * int pcs_full_vbo_offset_
 * 
 * int pcs_rsc_handle_
 * ObjectMatrices drw_matrix_buf[]
 * ubo ubo_view_matrices_
 */
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

	uint st_base_addr = (pcs_full_vbo_offset_ + VertID) * 3u;  
	ssbo_vbo_full_[st_base_addr+0] = vpos_ws.x; 
	ssbo_vbo_full_[st_base_addr+1] = vpos_ws.y;
	ssbo_vbo_full_[st_base_addr+2] = vpos_ws.z;
 }
#endif



#if defined(_KERNEL_MULTICOMPILE__COMPACT_EDGE_ADJ_IBO)
/*
 * uint IBO_BUF[]
 * uint pcs_meshbatch_edge_count_
 * int pcs_full_ibo_offset_
 * uint pcs_edge_count_
 * uint ssbo_edge_to_vert_[]
*/
 void main()
 {
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    const uint EdgeID = idx.x; 

    bool valid_thread = (EdgeID < pcs_edge_count_); 
	if (!valid_thread) return; 

	uvec4 vid; 
    load_and_decode_ibo__edge_adj(EdgeID, 4u, /*out*/vid); 

	uint st_base_addr = (pcs_full_ibo_offset_ + EdgeID) * 4u;  
	ssbo_edge_to_vert_[st_base_addr+0] = vid[0];
	ssbo_edge_to_vert_[st_base_addr+1] = vid[1]; 
	ssbo_edge_to_vert_[st_base_addr+2] = vid[2];
	ssbo_edge_to_vert_[st_base_addr+3] = vid[3];
 }

#endif

