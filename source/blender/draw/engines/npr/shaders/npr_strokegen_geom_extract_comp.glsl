


 


/* ---------------------------------- Compaction --------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
	#define GLOBAL_COMPACTION_COUNTER__MESH_EDGES ssbo_bnpr_mesh_pool_counters_.num_contour_edges
#endif

#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS)
shared uint LDS_digit_per_lane[32u]; 
shared uint LDS_offset_per_lane_slot[32u];  
shared uint LDS_hist_blk;    
shared uint LDS_scan_block_offset; 

uint compact(bool val, uint groupIdx)
{
	const uint wave_id = groupIdx >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */
    const uint lane_id = groupIdx % 32u;
    const uint num_waves = gl_WorkGroupSize.x >> 5u; 
    
	/* Clear LDS counters */
    if (wave_id == 0u) 
    { 
        if (lane_id == 0u) LDS_hist_blk = 0u; 
        LDS_digit_per_lane[lane_id] = 0u; 
    }
    barrier(); 
	/*  w0    w1    w2    w3         LDS_digit_per_lane
	 * 00:1  04:1  08:0  12:1  l0    0000
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0000
	 * 03:1  07:0  11:1  15:0  l3    0000
	 * ---------------------------
	 * LDS_hist_blk
	 *  0
	*/

    /* Mark 1/0 at bit #wave_id */ 
    uint compact_bitval = uint(val); 
    uint lds_compact_input = compact_bitval << wave_id; 
    atomicOr(LDS_digit_per_lane[lane_id], lds_compact_input); 
    barrier(); 
	/*  w0    w1    w2    w3         LDS_digit_per_lane
	 * 00:1  04:1  08:0  12:1  l0    1011
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0011
	 * 03:1  07:0  11:1  15:0  l3    0101
	*/

    /* Prefix sum on lane sums */
    uint lane_digit = LDS_digit_per_lane[lane_id]; 
    uint wave_mask = (~(0xffffffffu << wave_id));
	/*  w0    w1    w2    w3     
	 * 0000  0001  0011  0111        wave_mask
	 *
	 *  w0    w1    w2    w3         lane_digit
	 * 00:1  04:1  08:0  12:1  l0    1011
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0011
	 * 03:1  07:0  11:1  15:0  l3    0101
	*/
    uint lane_digit_masked = lane_digit & wave_mask; 
    uint num_1_bits_low = bitCount(lane_digit_masked);
    uint lane_offset = num_1_bits_low; 

    if (wave_id == 0u)
    {
        uint lane_digit_sum = bitCount(lane_digit); /* digit sum */ 
        LDS_offset_per_lane_slot[lane_id] = atomicAdd(LDS_hist_blk, lane_digit_sum);  
    }
    barrier(); 

    /* Add block sum to global counter. */
    if (groupIdx == gl_WorkGroupSize.x - 1u)
    {
        LDS_scan_block_offset = atomicAdd(
            GLOBAL_COMPACTION_COUNTER__MESH_EDGES, 
            LDS_hist_blk
        ); 
    }
    barrier(); 

    /* Compute final offset */
    uint local_offset = LDS_offset_per_lane_slot[lane_id] + lane_offset; 
    uint blk_offset   = LDS_scan_block_offset; 
    
    uint scanres = local_offset + blk_offset; 
	return scanres; 
}
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
	vec3 p1 = (v1+v2+v3) / 3.0f; /* bary center of Tri321 */
	vec3 p10 = p0 - p1;
	float p10_d_n3 = dot(normalize(p10), normalize(n3)); 
	bool concave_edge = p10_d_n3 < -.05f; 
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
		GLOBAL_COMPACTION_COUNTER__MESH_EDGES = 0; 
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

	uint compacted_idx = compact(is_contour, groupId); 

	if (is_contour)
	{
		/* transform to world space */
		for (uint i = 0; i < 4; ++i)
			v[i] = (model_to_world * vec4(v[i], 1.0f)).xyz; 
		v[1] -= edge_offset_dir * 0.001f; 
		v[2] -= edge_offset_dir * 0.001f; 
		
		/* write world pos to output buffer */
		uint base_addr = mesh_pool_addr__wpos(compacted_idx); 
		uint addr_p0 = rev_edge_dir ? base_addr + 3 : base_addr;  
		uint addr_p1 = rev_edge_dir ? base_addr : base_addr + 3; 
		buf_strokegen_mesh_pool[addr_p0+0]  = floatBitsToUint(v[1].x); 
		buf_strokegen_mesh_pool[addr_p0+1]  = floatBitsToUint(v[1].y); 
		buf_strokegen_mesh_pool[addr_p0+2]  = floatBitsToUint(v[1].z); 
		buf_strokegen_mesh_pool[addr_p1+0]  = floatBitsToUint(v[2].x); 
		buf_strokegen_mesh_pool[addr_p1+1]  = floatBitsToUint(v[2].y); 
		buf_strokegen_mesh_pool[addr_p1+2]  = floatBitsToUint(v[2].z); 
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
	return pos_ndc.xy; /* TODO: flip y or not? */
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
		uint base_addr = mesh_pool_addr__wpos(ContourEdgeIdx); 
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

		vpos_uv[0].xy = ndc_to_screen_uv(vpos_ndc[0]);
		vpos_uv[1].xy = ndc_to_screen_uv(vpos_ndc[1]);
		vec2 edge_dir = normalize(vpos_uv[1].xy - vpos_uv[0].xy); 
 
		uint addr_st = mesh_pool_addr__edgedir(ContourEdgeIdx, NumContourEdges); 
		buf_strokegen_mesh_pool[addr_st] = floatBitsToUint(edge_dir.x);
		buf_strokegen_mesh_pool[addr_st+1] = floatBitsToUint(edge_dir.y);
	}
}
#endif


