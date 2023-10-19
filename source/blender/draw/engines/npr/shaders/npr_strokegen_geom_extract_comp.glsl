

/*    v0
 *   /  \
 *  /    \                
 * v1----v2             
 *  \    /                
 *   \  /                     
 *    v3    winding 012, 321                
*/
bool is_contour_edge(vec3 v0, vec3 v1, vec3 v2, vec3 v3, vec3 cam_pos)
{ /* impl based on overlay_outline_prepass_vert_no_geom.glsl */
	vec3 v10 = v0 - v1;
   	vec3 v12 = v2 - v1;
	vec3 v13 = v3 - v1;

	vec3 n0 = cross(v12, v10);
	vec3 n3 = cross(v13, v12);

	vec3 view_dir = cam_pos - v1; 

	float face_orient_012 = dot(view_dir, n0);
  	float face_orient_321 = dot(view_dir, n3);
	bool is_contour = (sign(face_orient_012) != sign(face_orient_321)); 
	
	return is_contour; 
}


shared uint LDS_digit_per_lane[32u]; 
shared uint LDS_offset_per_lane_slot[32u]; 
shared uint LDS_hist_blk;  
shared uint LDS_scan_block_offset; 
#define GLOBAL_COMPACTION_COUNTER ssbo_bnpr_mesh_pool_acounters_.num_contour_edges

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
            GLOBAL_COMPACTION_COUNTER, 
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
		GLOBAL_COMPACTION_COUNTER = 0; 
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
	
	bool is_contour = is_contour_edge(v[0], v[1], v[2], v[3], cam_pos_loc); 
	if (false == valid_thread) is_contour = false; 

	uint compacted_idx = compact(is_contour, groupId); 

	uint base_addr = compacted_idx * 12u; 
	buf_strokegen_mesh_pool[base_addr+0]  = floatBitsToUint(v[0].x); 
	buf_strokegen_mesh_pool[base_addr+1]  = floatBitsToUint(v[0].y); 
	buf_strokegen_mesh_pool[base_addr+2]  = floatBitsToUint(v[0].z); 
	buf_strokegen_mesh_pool[base_addr+3]  = floatBitsToUint(v[1].x); 
	buf_strokegen_mesh_pool[base_addr+4]  = floatBitsToUint(v[1].y); 
	buf_strokegen_mesh_pool[base_addr+5]  = floatBitsToUint(v[1].z); 
	buf_strokegen_mesh_pool[base_addr+6]  = floatBitsToUint(v[2].x); 
	buf_strokegen_mesh_pool[base_addr+7]  = floatBitsToUint(v[2].y); 
	buf_strokegen_mesh_pool[base_addr+8]  = floatBitsToUint(v[2].z); 
	buf_strokegen_mesh_pool[base_addr+9]  = floatBitsToUint(v[3].x); 
	buf_strokegen_mesh_pool[base_addr+10] = floatBitsToUint(v[3].y); 
	buf_strokegen_mesh_pool[base_addr+11] = floatBitsToUint(v[3].z); 
}