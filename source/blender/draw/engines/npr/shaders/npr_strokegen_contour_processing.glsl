
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)


#if defined(_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION)
/*
uint ssbo_contour_edge_rank_in_/out_[]
uint ssbo_contour_edge_list_len_in_/out_[]
uint ssbo_contour_edge_list_head_in_/out_[]
uint ssbo_contour_edge_vpos_in_/out_[]

ssbo_bnpr_mesh_pool_counters_.num_contour_edges
*/
void main()
{
    const uint idx = gl_GlobalInvocationID.x; 
	const uint contour_id = idx; 

	const uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
	bool valid_thread = contour_id < num_contours; 

    uint head_contour_id = ssbo_contour_edge_list_head_in_[contour_id];
    uint list_len        = ssbo_contour_edge_list_len_in_[contour_id]; 
    uint rank            = ssbo_contour_edge_rank_in_[contour_id]; 
    
    uint output_addr = head_contour_id + rank; 
    if (valid_thread)
    {
        ssbo_contour_edge_rank_out_[output_addr]      = rank; 
        ssbo_contour_edge_list_len_out_[output_addr]  = list_len;
        ssbo_contour_edge_list_head_out_[output_addr] = head_contour_id;

        uvec3 vpos_0_enc, vpos_1_enc; 
        Load3(ssbo_contour_edge_vpos_in_,   contour_id*2u,     vpos_0_enc);
        Store3(ssbo_contour_edge_vpos_out_, output_addr*2u,    vpos_0_enc);
        Load3(ssbo_contour_edge_vpos_in_,   contour_id*2u+1u,  vpos_1_enc); 
        Store3(ssbo_contour_edge_vpos_out_, output_addr*2u+1u, vpos_1_enc);
    }
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
        uint head_contour_id = ssbo_contour_edge_list_head_[contour_id];
        uint list_len        = ssbo_contour_edge_list_len_[contour_id]; 
        uint rank            = ssbo_contour_edge_rank_[contour_id]; 
        uint addr_ld = head_contour_id + rank; 

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

	if (idx.x == 0) /* I dont care, this is cheap, just keep this up to date here */
	{
		uint num_contours = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
		
		ssbo_tree_scan_infos_contour_segmentation_.num_scan_items = num_contours; 
		ssbo_tree_scan_infos_contour_segmentation_.num_valid_scan_threads = compute_num_threads(
			num_contours, 2u
		);
      	ssbo_tree_scan_infos_contour_segmentation_.num_thread_groups = compute_num_groups(
			num_contours, GROUP_SIZE_BNPR_SCAN_SWEEP, 2u
      	);
	}

}
#endif