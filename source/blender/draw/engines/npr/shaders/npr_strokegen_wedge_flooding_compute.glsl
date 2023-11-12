#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)

uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING)
/*
 * uint pcs_edge_id_offset_; 
 * uint pcs_edge_count_; 
 * uint ssbo_wedge_flooding_pointers_in_[]; 
 * uint ssbo_wedge_flooding_pointers_out_[];
 * uint ssbo_edge_to_edges_[]; 
*/
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    /* Do not use EdgeID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < pcs_edge_count_); 
    
    const uint EdgeID = idx.x + pcs_edge_id_offset_; 


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER)
    WedgeFloodingPointer wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_in_[EdgeID]); 

    if (false == wfptr.is_seed)
    {
        /* Now try to find a seed */
        AdjWedgeInfo awi[4]; 
        WedgeFloodingPointer wfptr_bwedge[4]; 
        for (int ibwedge = 0; ibwedge < 4; ++ibwedge)
        {
            awi[ibwedge]          = decode_adj_wedge_info(
                ssbo_edge_to_edges_[EdgeID*4 + ibwedge]
            ); 
            wfptr_bwedge[ibwedge] = decode_wedge_flooding_pointer(
                ssbo_wedge_flooding_pointers_in_[awi[ibwedge].wedge_id]
            );

            if (wfptr_bwedge[ibwedge].is_seed)
            {
                wfptr.next_wedge_id = wfptr_bwedge[ibwedge].next_wedge_id; 
                wfptr.is_seed = true; 
                break; /* found a seed, quit */
            }
        }

        if (false == wfptr.is_seed)
        { /* no seed is found in adjacent wedges, jump with the long range ptr  */
            WedgeFloodingPointer wfptr_next = decode_wedge_flooding_pointer(
                ssbo_wedge_flooding_pointers_in_[wfptr.next_wedge_id]
            );
            if (wfptr_next.next_wedge_id != EdgeID)
            {
                wfptr.next_wedge_id = wfptr_next.next_wedge_id;
                wfptr.is_seed = wfptr_next.is_seed; 
            }else
            { /* long range ptr is pointing to this wedge */
                uint random_select = pcg(EdgeID) % 4u; 
                for (uint roll = 0; roll < 4u; ++ roll)
                {
                    uint ibwedge = (random_select + roll) % 4u; 
                    if (wfptr_bwedge[ibwedge].next_wedge_id != EdgeID)
                    {
                        wfptr.next_wedge_id = wfptr_bwedge[ibwedge].next_wedge_id;
                        wfptr.is_seed = wfptr_bwedge[ibwedge].is_seed;
                        break; 
                    }
                }
            }
        }
    }

    #if !defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER)
        if (valid_thread)
            ssbo_wedge_flooding_pointers_out_[EdgeID] = encode_wedge_flooding_pointer(wfptr); 
    #endif
 
    #if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER)
        /* Compact seed edges */
        bool is_seed_edge = wfptr.is_seed && valid_thread; 
        uint compacted_idx = compact_filtered_edge(is_seed_edge, groupId); 
        
        if (is_seed_edge) /* Store compacted seed edge */
            ssbo_filtered_edge_to_edge_[compacted_idx] = EdgeID;   
    #endif

#endif
}
#endif


#if defined(_KERNEL_MULTICOMPILE__COMPACT_FILTERED_VERTS)
/*
 * uint pcs_vert_id_offset_; 
 * uint pcs_vert_count_; 
 * uint ssbo_wedge_flooding_pointers_in_[];
 * uint ssbo_vert_to_edge_list_header_[]; 
 * uint ssbo_edge_to_edges_[]; 
 * uint ssbo_edge_to_vert_[]; 
 * uint ssbo_filtered_vert_to_vert_[]; 
 * compact_filtered_vert(bool valid_thread, uint groupId)
 * compact counter: uint ssbo_bnpr_mesh_pool_counters_.num_filtered_verts
*/
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < pcs_vert_count_); 
    
    const uint VertID = idx.x + pcs_vert_id_offset_; 
    
    bool vert_linked_to_seeded_wedge = false; 
    if (valid_thread)
    {
        uint wedge_id = ssbo_vert_to_edge_list_header_[VertID]; 
        uvec2 vids_cwedge = uvec2(
            ssbo_edge_to_vert_[wedge_id*4 + 1], 
            ssbo_edge_to_vert_[wedge_id*4 + 3]
        ); 

        /* Select how to fwd rotate the center wedge */
        uint ivert = (VertID == vids_cwedge[0]) ? 1u : 3u;
        uint iface_beg_rot_fwd = (ivert == 1u) ? 0u : 1u; 

        /* Rotate wedge around the vert */
        AdjWedgeInfo awi; 
        awi.iface_adj = iface_beg_rot_fwd; 
        awi.wedge_id  = wedge_id; 

        #define MAX_WEDGE_ROTATES 16u
        uint rotate_step = 0u; 
        do {
            WedgeFloodingPointer wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_in_[awi.wedge_id]);
            if (wfptr.is_seed) {
                vert_linked_to_seeded_wedge = true; 
                break; 
            }

            uint iwedge_next = mark__cwedge_rotate_next(awi.iface_adj); 
            awi = decode_adj_wedge_info(ssbo_edge_to_edges_[awi.wedge_id*4u + iwedge_next]);

            rotate_step++; 
        } while (
            rotate_step < MAX_WEDGE_ROTATES 
            && awi.wedge_id != wedge_id
        ); 
    }

    vert_linked_to_seeded_wedge = vert_linked_to_seeded_wedge && valid_thread; 
    /* Compact seed verts */
    uint compacted_idx = compact_filtered_vert(vert_linked_to_seeded_wedge, groupId); 
        
    if (vert_linked_to_seeded_wedge) /* Store compacted seed vert */
        ssbo_filtered_vert_to_vert_[compacted_idx] = VertID;   
}
#endif




#if defined(_KERNEL_MULTICOMPILE__WEDGE_QUADRICS)
void store_wedge_quadric(Quadric q)
{
    vec4 data_0, data_1; 
    vec2 data_2; 
    encode_quadric(q_e, /*out*/data_0, data_1, data_2); 

    ssbo_quadric_data_[FilteredEdgeID*10u + 0] = floatBitsToUint(data_0.x);
    ssbo_quadric_data_[FilteredEdgeID*10u + 1] = floatBitsToUint(data_0.y);
    ssbo_quadric_data_[FilteredEdgeID*10u + 2] = floatBitsToUint(data_0.z);
    ssbo_quadric_data_[FilteredEdgeID*10u + 3] = floatBitsToUint(data_0.w);
    ssbo_quadric_data_[FilteredEdgeID*10u + 4] = floatBitsToUint(data_1.x);
    ssbo_quadric_data_[FilteredEdgeID*10u + 5] = floatBitsToUint(data_1.y);
    ssbo_quadric_data_[FilteredEdgeID*10u + 6] = floatBitsToUint(data_1.z);
    ssbo_quadric_data_[FilteredEdgeID*10u + 7] = floatBitsToUint(data_1.w);
    ssbo_quadric_data_[FilteredEdgeID*10u + 8] = floatBitsToUint(data_2.x);
    ssbo_quadric_data_[FilteredEdgeID*10u + 9] = floatBitsToUint(data_2.y);
}
void load_edge_quadric(uint EdgeID, /*out*/out Quadric q)
{
    vec4 data_0, data_1; 
    vec2 data_2; 
    data_0.x = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 0]);
    data_0.y = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 1]);
    data_0.z = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 2]);
    data_0.w = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 3]);
    data_1.x = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 4]);
    data_1.y = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 5]);
    data_1.z = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 6]);
    data_1.w = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 7]);
    data_2.x = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 8]);
    data_2.y = uintBitsToFloat(ssbo_quadric_data_[EdgeID*10u + 9]);
    q = decode_quadric(data_0, data_1, data_2); 
}

/*
 * uint ssbo_bnpr_mesh_pool_counters_.num_filtered_edges 
 * uint ssbo_filtered_edge_to_edge_[]; 
 * uint ssbo_quadric_data_[]; 
 * ubo_view_matrices_
*/
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    

#if defined(_KERNEL_MULTICOMPILE__WEDGE_QUADRICS__CALC_EDGE)
    const uint FilteredEdgeID = idx.x; 
    const uint NumFilteredEdges = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
    bool valid_thread = (idx.x < NumFilteredEdges); 
    if (!valid_thread) return; /* quit if not valid thread */

    /* Do not use EdgeID here since it's offseted with current mesh batch */
    uint wedge_id = ssbo_filtered_edge_to_edge_[FilteredEdgeID]; 
    
    uvec4 vids_wedge;
    uint base_addr = wedge_id * 4u; 
    for (uint ivert = 0; ivert < 4u; ++ivert)
        vids_wedge[ivert] = ssbo_edge_to_vert_[base_addr+ivert];

    vec3 vpos_wedge[4]; 
	for (uint ivert = 0; ivert < 4; ++ivert)
		vpos_wedge[ivert] = ld_vbo(vids_wedge[ivert]);  

    /* transform matrices, see "common_view_lib.glsl" */ 
	mat4 view_to_world = ubo_view_matrices_.viewinv;
	bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

    Quadric q_e = compute_quadric(vpos_wedge[0], vpos_wedge[1], vpos_wedge[2], vpos_wedge[3], cam_pos_ws); 
    if (valid_thread)
       store_wedge_quadric(q_e); 
#endif


#if defined(_KERNEL_MULTICOMPILE__WEDGE_QUADRICS__CALC_VERT)
    const uint VertID = idx.x; 
#endif



}
#endif
