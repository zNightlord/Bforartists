#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)

uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}

struct VertFilteringInfo
{
    bool is_null_vert; 
    uint vert_id; 
    bool is_unstable; 
    bool is_border; 
};
uint encode_filtered_vert_info(VertFilteringInfo fvi)
{
    uint data = 0u; 
    data = (fvi.vert_id << 2u); 
    data |= ((fvi.is_unstable ? 1u : 0u) << 1u); 
    data |= (fvi.is_border   ? 1u : 0u); 

    if (fvi.is_null_vert)
        data = NULL_FILTERED_VERT; 
        
    return data; 
}
VertFilteringInfo decode_filtered_vert_info(uint data)
{
    VertFilteringInfo fvi; 
    fvi.is_null_vert = (data == NULL_FILTERED_VERT); 
    fvi.vert_id = (data >> 2u); 
    fvi.is_unstable = (data & 0x2u) != 0u; 
    fvi.is_border   = (data & 0x1u) != 0u; 
    return fvi; 
}

uint get_vert_count()
{
    return pcs_vert_count_ + ssbo_dyn_mesh_counters_out_.num_edges; 
}

uint get_edge_count()
{
    return pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_verts;
}

vec3 ld_vbo(uint GlobalVertID)
{
	uint base_addr = GlobalVertID * 3; 
	return vec3(ssbo_vbo_full_[base_addr], ssbo_vbo_full_[base_addr+1], ssbo_vbo_full_[base_addr+2]); 
}
void store_vbo(uint GlobalVertID, vec3 vpos)
{
    uint base_addr = GlobalVertID * 3; 
    ssbo_vbo_full_[base_addr]   = vpos.x;
    ssbo_vbo_full_[base_addr+1] = vpos.y;
    ssbo_vbo_full_[base_addr+2] = vpos.z;
}


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    /* Do not use EdgeID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < get_edge_count()); 
    
    const uint EdgeID = idx.x; 


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER)

    WedgeFloodingPointer wfptr; 

    #if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__INIT)
        uvec4 v; 
        v[0] = ssbo_edge_to_vert_[EdgeID*4u + 0];
        v[1] = ssbo_edge_to_vert_[EdgeID*4u + 1];
        v[2] = ssbo_edge_to_vert_[EdgeID*4u + 2];
        v[3] = ssbo_edge_to_vert_[EdgeID*4u + 3];
        bool is_border_wedge = wing_verts_is_border_edge(v); 
        
        vec3 vpos_wedge[4]; 
        for (uint ivert = 0; ivert < 4; ++ivert)
            vpos_wedge[ivert] = ld_vbo(v[ivert]);

        mat4 view_to_world = ubo_view_matrices_.viewinv;
        vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

        WedgeQuality wq = compute_wedge_quality(
            vpos_wedge[0], vpos_wedge[1], vpos_wedge[2], vpos_wedge[3], cam_pos_ws
        ); 
        
        wfptr.next_wedge_id = EdgeID; 
        wfptr.is_seed       = wq.unstable && !is_border_wedge; 
        wfptr.is_unstable   = wq.unstable_silouette && !is_border_wedge; 
        wfptr.is_border     = is_border_wedge; 

        if (valid_thread)
            ssbo_wedge_flooding_pointers_out_[EdgeID] = encode_wedge_flooding_pointer(wfptr); 

        return; /* quit after init */
    #else
        wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_in_[EdgeID]); 
    #endif
    
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

    if (valid_thread)
        ssbo_wedge_flooding_pointers_out_[EdgeID] = encode_wedge_flooding_pointer(wfptr); 
 
    #if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER)
        bool is_seed_edge = wfptr.is_seed && valid_thread; 

        #if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER__OUTPUT_FLAGS)
            EdgeFlags ef = load_edge_flags(EdgeID); 
            ef.selected = is_seed_edge; 
            if (valid_thread)
                store_edge_flags(EdgeID, ef);
        #endif

        #if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER__OUTPUT_CAMPACTION)
            /* Compact seed edges */
            uint compacted_idx = compact_filtered_edge(is_seed_edge, groupId); 
            
            /* Store edge <-> filtered_edge pointers */
            EdgeSelectionInfo fei; 
            fei.is_unstable = wfptr.is_unstable; 
            fei.is_border = wfptr.is_border;

            if (is_seed_edge && (0u < pcs_output_selected_edge_to_edge_))
            {
                fei.is_null_edge = false; 
                fei.edge_id = EdgeID; 
                ssbo_selected_edge_to_edge_[compacted_idx] = encode_edge_selection_info(fei);   
            }

            if (valid_thread && (0u < pcs_output_edge_to_selected_edge_))
            {
                fei.is_null_edge = (false == is_seed_edge); 
                fei.edge_id = compacted_idx;

                uint num_all_edges = get_edge_count(); 
                uint addr_edge_to_filtered_edge = ssbo_edge_to_edges_addr__edge_to_selected_edge(EdgeID, num_all_edges); 
                ssbo_edge_to_edges_[addr_edge_to_filtered_edge] = encode_edge_selection_info(fei);
            }
        #endif
    #endif

#endif
}
#endif




#if defined(_KERNEL_MULTICOMPILE__COMPACT_FILTERED_VERTS)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < get_vert_count()); 
    
    const uint VertID = idx.x; 
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[VertID]);
    uint wedge_id = vwlh.wedge_id; 
    uint vid_ivert1 = ssbo_edge_to_vert_[wedge_id*4 + 1]; 
    /* Select how to fwd rotate the center wedge */
    uint ivert = (VertID == vid_ivert1) ? 1u : 3u;
    uint iface_beg_rot_fwd = (ivert == 1u) ? 0u : 1u; 

    /* Traverse adj edges */
    VertFilteringInfo vfi;
    vfi.is_null_vert = false; 
    vfi.vert_id = VertID; 
    vfi.is_unstable = false; 
    vfi.is_border = false; 
    /* only select vert with 1-ring neighbors of full filtered wedges, 
     * since other non-filtered wedges do not have valid filtered data */
    bool vert_linked_to_seeded_wedge = true; 

    if (valid_thread)
    {
        /* Rotate wedge around the vert */
        AdjWedgeInfo awi; 
        awi.iface_adj = iface_beg_rot_fwd; 
        awi.wedge_id  = wedge_id; 

        #define MAX_WEDGE_ROTATES 16u
        uint rotate_step = 0u; 
        do {
            WedgeFloodingPointer wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_in_[awi.wedge_id]);
            if (!wfptr.is_seed) {
                vert_linked_to_seeded_wedge = false; 
                break; 
            }
            vfi.is_unstable = vfi.is_unstable || wfptr.is_unstable; 
            vfi.is_border   = vfi.is_border   || wfptr.is_border; 

            uint iwedge_next = mark__cwedge_rotate_next(awi.iface_adj); 
            awi = decode_adj_wedge_info(ssbo_edge_to_edges_[awi.wedge_id*4u + iwedge_next]);

            rotate_step++; 
        } while (
            rotate_step < MAX_WEDGE_ROTATES 
            && awi.wedge_id != wedge_id
        ); 
    }

    vert_linked_to_seeded_wedge = vert_linked_to_seeded_wedge && valid_thread; 

    uint vert_id_merged = ssbo_vert_merged_id_[VertID];
    vert_linked_to_seeded_wedge = vert_linked_to_seeded_wedge && (vert_id_merged == VertID); 

    /* Compact seed verts */
    uint compacted_idx = compact_filtered_vert(vert_linked_to_seeded_wedge, groupId); 
            
    if (vert_linked_to_seeded_wedge) /* Store compacted seed vert */
    {
        vfi.is_null_vert = false; 
        vfi.vert_id = VertID; 
        ssbo_selected_vert_to_vert_[compacted_idx] = encode_filtered_vert_info(vfi);   
    }

    if (valid_thread)
    {        
        VertWedgeListHeader vwlh;
        vwlh.wedge_id = wedge_id; 
        vwlh.ivert = ivert; 
        ssbo_vert_to_edge_list_header_[VertID] = encode_vert_wedge_list_header(vwlh);

        vfi.is_null_vert = (false == vert_linked_to_seeded_wedge);
        vfi.vert_id = compacted_idx;
        uint num_all_verts = get_vert_count();
        uint addr_vert_to_filtered_vert = ssbo_vert_to_edge_list_header_addr__vert_to_selected_vert(VertID, num_all_verts);
        ssbo_vert_to_edge_list_header_[addr_vert_to_filtered_vert] = encode_filtered_vert_info(vfi);
    }
}
#endif









#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING_)
void store_wedge_quadric(uint FilteredEdgeID, Quadric q)
{
    vec4 data_0, data_1; 
    vec3 data_2; 
    encode_quadric(q, /*out*/data_0, data_1, data_2); 

    uint st_addr = FilteredEdgeID * 11u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_edge_quadric_data_[st_addr + i] = floatBitsToUint(data_0[i]);
    st_addr += 4u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_edge_quadric_data_[st_addr + i] = floatBitsToUint(data_1[i]);
    st_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        ssbo_edge_quadric_data_[st_addr + i] = floatBitsToUint(data_2[i]);
}
Quadric load_wedge_quadric(uint FilteredEdgeID)
{
    Quadric q; 

    uint ld_addr = FilteredEdgeID * 11u; 
    vec4 data_0, data_1;
    vec3 data_2;
    for (uint i = 0; i < 4; ++i)
        data_0[i] = uintBitsToFloat(ssbo_edge_quadric_data_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 4; ++i)
        data_1[i] = uintBitsToFloat(ssbo_edge_quadric_data_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        data_2[i] = uintBitsToFloat(ssbo_edge_quadric_data_[ld_addr + i]);

    q = decode_quadric(data_0, data_1, data_2); 

    return q; 
}

void store_vert_quadric(uint FilteredVertID, Quadric q)
{
    vec4 data_0, data_1; 
    vec3 data_2; 
    encode_quadric(q, /*out*/data_0, data_1, data_2); 

    uint st_addr = FilteredVertID * 11u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_0[i]);
    st_addr += 4u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_1[i]);
    st_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_2[i]);
}

Quadric load_vert_quadric(uint FilteredVertID)
{
    Quadric q; 

    uint ld_addr = FilteredVertID * 11u; 
    vec4 data_0, data_1;
    vec3 data_2;
    for (uint i = 0; i < 4; ++i)
        data_0[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 4; ++i)
        data_1[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        data_2[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);

    q = decode_quadric(data_0, data_1, data_2); 

    return q; 
}


void store_filtered_edge_normal(uint FilteredEdgeID, vec3 normal)
{
    ssbo_filtered_normal_edge_out_[FilteredEdgeID*3u + 0u] = normal.x;
    ssbo_filtered_normal_edge_out_[FilteredEdgeID*3u + 1u] = normal.y;
    ssbo_filtered_normal_edge_out_[FilteredEdgeID*3u + 2u] = normal.z;
}
vec3 load_filtered_edge_normal(uint FilteredEdgeID)
{
    return vec3(
        ssbo_filtered_normal_edge_in_[FilteredEdgeID*3u + 0u],
        ssbo_filtered_normal_edge_in_[FilteredEdgeID*3u + 1u],
        ssbo_filtered_normal_edge_in_[FilteredEdgeID*3u + 2u]
    );
}

void store_filtered_vert_normal(uint FilteredVertID, vec3 normal)
{
    ssbo_filtered_normal_vert_[FilteredVertID*3u + 0u] = normal.x;
    ssbo_filtered_normal_vert_[FilteredVertID*3u + 1u] = normal.y;
    ssbo_filtered_normal_vert_[FilteredVertID*3u + 2u] = normal.z;
}
vec3 load_filtered_vert_normal(uint FilteredVertID)
{
    return vec3(
        ssbo_filtered_normal_vert_[FilteredVertID*3u + 0u],
        ssbo_filtered_normal_vert_[FilteredVertID*3u + 1u],
        ssbo_filtered_normal_vert_[FilteredVertID*3u + 2u]
    );
}

uint addr_filtered_vert_pos(uint FilteredVertID, uint num_filtered_verts)
{
    uint subbuff_start = num_filtered_verts * 11u; 
    // subbuff_start = (((subbuff_start + 3u) >> 2u) << 2u); /* align to 4 */
    return subbuff_start + FilteredVertID * 3u;
}
void store_filtered_vert_pos(uint FilteredVertID, uint num_filtered_verts, vec3 pos)
{
    uint base_addr = addr_filtered_vert_pos(FilteredVertID, num_filtered_verts); 
    ssbo_vert_quadric_data_out_[base_addr + 0u] = floatBitsToUint(pos.x);
    ssbo_vert_quadric_data_out_[base_addr + 1u] = floatBitsToUint(pos.y);
    ssbo_vert_quadric_data_out_[base_addr + 2u] = floatBitsToUint(pos.z);
}
vec3 load_filtered_vert_pos(uint FilteredVertID, uint num_filtered_verts)
{
    uint base_addr = addr_filtered_vert_pos(FilteredVertID, num_filtered_verts); 
    return vec3(
        uintBitsToFloat(ssbo_vert_quadric_data_in_[base_addr + 0u]),
        uintBitsToFloat(ssbo_vert_quadric_data_in_[base_addr + 1u]),
        uintBitsToFloat(ssbo_vert_quadric_data_in_[base_addr + 2u])
    );
}

/*
 * uint ssbo_bnpr_mesh_pool_counters_.num_filtered_edges/verts
 * uint ssbo_edge_to_edges_[];
 * uint ssbo_vert_to_edge_list_header_[]; 
 * uint ssbo_selected_edge_to_edge_[]; 
 * uint ssbo_selected_vert_to_vert_[]; 
 * uint ssbo_vert_quadric_data_in_[];  <- (reuse) ssbo_bnpr_mesh_pool_
 * uint ssbo_vert_quadric_data_out_[]; <- (reuse) ssbo_mesh_buffer_reuse_0_
 * uint ssbo_edge_quadric_data_[];     <- (reuse) ssbo_vert_merged_id_
 * float ssbo_vbo_full_[]; rw
 * float ssbo_filtered_normal_vert_[];  <- (reuse) ssbo_mesh_buffer_reuse_1_
 * float ssbo_filtered_normal_edge_[];  <- (reuse) ssbo_vert_quadric_data_in_ at iter#0
 * ubo_view_matrices_
*/
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    
#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__EDGE_NORMAL)
    const uint FilteredEdgeID = idx.x; 
    const uint NumFilteredEdges = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
    uint num_all_edges = get_edge_count(); 

    bool valid_thread = (idx.x < NumFilteredEdges); 
    if (!valid_thread) return; /* quit if not valid thread */

    EdgeSelectionInfo fei = decode_edge_selection_info(ssbo_selected_edge_to_edge_[FilteredEdgeID]);
    uint wedge_id = fei.edge_id; 
    
    uvec4 vids_wedge;
    uint base_addr = wedge_id * 4u; 
    for (uint ivert = 0; ivert < 4u; ++ivert)
        vids_wedge[ivert] = ssbo_edge_to_vert_[base_addr+ivert];

    vec3 vpos_curr[4]; 
	for (uint ivert = 0; ivert < 4; ++ivert)
		vpos_curr[ivert] = ld_vbo(vids_wedge[ivert]);  
    
    /* transform matrices, see "common_view_lib.glsl" */ 
	mat4 view_to_world = ubo_view_matrices_.viewinv;
	bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

    vec3 normal_curr; 
    if (pcs_first_iter_ > 0)
    {
        normal_curr = compute_wedge_normal(vpos_curr[0], vpos_curr[1], vpos_curr[2], vpos_curr[3]); 
    }else{
        normal_curr = load_filtered_edge_normal(FilteredEdgeID); 
    }

    float weight_sum = 0.01f; 
    vec3 normal_filtered = normal_curr; 
    for (uint ibwedge = 0u; ibwedge < 4u; ++ibwedge)
    {
        AdjWedgeInfo awi = decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + ibwedge]);  
        uint ivert_adj = mark__center_wedge_to_oppo_vert__at_face(awi.iface_adj);
        
        uint vid_adj = ssbo_edge_to_vert_[awi.wedge_id*4u + ivert_adj];
        vec3 vpos_not_overlap = ld_vbo(vid_adj); 

        uint iface_overlap = mark__bwedge_to_face(ibwedge); 
        uvec3 iverts_overlap;
        iverts_overlap[0] = mark__border_wedge_to_oppo_vert(ibwedge);  
        iverts_overlap[1] = mark__vert_to_next_vert(iface_overlap, iverts_overlap[0]); 
        iverts_overlap[2] = mark__vert_to_next_vert(iface_overlap, iverts_overlap[1]); 

        /* Arrange verts to fit their order in the adjacent wedge */
        vec3 vpos_adj[4]; 
        if (ivert_adj == 0u)
        {
            vpos_adj[0] = vpos_not_overlap; 
            vpos_adj[1] = vpos_curr[iverts_overlap[2]]; 
            vpos_adj[2] = vpos_curr[iverts_overlap[0]]; 
            vpos_adj[3] = vpos_curr[iverts_overlap[1]]; 
        }else{ /*ivert_adj == 2u*/
            vpos_adj[0] = vpos_curr[iverts_overlap[0]];
            vpos_adj[1] = vpos_curr[iverts_overlap[1]];
            vpos_adj[3] = vpos_not_overlap; 
            vpos_adj[2] = vpos_curr[iverts_overlap[2]];
        }

        uint addr = ssbo_edge_to_edges_addr__edge_to_selected_edge(awi.wedge_id, num_all_edges); 
        EdgeSelectionInfo fei = decode_edge_selection_info(ssbo_edge_to_edges_[addr]);
        uint filtered_adj_wedge_id = fei.edge_id; 
        
        vec3 normal_bwedge;
        if ((pcs_first_iter_ > 0) || (fei.is_null_edge)) 
            normal_bwedge = compute_wedge_normal(vpos_adj[0], vpos_adj[1], vpos_adj[2], vpos_adj[3]); 
        else 
            normal_bwedge = load_filtered_edge_normal(filtered_adj_wedge_id); 
        

        vec3 adjusted_normal; 
        float weight; 
        bilateral_filter_wedge_normal(
            vpos_curr[0], vpos_curr[1], vpos_curr[2], vpos_curr[3], normal_curr, 
            vpos_adj[0],  vpos_adj[1],  vpos_adj[2],  vpos_adj[3],  normal_bwedge, 
            cam_pos_ws, 
            /* out */ adjusted_normal, weight
        ); 

        normal_filtered += weight * adjusted_normal; 
        weight_sum += weight; 
    }

    normal_filtered = normalize(normal_filtered / weight_sum); 
    if (pcs_use_normal_filtering_ == 0)
        normal_filtered = normal_curr; 

    if (valid_thread)
        store_filtered_edge_normal(FilteredEdgeID, normal_filtered); 

#endif

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__EDGE_QUADRIC)
    const uint FilteredEdgeID = idx.x; 
    const uint NumFilteredEdges = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
    bool valid_thread = (idx.x < NumFilteredEdges); 
    if (!valid_thread) return; /* quit if not valid thread */

    EdgeSelectionInfo fei = decode_edge_selection_info(ssbo_selected_edge_to_edge_[FilteredEdgeID]); 
    uint wedge_id = fei.edge_id; 
    
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

    vec3 wedge_normal = load_filtered_edge_normal(FilteredEdgeID); 
    Quadric q_e = compute_wedge_quadric(
        vpos_wedge[0], vpos_wedge[1], vpos_wedge[2], vpos_wedge[3], 
        cam_pos_ws, 
        0.05f, 
        wedge_normal
    ); 
    if (valid_thread)
        store_wedge_quadric(FilteredEdgeID, q_e); 
#endif


#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC)
    const uint FilteredVertID = idx.x; 
    const uint NumFilteredVerts = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < NumFilteredVerts); 
    if (!valid_thread) return; /* quit if not valid thread */
    
    uint num_all_edges = get_edge_count(); 
    uint num_all_verts = get_vert_count(); 

    VertFilteringInfo vfi = decode_filtered_vert_info(ssbo_selected_vert_to_vert_[FilteredVertID]); 
    uint vert_id_global = vfi.vert_id; 
    vec3 vpos = ld_vbo(vert_id_global);


    Quadric q_v; 
    Quadric q_v_filtered; 
    float weight_sum = .0f; 

    vec3 vpos_filtered; 
    vec3 delta_bilateral = vec3(.0f, .0f, .0f);
    float bilateral_weight_sum = .0f;  

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
    /* note: this inits at diagonals with param, other places cleared to .0 */
    q_v.quadric = mat4(0.0); 
    q_v.area = .0f; 
    q_v_filtered = q_v; 
    
    vec3 filtered_vert_normal = vec3(.0f, .0f, .0f); 
    float filtered_vert_normal_weight = .0f; 

    vpos_filtered = vpos; 
#else
    q_v = load_vert_quadric(FilteredVertID); 
    float dev_q = pcs_quadric_deviation_; 
    float dev_g = pcs_geodist_deviation_; 
    weight_sum = 1.0f; // compute_vert_quadric_weight(vpos, q_v, vpos, q_v, dev_q, dev_g); 
    q_v_filtered.quadric = weight_sum * q_v.quadric;
    q_v_filtered.area = q_v.area;

    vpos_filtered = load_filtered_vert_pos(FilteredVertID, NumFilteredVerts);
#endif


    { /* Rotate wedge around the vert */
        VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id_global]); 
        uint ivert_oppo = (vwlh.ivert == 1u) ? 3 : 1; 
        
        AdjWedgeInfo awi; 
        awi.iface_adj = (vwlh.ivert == 1u) ? 0 : 1; /* which face are we begin to rotate */
        awi.wedge_id  = vwlh.wedge_id; 

    #define MAX_WEDGE_ROTATES 16u
        uint rotate_step = 0u; 
        do {
            uint vid_j = ssbo_edge_to_vert_[awi.wedge_id*4u + ivert_oppo]; 
            vec3 vpos_j = ld_vbo(vid_j);

            uint wid_addr = ssbo_edge_to_edges_addr__edge_to_selected_edge(awi.wedge_id, num_all_edges); 
            EdgeSelectionInfo fei = decode_edge_selection_info(ssbo_edge_to_edges_[wid_addr]); 
            uint filtered_wedge_id = fei.edge_id; 

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
            Quadric q_e = load_wedge_quadric(filtered_wedge_id);
            q_v_filtered.area += q_e.area; 

            float w = compute_edge_quadric_weight(vpos, (vpos_j + vpos) * .5f, q_e);
            q_v_filtered.quadric += w * q_e.quadric;

            filtered_vert_normal += q_e.area * load_filtered_edge_normal(filtered_wedge_id); 
            filtered_vert_normal_weight += q_e.area; 
#else
            uint vid_addr = ssbo_vert_to_edge_list_header_addr__vert_to_selected_vert(vid_j, num_all_verts); 
            VertFilteringInfo vfi = decode_filtered_vert_info(ssbo_vert_to_edge_list_header_[vid_addr]); 
            uint filtered_vid_j = vfi.vert_id; 

            /* Vertex Pos Filtering */
            vec3 filtered_vpos_j;
            if (false == vfi.is_null_vert) 
                filtered_vpos_j = load_filtered_vert_pos(filtered_vid_j, NumFilteredVerts); 
            else
                filtered_vpos_j = .5f * (vpos + vpos_j); 

            vec3 wedge_normal_filtered = load_filtered_edge_normal(filtered_wedge_id);
            vec3 wedge_pos_filtered = .5f * (vpos_filtered + filtered_vpos_j); 
            delta_bilateral += wedge_normal_filtered * dot(wedge_normal_filtered, wedge_pos_filtered - vpos_filtered);
            bilateral_weight_sum += 1.0f; 

            /* Quadric Diffusion */
            float w = .0f; 
            Quadric q_oppo; 
            vec3 q_pos = vpos_j; 
            if (false == vfi.is_null_vert)
                q_oppo = load_vert_quadric(filtered_vid_j); 
            else{
                q_oppo = load_wedge_quadric(filtered_wedge_id);
                q_pos = .5f * (vpos + vpos_j); 
            }
            w = compute_vert_quadric_weight(vpos, q_v, q_pos, q_oppo, dev_q, dev_g);
            q_v_filtered.quadric += w * q_oppo.quadric;

#endif
            weight_sum += w; 


            uint iwedge_next = mark__cwedge_rotate_next(awi.iface_adj); 
            awi = decode_adj_wedge_info(ssbo_edge_to_edges_[awi.wedge_id*4u + iwedge_next]);
            ivert_oppo = mark__cwedge_to_beg_vert(awi.iface_adj); /* we are rotating around the end vert */

            rotate_step++; 
        } while (
            rotate_step < MAX_WEDGE_ROTATES 
            && awi.wedge_id != vwlh.wedge_id
        ); 
    }

    if (weight_sum > .0f)
        q_v_filtered.quadric = q_v_filtered.quadric / weight_sum; 
    if (valid_thread)
        store_vert_quadric(FilteredVertID, q_v_filtered); 

    if (bilateral_weight_sum > .0f)
    {
        // bilateral_weight_sum *= 3.0f; 
        delta_bilateral = delta_bilateral / bilateral_weight_sum;
        vpos_filtered += delta_bilateral; 
    }
    if (valid_thread)
        store_filtered_vert_pos(FilteredVertID, NumFilteredVerts, vpos_filtered);
    

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
    /*TODO: optimize normal storage*/
    filtered_vert_normal = normalize(filtered_vert_normal / filtered_vert_normal_weight); 
    if (valid_thread)
        store_filtered_vert_normal(FilteredVertID, filtered_vert_normal);
#endif

#endif 


#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__MOVE_VERTS)
    const uint FilteredVertID = idx.x; 
    const uint NumFilteredVerts = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
    uint num_all_verts = get_vert_count(); 

    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < NumFilteredVerts); 
    if (!valid_thread) return; /* quit if not valid thread */
    
    VertFilteringInfo vfi = decode_filtered_vert_info(ssbo_selected_vert_to_vert_[FilteredVertID]); 
    uint vert_id_global = vfi.vert_id; 
    if ((!vfi.is_unstable) || vfi.is_border) return; /* we only move unstable,non-border verts */

    /* vert pos from last filtering iter */
    vec3 vpos = ld_vbo(vert_id_global);
    vec3 vpos_orig = vpos; 
    /* bilateral filtered vert pos */
    vec3 vpos_filtered = load_filtered_vert_pos(FilteredVertID, NumFilteredVerts);

    /* Project vert pos to quadric surface */
    Quadric q_v = load_vert_quadric(FilteredVertID); 
    vec3 nv = load_filtered_vert_normal(FilteredVertID); 

    #if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__MOVE_VERTS_CONSTRAINED_SOLVE)
        mat3 A = mat3(q_v.quadric); /* slice upper 3x3 */
        vec3 b = vec3(q_v.quadric[3][0], q_v.quadric[3][1], q_v.quadric[3][2]); 
        float lambda = -(dot(vpos, A * vpos) + dot(b, nv)) / dot(nv, A * nv); 
        vpos += nv * lambda; 
    #else /* Directly solve by minimizing the quadrics */
        mat4 q_d = mat4(q_v.quadric); 
 
        mat4 q_reg = compute_plane_quadric(nv, vpos_orig); 
        q_d = q_d + pcs_positiion_regularization_scale_ * q_reg; /*damping*/

        dmat3 A = dmat3(q_d); 
        dvec3 b = dvec3(q_d[3][0], q_d[3][1], q_d[3][2]); 
        vpos = vec3(-inverse(A)*b); 
    #endif

    if (any(isnan(vpos)) || any(isinf(vpos)))
        vpos = vpos_orig; 

    if (valid_thread)
        store_vbo(vert_id_global, vpos);
#endif

}


#endif



