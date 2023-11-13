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
    uint wedge_id = ssbo_vert_to_edge_list_header_[VertID]; 
    uint vid_ivert1 = ssbo_edge_to_vert_[wedge_id*4 + 1]; 
    /* Select how to fwd rotate the center wedge */
    uint ivert = (VertID == vid_ivert1) ? 1u : 3u;
    uint iface_beg_rot_fwd = (ivert == 1u) ? 0u : 1u; 

    bool vert_linked_to_seeded_wedge = false; 
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

    if (valid_thread)
    {
        VertWedgeListHeader vwlh;
        vwlh.wedge_id = wedge_id; 
        vwlh.ivert = ivert; 
        ssbo_vert_to_edge_list_header_[VertID] = encode_vert_wedge_list_header(vwlh);
    }
}
#endif




#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING_)
void store_wedge_quadric(uint FilteredEdgeID, Quadric q)
{
    vec4 data_0, data_1; 
    vec3 data_2; 
    encode_quadric(q_e, /*out*/data_0, data_1, data_2); 

    uint st_addr = FilteredEdgeID * 11u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_wedge_quadric_data_[st_addr + i] = floatBitsToUint(data_0[i]);
    st_addr += 4u; 
    for (uint i = 0; i < 4; ++i)
        ssbo_wedge_quadric_data_[st_addr + i] = floatBitsToUint(data_1[i]);
    st_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        ssbo_wedge_quadric_data_[st_addr + i] = floatBitsToUint(data_2[i]);
}
Quadric load_wedge_quadric(uint FilteredEdgeID)
{
    Quadric q; 

    uint ld_addr = FilteredEdgeID * 11u; 
    vec4 data_0, data_1;
    vec3 data_2;
    for (uint i = 0; i < 4; ++i)
        data_0[i] = uintBitsToFloat(ssbo_wedge_quadric_data_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 4; ++i)
        data_1[i] = uintBitsToFloat(ssbo_wedge_quadric_data_[ld_addr + i]);
    ld_addr += 4u;
    for (uint i = 0; i < 3; ++i)
        data_2[i] = uintBitsToFloat(ssbo_wedge_quadric_data_[ld_addr + i]);

    q = decode_quadric(data_0, data_1, data_2); 

    return q; 
}

void store_vert_quadric(uint FilteredVertID, Quadric q)
{
    vec4 data_0, data_1; 
    vec3 data_2; 
    encode_quadric(q_e, /*out*/data_0, data_1, data_2); 

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
    ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 0u] = normal.x;
    ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 1u] = normal.y;
    ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 2u] = normal.z;
}
vec3 load_filtered_edge_normal(uint FilteredEdgeID)
{
    return vec3(
        ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 0u],
        ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 1u],
        ssbo_filtered_normal_edge_[FilteredEdgeID*3u + 2u]
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

vec3 ld_vbo(uint GlobalVertID)
{
	uint base_addr = GlobalVertID * 3; 
	return vec3(ssbo_vbo_full_[base_addr], ssbo_vbo_full_[base_addr+1], ssbo_vbo_full_[base_addr+2]); 
}
vec3 store_vbo(uint GlobalVertID, vec3 vpos)
{
    uint base_addr = GlobalVertID * 3; 
    ssbo_vbo_full_[base_addr]   = vpos.x;
    ssbo_vbo_full_[base_addr+1] = vpos.y;
    ssbo_vbo_full_[base_addr+2] = vpos.z;
}

/*
 * uint ssbo_bnpr_mesh_pool_counters_.num_filtered_edges/verts
 * uint ssbo_edge_to_edges_[];
 * uint ssbo_vert_to_edge_list_header_[]; 
 * uint ssbo_filtered_edge_to_edge_[]; 
 * uint ssbo_filtered_vert_to_vert_[]; 
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
    

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__EDGE_QUADRIC)
    const uint FilteredEdgeID = idx.x; 
    const uint NumFilteredEdges = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges; 
    /* Do not use EdgeID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < NumFilteredEdges); 
    if (!valid_thread) return; /* quit if not valid thread */

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

    vec3 wedge_normal; 
    Quadric q_e = compute_wedge_quadric(
        vpos_wedge[0], vpos_wedge[1], vpos_wedge[2], vpos_wedge[3], cam_pos_ws, 
        /*out*/ wedge_normal
    ); 
    if (valid_thread)
    {
        store_wedge_quadric(FilteredEdgeID, q_e); 
        store_filtered_edge_normal(FilteredEdgeID, wedge_normal)
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC)
    const uint FilteredVertID = idx.x; 
    const uint NumFilteredVerts = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < NumFilteredVerts); 
    if (!valid_thread) return; /* quit if not valid thread */
    
    uint vert_id_global = ssbo_filtered_vert_to_vert_[FilteredVertID]; 
    vec3 vpos = ld_vbo(vert_id_global);

    Quadric q_v; 
    float weight_sum = .0f; 
#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
    /* note: this inits at diagonals with param, other places cleared to .0 */
    q_v.quadric = mat4(0.0); 
    q_v.area = .0f; 
    
    vec3 filtered_vert_normal = vec3(.0f, .0f, .0f); 
    float filtered_vert_normal_weight = .0f; 
#else
    q_v = load_vert_quadric(FilteredVertID); 
    weight_sum += compute_vert_quadric_weight(vpos, q_v, vpos, q_v); 
#endif
    Quadric q_v_filtered = q_v; 

    
    { /* Rotate wedge around the vert */
        VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id_global]); 
        uint ivert_oppo = (vmlh.ivert == 1u) ? 3 : 1; 
        
        AdjWedgeInfo awi; 
        awi.iface_adj = (vmlh.ivert == 1u) ? 0 : 1; /* which face are we begin to rotate */
        awi.wedge_id  = vwlh.wedge_id; 

    #define MAX_WEDGE_ROTATES 16u
        uint rotate_step = 0u; 
        do {
            vec3 vpos_oppo = ld_vbo(ssbo_edge_to_vert_[awi.wedge_id*4u + ivert_oppo]);

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
            Quadric q_e = load_wedge_quadric(awi.wedge_id);
            q_v_filtered.area += q_e.area; 

            float w = compute_edge_quadric_weight(vpos, (vpos_oppo + vpos) * .5f, q_e);
            q_v_filtered.quadric += w * q_e.quadric;

            filtered_vert_normal += q_e.area * load_filtered_edge_normal(awi.wedge_id); 
            filtered_vert_normal_weight += q_e.area; 
#else
            Quadric q_oppo = load_vert_quadric(awi.wedge_id); 
            float w = compute_vert_quadric_weight(vpos, q_v, vpos_oppo, q_oppo); 
            q_v_filtered.quadric += w * q_oppo.quadric; 
#endif
            weight_sum += w; 


            uint iwedge_next = mark__cwedge_rotate_next(awi.iface_adj); 
            awi = decode_adj_wedge_info(ssbo_edge_to_edges_[awi.wedge_id*4u + iwedge_next]);
            ivert_oppo = mark__cwedge_to_beg_vert(awi.iface_adj); /* we are rotating around the end vert */

            rotate_step++; 
        } while (
            rotate_step < MAX_WEDGE_ROTATES 
            && awi.wedge_id != wedge_id
        ); 
    }

    if (weight_sum > .0f)
        q_v_filtered.quadric = q_v_filtered.quadric / weight_sum; 
    if (valid_thread)
        store_vert_quadric(FilteredVertID, q_v_filtered); 

#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT)
    if (valid_thread)
        store_filtered_vert_normal(FilteredVertID, filtered_vert_normal / filtered_vert_normal_weight);
#endif

#endif


#if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__MOVE_VERTS)
    const uint FilteredVertID = idx.x; 
    const uint NumFilteredVerts = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts; 
    /* Do not use VertID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < NumFilteredVerts); 
    if (!valid_thread) return; /* quit if not valid thread */
    
    Quadric q_v = load_vert_quadric(FilteredVertID); 
    vec4 nv = load_filtered_vert_normal(FilteredVertID); 
    
    uint vert_id_global = ssbo_filtered_vert_to_vert_[FilteredVertID]; 
    vec3 vpos = ld_vbo(vert_id_global);

    #if defined(_KERNEL_MULTICOMPILE__MESH_FILTERING__MOVE_VERTS_CONSTRAINED_SOLVE)
        mat3 A = mat3(q_v.quadric); /* slice upper 3x3 */
        vec3 b = vec3(q_v.quadric[3][0], q_v.quadric[3][1], q_v.quadric[3][2]); 
        float lambda = -(dot(vpos, A * vpos) + dot(b, nv)) / dot(nv, A * nv); 
        vpos += nv * lambda; 
    #else /* Directly solve by minimizing the quadrics */
        dmat4 q_d = dmat4(q_v.quadric); 
        dmat3 A = dmat3(q_d); 
        dvec3 b = dvec3(q_d.quadric[3][0], q_d.quadric[3][1], q_d.quadric[3][2]); 
        vpos = vec3(-inverse(A)*b); 
    #endif

    if (valid_thread)
        store_vbo(vert_id_global, vpos);
#endif



}
#endif
