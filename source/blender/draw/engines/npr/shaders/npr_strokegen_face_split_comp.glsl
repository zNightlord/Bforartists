
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_allocation_lib.glsl)

// pcs_split_iter_
// SSBOData_StrokeGenFaceSplitCounters   ssbo_face_split_counters_[]
// SSBOData_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_in_[]
// SSBOData_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_out_[]
// uint ssbo_per_face_split_info_[]

vec3 ld_vpos(uint vtx_id)
{
	vec3 vpos; 
    Load3(ssbo_vbo_full_, vtx_id, vpos);

    return vpos; 
}
void st_vpos(uint vtx_id, vec3 vpos)
{
    Store3(ssbo_vbo_full_, vtx_id, vpos);
}

struct FaceSplitInfo
{
    uint wedge_id;
    uint iface;  
}; 
uint encode_face_split_info(FaceSplitInfo fsi)
{
    return ((fsi.wedge_id << 1u) | (fsi.iface)); 
}
FaceSplitInfo decode_face_split_info(uint enc)
{
    FaceSplitInfo fsi; 
    fsi.wedge_id = enc >> 1u; 
    fsi.iface = enc & 1u; 
    return fsi; 
}
void store_per_split_face_info(uint split_face_id, FaceSplitInfo psfi)
{
    ssbo_per_face_split_info_[split_face_id] = encode_face_split_info(psfi);
}
FaceSplitInfo load_per_split_face_info(uint split_face_id)
{
    return decode_face_split_info(ssbo_per_face_split_info_[split_face_id]); 
} 


#define SPLIT_FACE_COUNTER (ssbo_face_split_counters_[pcs_split_iter_].num_split_faces)


#if defined(_KERNEL_MULTICOMPILE__FACE_SPLIT_INIT)
void main()
{ /* Dispatched at beginning of each split sequence */
    if (gl_GlobalInvocationID.x == 0)
        ssbo_face_split_counters_[0].num_split_faces = 0u;
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FACE_SPLIT_WORK_GEN)


// #define AC_TAG split_face

// DECL_LDS_COUNTERS_PER_LANE(AC_TAG)
// #define LDS_COUNTERS_PER_LANE CAT(LDS_counters_per_lane_, AC_TAG)
// DECL_LDS_COUNTERS_PER_LANE_SUM(AC_TAG)
// #define LDS_COUNTERS_PER_LANE_SUM CAT(LDS_counters_per_lane_sum_, AC_TAG)
// DECL_LDS_ALLOC_BLK_COUNTER(AC_TAG)
// #define LDS_ALLOC_BLK_COUNTER CAT(LDS_alloc_blk_counter_, AC_TAG)
// DECL_LDS_ALLOC_BLOCK_OFFSET(AC_TAG)
// #define LDS_ALLOC_BLOCK_OFFSET CAT(LDS_alloc_blk_offset_, AC_TAG)

// DECL_ALLOCATION_FUNC(AC_TAG, ssbo_face_split_counters_[pcs_split_iter_].num_split_faces)


bool should_wedge_generate_iface(uint wedge_id, uint iface, AdjWedgeInfo w[4], EdgeFlags efs[4])
{ // Returns true if the wedge should generate the face
    bool gen_face = true;

    uvec2 iwedges_fi = (iface == 0u) ? uvec2(1, 2) : uvec2(3, 0); 
    for (uint iiw = 0u; iiw < 2u; ++iiw)
    { // Only wedge with the highest id generates the face 
        uint iwedge = iwedges_fi[iiw]; 

        if (wedge_id < w[iwedge].wedge_id)
            gen_face = false;

        if (!efs[iwedge].selected)
            gen_face = false; 
    }

    return gen_face; 
}

void main()
{
    const uint groupIdx = gl_LocalInvocationID.x; 
    
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 
    
    if (gl_GlobalInvocationID.x == 0u)
    {
        ssbo_dyn_mesh_counters_in_.num_edges = ssbo_dyn_mesh_counters_out_.num_edges; 
        ssbo_dyn_mesh_counters_in_.num_verts = ssbo_dyn_mesh_counters_out_.num_verts; 
    }

    AdjWedgeInfo w[4]; 
    w[0] = decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 0u]);
    w[1] = decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 1u]);
    w[2] = decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 2u]);
    w[3] = decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 3u]);

    EdgeFlags ef = load_edge_flags(wedge_id); 
    bool is_split_ok = valid_thread 
        && (ef.selected)
        && (!ef.sel_border)
        && (!ef.dupli)
        && (!ef.border) 
        && (!ef.del_by_split)
        && (!ef.del_by_collapse); 

    EdgeFlags efs[4] = {
        load_edge_flags(w[0].wedge_id), 
        load_edge_flags(w[1].wedge_id), 
        load_edge_flags(w[2].wedge_id), 
        load_edge_flags(w[3].wedge_id)
    };
    if (valid_thread)
        update_edge_flags__reset_face_split_new_edge(wedge_id, ef); // reset to false for sqrt-3 subdiv

    bvec2 gen_face; 
    for (uint iface = 0; iface < 2u; ++iface)
    {
        gen_face[iface] = valid_thread && is_split_ok && should_wedge_generate_iface(wedge_id, iface, w, efs); 
    }
    uint num_gen_faces = uint(gen_face[0]) + uint(gen_face[1]);  


    uint alloc_offset = 0u; 
    alloc_offset = alloc_split_face(groupIdx, num_gen_faces);  


    // Store per-face split info
    uint split_face_id = alloc_offset; 
    for (uint iface = 0u; iface < 2u; ++iface)
        if (gen_face[iface] && valid_thread)
        {
            FaceSplitInfo pfsi;
            pfsi.wedge_id = wedge_id; 
            pfsi.iface = iface;

            store_per_split_face_info(split_face_id, pfsi);
            split_face_id++; 
        }
}
#endif



#if defined(_KERNEL_MULTICOMPILE__FACE_SPLIT_EXECUTE)
void Store4_EEAdj(uint w, uvec4 adj_wedges, uvec4 adj_iface_adjs)
{
    uvec4 data;
    data[0] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[0], adj_iface_adjs[0])); 
    data[1] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[1], adj_iface_adjs[1]));
    data[2] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[2], adj_iface_adjs[2]));
    data[3] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[3], adj_iface_adjs[3]));

    Store4(ssbo_edge_to_edges_, w, data); 
}

void main()
{
    /* Update dynamic mesh counters */
    if (gl_GlobalInvocationID.x == 0u) 
    {
        /* accumulate global mesh counters */
        ssbo_dyn_mesh_counters_out_.num_edges = 
            ssbo_dyn_mesh_counters_in_.num_edges + 3u * SPLIT_FACE_COUNTER;
        ssbo_dyn_mesh_counters_out_.num_verts =
            ssbo_dyn_mesh_counters_in_.num_verts + SPLIT_FACE_COUNTER;
        
        /* prep local counters for next split iteration */
        ssbo_face_split_counters_[pcs_split_iter_ + 1].num_split_faces = 0u;
    } 


    const uint groupId = gl_LocalInvocationID.x; 

    uint split_face_id = gl_GlobalInvocationID.x;
    uint num_split_faces = SPLIT_FACE_COUNTER; 
    bool valid_thread = split_face_id < num_split_faces;

    uint num_existing_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_in_.num_edges; 
    uint num_existing_verts = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts; 

    if (!valid_thread) return; 


    /*                                                  
     *         P2           New edges: E0,1,2                
     *        /| \             
     *       / |  \            
     *      /  E2  \                                    
     *    L2   |    L1      New Vertex: Vn                            
     *    /   _Vn_   \                                                             
     *   /  E0    E1  \     Old edges: L0,1,2                                               
     *  / _/  --->  \_ \    
     * P0 ---- L0 ---- P1                                          
    */
    FaceSplitInfo pfsi = load_per_split_face_info(split_face_id); 

    /* Existing Edges */ 
    uvec2 iwedges_fi = (pfsi.iface == 0u) ? uvec2(1, 2) : uvec2(3, 0); 
    AdjWedgeInfo L[3] = {
        AdjWedgeInfo(pfsi.wedge_id, pfsi.iface == 0u ? 1u : 0u),
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfsi.wedge_id*4u + iwedges_fi[0u]]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfsi.wedge_id*4u + iwedges_fi[1u]])
    };
    /* Existing Verts */ 
    uvec3 P;
    P[0] = ssbo_edge_to_vert_[pfsi.wedge_id*4u + mark__cwedge_to_beg_vert(pfsi.iface)];
    P[1] = ssbo_edge_to_vert_[pfsi.wedge_id*4u + mark__cwedge_to_end_vert(pfsi.iface)];
    P[2] = ssbo_edge_to_vert_[pfsi.wedge_id*4u + mark__center_wedge_to_oppo_vert__at_face(pfsi.iface)];

    /* New Edges  x3 */
    uvec3 E = (num_existing_edges + split_face_id * 3u).xxx + uvec3(0u, 1u, 2u);
    /* New Vertex x1 */
    uint Vn = num_existing_verts + split_face_id; 


    /* Store new edges */
    /*           E0                        E1                        E2
     *           P2                        P2                        P2                        
     *          / |                        | \                      /| \                
     *         /  |                        |  \                    / |  \               
     *        /   E2                       E2  \                  /  E2  \              
     *      L2 f1 |                        | f1 L1              L2 f1| f0 L1            
     *      /   _Vn_                      _Vn_   \              /   _Vn_   \            
     *     /  E0 f0 E1                  E0 f0 E1  \            /  E0    E1  \           
     *    / _/  --->  \_              _/  --->  \_ \          / _/        \_ \          
     *   P0 ---- L0 ---- P1        P0 ---- L0 ---- P1        P0              P1                                
     *   ==================        ==================        ==================       
     *   { p2, p0, p1, vn }        { p2, vn, p0, p1 }        { p0, vn, p1, p2 }   { v0~3 }       
     *   { l2, l0, e1, e2 }        { e2, e0, l0, l1 }        { e0, e1, l1, l2 }   { w0~3 }
     *   { ~~, ~~, f1, f0 }        { f1, f1, ~~, ~~ }        { f0, f0, ~~, ~~ }   { w0~3.iface_adj }
    */
    Store4(ssbo_edge_to_vert_, E[0], uvec4(P[2], P[0], P[1], Vn));
    Store4_EEAdj(E[0], uvec4(L[2].wedge_id, L[0].wedge_id, E[1], E[2]), uvec4(L[2].iface_adj, L[0].iface_adj, 1, 0)); 
    Store4(ssbo_edge_to_vert_, E[1], uvec4(P[2], Vn, P[0], P[1])); 
    Store4_EEAdj(E[1], uvec4(E[2], E[0], L[0].wedge_id, L[1].wedge_id), uvec4(1, 1, L[0].iface_adj, L[1].iface_adj)); 
    Store4(ssbo_edge_to_vert_, E[2], uvec4(P[0], Vn, P[1], P[2]));
    Store4_EEAdj(E[2], uvec4(E[0], E[1], L[1].wedge_id, L[2].wedge_id), uvec4(0, 0, L[1].iface_adj, L[2].iface_adj));

    EdgeFlags ef_new = init_edge_flags__new_face_split_edge();
    store_edge_flags(E[0], ef_new);
    store_edge_flags(E[1], ef_new);
    store_edge_flags(E[2], ef_new);

    /* Adjust old edges */ 
    /*                                                  
     *         P2           
     *        /| \          
     *       / |  \         
     *      /  E2  \        
     *    L2   |    L1      
     *    /   _Vn_   \                                                             
     *   /  E0    E1  \     Old edges:                                                      
     *  / _/  --->  \_ \    Adjust ev-links: l0/l1/l2 points to Vn
     * P0 ---- L0 ---- P1   Adjust ee-links: l0/l1/l2 points to (e1f1,e0f1)/(e2f1,e1f0)/(e0f0,e2f0)                                                
    */
    AdjWedgeInfo ee_links[6] = {
        /* prev */             /* next */
        AdjWedgeInfo(E[0], 1), AdjWedgeInfo(E[1], 1), // L0
        AdjWedgeInfo(E[1], 0), AdjWedgeInfo(E[2], 1), // L1
        AdjWedgeInfo(E[2], 0), AdjWedgeInfo(E[0], 0)  // L2
    }; 
    for (uint k = 0; k < 3; ++k)
    { 
        /* Adjust ev-links */
        uint iface_overlap = L[k].iface_adj == 1u ? 0u : 1u;
        uint ivert_oppo = mark__center_wedge_to_oppo_vert__at_face(iface_overlap); 
        ssbo_edge_to_vert_[L[k].wedge_id*4u + ivert_oppo] = Vn;

        /* Adjust ee-links */
        uvec2 ibwedges = {
            mark__cwedge_rotate_back(iface_overlap), 
            mark__cwedge_rotate_next(iface_overlap)
        }; 
        ssbo_edge_to_edges_[L[k].wedge_id*4u + ibwedges[0]] = encode_adj_wedge_info(ee_links[2*k]); 
        ssbo_edge_to_edges_[L[k].wedge_id*4u + ibwedges[1]] = encode_adj_wedge_info(ee_links[2*k+1]); 
    }

    /* Store new vertex */
    VertWedgeListHeader vwlh; 
    vwlh.wedge_id = E[1]; 
    vwlh.ivert = 1u; /* check the graph for E1 above */
    ssbo_vert_to_edge_list_header_[Vn] = encode_vert_wedge_list_header(vwlh); 

    vec3 vpos[3] = {ld_vpos(P[0]), ld_vpos(P[1]), ld_vpos(P[2])}; 
    vec3 vpos_new = (vpos[0] + vpos[1] + vpos[2]) / 3.0f; 
    st_vpos(Vn, vpos_new); 

    VertFlags vf = init_vert_flags__new_split_face(); 
    store_vert_flags(Vn, vf); 


}
#endif
