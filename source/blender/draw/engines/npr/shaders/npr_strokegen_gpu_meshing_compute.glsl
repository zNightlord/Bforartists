
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)


/* Hash Funcation Primitives ------------------------------------- 
 * code from https://jcgt.org/published/0009/03/02/supplementary.pdf
*/
/* https://www.pcg-random.org/ */
uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
/* http://www.jcgt.org/published/0009/03/02/ */
uvec3 pcg3d(uvec3 v) {

    v = v * 1664525u + 1013904223u;

    v.x += v.y*v.z;
    v.y += v.z*v.x;
    v.z += v.x*v.y;

    v ^= v >> 16u;

    v.x += v.y*v.z;
    v.y += v.z*v.x;
    v.z += v.x*v.y;

    return v;
}
/* Hashing 3d vert position
 * Converting 1D hashing func into ND via nesting.
 * this is better than linear combination of N hashes, 
 * see Ch 4.2 of https://jcgt.org/published/0009/03/02/paper.pdf */
uint pcg_nested_3d(vec3 v)
{
    uvec3 v3 = floatBitsToUint(v);
    return pcg(v3.z + pcg(v3.y + pcg(v3.x))); 
}


/* Setup Vertex Hash ------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE)

#define HASH_MAP_SIZE pcs_hash_map_size_
#define SSBO_HASH_MAP_HEADERS ssbo_vert_spatial_map_headers_
#define HASH_KEY_TYPE vec3

#define NOT_FOUND 0xffffffffu

#define CHECKSUM_EMPTY 0u
#define MAX_CHECKSUM_ITER 64u
#define MAX_PROBE_STEPS 32u

uint HASH_FUNC(vec3 v)
{
    return (((pcg_nested_3d(v) % pcs_hash_map_size_))); 
}
uint HASH_PROBE_NEXT(uint hid)      
{
    return (((pcg(hid)) % pcs_hash_map_size_)); 
}
uint CHECKSUM_FUNC(vec3 v)
{
    uvec3 v3 = floatBitsToUint(v);
    v3 = pcg3d(v3); 
    uint c = pcg(v3.z + pcg(v3.y + pcg(v3.x)));
    
    if (c == CHECKSUM_EMPTY) 
    { /* try to roll out a valid checksum */
        for (uint i = 0u; i < MAX_CHECKSUM_ITER; i++)
        {
            c = pcg(c); 
            if (c != CHECKSUM_EMPTY)
                return c; 
        }
    }

    return c; 
}
#endif



/* Setup Edge Hash ------------------------------------------------- */
#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY)

#define HASH_MAP_SIZE pcs_hash_map_size_
#define SSBO_HASH_MAP_HEADERS ssbo_edge_index_map_headers_
#define HASH_KEY_TYPE uvec2

#define NOT_FOUND 0xffffffffu

#define CHECKSUM_EMPTY 0u
#define MAX_CHECKSUM_ITER 64u
#define MAX_PROBE_STEPS 24u

uvec2 MAKE_HASH_KEY(uint v0, uint v1)
{ /* note: ensure a deterministic order */
    return uvec2(min(v0, v1), max(v0, v1)); 
}
uint HASH_FUNC(uvec2 v)
{
    v = MAKE_HASH_KEY(v.x, v.y); 
    return pcg(v.x + pcg(v.y)) % pcs_hash_map_size_; 
}
uint HASH_PROBE_NEXT(uint hid)      
{
    return (((pcg(hid)) % pcs_hash_map_size_)); 
}
uint CHECKSUM_FUNC(uvec2 v)
{
    v = MAKE_HASH_KEY(v.x, v.y); 
    uint c = pcg(v.y + pcg(v.x + 13u));
    
    if (c == CHECKSUM_EMPTY) 
    { /* try to roll out a valid checksum */
        for (uint i = 0u; i < MAX_CHECKSUM_ITER; i++)
        {
            c = pcg(c); 
            if (c != CHECKSUM_EMPTY)
                return c; 
        }
    }

    return c; 
}
#endif



#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE) || defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY)
uint FUNC_HASHMAP_INSERT(HASH_KEY_TYPE vpos, out bool inserted)
{
    /* Entry allocation */
    uint hash_id = HASH_FUNC(vpos);
    uint checksum = CHECKSUM_FUNC(vpos);
    inserted = false; 

    /* Probing the hash table. Allow limited number of searches */
    uint probe_step = 0u;
    while (probe_step < MAX_PROBE_STEPS)
    {
        uint stored_checksum = atomicCompSwap(
            SSBO_HASH_MAP_HEADERS[hash_id], 
            CHECKSUM_EMPTY, 
            checksum
        );

        if (stored_checksum == CHECKSUM_EMPTY || stored_checksum == checksum)
        { /* The entry at that index was empty or already occupied */
            inserted = (stored_checksum == CHECKSUM_EMPTY); 
            return hash_id; 
        }

        if (stored_checksum != checksum)
        { /* collision, compute another index */
            hash_id = HASH_PROBE_NEXT(hash_id); 
            probe_step++; 
        }
    }

    return NOT_FOUND; 
}

uint FUNC_HASHMAP_SEARCH(HASH_KEY_TYPE pos)
{
    /* Entry Alloc */
    uint hash_id = HASH_FUNC(pos); 
    uint checksum = CHECKSUM_FUNC(pos);

    /* Probing the hash table. Allow limited number of searches */
    uint probe_step = 0u;
    while (probe_step < MAX_PROBE_STEPS)
    {
        uint stored_checksum = SSBO_HASH_MAP_HEADERS[hash_id];

        if (stored_checksum == checksum)
            return hash_id; 

        /* collision, compute another index */
        hash_id = HASH_PROBE_NEXT(hash_id); 
        probe_step++; 
    }

    return NOT_FOUND; 
}
#endif



#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE) || defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY)
vec3 ld_vbo(uint vert)
{
	uint base_addr = vert * 3; 
	return vec3(ssbo_vbo_full_[base_addr], ssbo_vbo_full_[base_addr+1], ssbo_vbo_full_[base_addr+2]); 
}
#endif



#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE)

#define GLOBAL_COMPACTION_COUNTER__MESH_VERTS ssbo_bnpr_mesh_pool_counters_.num_verts
#define GLOBAL_COMPACTION_COUNTER__DBG_VERTS ssbo_bnpr_mesh_pool_counters_.num_edges

void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    const uint VertID = idx.x; 
    bool valid_thread = (idx.x < pcs_vert_count_); /* Do not use VertID here since it's offseted with current mesh batch  */

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRAP)
    const uint st_hash_addr = idx.x; /* TODO: use GL buffer copy function rather than this */
    SSBO_HASH_MAP_HEADERS[st_hash_addr] = CHECKSUM_EMPTY; 
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BUILD_HASHMAP)

    vec3 vpos = ld_vbo(VertID); 
    bool inserted = false; 
    uint hash_id = NOT_FOUND;
    if (valid_thread)
        hash_id = FUNC_HASHMAP_INSERT(vpos, inserted); 

    barrier(); /* must be outside of any dynamic control divergence */

    if (inserted) /* store merged vert id. safe since we ensure <=1 vert doing this */
    {
        ssbo_vert_spatial_map_payloads_[hash_id] = VertID; 
        /* atomicAdd(GLOBAL_COMPACTION_COUNTER__MESH_VERTS, 1); */
    }

    if (valid_thread) 
    { 
        /* avoid hash search for non-duplicated verts */
        ssbo_vert_merged_id_[VertID] = inserted ? VertID : NOT_FOUND;

        /* Init vert flags */
        VertFlags vf = init_vert_flags(!inserted/*duplicated?*/);
        store_vert_flags(VertID, vf);  
    }
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_DEDUPLICATE)
    uint merged_vert_id = ssbo_vert_merged_id_[VertID];
    bool is_depli_vert = merged_vert_id == NOT_FOUND; 

    /* Only apply expensive hash search for duplicated verts */
    uint hash_id = NOT_FOUND; 
    if (is_depli_vert && valid_thread)
    { 
        vec3 pos = ld_vbo(VertID); 
        hash_id = FUNC_HASHMAP_SEARCH(pos); 
    }
    barrier(); /* must be outside of any dynamic control divergence */ 
    if (is_depli_vert && valid_thread)
        merged_vert_id = ssbo_vert_spatial_map_payloads_[hash_id]; 
    
    
    if (is_depli_vert && valid_thread)
    {
        merged_vert_id = (merged_vert_id != NOT_FOUND) ? merged_vert_id 
            : VertID /*fucked up but at least we can keep it sane here */; 
        ssbo_vert_merged_id_[VertID] = merged_vert_id; 
        
        /* Debug snippet: Find if different position hashed to the same position with the same checksum. 
         * if (merged_vert_id != NOT_FOUND) {
         *  vec3 merged_pos = ld_vbo(merged_vert_id); 
         *  if (any(notEqual(merged_pos, pos)))
         *  atomicAdd(GLOBAL_COMPACTION_COUNTER__DBG_VERTS, 1);  
         * }
         */
    }
#endif
#endif
#endif /* ! _KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRAP */
}

#endif /* _KERNEL_MULTICOMPILE__VERT_MERGE */



#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY)

#define GLOBAL_COMPACTION_COUNTER__NUM_EDGES ssbo_bnpr_mesh_pool_counters_.num_edges

void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    bool valid_thread = (idx.x < pcs_edge_count_); /* Do not use EdgeID here since it's offseted with current mesh batch */

    const uint EdgeID = idx.x; 


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_BOOTSTRAP)
    const uint st_hash_addr = idx.x; /* TODO: use GL buffer copy function rather than this */
    /* Init hash map */
    SSBO_HASH_MAP_HEADERS[st_hash_addr] = CHECKSUM_EMPTY;
    /* Prep for atomic min op on hash payloads */
    ssbo_edge_spatial_map_payloads_[st_hash_addr] = 0xffffffffu; 
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_HASHING)
    uvec4 vid; 
    uint base_addr = EdgeID * 4u; 
    vid[0] = ssbo_edge_to_vert_[base_addr+0];
    vid[1] = ssbo_edge_to_vert_[base_addr+1]; 
    vid[2] = ssbo_edge_to_vert_[base_addr+2];
    vid[3] = ssbo_edge_to_vert_[base_addr+3];
    bool is_border = line_adj_is_border_edge(vid);
   
    /* fetch merged vtx id */
    for (uint i = 0u; i < 4u; i++) {
        uint vert_id = vid[i]; 
        uint merged_id = ssbo_vert_merged_id_[vert_id]; 
        vid[i] = merged_id; 
    }
    uvec2 edge_verts = uvec2(vid[1], vid[2]); /*v1, 2 is on the edge*/

    /* Store e-v linkage, also transform to wedge index order */
    uvec4 wing_verts = line_adj_to_wing_verts(vid); 
    if (valid_thread) {
        ssbo_edge_to_vert_[base_addr+0] = wing_verts[0];
        ssbo_edge_to_vert_[base_addr+1] = wing_verts[1]; 
        ssbo_edge_to_vert_[base_addr+2] = wing_verts[2];
        ssbo_edge_to_vert_[base_addr+3] = wing_verts[3];
    }

    /* Hash Edge using 2 deduplicated vertex ids*/
    bool inserted = false; 
    uint hash_id = NOT_FOUND;
    if (valid_thread)
        hash_id = FUNC_HASHMAP_INSERT(edge_verts, inserted); 

    barrier(); /* must be outside of any dynamic control divergence */

    if (valid_thread && (hash_id != NOT_FOUND)) 
    { /* Use the minimum edge id as merged id, so that it's deterministic.
       * This is required for temporal coherency.
      */
        atomicMin(ssbo_edge_spatial_map_payloads_[hash_id], EdgeID); 
        /* Cache temporary data */
        ssbo_edge_to_edges_[EdgeID*4u + 0u] = hash_id;
        ssbo_edge_to_edges_[EdgeID*4u + 1u] = (is_border ? 1u : 0u);
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_HASHING_FINISH)
    /* Load temp data */
    uvec2 temp_data = uvec2(
        ssbo_edge_to_edges_[EdgeID*4u + 0u], 
        ssbo_edge_to_edges_[EdgeID*4u + 1u]
    );
    const uint hash_id   = temp_data.x;
    const bool valid_hash = (NOT_FOUND != hash_id); 
    const bool is_border = ((temp_data.y & 1u) == 1u);

    /* Select minimum id among edges with the same hash */
    uint deduplicated_edge_id = ssbo_edge_spatial_map_payloads_[hash_id]; 
    bool inserted = valid_hash && (deduplicated_edge_id == EdgeID); 
    if (valid_thread)
    { 
        /* temp cache inserted flag here to avoid hash search for deduplication */
        ssbo_edge_to_edges_[EdgeID*4u + 0u] = inserted ? 1u : 0u;

        /* Init edge states */
        EdgeFlags ef = init_edge_flags(!inserted, is_border);
        store_edge_flags(EdgeID, ef); 
        
        /* Setup v-e link */
        uvec4 wing_verts; 
        uint base_addr = EdgeID * 4u; 
        wing_verts[0] = ssbo_edge_to_vert_[base_addr+0];
        wing_verts[1] = ssbo_edge_to_vert_[base_addr+1]; 
        wing_verts[2] = ssbo_edge_to_vert_[base_addr+2];
        wing_verts[3] = ssbo_edge_to_vert_[base_addr+3];
        uvec2 iverts_cwedge = mark__wedge_to_verts(4u); 
        if (inserted)
            for (uint iiverts = 0u; iiverts < 2u; iiverts++)
            { /* Put an arbitrary neighbor wedge id for each vert */
                uint ivert = iverts_cwedge[iiverts]; 
                uint vert_id = wing_verts[ivert];
                ssbo_vert_to_edge_list_header_[vert_id] = 
                    encode_vert_wedge_list_header(VertWedgeListHeader(EdgeID, ivert)); 
            }
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_FIND_ADJ)
    uint WedgeId = EdgeID; 
    uint temp_cached_hash_info = ssbo_edge_to_edges_[WedgeId*4u + 0u];
    bool is_wedge_unique_hashed = (temp_cached_hash_info == 1u);  
    bool valid_edge = (is_wedge_unique_hashed); 

    uvec4 vids_wedge; /*bwedge:4 border wedges*/
    uint base_addr = WedgeId * 4u; 
    { /* load edge adjacency */
        vids_wedge[0] = ssbo_edge_to_vert_[base_addr+0];
        vids_wedge[1] = ssbo_edge_to_vert_[base_addr+1];
        vids_wedge[2] = ssbo_edge_to_vert_[base_addr+2];
        vids_wedge[3] = ssbo_edge_to_vert_[base_addr+3];
    }
    bool is_border_wedge = wing_verts_is_border_edge(vids_wedge); 

    /* Search for adj wedges */
    uvec2 vids_bwedge[4]; 
    for (uint iwedge = 0; iwedge < 4; ++iwedge) 
    {
        uvec2 iverts_iwedge = mark__wedge_to_verts(iwedge); 
        vids_bwedge[iwedge] = uvec2(vids_wedge[iverts_iwedge[0]], vids_wedge[iverts_iwedge[1]]); 
    }
    
    uvec4 hashmap_index = uvec4(NOT_FOUND, NOT_FOUND, NOT_FOUND, NOT_FOUND); 
    for (uint i = 0u; i < 4u; i++)
        if (valid_thread)
            hashmap_index[i] = FUNC_HASHMAP_SEARCH(vids_bwedge[i]); 
    
    barrier(); /* must be outside of any dynamic control divergence */ 

    uvec4 wedge_id;
    uint st_base_addr = EdgeID * 4u; 
    if (valid_thread)
        for (uint iwedge = 0u; iwedge < 4u; iwedge++)
        {
            /* load deduplicated id of adjacent wedges  */
            if (valid_edge && (hashmap_index[iwedge] != NOT_FOUND))
                wedge_id[iwedge] = ssbo_edge_spatial_map_payloads_[hashmap_index[iwedge]]; 
            else
                wedge_id[iwedge] = NULL_EDGE; 
            
            /* find non-overlapping iface in adj. wedge */
            uint iface_adj_non_overlapping = 0u; 
            if (wedge_id[iwedge] != NULL_EDGE)
            { /* one opposite vert of the adj. wedge is inside curr edge, while another vert is outside */
                uint oppo_vert_id = vids_wedge[mark__border_wedge_to_oppo_vert(iwedge)];
                
                uint ld_addr_adj_wedge_data = wedge_id[iwedge] * 4u; 
                uvec2 vids_oppo_adj_wedge = uvec2(
                    ssbo_edge_to_vert_[ld_addr_adj_wedge_data + mark__center_wedge_to_oppo_vert__at_face(0)], 
                    ssbo_edge_to_vert_[ld_addr_adj_wedge_data + mark__center_wedge_to_oppo_vert__at_face(1)]
                );

                if (oppo_vert_id != vids_oppo_adj_wedge[1])
                    iface_adj_non_overlapping = 1u;
            } 
            
            AdjWedgeInfo awi; 
            awi.wedge_id = wedge_id[iwedge]; 
            awi.iface_adj = iface_adj_non_overlapping; 
            ssbo_edge_to_edges_[EdgeID*4u + iwedge] = encode_adj_wedge_info(awi); 
            /* note: == NOT_FOUND when this edge is a duplicated one */
        }


    /* Adjust Vert-to-Wedge link to ensure border edges are always prefered */
    if (valid_edge && valid_thread && is_border_wedge)
    { /* end vert points the border edge */
        uint ivert_end = mark__cwedge_to_end_vert(mark__border_iface_mainfold()); 
        uint vtxid_end = vids_wedge[ivert_end]; 
        ssbo_vert_to_edge_list_header_[vtxid_end] = encode_vert_wedge_list_header(VertWedgeListHeader(EdgeID, ivert_end)); 
    }

    /* Initialize subdiv tree roots */
    if (valid_thread)
    { 
        LoopSubdEdgeTreeNode edge_tree_ptr = init_loop_subd_tree_root(EdgeID, wedge_id); 
        ssbo_subd_edge_tree_node_[EdgeID] = encode_loop_subd_tree_node(edge_tree_ptr);  
    }

    /* initialze dynamic topology counters */
    if (WedgeId == 0u)
    {
        ssbo_dyn_mesh_counters_in_.num_edges = 0u;
        ssbo_dyn_mesh_counters_in_.num_verts = 0u;

        ssbo_dyn_mesh_counters_out_.num_edges = 0u;
        ssbo_dyn_mesh_counters_out_.num_verts = 0u;

        ssbo_edge_split_counters_[0].num_split_edges_pass_1 = 0u;
        ssbo_edge_split_counters_[0].num_split_edges = 0u; 
        ssbo_edge_split_counters_[1].num_split_edges_pass_1 = 0u;
        ssbo_edge_split_counters_[1].num_split_edges = 0u; 

        ssbo_bnpr_vert_debug_draw_args_.vertex_len = 0u; 
        ssbo_bnpr_vert_debug_draw_args_.instance_len = 1u;
        ssbo_bnpr_vert_debug_draw_args_.vertex_first = 0u;
        ssbo_bnpr_vert_debug_draw_args_.base_index = 0u; 
        ssbo_bnpr_vert_debug_draw_args_.instance_first_indexed = 0u; 
    }

#endif

}

#endif