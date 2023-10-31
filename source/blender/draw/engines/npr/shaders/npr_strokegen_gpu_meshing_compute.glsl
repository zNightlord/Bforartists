

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



#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE)

#define GLOBAL_COMPACTION_COUNTER__MESH_VERTS ssbo_bnpr_mesh_pool_counters_.num_verts
#define GLOBAL_COMPACTION_COUNTER__DBG_VERTS ssbo_bnpr_mesh_pool_counters_.num_edges
vec3 ld_vbo(uint vert)
{
	uint base_addr = vert * 3; 
	return vec3(ssbo_vbo_full_[base_addr], ssbo_vbo_full_[base_addr+1], ssbo_vbo_full_[base_addr+2]); 
}


void main()
{
	const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    const uint VertID = idx.x; 

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRAP)
    const uint st_hash_addr = idx.x; /* TODO: use GL buffer copy function rather than this */
    SSBO_HASH_MAP_HEADERS[st_hash_addr] = CHECKSUM_EMPTY; 
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BUILD_HASHMAP)
    bool valid_thread = (VertID < pcs_vert_count_);

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
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_DEDUPLICATE)
    bool valid_thread = (VertID < pcs_vert_count_);

    vec3 pos = ld_vbo(VertID); 
    uint hash_id = NOT_FOUND; 
    if (valid_thread)
        hash_id = FUNC_HASHMAP_SEARCH(pos); 
    
    barrier(); /* must be outside of any dynamic control divergence */ 

    uint merged_vert_id = ssbo_vert_spatial_map_payloads_[hash_id]; 
    if (valid_thread)
    {
        merged_vert_id = (merged_vert_id != NOT_FOUND) ? merged_vert_id 
            : VertID/*fucked up but at least we can keep it sane here */; 
        ssbo_vert_merged_id_[VertID] = merged_vert_id; 
        
        /* Find if different position hashed to the same position with the same checksum. 
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
    const uint EdgeID = idx.x; 

    bool valid_thread = (EdgeID < pcs_edge_count_); 


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_BOOTSTRAP)
    const uint st_hash_addr = idx.x; /* TODO: use GL buffer copy function rather than this */
    SSBO_HASH_MAP_HEADERS[st_hash_addr] = CHECKSUM_EMPTY;
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_HASHING)
    uvec4 vid; 
    uint base_addr = EdgeID * 4u; 
    vid[0] = ssbo_edge_to_vert_[base_addr+0];
    vid[1] = ssbo_edge_to_vert_[base_addr+1]; 
    vid[2] = ssbo_edge_to_vert_[base_addr+2];
    vid[3] = ssbo_edge_to_vert_[base_addr+3];
   
    /* fetch merged id */
    for (uint i = 0u; i < 4u; i++) {
        uint vert_id = vid[i]; 
        uint merged_id = ssbo_vert_merged_id_[vert_id]; 
        vid[i] = merged_id; 
    }
    uvec2 edge_verts = uvec2(vid[1], vid[2]); /*v1, 2 is on the edge*/

    /* store merged id */
    if (valid_thread) {
        ssbo_edge_to_vert_[base_addr+0] = vid[0];
        ssbo_edge_to_vert_[base_addr+1] = vid[1]; 
        ssbo_edge_to_vert_[base_addr+2] = vid[2];
        ssbo_edge_to_vert_[base_addr+3] = vid[3];
    }

    /* Hash */
    bool inserted = false; 
    uint hash_id = NOT_FOUND;
    if (valid_thread)
        hash_id = FUNC_HASHMAP_INSERT(edge_verts, inserted); 

    barrier(); /* must be outside of any dynamic control divergence */

    if (inserted) 
    { /* store edge id. safe since we ensure <=1 vert doing this */
        ssbo_edge_spatial_map_payloads_[hash_id] = EdgeID; 
        /*atomicAdd(GLOBAL_COMPACTION_COUNTER__NUM_EDGES, 1);*/ /*debug only*/
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_FIND_ADJ)
    uvec4 edge_adj_verts;  
    uint base_addr = EdgeID * 4u; 
    { /* load edge adjacency */
        edge_adj_verts[0] = ssbo_edge_to_vert_[base_addr+0];
        edge_adj_verts[1] = ssbo_edge_to_vert_[base_addr+1];
        edge_adj_verts[2] = ssbo_edge_to_vert_[base_addr+2];
        edge_adj_verts[3] = ssbo_edge_to_vert_[base_addr+3];
    }
    /*    v0
     *   /  \
     *  /    \                
     * v1====v2 <-- curr edge            
     *  \    /                
     *   \  /                     
     *    v3    winding 012, 321                
    */

    /* Deduplicate */
    uvec2 curr_edge = uvec2(edge_adj_verts[1], edge_adj_verts[2]); 
    uint curr_edge_hash_id = NOT_FOUND; 
    uint curr_edge_dedup_id = NOT_FOUND; 
    if (valid_thread)
    {
        curr_edge_hash_id = FUNC_HASHMAP_SEARCH(curr_edge);
        curr_edge_dedup_id = ssbo_edge_spatial_map_payloads_[curr_edge_hash_id];
    }
    bool valid_edge = (curr_edge_hash_id != NOT_FOUND)
        && (curr_edge_dedup_id == EdgeID); 

    /* Search for adj edges */
    uvec2 adj_edges[4] = 
    {
        uvec2(edge_adj_verts[0], edge_adj_verts[1]), 
        uvec2(edge_adj_verts[2], edge_adj_verts[0]), 
        uvec2(edge_adj_verts[3], edge_adj_verts[2]),
        uvec2(edge_adj_verts[1], edge_adj_verts[3])
    };
    uvec4 hashmap_index = uvec4(NOT_FOUND, NOT_FOUND, NOT_FOUND, NOT_FOUND); 
    for (uint i = 0u; i < 4u; i++)
        if (valid_thread)
            hashmap_index[i] = FUNC_HASHMAP_SEARCH(adj_edges[i]); 
    
    barrier(); /* must be outside of any dynamic control divergence */ 

    uvec4 edge_id = uvec4(NOT_FOUND, NOT_FOUND, NOT_FOUND, NOT_FOUND);  
    uint st_base_addr = EdgeID * 4u; 
    if (valid_thread)
        for (uint i = 0u; i < 4u; i++)
        {
            if (valid_edge)
                edge_id[i] = ssbo_edge_spatial_map_payloads_[hashmap_index[i]]; 
            ssbo_edge_to_edges_[EdgeID+0] = edge_id[i]; 
            /* note: == NOT_FOUND when this edge is a duplicated one */
        }
#endif

}

#endif
