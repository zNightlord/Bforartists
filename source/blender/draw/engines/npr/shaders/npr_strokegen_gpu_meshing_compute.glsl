




#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE)
#define GLOBAL_COMPACTION_COUNTER__MESH_VERTS ssbo_bnpr_mesh_pool_counters_.num_verts
#define GLOBAL_COMPACTION_COUNTER__DBG_VERTS ssbo_bnpr_mesh_pool_counters_.num_edges

/* Hash Function for 3D Vert Pos ------------------------------------- 
 * code from https://jcgt.org/published/0009/03/02/supplementary.pdf
*/
/* https://www.pcg-random.org/ */
uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
uint iqint1(uint n)
{
    n = (n << 13u) ^ n;
    n = n * (n * n * 15731u + 789221u) + 1376312589u;
    return n;
}
uint xxhash32(uvec3 p)
{
    const uint PRIME32_2 = 2246822519u, PRIME32_3 = 3266489917u;
    const uint PRIME32_4 = 668265263u, PRIME32_5 = 374761393u;
    uint h32 = p.z + PRIME32_5 + p.x*PRIME32_3;
    h32 = PRIME32_4 * ((h32 << 17) | (h32 >> (32 - 17)));
    h32 += p.y * PRIME32_3;
    h32 = PRIME32_4 * ((h32 << 17) | (h32 >> (32 - 17)));
    h32 = PRIME32_2 * (h32 ^ (h32 >> 15));
    h32 = PRIME32_3 * (h32 ^ (h32 >> 13));
    return h32 ^ (h32 >> 16);
}


/* 3d position hash
 * Converting 1D hashing func into ND via nesting.
 * this is better than linear combination of N hashes, 
 * see Ch 4.2 of https://jcgt.org/published/0009/03/02/paper.pdf */
uint pcg_nested_3d(vec3 v)
{
    uvec3 v3 = floatBitsToUint(v);
    return pcg(v3.z + pcg(v3.y + pcg(v3.x))); 
}
uint iqint1_nested_3d(vec3 v)
{
    uvec3 v3 = floatBitsToUint(v);
    return iqint1(v3.z + iqint1(v3.y + iqint1(v3.x))); 
}
uint xxhash32_3d(vec3 p)
{
    uvec3 v3 = floatBitsToUint(p);
    return xxhash32(v3); 
}



#define NOT_FOUND 0xffffffffu
#define MAX_PROBE_STEPS 32u

uint VERT_HASH(vec3 v)
{
    return (((pcg_nested_3d(v) % pcs_hash_map_size_))); 
}
uint VERT_PROBE_NEXT(uint hid)      
{
    return (((pcg(hid)) % pcs_hash_map_size_)); 
}

#define CHECKSUM_EMPTY 0u
#define MAX_CHECKSUM_ITER 64u
uint VERT_CHECKSUM(vec3 v)
{
    uint c = xxhash32_3d(v.zyx);
    
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

uint spatial_map_insert_vtx(vec3 vpos, out bool inserted)
{
    /* Entry allocation */
    uint hash_id = VERT_HASH(vpos);
    uint checksum = VERT_CHECKSUM(vpos);
    inserted = false; 

    /* Probing the hash table. Allow limited number of searches */
    uint probe_step = 0u;
    while (probe_step < MAX_PROBE_STEPS)
    {
        uint stored_checksum = atomicCompSwap(
            ssbo_vert_spatial_map_headers_[hash_id], 
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
            hash_id = VERT_PROBE_NEXT(hash_id); 
            probe_step++; 
        }
    }

    return NOT_FOUND; 
}

uint spatial_map_search_vtx(vec3 pos)
{
    /* Entry Alloc */
    uint hash_id = VERT_HASH(pos); 
    uint checksum = VERT_CHECKSUM(pos);

    /* Probing the hash table. Allow limited number of searches */
    uint probe_step = 0u;
    while (probe_step < MAX_PROBE_STEPS)
    {
        uint stored_checksum = ssbo_vert_spatial_map_headers_[hash_id];

        if (stored_checksum == checksum)
            return hash_id; 

        /* collision, compute another index */
        hash_id = VERT_PROBE_NEXT(hash_id); 
        probe_step++; 
    }

    return NOT_FOUND; 
}



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

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRP)
    const uint st_hash_addr = idx.x; /* TODO: use GL buffer copy function rather than this */
    ssbo_vert_spatial_map_headers_[st_hash_addr] = CHECKSUM_EMPTY; 
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_BUILD_HASHMAP)
    bool valid_thread = (VertID < pcs_vert_count_);

    vec3 vpos = ld_vbo(VertID); 
    bool inserted = false; 
    uint hash_id = NOT_FOUND;
    if (valid_thread)
        hash_id = spatial_map_insert_vtx(vpos, inserted); 

    barrier(); /* must be outside of any dynamic control divergence */

    if (inserted) /* store merged vert id. safe since we ensure <=1 vert doing this */
    {
        ssbo_vert_spatial_map_payloads_[hash_id] = VertID; 
        atomicAdd(GLOBAL_COMPACTION_COUNTER__MESH_VERTS, 1); 
    }
#else

#if defined(_KERNEL_MULTICOMPILE__VERT_MERGE_DEDUPLICATE)
    bool valid_thread = (VertID < pcs_vert_count_);

    vec3 pos = ld_vbo(VertID); 
    uint hash_id = NOT_FOUND; 
    if (valid_thread)
        spatial_map_search_vtx(pos); 
    
    barrier(); /* must be outside of any dynamic control divergence */ 

    uint merged_vert_id = ssbo_vert_spatial_map_payloads_[hash_id]; 
    if (valid_thread && merged_vert_id != NOT_FOUND)
    {
        ssbo_vert_merged_id_[VertID] = merged_vert_id; 
        
        /* Find if different position hashed to the same position with the same checksum. 
            * vec3 merged_pos = ld_vbo(merged_vert_id); 
            * if (any(notEqual(merged_pos, pos)))
            atomicAdd(GLOBAL_COMPACTION_COUNTER__DBG_VERTS, 1);  */
    }
#endif
#endif
#endif /* ! _KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRP */
}

#endif /* _KERNEL_MULTICOMPILE__VERT_MERGE */



#if defined(_KERNEL_MULTICOMPILE__EDGE_LINK)

void main()
{
    
}

#endif
