
/* Inputs: 
float ssbo_vbo_full_[]

.define("EDGE_FLAGS_INCLUDED", "1")
.define("EDGE_FLAGS_INCLUDED", "1")
uint ssbo_vert_flags_[]
uint ssbo_edge_flags_[]

uint ssbo_per_edge_split_info_[]
uint ssbo_per_split_edge_info_[]

uint pcs_split_iter_; 
uint ssbo_dyn_edge_count_[]; 
uint ssbo_dyn_vert_count_[]; 
*/

vec3 ld_vpos(uint vtx_id)
{
	uint base_addr = vtx_id * 3; 
	return vec3(
        ssbo_vbo_full_[base_addr], 
        ssbo_vbo_full_[base_addr+1], 
        ssbo_vbo_full_[base_addr+2]
    ); 
}
void st_vpos(uint vtx_id, vec3 vpos)
{
    uint base_addr = vtx_id * 3; 
    ssbo_vbo_full_[base_addr]   = vpos.x;
    ssbo_vbo_full_[base_addr+1] = vpos.y;
    ssbo_vbo_full_[base_addr+2] = vpos.z;
}


struct EdgeSplitInfo
{
    bool should_split; 
    float edge_len; 
    uint id; 
}; 
uvec2 encode_edge_split_info(EdgeSplitInfo esi)
{
    uvec2 esi_enc = 0; 
    esi_enc.x = floatBitsToUint(esi.edge_len); 
    esi_enc.y = ((id << 1u) | esi.should_split);

    return esi_enc;  
}
EdgeSplitInfo decode_edge_split_info(uvec2 esi_enc)
{
    EdgeSplitInfo esi; 
    esi.edge_len = uintBitsToFloat(esi_enc.x);
    esi.should_split = (esi_enc.y & 1u);
    esi_enc.y >>= 1u; 
    esi.id = esi_enc.y;  

    return esi; 
}

void store_per_split_edge_info(uint split_edge_id, EdgeSplitInfo psei)
{
    uvec2 psei_enc = encode_edge_split_info(psei); 
    ssbo_per_split_edge_info_[split_edge_id*2] = psei_enc[0]; 
    ssbo_per_split_edge_info_[split_edge_id*2 + 1] = psei_enc[1]; 
}
EdgeSplitInfo load_per_split_edge_info(uint split_edge_id)
{
    uvec2 psei_enc; 
    psei_enc[0] = ssbo_per_split_edge_info_[split_edge_id*2]; 
    psei_enc[1] = ssbo_per_split_edge_info_[split_edge_id*2 + 1]; 
    return decode_edge_split_info(psei_enc); 
}
void store_per_edge_split_info(uint wedge_id, EdgeSplitInfo pesi)
{
    uvec2 pesi_enc = encode_edge_split_info(pesi); 
    ssbo_per_edge_split_info_[wedge_id*2] = pesi_enc[0]; 
    ssbo_per_edge_split_info_[wedge_id*2 + 1] = pesi_enc[1]; 
}
void load_per_edge_split_info(uint wedge_id)
{
    uvec2 pesi_enc; 
    pesi_enc[0] = ssbo_per_edge_split_info_[wedge_id*2]; 
    pesi_enc[1] = ssbo_per_edge_split_info_[wedge_id*2 + 1]; 
    return decode_edge_split_info(pesi_enc); 
}



bool split_score_larger(EdgeSplitInfo esi_0, EdgeSplitInfo esi_1)
{
}





#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_COMPACT)
void main()
{




}
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT)
void main()
{
/* Parallel Edge Split
 *          v0                   v0          
 *         /  \                 /||\         
 *        /    \               / || \        
 *       W      3             W  ||  3       
 *      0        W           0   E0   W      
 *     /          \         /    ||    \     
 *    /            \       /     ||     \    
 *  v1 ---- W4 ---- v3   v1 =E1= v4 =E3= v3  
 *    \            /       \     ||     /    
 *     \          /         \    ||    /     
 *      W        2           W   E2   2      
 *       1      W             1  ||  W       
 *        \    /               \ || /        
 *         \  /                 \||/         
 *          v2                   v2          
 */
#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_RESOLVE_COLLISION)
#endif

#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXECUTE)
#endif

 }
#endif