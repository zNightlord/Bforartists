
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)

/* Inputs: 
uint pcs_flip_iter_;
uint pcs_flip_opti_goal_type_; 
uint pcs_vert_count_; 
uint ssbo_dyn_edge_count_[];
ssbo_dyn_mesh_counters_in_; 
uint ssbo_per_edge_flip_info_[]; 
uint ssbo_per_flip_edge_info_[]; 
float ssbo_vbo_full_[]; 

compact_flip_edge

#ifdef _KERNEL_MULTICOMPILE__EDGE_FLIP_COMPACT
.define("GLOBAL_COUNTER", "ssbo_edge_flip_counters_[pcs_flip_iter_].num_flip_edges_pass_1")
.define("CP_TAG", "flip_inital_selection")
#endif
*/


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



float calc_dihedral_angle(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{ /* TODO: not verified, can be buggy */
	vec3 v10 = v0 - v1;
	vec3 v13 = v3 - v1;
   	vec3 v12 = v2 - v1;

	vec3 n0 = normalize(cross(v13, v10));
	vec3 n2 = normalize(cross(v12, v13));
	
    /* https://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleMeshDerivativesCheatSheet.pdf */
    float dihedral = atan(dot(v13, cross(n0, n2)), dot(n0, n2));
    return dihedral; 
}


/* Edge flipping for different opti goals.
 * See "Automatic and High-quality Surface Mesh Generation for CAD Models" 
 */
#define PI 3.1415926535614f

#define EDGE_SPLIT_OPTI_VALENCE 0u
#define EDGE_SPLIT_OPTI_DELAUNAY 1u

bool should_edge_flip_common(vec3 p0, vec3 p1, vec3 p2, vec3 p3)
{
    /* only flip when edge is flat */
    float dihedral_angle = calc_dihedral_angle(p0, p1, p2, p3); 
    if (abs(dihedral_angle - PI) > PI * .07f /*12.6 degree*/)
    { 
        return false; 
    }
    return true; 
}

/* larger score == higher priority */
bool flip_score_larger(EdgeFlipInfo efi_0, uint wedge_id_0, EdgeFlipInfo efi_1, uint wedge_id_1)
{
    if (!efi_0.is_flip_ok) return false; 
    if (!efi_1.is_flip_ok) return true; 
    if (efi_0.score == efi_1.score)
        return wedge_id_0 > wedge_id_1; 
    return efi_0.score > efi_1.score; 
}
float comp_edge_flip_valence_score(
    vec3 p0, vec3 p1, vec3 p2, vec3 p3, 
    uvec4 valences, bvec4 border_edge
){ 
    if (!should_edge_flip_common(p0, p1, p2, p3))
        return -10000000.0f; /* no, fuck off */

    /* only flip when gets closer to optimal valence */
    vec4 valences_before = vec4(valences); 
    vec4 valences_after = valences_before + vec4(+1.0f, -1.0f, +1.0f, -1.0f); 

    vec4 opti_valence = border_edge ? 4.0f : 6.0f; 
    
    vec4 vec_to_opti_before = opti_valence - valences_before; 
    float dist_before = dot(vec_to_opti_before, vec_to_opti_before); 
    vec4 vec_to_opti_after = opti_valence - valences_after; 
    float dist_after = dot(vec_to_opti_after, vec_to_opti_after); 

    return dist_before - dist_after; 
}
float comp_edge_flip_delaunay_score(
    vec3 p0, vec3 p1, vec3 p2, vec3 p3, 
    uvec4 valences, bvec4 border_edge
){
    if (!should_edge_flip_common(p0, p1, p2, p3))
        return -10000000.0f; 

    /* Delaunay edge flip condition: two opposite angles sum larger than Pi */
    vec3 v21 = p2 - p1; 
    vec3 v23 = p2 - p3; 
    float angle_123 = acos(dot(v21, v23)); 

    vec3 v01 = p0 - p1; 
    vec3 v03 = p0 - p3; 
    float angle_103 = acos(dot(v01, v03)); 
    if (PI < (angle_123 + angle_103))
        return 1.0f; 
    return .0f; 
}






struct EdgeFlipInfo
{
    bool is_flip_ok; 
    float score; 
    uint id; 
}; 
uvec2 encode_edge_flip_info(EdgeFlipInfo efi)
{
    uvec2 efi_enc = uvec2(0, 0); 
    efi_enc.x = floatBitsToUint(efi.edge_len); 
    efi_enc.y = ((efi.id << 1u) | uint(efi.is_flip_ok));

    return efi_enc;  
}
EdgeFlipInfo decode_edge_flip_info(uvec2 efi_enc)
{
    EdgeFlipInfo efi; 
    efi.edge_len = uintBitsToFloat(efi_enc.x);
    efi.is_flip_ok = (0u != (efi_enc.y & 1u));
    efi_enc.y >>= 1u; 
    efi.id = efi_enc.y;  

    return efi; 
}
void store_per_flip_edge_info(uint flip_edge_id, EdgeFlipInfo psei)
{
    uvec2 psei_enc = encode_edge_flip_info(psei); 
    Store2(ssbo_per_flip_edge_info_, flip_edge_id, psei_enc);
}
EdgeFlipInfo load_per_flip_edge_info(uint flip_edge_id)
{
    uvec2 psei_enc; 
    Load2(ssbo_per_flip_edge_info_, flip_edge_id, psei_enc); 
    return decode_edge_flip_info(psei_enc); 
}
void store_per_edge_flip_info(uint wedge_id, EdgeFlipInfo pefi)
{
    uvec2 pefi_enc = encode_edge_flip_info(pefi); 
    Store2(ssbo_per_edge_flip_info_, wedge_id, pefi_enc); 
}
EdgeFlipInfo load_per_edge_flip_info(uint wedge_id)
{
    uvec2 pefi_enc;
    Load2(ssbo_per_edge_flip_info_, wedge_id, pefi_enc); 
    return decode_edge_flip_info(pefi_enc); 
}






#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_INIT)
struct FlipPrepContext
{
    uint vert_valence; 
}; 
bool vertex_1ring_visitor_prep_flip(CirculatorIterData iter, inout FlipPrepContext ctx)
{
    
}
void main()
{
    /* Dispatched at beginning of each flip operation */
    if (gl_GlobalInvocationID.x == 0u)
        ssbo_edge_flip_counters_[0].num_flip_edges_pass_1 = 0u;

    if (pcs_flip_opti_goal_type_ == EDGE_SPLIT_OPTI_VALENCE)
    { /* compute vertex valence, so dispatch per vertex */
        uint vert_id = gl_GlobalInvocationID.x; 
        uint num_verts  = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts; 
        
        VertWedgeListHeader vwlh = ssbo_vert_wedge_list_headers_[vert_id]; 
        FlipPrepContext ctx; 
        ctx.vert_valence = 0u; 


    }
}
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_COMPACT)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
    uint wedge_id = gl_GlobalInvocationID.x; 

    EdgeFlags ef = load_edge_flags(wedge_id); 
    uint v0 = ssbo_edge_to_vert_[wedge_id*4u + 0]; 
    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + 1]; 
    uint v2 = ssbo_edge_to_vert_[wedge_id*4u + 2]; 
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + 3]; 

    vec3 vpos[4] = { ld_vpos(v0), ld_vpos(v1), ld_vpos(v2), ld_vpos(v3) }; 

    float score = .0f; 
    if (pcs_flip_opti_goal_type_ == EDGE_SPLIT_OPTI_VALENCE)
    { 
        
        score = should_edge_flip_common(vpos[0], vpos[1], vpos[2], vpos[3]) ? 1.0f : -1.0f; 
    }
    
}
#endif



#if define(_KERNEL_MULTICOMPILE__EDGE_FLIP)
void main()
{
    uint flip_edge_id = gl_GlobalInvocationID.x;
    uint flip_edges = flip_counters_[flip_iter_].flip_edges_pass_1; 
    bool valid_thread = flip_edge_id < flip_edges;

}
#endif


