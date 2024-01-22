
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
compact_flip_edge_I

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
struct EdgeFlipInfo
{
    bool is_flip_ok; 
    float score; 
    uint id; 
}; 
uvec2 encode_edge_flip_info(EdgeFlipInfo efi)
{
    uvec2 efi_enc = uvec2(0, 0); 
    efi_enc.x = floatBitsToUint(efi.score); 
    efi_enc.y = ((efi.id << 1u) | uint(efi.is_flip_ok));

    return efi_enc;  
}
EdgeFlipInfo decode_edge_flip_info(uvec2 efi_enc)
{
    EdgeFlipInfo efi; 
    efi.score = uintBitsToFloat(efi_enc.x);
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


struct VertexEdgeFlipInfo
{
    uint valence; 
}; 
void store_ssbo_vertex_edge_flip_info(uint vert_id, VertexEdgeFlipInfo vefi)
{
    uint vefi_enc = vefi.valence; 
    ssbo_vertex_edge_flip_info_[vert_id] = vefi_enc;
}
VertexEdgeFlipInfo load_ssbo_vertex_edge_flip_info(uint vert_id)
{
    uint vefi_enc = ssbo_vertex_edge_flip_info_[vert_id]; 
    
    VertexEdgeFlipInfo vefi; 
    vefi.valence = vefi_enc; 
    return vefi; 
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

#define EDGE_FLIP_OPTI_VALENCE 0u
#define EDGE_FLIP_OPTI_DELAUNAY 1u

bool should_edge_flip_common(vec3 p0, vec3 p1, vec3 p2, vec3 p3)
{
    /* only flip when edge is flat */
    // float dihedral_angle = calc_dihedral_angle(p0, p1, p2, p3); // buggy
    {
        vec3 v10 = p0 - p1;
        vec3 v13 = p3 - p1;
        vec3 v12 = p2 - p1;

        vec3 n0 = normalize(cross(v13, v10));
        vec3 n2 = normalize(cross(v12, v13));

        if (acos(dot(n0, n2)) > PI * .07f /*12.6 degree*/)
            return false; 
    }

    /* wedge quad must be convex */
    vec3 v01 = p0 - p1;
    vec3 v21 = p2 - p1;
    vec3 v31 = p3 - p1;
    float angle_013 = acos(dot(normalize(v01), normalize(v31)));
    float angle_213 = acos(dot(normalize(v21), normalize(v31)));
    float angle_012 = angle_013 + angle_213; // triangle inner angle always < pi, but not for quad angle
    if (angle_012 > PI * .6f) // < .6Pi, this also avoids sliver triangle after flipping
    {
        return false;
    }
    vec3 v03 = p0 - p3;
    vec3 v13 = -v31; 
    vec3 v23 = p2 - p3;
    float angle_031 = acos(dot(normalize(v03), normalize(v13)));
    float angle_231 = acos(dot(normalize(v23), normalize(v13)));
    float angle_032 = angle_031 + angle_231;
    if (angle_032 > PI * .6f)
    {
        return false;
    }

    return true; 
}


uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
/* larger score == higher priority */
bool flip_score_larger(EdgeFlipInfo efi_0, uint wedge_id_0, EdgeFlipInfo efi_1, uint wedge_id_1)
{
    if (!efi_0.is_flip_ok) return false; 
    if (!efi_1.is_flip_ok) return true; 
    if (efi_0.score == efi_1.score)
        return pcg(wedge_id_0) > pcg(wedge_id_1); 
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

    vec4 opti_valence; 
    opti_valence.x = border_edge.x ? 4.0f : 6.0f; 
    opti_valence.y = border_edge.y ? 4.0f : 6.0f; 
    opti_valence.z = border_edge.z ? 4.0f : 6.0f; 
    opti_valence.w = border_edge.w ? 4.0f : 6.0f; 
    
    vec4 vec_to_opti_before = opti_valence - valences_before; 
    float dist_before = dot(vec_to_opti_before, vec_to_opti_before); 
    vec4 vec_to_opti_after = opti_valence - valences_after; 
    float dist_after = dot(vec_to_opti_after, vec_to_opti_after); 

    return dist_before - dist_after; 
}
float comp_edge_flip_delaunay_score(
    vec3 p0, vec3 p1, vec3 p2, vec3 p3
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
    return (angle_123 + angle_103) - PI; 
}










#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_INIT)
struct FlipValidateContext
{
    uint vert_valence; 
}; 
bool vertex_1ring_visitor_prep_flip(CirculatorIterData iter, inout FlipValidateContext ctx)
{
    ctx.vert_valence = ctx.vert_valence + 1u; 
    return true; 
}
void main()
{
    /* Dispatched at beginning of each flip operation */
    if (gl_GlobalInvocationID.x == 0u)
        ssbo_edge_flip_counters_[pcs_flip_iter_].num_flip_edges_pass_1 = 0u;

    if (pcs_flip_opti_goal_type_ == EDGE_FLIP_OPTI_VALENCE)
    { /* compute vertex valence, so dispatch per vertex */
        uint vert_id = gl_GlobalInvocationID.x; 
        uint num_verts  = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts;
        bool valid_thread = vert_id < num_verts;  
        
        if (valid_thread)
        {
            bool rot_fwd = true; 
            VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
            FlipValidateContext ctx; 
            ctx.vert_valence = 0u; 

            VE_CIRCULATOR(vwlh, vertex_1ring_visitor_prep_flip, ctx, rot_fwd); 

            store_ssbo_vertex_edge_flip_info(vert_id, VertexEdgeFlipInfo(ctx.vert_valence));
        }
    }
}
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_COMPACT)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 

    EdgeFlags ef = load_edge_flags(wedge_id); 

    uint v0 = ssbo_edge_to_vert_[wedge_id*4u + 0]; 
    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + 1]; 
    uint v2 = ssbo_edge_to_vert_[wedge_id*4u + 2]; 
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + 3]; 

    vec3 vpos[4] = { ld_vpos(v0), ld_vpos(v1), ld_vpos(v2), ld_vpos(v3) }; 

    
    bool is_flip_ok = (v0 != v2) 
        && (ef.selected)
        && (!ef.dupli)
        && (!ef.border) 
        && (!ef.del_by_split)
        && (!ef.del_by_collapse)
        && valid_thread; 
    float score = .0f; 
    /* Edge flip for optimizing valence */
    if (pcs_flip_opti_goal_type_ == EDGE_FLIP_OPTI_VALENCE)
    { 
        VertexEdgeFlipInfo vefi[4] = {
            load_ssbo_vertex_edge_flip_info(v0), 
            load_ssbo_vertex_edge_flip_info(v1),
            load_ssbo_vertex_edge_flip_info(v2),
            load_ssbo_vertex_edge_flip_info(v3)
        };
        AdjWedgeInfo w[4] = {
            decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 0u]), 
            decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 1u]),
            decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 2u]),
            decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 3u])
        }; 
        EdgeFlags ef[4] = {
            load_edge_flags(w[0].wedge_id), 
            load_edge_flags(w[1].wedge_id), 
            load_edge_flags(w[2].wedge_id), 
            load_edge_flags(w[3].wedge_id)
        }; 

        score = comp_edge_flip_valence_score(
            vpos[0], vpos[1], vpos[2], vpos[3], 
            uvec4(vefi[0].valence, vefi[1].valence, vefi[2].valence, vefi[3].valence), 
            bvec4(ef[0].border,    ef[1].border,    ef[2].border,    ef[3].border)
        );
    }
    /* Edge flip for optimizing delaunay property */
    if (pcs_flip_iter_ == EDGE_FLIP_OPTI_DELAUNAY)
    {
        score = comp_edge_flip_delaunay_score(
            vpos[0], vpos[1], vpos[2], vpos[3]
        );
    }
    
    is_flip_ok = is_flip_ok && (score > 0.0f); 
    uint flip_edge_id = compact_flip_inital_selection(is_flip_ok, groupId); 

    if (is_flip_ok)
    {
        EdgeFlipInfo pfei; 
        pfei.is_flip_ok = true; 
        pfei.score = score; 
        pfei.id = wedge_id; 
        store_per_flip_edge_info(flip_edge_id, pfei); 
    }

    if (valid_thread)
    {
        EdgeFlipInfo pefi; 
        pefi.is_flip_ok = is_flip_ok; 
        pefi.score = score; 
        pefi.id = is_flip_ok ? flip_edge_id : NULL_EDGE; 
        store_per_edge_flip_info(wedge_id, pefi);
    }
}
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP)

#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_VALIDATE)
struct FlipValidateContext
{ /* rotate around v0, cancel flip if connected to v2 (causes topo inconsistency) */
    uint v2; 
    EdgeFlipInfo pfei; 
    bool is_flip_valid;  
}; 
bool vertex_1ring_visitor_validate_flip(CirculatorIterData iter, inout FlipValidateContext ctx)
{
    uint wi = iter.awi.wedge_id; 
    uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
    uint vi = ssbo_edge_to_vert_[wi + ivert_vi]; 
    
    ctx.is_flip_valid = ctx.is_flip_valid && (vi != ctx.v2); 

    return ctx.is_flip_valid; 
}
#endif

#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_EXECUTE)
void Store4_EVAdj(uint wedge, uvec4 verts)
{
    Store4(ssbo_edge_to_vert_, wedge, verts); 
}
void Store4_EEAdj(uint wedge, AdjWedgeInfo awi[4])
{
    uvec4 data;
    data[0] = encode_adj_wedge_info(awi[0]); 
    data[1] = encode_adj_wedge_info(awi[1]);
    data[2] = encode_adj_wedge_info(awi[2]);
    data[3] = encode_adj_wedge_info(awi[3]);

    Store4(ssbo_edge_to_edges_, wedge, data); 
}
#endif

void main()
{
    uint flip_edge_id = gl_GlobalInvocationID.x;
    uint num_flip_edges = ssbo_edge_flip_counters_[pcs_flip_iter_].num_flip_edges_pass_1; 
    bool valid_thread = flip_edge_id < num_flip_edges;

    EdgeFlipInfo pfei = load_per_flip_edge_info(flip_edge_id); 

#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_VALIDATE)
    uint v0, v2; 
    v0 = ssbo_edge_to_vert_[pfei.id*4u + 0u];
    v2 = ssbo_edge_to_vert_[pfei.id*4u + 2u];

    bool ve_rot_fwd = true; 
    FlipValidateContext ctx; 
    ctx.v2 = v2; 
    ctx.pfei = pfei; 
    ctx.is_flip_valid = true; 
    VE_CIRCULATOR(VertWedgeListHeader(v0, 0u), vertex_1ring_visitor_validate_flip, ctx, ve_rot_fwd); 

    if (valid_thread && !ctx.is_flip_valid)
    { /* flip canceled */
        pfei.is_flip_ok = false;
        store_per_flip_edge_info(flip_edge_id, pfei);

        EdgeFlipInfo pefi; 
        pefi.is_flip_ok = false; 
        pefi.score = pfei.score;
        pefi.id = flip_edge_id;
        store_per_edge_flip_info(pfei.id, pefi); 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_RESOLVE_CONFLICT)
    AdjWedgeInfo w[4] = {
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfei.id*4u + 0u]), 
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfei.id*4u + 1u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfei.id*4u + 2u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[pfei.id*4u + 3u])
    }; 
    EdgeFlipInfo pefi[4] = {
        load_per_edge_flip_info(w[0].wedge_id), 
        load_per_edge_flip_info(w[1].wedge_id), 
        load_per_edge_flip_info(w[2].wedge_id), 
        load_per_edge_flip_info(w[3].wedge_id)
    }; 

    EdgeFlipInfo pfei_new = pfei; 
    for (uint i = 0u; i < 4u; ++i)
        if (flip_score_larger(pefi[i], w[i].wedge_id, pfei, pfei.id))
        {
            pfei_new.is_flip_ok = false;
            break; 
        }

    if (valid_thread && !pfei_new.is_flip_ok)
    { /* flip canceled */
        pfei.is_flip_ok = false;
        store_per_flip_edge_info(flip_edge_id, pfei);
        // not needed now
        // store_per_edge_flip_info(pfei.id, pefi); 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__EDGE_FLIP_EXECUTE)
    if (!valid_thread) return; 
    if (!pfei.is_flip_ok)
    { /* edge flip canceled */
        // EdgeFlags ef = load_edge_flags(pfei.id); 
        // no need for rollback any flag(s), for now...
        return; 
    }
    
    uint wedge_id = pfei.id; 
    uint v0 = ssbo_edge_to_vert_[wedge_id*4u + 0]; 
    uint v1 = ssbo_edge_to_vert_[wedge_id*4u + 1]; 
    uint v2 = ssbo_edge_to_vert_[wedge_id*4u + 2]; 
    uint v3 = ssbo_edge_to_vert_[wedge_id*4u + 3]; 

    AdjWedgeInfo w[4] = {
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 0u]), 
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 1u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 2u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 3u])
    }; 
#define w0 w[0].wedge_id
#define w1 w[1].wedge_id
#define w2 w[2].wedge_id
#define w3 w[3].wedge_id
#define w4 ((pfei.id))
   /* ------------- Edge Flip -------------
    *         v0                  v0                                                             
    *        /  \                /|\                                                            
    *       /    \              / | \                                                           
    *      W      3            W  |  3                                                          
    *     0   f1   W          0   |   W                                                         
    *    /          \        /  . |    \                                                        
    *   /  ------->  \      /  /| | |   \                                                       
    * v1 ---- W4 ---- v3  v1  f1| W4|f0  v3                                                     
    *   \  <-------  /      \   | | |/  /                                                  
    *    \          /        \    | '  /                                                        
    *     W   f0   2          W   |   2                                                         
    *      1      W            1  |  W                                                          
    *       \    /              \ | /                                                           
    *        \  /                \|/                                                            
    *         v2                  v2
    * 
    * Affected Topology:
    * e-v[w4] = { v1, v2, v3, v0 }  
    * e-e[w4] = { w1, w2, w3, w0 } */
    Store4_EVAdj(w4, uvec4(v1, v2, v3, v0)); 
    AdjWedgeInfo edges_shuffled[4] = { w[1], w[2], w[3], w[0] };
    Store4_EEAdj(w4, edges_shuffled); 
    /*
    * e-e[w0]: { w4, w3 } |-> { w1, {w4, 0} }
    * e-v[w0]: oppo_vert(v3) |-> v2
    * 
    * e-e[w1]: { w2, w4 } |-> { {w4, 0}, w0 }                       
    * e-v[w1]: oppo_vert(v3) |-> v0
    * 
    * e-e[w2]: { w4, w1 } |-> { w3, {w4, 1} }                       
    * e-v[w2]: oppo_vert(v1) |-> v0
    *
    * e-e[w3]: { w0, w4 } |-> { {w4, 1}, w2 }
    * e-v[w3]: oppo_vert(v1) |-> v2 */
    AdjWedgeInfo awi_updates[8] = 
    {
        w[1], AdjWedgeInfo(w4, 0u), // update for w0
        AdjWedgeInfo(w4, 0u), w[0], // w1
        w[3], AdjWedgeInfo(w4, 1u), // w2
        AdjWedgeInfo(w4, 1u), w[2]  // w3
    }; 
    uint voppo_updates[4] = uint[4](v2, v0, v0, v2); 
    for (uint ibwedge = 0; ibwedge < 4; ++ibwedge)
    {
        uint wi = w[ibwedge].wedge_id; 
        uint iface_overlap = w[ibwedge].iface_adj == 0u ? 1u : 0u; 
        /* Update edge links */
        ssbo_edge_to_edges_[wi*4u + mark__he_to_wedge(iface_overlap, 1u)] = encode_adj_wedge_info(awi_updates[2*ibwedge + 0]); 
        ssbo_edge_to_edges_[wi*4u + mark__he_to_wedge(iface_overlap, 2u)] = encode_adj_wedge_info(awi_updates[2*ibwedge + 1]); 
        /* Update opposite vertex id */
        ssbo_edge_to_vert_[wi*4u + mark__center_wedge_to_oppo_vert__at_face(iface_overlap)] = voppo_updates[ibwedge]; 
    }

    /* v-e: v1, v3 if v-e[vi] linked to w4 */
    uint ivert_v1_w1 = mark__cwedge_to_end_vert(w[1].iface_adj); 
    try_update_ve_link(v1, w4, VertWedgeListHeader(w1, ivert_v1_w1)); 
    uint ivert_v3_w3 = mark__cwedge_to_end_vert(w[3].iface_adj); 
    try_update_ve_link(v3, w4, VertWedgeListHeader(w3, ivert_v3_w3)); 

    /* update edge flags */
    update_edge_flags__flipped(w4); 

#undef w0
#undef w1
#undef w2
#undef w3

#endif

}
#endif


