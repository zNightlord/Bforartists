
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)


/* Inputs: 
int pcs_collapse_iter_

!!! Remember to zero out these counters !!!
* SSBO_StrokeGenEdgeCollapseCounters    ssbo_edge_collapse_counters_[]
* SSBO_StrokeGenDynamicMeshCounters     ssbo_dyn_mesh_counters_in_
* SSBO_StrokeGenDynamicMeshCounters     ssbo_dyn_mesh_counters_out_

#ifdef _KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT
.define("GLOBAL_COUNTER", "ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapsed_edges_pass_1")
.define("CP_TAG", "collapse_inital_selection")
#endif

 float ssbo_vbo_full_[]
 float pcs_remesh_edge_len_

 uint ssbo_per_collapse_edge_info_[]
 uint ssbo_per_edge_collapse_info_in_[]
 uint ssbo_per_edge_collapse_info_out_[]
 uint ssbo_per_vert_collapse_wedge_id_[] // !!!clean this buffer to NULL_EDGEs
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

float get_collapse_edge_len_max()
{
    return calc_remesh_edge_len_min(pcs_remesh_edge_len_); 
}

struct EdgeCollapseInfo
{
    bool is_collapse_ok; 
    float edge_len; 
    uint id; 
}; 
uvec2 encode_edge_collapse_info(EdgeCollapseInfo eci)
{
    uvec2 eci_enc = uvec2(0, 0); 
    eci_enc.x = floatBitsToUint(eci.edge_len); 
    eci_enc.y = ((eci.id << 1u) | uint(eci.is_collapse_ok));

    return eci_enc;  
}
EdgeCollapseInfo decode_edge_collapse_info(uvec2 eci_enc)
{
    EdgeCollapseInfo eci; 
    eci.edge_len = uintBitsToFloat(eci_enc.x);
    eci.is_collapse_ok = (0u != (eci_enc.y & 1u));
    eci_enc.y >>= 1u; 
    eci.id = eci_enc.y;  

    return eci; 
}
void store_per_collapse_edge_info(uint collapse_edge_id, EdgeCollapseInfo pcei)
{
    uvec2 pcei_enc = encode_edge_collapse_info(pcei); 
    Store2(ssbo_per_collapse_edge_info_, collapse_edge_id, pcei_enc);
}
EdgeCollapseInfo load_per_collapse_edge_info(uint collapse_edge_id)
{
    uvec2 pcei_enc; 
    Load2(ssbo_per_collapse_edge_info_, collapse_edge_id, pcei_enc); 
    return decode_edge_collapse_info(pcei_enc); 
}

void store_per_edge_collapse_info(uint wedge_id, EdgeCollapseInfo peci)
{
    uvec2 peci_enc = encode_edge_collapse_info(peci); 
    Store2(ssbo_per_edge_collapse_info_out_, wedge_id, peci_enc); 
}
void store_per_edge_collapse_info_two_buffers(uint wedge_id, EdgeCollapseInfo peci)
{
    uvec2 peci_enc = encode_edge_collapse_info(peci); 
    Store2(ssbo_per_edge_collapse_info_in_, wedge_id, peci_enc); 
    Store2(ssbo_per_edge_collapse_info_out_, wedge_id, peci_enc); 
}
void store_per_edge_collapse_info__id_and_flag(uint wedge_id, EdgeCollapseInfo peci)
{
    uvec2 peci_enc = encode_edge_collapse_info(peci); 
    ssbo_per_edge_collapse_info_out_[wedge_id*2u + 1u] = peci_enc.y;
}
EdgeCollapseInfo load_per_edge_collapse_info(uint wedge_id)
{
    uvec2 peci_enc;
    Load2(ssbo_per_edge_collapse_info_in_, wedge_id, peci_enc); 
    return decode_edge_collapse_info(peci_enc); 
}
bool load_per_edge_collapse_info__is_collapse_ok(uint wedge_id)
{
    uvec2 peci_enc;
    peci_enc.x = 0; // does not matter
    peci_enc.y = ssbo_per_edge_collapse_info_in_[2*wedge_id + 1u]; 
    EdgeCollapseInfo eci = decode_edge_collapse_info(peci_enc); 
    return eci.is_collapse_ok; 
}





void store_per_vert_collapse_wedge_id(uint vid, uint wedge_id)
{
    ssbo_per_vert_collapse_wedge_id_[vid] = wedge_id; 
}
uint load_per_vert_collapse_wedge_id(uint vid)
{
    return ssbo_per_vert_collapse_wedge_id_[vid]; 
}





/* Prefer collapsing on shorter edges */
bool collapse_score_larger(EdgeCollapseInfo eci_0, uint wedge_id_0, EdgeCollapseInfo eci_1, uint wedge_id_1)
{
    if (!eci_0.is_collapse_ok) return false;
    if (!eci_1.is_collapse_ok) return true;
    if (eci_0.edge_len == eci_1.edge_len)
        return wedge_id_0 < wedge_id_1; 
    return eci_0.edge_len < eci_1.edge_len; 
}



#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_INIT)
void main()
{ 
    if (pcs_collapse_iter_ == 0 && gl_GlobalInvocationID.x == 0u)
        ssbo_edge_collapse_counters_[0].num_collapsed_edges_pass_1 = 0u;

    /* initialze vertex marks */
    uint vid = gl_GlobalInvocationID.x; 
    uint num_verts  = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts;
    bool valid_thread = vid < num_verts;  
    
    if (valid_thread)
        store_per_vert_collapse_wedge_id(vid, NULL_EDGE); 
}
#endif







#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT)
void main()
{
    uint wedge_id = gl_GlobalInvocationID.x; 
    const uint groupId = gl_LocalInvocationID.x; 

    bool valid_thread = wedge_id < (pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_edges); 
    if (wedge_id == 0u)
    {
        ssbo_dyn_mesh_counters_in_.num_edges = ssbo_dyn_mesh_counters_out_.num_edges; 
        ssbo_dyn_mesh_counters_in_.num_verts = ssbo_dyn_mesh_counters_out_.num_verts; 
    }

    /* Calc edge length */
    uvec2 iverts_cwedge = mark__wedge_to_verts(4/*mark of cwedge*/); 
    uvec2 verts_cwedge = uvec2(
        ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge.x],
        ssbo_edge_to_vert_[wedge_id*4u + iverts_cwedge.y]
    ); 
    vec3 vpos_cwedge[2]; 
    vpos_cwedge[0] = ld_vpos(verts_cwedge.x);
    vpos_cwedge[1] = ld_vpos(verts_cwedge.y);
    float edge_len = length(vpos_cwedge[0] - vpos_cwedge[1]);

    EdgeFlags ef = load_edge_flags(wedge_id);

    /* Filter collapse edges */
    bool is_collapse_ok = (edge_len < get_collapse_edge_len_max()) 
        && (ef.selected)
        && (!ef.dupli)
        && (!ef.border)
        && (!ef.del_by_collapse)
        && (!ef.del_by_split)
        && valid_thread; 

    uint collapse_edge_id = compact_collapse_select_short_edges(is_collapse_ok, groupId);
    if (is_collapse_ok)
    {
        EdgeCollapseInfo pcei; 
        pcei.is_collapse_ok = true; 
        pcei.edge_len = edge_len;
        pcei.id = wedge_id;
        store_per_collapse_edge_info(collapse_edge_id, pcei);
    }

    if (valid_thread)
    {
        EdgeCollapseInfo peci; 
        peci.is_collapse_ok = is_collapse_ok;
        peci.edge_len = edge_len;
        peci.id = is_collapse_ok ? collapse_edge_id : NULL_EDGE;
        store_per_edge_collapse_info_two_buffers(wedge_id, peci); 
    }
}
#endif





#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_0)
    struct CollapseValidateContext
    {
        bool is_collapse_ok; 
        
        vec3 collapse_pos; 
        /*    
        *     vp       
        *    /  \      
        *   wp   \     
        *  /      \    
        * v --wi--vi   
        *  \<-----/    
        *  wn fi wo    
        *    \  /      
        *     vn        
        */
        uint v; // center vtx of the circulation
        vec3 vpos; 
        uint vp; // vtx of the previous visited wedge
        vec3 vpos_p; 
        uint vi; // vtx of the currently visited wedge
        vec3 vpos_i; 

        uint v0; 
        uint v1; 
        uint v2; 
        uint v3; 
        bool check_topo; // wether to check topological consistency
    }; 
    CollapseValidateContext init_collapse_validate_context(
        bool rotate_v1, vec3 vpos[4], vec3 collapse_pos, 
        uvec4 v, uvec4 w
    )
    {
        CollapseValidateContext ctx; 
        ctx.is_collapse_ok = true; 
        
        ctx.collapse_pos = collapse_pos; 
            
        uint ivert = rotate_v1 ? 1 : 3; 
        ctx.v = v[ivert]; 
        ctx.vpos = vpos[ivert]; 
        
        uint ivert_p = rotate_v1 ? 0 : 2; 
        ctx.vp = v[ivert_p]; 
        ctx.vpos_p = vpos[ivert_p]; 
        
        uint ivert_i = rotate_v1 ? 3 : 1; 
        ctx.vi = v[ivert_i]; 
        ctx.vpos_i = vpos[ivert_i]; 
        
        ctx.v0 = v[0]; 
        ctx.v1 = v[1]; 
        ctx.v2 = v[2]; 
        ctx.v3 = v[3]; 
        ctx.check_topo = rotate_v1; 
        // need to check for both sides

        return ctx; 
    }

    /*
    *        v0  ...  vp                         v_oppo --- vn  ...  v0             
    *       /  \     /  \                            \     /  \     /  \            
    *      /    \   wp   \                            \  wo fi wn  /    \           
    *     / ---> \ /      \                            \ / ---> \ / ---> \          
    *   v1------- v--wi-- vi                           vi --wi-- v ------ v3        
    *     \ <--- / \<-----/ \        fi:awi.iface_adj    \      / \ <--- /          
    *      \    /  wn fi wo  \       wi:awi.wedge_id      \   wp   \    /           
    *       \  /     \  /     \                            \  /     \  /            
    *        v2  ...  vn ----- v_oppo                       vp  ...  v2             
    */    
    bool collapse_validate(
        CirculatorIterData iter, 
        inout CollapseValidateContext ctx
    )      
    {
        uint wi = iter.awi.wedge_id; 

        uint iwedge_wo = mark__cwedge_rotate_back(iter.awi.iface_adj); // outside edge
        AdjWedgeInfo awi_wo = decode_adj_wedge_info(ssbo_edge_to_edges_[iter.awi.wedge_id*4u + iwedge_wo]); 
        uint ivert_voppo = mark__center_wedge_to_oppo_vert__at_face(awi_wo.iface_adj);
        uint v_oppo = ssbo_edge_to_vert_[awi_wo.wedge_id*4u + ivert_voppo]; 

        uint ivert_vn = mark__center_wedge_to_oppo_vert__at_face(iter.awi.iface_adj);  
        uint vn = ssbo_edge_to_vert_[iter.awi.wedge_id*4u + ivert_vn]; 
        vec3 vpos_n = ld_vpos(vn); 


        /* Check topology consistency */
        /* must meet the link condition, 
        * see "SurfaceMesh::is_collapse_ok" in the pmp library 
        * or "Generic Memoryless Polygonal Simplification" 
        * For further theory see "Topology preserving edge contraction" */
        uvec3 valid_link_verts = (ctx.v == ctx.v1) ? uvec3(ctx.v2, ctx.v3, ctx.v0) : uvec3(ctx.v0, ctx.v1, ctx.v2);
        uint v_oppo_invalid = (ctx.v == ctx.v1) ? ctx.v3 : ctx.v1;
        
        /* The outer edge is linked to another vert on the collapsed wedge, 
        * <=> 1-ring of v1 & v3 has a common edge (excluding the collapsed faces v013 and v123) */
        if (ctx.check_topo && all(ctx.vi.xxx != valid_link_verts) && v_oppo == v_oppo_invalid)
            ctx.is_collapse_ok = false; 

        
        /* Check geometry consistency */
        /* Foldover avoidance, see Ch.3.7 in https://www.cs.cmu.edu/~garland/thesis/thesis-onscreen.pdf */
        vec3 face_normal = normalize(cross(vpos_n-ctx.vpos, ctx.vpos_i-ctx.vpos)); // cross dot in CCW order
        vec3 clip_plane_normal = normalize(cross(face_normal, ctx.vpos_i - vpos_n)); // cross dot in CCW order

        if (dot(ctx.collapse_pos - vpos_n, clip_plane_normal) < .0f)
            ctx.is_collapse_ok = false; /* falls outside the clipping plane */

        /* Update context */
        ctx.vp = ctx.vi; 
        ctx.vpos_p = ctx.vpos_i; 
        ctx.vi = vn; 
        ctx.vpos_i = vpos_n; 

        return ctx.is_collapse_ok; // cont'/break the loop
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_1)
    struct CollapseCollisionContext
    {
        EdgeCollapseInfo eci;
        bool is_collapse_ok; 
        uint wedge_id; 
    }; 
    bool collapse_resolve_collision_I(
        CirculatorIterData iter, 
        inout CollapseCollisionContext ctx
    ){
        EdgeCollapseInfo peci = load_per_edge_collapse_info(iter.awi.wedge_id); 
        if (peci.is_collapse_ok && ctx.is_collapse_ok && ctx.wedge_id != iter.awi.wedge_id)
        { /* found a better collapse alternative in 1-ring */
            if (collapse_score_larger(peci, iter.awi.wedge_id, ctx.eci, ctx.wedge_id))
                ctx.is_collapse_ok = false; 
        }

        return ctx.is_collapse_ok; 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_2)
    void SolveConflictionAtVertex(uint vi, uint wedge_id, inout EdgeCollapseInfo eci)
    {
        bool lose_collapse = false; 
        uint wk = load_per_vert_collapse_wedge_id(vi); 
        bool wk_collapsed = (wk != NULL_EDGE); 
        EdgeCollapseInfo peci_wk = EdgeCollapseInfo(false, 10000.0f, NULL_EDGE); 
        if (wk_collapsed)
        {
            peci_wk = load_per_edge_collapse_info(wk); 
            lose_collapse = collapse_score_larger(
                peci_wk, wk, 
                eci, wedge_id
            ); 
            /* cancel collapsing for the loser */
            if (lose_collapse) 
                eci.is_collapse_ok = false; 
            else 
            { 
                peci_wk.is_collapse_ok = false;
                store_per_edge_collapse_info__id_and_flag(wk, peci_wk);
            }
        }
    }

    struct CollapseCollisionContext
    {
        EdgeCollapseInfo eci;
        bool is_collapse_ok; 
        uint wedge_id; 
    }; 
    bool collapse_resolve_collision_II(
        CirculatorIterData iter, 
        inout CollapseCollisionContext ctx
    ){
        uint wi = iter.awi.wedge_id; 
        if (wi == ctx.wedge_id)
            return true;

        uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
        uint vi = ssbo_edge_to_vert_[iter.awi.wedge_id*4u + ivert_vi];
        SolveConflictionAtVertex(vi, ctx.wedge_id, /*inout*/ctx.eci);
        return true; 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE)
    struct CollapseAdjustEVLinkContext
    {
        uint v1; // we are removing v3
        uint cwedge; 
        uint w2; 
        uint w3; 
        
        AdjWedgeInfo e30; 
        AdjWedgeInfo e31; 
        AdjWedgeInfo e20; 
        AdjWedgeInfo e21; 

        uint q2; 
        uint q3; 
    }; 

    bool collapse_adjust_ev_links_to_v1(
        CirculatorIterData iter, 
        inout CollapseAdjustEVLinkContext ctx
    ){ /* fwd rotate around v3 */
        uint ivert_wo = mark__cwedge_rotate_back(iter.awi.iface_adj); 
        AdjWedgeInfo wo = decode_adj_wedge_info(ssbo_edge_to_edges_[iter.awi.wedge_id*4u + ivert_wo]); 

        if (all(bvec2(iter.awi.wedge_id.xx != uvec2(ctx.cwedge, ctx.w2))))
        { /* adjust v-e link for wi and wo */
            uint ivert_v3; 
            /* wi */
            if (iter.awi.wedge_id != ctx.w3)
            { /* dont do this when wi==w3*/
                ivert_v3 = mark__cwedge_to_end_vert(iter.awi.iface_adj); 
                ssbo_edge_to_vert_[iter.awi.wedge_id*4u + ivert_v3] = ctx.v1; /* no race cond' since these edges already locked */
            } 
            /* wo */
            ivert_v3 = mark__center_wedge_to_oppo_vert__at_face(wo.iface_adj == 1 ? 0 : 1); 
            ssbo_edge_to_vert_[wo.wedge_id*4u + ivert_v3] = ctx.v1; 
        }
        
        /* Cache pointers */
        if (iter.awi.wedge_id == ctx.w3)
        {
            ctx.e31 = wo; 
            ctx.e30 = iter.awi_next; // wn

            uint ivert_q3 = mark__cwedge_to_beg_vert(ctx.e30.iface_adj);
            ctx.q3 = ssbo_edge_to_vert_[ctx.e30.wedge_id*4u + ivert_q3];
        }

        if (iter.awi_next.wedge_id == ctx.w2)
        {
            ctx.e20 = wo; 
            ctx.e21 = iter.awi; // wi

            uint ivert_q2 = mark__cwedge_to_beg_vert(ctx.e21.iface_adj);
            ctx.q2 = ssbo_edge_to_vert_[ctx.e21.wedge_id*4u + ivert_q2]; 

            // flip the iface to change e21 from e21 to w2's ee adjacency
            ctx.e21.iface_adj = ctx.e21.iface_adj == 1u ? 0u : 1u; 
        }


        return true; 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 

    uint collapse_edge_id = gl_GlobalInvocationID.x; 
    uint num_selected_edges = ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapsed_edges_pass_1;
    bool valid_thread = collapse_edge_id < num_selected_edges; 

    EdgeCollapseInfo pcei = load_per_collapse_edge_info(collapse_edge_id);

    AdjWedgeInfo w[4]; 
    w[0] = decode_adj_wedge_info(ssbo_edge_to_edges_[pcei.id*4u + 0u]);
    w[1] = decode_adj_wedge_info(ssbo_edge_to_edges_[pcei.id*4u + 1u]);
    w[2] = decode_adj_wedge_info(ssbo_edge_to_edges_[pcei.id*4u + 2u]);
    w[3] = decode_adj_wedge_info(ssbo_edge_to_edges_[pcei.id*4u + 3u]);
#define w0 ((w[0].wedge_id))
#define w1 ((w[1].wedge_id))
#define w2 ((w[2].wedge_id))
#define w3 ((w[3].wedge_id))
#define w4 ((pcei.id))

    uint v0, v1, v2, v3; 
    v0 = ssbo_edge_to_vert_[pcei.id*4u + 0u];
    v1 = ssbo_edge_to_vert_[pcei.id*4u + 1u];
    v2 = ssbo_edge_to_vert_[pcei.id*4u + 2u];
    v3 = ssbo_edge_to_vert_[pcei.id*4u + 3u];


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_0)
    /* Validate collapse with geometry/topology consistency */
    EdgeFlags ef_w0 = load_edge_flags(w0); 
    EdgeFlags ef_w1 = load_edge_flags(w1); 
    EdgeFlags ef_w2 = load_edge_flags(w2); 
    EdgeFlags ef_w3 = load_edge_flags(w3); 
    if (any(bvec4(ef_w0.border, ef_w1.border, ef_w2.border, ef_w3.border)))
        pcei.is_collapse_ok = false; // has a border edge in 1-ring neighborhood
        // if need to collapse on a border edge, the affected v-e adj. needs to updated to border edge

    /* Validate collapse with geometry/topology consistency 
     * by visiting 1-ring neighborhoods around 2 edge verts */
    if (valid_thread && pcei.is_collapse_ok)
    {
        vec3 vpos[4] = { ld_vpos(v0), ld_vpos(v1), ld_vpos(v2), ld_vpos(v3) };
        vec3 collapse_pos = .5f * (vpos[1] + vpos[3]); 
        bool rotate_fwd = true; 
        
        { /* check around the vert */
            VertWedgeListHeader vwlh_v1 = VertWedgeListHeader(w4, 1); /* rotate around v1 */
            CollapseValidateContext ctx_v1 = init_collapse_validate_context(
                true, vpos, collapse_pos, 
                uvec4(v0, v1, v2, v3), uvec4(w0, w1, w2, w3) 
            ); 

            VE_CIRCULATOR(vwlh_v1, collapse_validate, ctx_v1, rotate_fwd); 
            pcei.is_collapse_ok = ctx_v1.is_collapse_ok; 
        }        

        if (pcei.is_collapse_ok)
        { 
            VertWedgeListHeader vwlh_v3 = VertWedgeListHeader(w4, 3); /* rotate around v3 */
            CollapseValidateContext ctx_v3 = init_collapse_validate_context(
                false, vpos, collapse_pos, 
                uvec4(v0, v1, v2, v3), uvec4(w0, w1, w2, w3) 
            ); 
            
            VE_CIRCULATOR(vwlh_v3, collapse_validate, ctx_v3, rotate_fwd)   
            pcei.is_collapse_ok = ctx_v3.is_collapse_ok; 
        }
    }

    if (valid_thread)
    {
        store_per_collapse_edge_info(collapse_edge_id, pcei); 

        /* Ensure two ping-pong buffers are coherent */
        EdgeCollapseInfo peci;
        peci.edge_len = pcei.edge_len; 
        peci.id = collapse_edge_id;  
        peci.is_collapse_ok = load_per_edge_collapse_info__is_collapse_ok(pcei.id); 
        if (false == pcei.is_collapse_ok)
        { /* update per edge data */
            peci.is_collapse_ok = false; 
        }
        store_per_edge_collapse_info__id_and_flag(pcei.id, peci); 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_1)
    /* Resolve collision pass I 
     * for each edge, traverse 1-ring of its 2 verts */
    EdgeCollapseInfo pcei_new = pcei; 
    bool rotate_fwd = true; 
    if (pcei_new.is_collapse_ok)
    {
        {
            VertWedgeListHeader vwlh_v1 = VertWedgeListHeader(w4, 1); 
            CollapseCollisionContext ctx; 
            ctx.eci = pcei; 
            ctx.is_collapse_ok = true;
            ctx.wedge_id = pcei.id; 
            
            VE_CIRCULATOR(vwlh_v1, collapse_resolve_collision_I, ctx, rotate_fwd); 
            pcei_new.is_collapse_ok = ctx.is_collapse_ok; 
        }
        
        if (pcei_new.is_collapse_ok)
        {
            VertWedgeListHeader vwlh_v3 = VertWedgeListHeader(w4, 3); 
            CollapseCollisionContext ctx; 
            ctx.eci = pcei; 
            ctx.is_collapse_ok = true; 
            ctx.wedge_id = pcei.id; 
            
            VE_CIRCULATOR(vwlh_v3, collapse_resolve_collision_I, ctx, rotate_fwd); 
            pcei_new.is_collapse_ok = ctx.is_collapse_ok; 
        }
    }


    EdgeCollapseInfo peci; 
    peci.edge_len = pcei.edge_len; 
    peci.id = collapse_edge_id; 
    peci.is_collapse_ok = load_per_edge_collapse_info__is_collapse_ok(pcei.id); 
    if (valid_thread && pcei.is_collapse_ok && !pcei_new.is_collapse_ok)
    { /* Update info when the collapse gets canceled */
        store_per_collapse_edge_info(collapse_edge_id, pcei_new); 
        peci.is_collapse_ok = false; 
    }
    if (valid_thread)
        store_per_edge_collapse_info__id_and_flag(pcei.id, peci); 


    if (valid_thread && pcei_new.is_collapse_ok)
    { /* mark edge score for 2 edge verts */
        store_per_vert_collapse_wedge_id(v1, pcei.id); 
        store_per_vert_collapse_wedge_id(v3, pcei.id); 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_2)
    EdgeCollapseInfo peci;
    peci.id = collapse_edge_id; 
    peci.edge_len = pcei.edge_len; 
    peci.is_collapse_ok = load_per_edge_collapse_info__is_collapse_ok(pcei.id); 

    EdgeCollapseInfo pcei_new = pcei; 
    if (valid_thread/*important since we store data here*/ 
        && pcei.is_collapse_ok)
    {
        if (pcei_new.is_collapse_ok){
            if (pcei_new.is_collapse_ok)
            {
                bool rotate_fwd = true; 

                VertWedgeListHeader vwlh_v1 = VertWedgeListHeader(pcei.id, 1);
                CollapseCollisionContext ctx;
                ctx.eci            = pcei;
                ctx.is_collapse_ok = true; 
                ctx.wedge_id       = pcei.id; 
                
                VE_CIRCULATOR(vwlh_v1, collapse_resolve_collision_II, ctx, rotate_fwd);

                pcei_new.is_collapse_ok = ctx.is_collapse_ok;
                if (pcei_new.is_collapse_ok)
                {
                    VertWedgeListHeader vwlh_v3 = VertWedgeListHeader(pcei.id, 3);
                    CollapseCollisionContext ctx;
                    ctx.eci            = pcei;
                    ctx.is_collapse_ok = true; 
                    ctx.wedge_id       = pcei.id; 
                    
                    VE_CIRCULATOR(vwlh_v3, collapse_resolve_collision_II, ctx, rotate_fwd);

                    pcei_new.is_collapse_ok = ctx.is_collapse_ok; 
                }
            }
        }


        if (valid_thread && pcei.is_collapse_ok && !pcei_new.is_collapse_ok)
        {
            store_per_collapse_edge_info(collapse_edge_id, pcei_new); 
            peci.is_collapse_ok = false; 
        }
    }

    if (valid_thread && peci.is_collapse_ok == false)
        store_per_edge_collapse_info__id_and_flag(pcei.id, peci); 
#endif


/* Apply the edge collapse */
#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE)
    if (!valid_thread)
        return;
    
    if (gl_GlobalInvocationID.x == 0u) 
        /* prep local counters for next split iteration, have to do this before early exit */
        ssbo_edge_collapse_counters_[pcs_collapse_iter_ + 1].num_collapsed_edges_pass_1 = 0u;

    if (!pcei.is_collapse_ok)
        return;

    EdgeCollapseInfo peci = load_per_edge_collapse_info(pcei.id);
    if (!peci.is_collapse_ok) /* canceled by conflit resolve pass #2*/
        return; 

    /* Parallel Edge Collapse * 
    * (Following graph are only to clearly show how topology changes, actual position of vertices are not accurate)
    *                                  
    *             v0 -----E31---- q3         v0 -----E31---- q3                .v0/q2    Note: in this case                             
    *            /  \  <-------  /             \  <-------  /                _/ /|       E30==W2, E21==W3.                                
    *           /    \          0               \          0               _/  / |                               
    *          W      3        3                 W  W     3              W/   3  |                               
    *         0   f1   W      E                   0<=3   E             _0    W   |                               
    *        /          \    /                     \    /            _/ f1  E21  |                               
    *       /  ------->  \  /                       \  /           _/ ---> /     E31                             
    *     v1 ---- W4 ---- v3      ==collapse=>       v1 <= v3     v1- W4 -v3     |                             
    *       \  <-------  /  \                       /  \           \_ <--- \     E20                              
    *        \          /    \                     /    \            \_ f0  E30  |                                
    *         W   f0   2      1                   W  W   1             \_    2   |                                
    *          1      W        2                 1<=2     2              W_   W  |                                
    *           \    /          E               /          E              1\_  \ |                               
    *            \  /  ------->  \             /  ------->  \                \_ \|                               
    *             v2 ----E20----- q2         v2 ----E20----- q2                `v2/q3                            
    *    Changes: 
    *    1. del v3, replace by v1; 
    *    2. del W2, replace by w1;      
    *    3. del W3, replace by w0; 
    *    4. del E4; 
    *    Affected elems: 
    *    by v3: e-v link for all edges in v3's 1-ring closure, (special)e-v link for W1-q2, e-v link for W0-q3
    *    by W2: v-e link for v2, e-e link for W1,E20,E21, 
    *    by W3: v-e link for v0, e-e link for W0,E30,E31, 
    *    by W4: v-e link for v1, e-e link for W0,W1         
    * 
    *    Implementation: 
    *                                  
    *    mark v3, w2, w3, w4 as deleted    */
    update_vert_flags__del_by_collapse(v3); 
    update_edge_flags__del_by_collapse(w2); 
    update_edge_flags__del_by_collapse(w3); 
    update_edge_flags__del_by_collapse(w4); 
    /*                                  
    *    for all edges in v3's 1-ring closure:
    *      change e-v link from v3 to v1
    */
    bool rot_fwd = true; 
    VertWedgeListHeader vwlh_v3 = VertWedgeListHeader(w4, 3); 
    AdjWedgeInfo awi_null = AdjWedgeInfo(NULL_EDGE, 0); 
    CollapseAdjustEVLinkContext ctx = CollapseAdjustEVLinkContext(
        v1, pcei.id, w2, w3, 
        awi_null, awi_null, awi_null, awi_null, 0u, 0u
    ); 
    VE_CIRCULATOR(vwlh_v3, collapse_adjust_ev_links_to_v1, ctx, rot_fwd)
    /* Specially, link q2 to w1, q3 to w0 */
    ssbo_edge_to_vert_[w1*4u + mark__center_wedge_to_oppo_vert__at_face(w[1].iface_adj == 0u ? 1u : 0u)] = ctx.q2;
    ssbo_edge_to_vert_[w0*4u + mark__center_wedge_to_oppo_vert__at_face(w[0].iface_adj == 0u ? 1u : 0u)] = ctx.q3;  
    /*                                  
    *    if v-e[v2]==w2: 
    *      v-e[v2]=w1
    *    if v-e[v0]==w3:
    *      v-e[v0]=w0
    *    if v-e[v1]==w4:
    *      v-e[v0]=w0 // todo: border not handled 
    */
    try_update_ve_link(v2, w2, VertWedgeListHeader(w1, mark__cwedge_to_beg_vert(w[1].iface_adj))); 
    try_update_ve_link(v0, w3, VertWedgeListHeader(w0, mark__cwedge_to_end_vert(w[0].iface_adj))); 
    try_update_ve_link(v1, w4, VertWedgeListHeader(w0, mark__cwedge_to_beg_vert(w[0].iface_adj))); 
    /*    
    *    link  w1.next to e20: e-e[ w1*4 + mark__wedge_to_next_wedge(w1.iface_adj == 1 ? 0 : 1)] = e20
    *    link e20.prev to  w1: e-e[e20*4 + mark__wedge_to_prev_wedge(e20.iface_adj == 1 ? 0 : 1)] = w1
    *    link  w1.prev to e21: e-e[ w1*4 + mark__wedge_to_prev_wedge(w1.iface_adj == 1 ? 0 : 1)] = e21
    *    link e21.next to  w1: e-e[e21*4 + mark__wedge_to_next_wedge(e21.iface_adj == 1 ? 0 : 1)] = w1
    *    
    *    link  w0.next to e30: e-e[ w0*4 + mark__wedge_to_next_wedge(w0.iface_adj == 1 ? 0 : 1)] = e30
    *    link e30.prev to  w0: e-e[e30*4 + mark__wedge_to_prev_wedge(e30.iface_adj == 1 ? 0 : 1)] = w0
    *    link  w0.prev to e31: e-e[ w0*4 + mark__wedge_to_prev_wedge(w0.iface_adj == 1 ? 0 : 1)] = e31
    *    link e31.next to  w0: e-e[e31*4 + mark__wedge_to_next_wedge(e31.iface_adj == 1 ? 0 : 1)] = w0
    */
    
    if (ctx.e21.wedge_id == w3) ctx.e21 = w[0]; /* Special case, see main graph above */
    if (ctx.e30.wedge_id == w2) ctx.e30 = w[1]; 
    uint e20 = ctx.e20.wedge_id;
    uint e21 = ctx.e21.wedge_id; 
    uint e30 = ctx.e30.wedge_id; 
    uint e31 = ctx.e31.wedge_id; 
    ssbo_edge_to_edges_[ w1*4 + mark__wedge_to_next_wedge(w[1].iface_adj    == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(ctx.e20); 
    ssbo_edge_to_edges_[e20*4 + mark__wedge_to_prev_wedge(ctx.e20.iface_adj == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(w[1]); 
    ssbo_edge_to_edges_[ w1*4 + mark__wedge_to_prev_wedge(w[1].iface_adj    == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(ctx.e21); 
    ssbo_edge_to_edges_[e21*4 + mark__wedge_to_next_wedge(ctx.e21.iface_adj == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(w[1]); 
    ssbo_edge_to_edges_[ w0*4 + mark__wedge_to_next_wedge(w[0].iface_adj    == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(ctx.e30); 
    ssbo_edge_to_edges_[e30*4 + mark__wedge_to_prev_wedge(ctx.e30.iface_adj == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(w[0]); 
    ssbo_edge_to_edges_[ w0*4 + mark__wedge_to_prev_wedge(w[0].iface_adj    == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(ctx.e31); 
    ssbo_edge_to_edges_[e31*4 + mark__wedge_to_next_wedge(ctx.e31.iface_adj == 1 ? 0 : 1, 4u)] = encode_adj_wedge_info(w[0]); 

    /* Update v1 pos */
    vec3 vpos = .5f * (ld_vpos(v1) + ld_vpos(v3)); 
    st_vpos(v1, vpos); 

#endif

#undef w0
#undef w1
#undef w2
#undef w3

}
#endif