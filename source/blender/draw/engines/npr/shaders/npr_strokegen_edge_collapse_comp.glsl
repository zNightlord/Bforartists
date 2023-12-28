
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
.define("GLOBAL_COUNTER", "ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapse_edges_pass_1")
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
EdgeCollapseInfo load_per_edge_collapse_info(uint wedge_id)
{
    uvec2 peci_enc;
    Load2(ssbo_per_edge_collapse_info_in_, wedge_id, peci_enc); 
    return decode_edge_collapse_info(peci_enc); 
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
    if (eci_0.edge_len == eci_0.edge_len)
        return wedge_id_0 < wedge_id_1; 
    return eci_0.edge_len < eci_0.edge_len; 
}



#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_INIT)
void main()
{ 
    uint vid = gl_GlobalInvocationID.x; 
    if (pcs_collapse_iter_ == 0 && gl_GlobalInvocationID.x == 0u)
    { 
        ssbo_edge_collapse_counters_[0].num_collapse_edges_pass_1 = 0u;
    }
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
        && (!ef.dupli)
        && (!ef.border)
        && (!ef.del_by_collapse)
        && (!ef.del_by_split)
        && valid_thread;; 

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
        store_per_edge_collapse_info(wedge_id, peci); 
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
        bool rotate_v1, vec3 vpos[4], 
        vec3 collapse_pos, 
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
        ctx.check_topo = rotate_v1; // costly check, only do this once

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
        uint vi = iter.awi.wedge_id; 

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
        if (ctx.check_topo && all(vi.xxx != valid_link_verts) && v_oppo == v_oppo_invalid)
            ctx.is_collapse_ok = false; 

        
        /* Check geometry consistency */
        /* Foldover avoidance, see Ch.3.7 in https://www.cs.cmu.edu/~garland/thesis/thesis-onscreen.pdf */
        vec3 face_normal = normalize(cross(vpos_n-ctx.vpos, ctx.vpos_i-ctx.vpos)); // cross dot in CCW order
        vec3 clip_plane_normal = normalize(cross(face_normal, ctx.vpos_i - vpos_n)); // cross dot in CCW order

        if (dot(ctx.collapse_pos - vpos_n, clip_plane_normal) < .0f)
            ctx.is_collapse_ok = false; /* falls outside the clipping plane */

        /* Update context */
        ctx.vp = iter.awi.wedge_id; 
        ctx.vpos_p = ctx.vpos_i; 
        ctx.vi = iter.awi_next.wedge_id; 
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
        uint pcs_collapse_iter; 
    }; 
    bool collapse_resolve_collision_I(
        CirculatorIterData iter, 
        inout CollapseCollisionContext ctx
    ){
        EdgeCollapseInfo peci = load_per_edge_collapse_info(iter.awi.wedge_id, ctx.pcs_collapse_iter)
        if (peci.is_collapse_ok && ctx.is_collapse_ok)
        { /* found a better collapse alternative in 1-ring */
            if (collapse_score_larger(peci, iter.awi.wedge_id, eci, ctx.wedge_id))
                ctx.is_collapse_ok = false; 
        }

        return ctx.is_collapse_ok; 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 

    uint collapse_edge_id = gl_GlobalInvocationID.x; 
    uint num_selected_edges = ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapse_edges_pass_1;
    bool valid_thread = collapse_edge_id < num_selected_edges; 

    EdgeCollapseInfo pcei = load_per_collapse_edge_info(collapse_edge_id);

    AdjWedgeInfo w[4]; 
    w[0] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 0u]);
    w[1] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 1u]);
    w[2] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 2u]);
    w[3] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 3u]);
#define w0 ((w[0].wedge_id))
#define w1 ((w[1].wedge_id))
#define w2 ((w[2].wedge_id))
#define w3 ((w[3].wedge_id))

    uint v0, v1, v2, v3; 
    v0 = ssbo_edge_to_vert_[psei_curr.id*4u + 0u];
    v1 = ssbo_edge_to_vert_[psei_curr.id*4u + 1u];
    v2 = ssbo_edge_to_vert_[psei_curr.id*4u + 2u];
    v3 = ssbo_edge_to_vert_[psei_curr.id*4u + 3u];


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
    vec3 collapse_pos = vec3(.0f); 
    if (validThread && pcei.is_collapse_ok)
    {
        collapse_pos = .5f * (ld_vpos(v1) + ld_vpos(v3)); 
        CollapseValidateContext ctx_init = init_collapse_validate_context(
            true/*rotate_v1*/, vpos, collapse_pos, 
            uvec4(v0, v1, v2, v3), uvec4(w0, w1, w2, w3) 
        ); 

        CollapseValidateContext ctx_v1 = ctx_init; 
        ctx_v1.rotate_v1 = true;
        VE_CIRCULATOR_FWD(v1, collapse_validate, ctx_v1); 
        pcei.is_collapse_ok = ctx_v1.is_collapse_ok; 

        if (pcei.is_collapse_ok)
        {
            CollapseValidateContext ctx_v3 = ctx_init; 
            ctx_v3.rotate_v1 = false; 
            VE_CIRCULATOR_FWD(v3, collapse_validate, ctx_v3)   
            pcei.is_collapse_ok = ctx_v3.is_collapse_ok; 
        }
    }

    if (valid_thread)
    {
        store_per_collapse_edge_info(collapse_edge_id, pcei); 
        if (false == pcei.is_collapse_ok)
        { /* update per edge data */
            EdgeCollapseInfo peci; 
            peci.is_collapse_ok = false; 
            peci.edge_len = .0f; 
            peci.id = collapse_edge_id; /* keep link */
            store_per_edge_collapse_info(pcei.id, peci); 
        }

        /* initialze vertex marks */
        if (valid_thread)
        {
            store_per_vert_collapse_wedge_id(v0, NULL_EDGE); 
            store_per_vert_collapse_wedge_id(v1, NULL_EDGE); 
            store_per_vert_collapse_wedge_id(v2, NULL_EDGE); 
            store_per_vert_collapse_wedge_id(v3, NULL_EDGE);
        }
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_1)
    /* Resolve collision pass I 
    * for each edge, traverse 1-ring of its 2 verts */
    EdgeCollapseInfo pcei_new = pcei; 
    if (pcei_new.is_collapse_ok)
    {
        CollapseCollisionContext ctx; 
        ctx.is_collapse_ok = true; 
        ctx.edge_len = pcei.edge_len; 
        ctx.wedge_id = pcei.id; 
        
        VE_CIRCULATOR_FWD(v1, ctx, collapse_resolve_collision_I); 
        
        pcei_new.is_collapse_ok = ctx.is_collapse_ok; 
        if (pcei_new.is_collapse_ok)
        {
            ctx.is_collapse_ok = true; 
            ctx.edge_len = pcei.edge_len; 
            ctx.wedge_id = pcei.id; 
            
            VE_CIRCULATOR_FWD(v3, ctx, collapse_resolve_collision_I); 
            
            pcei_new.is_collapse_ok = ctx.is_collapse_ok; 
        }
    }

    if (valid_thred && pcei.is_collapse_ok && !pcei_new.is_collapse_ok)
    { /* Update info when the collapse gets canceled */
        store_per_collapse_edge_info(collapse_edge_id, pcei_new); 
        
        EdgeCollapseInfo peci; 
        peci.is_collapse_ok = false; 
        peci.edge_len = .0f; 
        peci.id = collapse_edge_id; 
        store_per_edge_collapse_info(pcei.id, peci); 
    }

    if (valid_thread && pcei_new.is_collapse_ok)
    { /* mark edge score for 2 edge verts */
        store_per_vert_collapse_wedge_id(v1, pcei.id); 
        store_per_vert_collapse_wedge_id(v3, pcei.id); 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_2)
    /* For each edge, find the min collapse score in its 4 adj verts */
    EdgeCollapseInfo pcei_new = pcei; 
    if (valid_thread && pcei.is_collapse_ok)
    {
        uvec4 collapse_links = uvec4(
            load_per_vert_collapse_wedge_id(v0), 
            load_per_vert_collapse_wedge_id(v1), 
            load_per_vert_collapse_wedge_id(v2), 
            load_per_vert_collapse_wedge_id(v3)
        )
        for (uint ivert = 0; ivert < 4; ++ivert)
        {
            uint wedge_collapsed = collapse_links[ivert]; 
            if (wedge_collapsed != NULL_EDGE)
            {
                EdgeCollapseInfo peci_adj = load_per_edge_collapse_info(wedge_collapsed, pcs_collapse_iter_); 
                if (collapse_score_larger(peci_adj, wedge_collapsed, pcei, pcei.id))
                {
                    pcei_new.is_collapse_ok = false; 
                    break; 
                }
            }
        }

        // Not needed
        // if (valid_thread && pcei.is_collapse_ok && !pcei_new.is_collapse_ok)
        // {
        //     // EdgeCollapseInfo peci; 
        //     // peci.is_collapse_ok = false; 
        //     // peci.score = .0f; 
        //     // peci.id = collapse_edge_id; 
        //     // store_per_edge_collapse_info(pcei.id, peci, pcs_collapse_iter_); 
        // }
    }
#endif


/* Apply the edge collapse */
#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE)
    if (!valid_thread)
        return;
    if (valid_thread && !psei_curr.is_collapse_ok)
    {
        EdgeFlags ef = load_edge_flags(psei_curr.id); 
        ef.del_by_collapse = false; /* cancel the collapse and quite */
        store_edge_flags(psei_curr.id, ef); 
        return; 
    }

    uint num_existing_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_in_.num_edges; 
    uint num_existing_verts = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts; 

#endif

#undef w0
#undef w1
#undef w2
#undef w3

}
#endif