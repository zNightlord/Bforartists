
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)


/* Inputs: 
uint pcs_split_iter_; 
float pcs_remesh_edge_len_; // TODO: this will be a ssbo in the future
int pcs_edge_count_; 
int pcs_vert_count_; 

!!! Remember to zero out these counters !!!
* SSBO_StrokeGenDynMeshCounters     ssbo_edge_split_counters_[]
* SSBO_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_in_

#ifdef _KERNEL_MULTICOMPILE__EDGE_SPLIT_COMPACT
.define("GLOBAL_COUNTER", "ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges_pass_1")
.define("CP_TAG", "split_select_long_edges")
#endif

#ifdef _KERNEL_MULTICOMPILE__EDGE_SPLIT_RESOLVE_CONFLICT
.define("GLOBAL_COUNTER", "ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges")
.define("CP_TAG", "split_resolve_conflict")
#endif

* float ssbo_vbo_full_[]
* uint ssbo_edge_to_vert_[]
* uint ssbo_edge_to_edges_[]
* uint ssbo_vert_to_edge_list_header_[]

* .define("EDGE_FLAGS_INCLUDED", "1")
* .define("EDGE_FLAGS_INCLUDED", "1")
* uint ssbo_vert_flags_[]
* uint ssbo_edge_flags_[]

* uint ssbo_per_edge_split_info_[]
* uint ssbo_per_split_edge_info_[]
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

float get_split_edge_len_min()
{
    return calc_remesh_edge_len_max(pcs_remesh_edge_len_); 
}


struct EdgeSplitInfo
{
    bool is_split_ok; 
    float edge_len; 
    uint id; 
}; 
uvec2 encode_edge_split_info(EdgeSplitInfo esi)
{
    uvec2 esi_enc = uvec2(0, 0); 
    esi_enc.x = floatBitsToUint(esi.edge_len); 
    esi_enc.y = ((esi.id << 1u) | uint(esi.is_split_ok));

    return esi_enc;  
}
EdgeSplitInfo decode_edge_split_info(uvec2 esi_enc)
{
    EdgeSplitInfo esi; 
    esi.edge_len = uintBitsToFloat(esi_enc.x);
    esi.is_split_ok = (0u != (esi_enc.y & 1u));
    esi_enc.y >>= 1u; 
    esi.id = esi_enc.y;  

    return esi; 
}



void store_per_split_edge_info(uint split_edge_id, EdgeSplitInfo psei)
{
    uvec2 psei_enc = encode_edge_split_info(psei); 
    Store2(ssbo_per_split_edge_info_, split_edge_id, psei_enc);
}
EdgeSplitInfo load_per_split_edge_info(uint split_edge_id)
{
    uvec2 psei_enc; 
    Load2(ssbo_per_split_edge_info_, split_edge_id, psei_enc); 
    return decode_edge_split_info(psei_enc); 
}
void store_per_edge_split_info(uint wedge_id, EdgeSplitInfo pesi)
{
    uvec2 pesi_enc = encode_edge_split_info(pesi); 
    Store2(ssbo_per_edge_split_info_, wedge_id, pesi_enc); 
}
EdgeSplitInfo load_per_edge_split_info(uint wedge_id)
{
    uvec2 pesi_enc;
    Load2(ssbo_per_edge_split_info_, wedge_id, pesi_enc); 
    return decode_edge_split_info(pesi_enc); 
}


/* true/false if esi_0 hase higher/lower priority than esi_1  */
struct EdgeSplitPriorityContext
{
    EdgeSplitInfo pesi; 
    uint wedge_id; 
}; 
bool split_priority_higher(EdgeSplitPriorityContext ctx_0, EdgeSplitPriorityContext ctx_1)
{
    if (!ctx_0.pesi.is_split_ok)
        return false;
    if (!ctx_1.pesi.is_split_ok)
        return true;

    if (ctx_0.pesi.edge_len == ctx_1.pesi.edge_len)
        return ctx_0.wedge_id < ctx_1.wedge_id; /* use lower id */

    return ctx_0.pesi.edge_len > ctx_1.pesi.edge_len; /* prefer longer edge */
}


#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_INIT)
void main()
{ /* Dispatched at beginning of each split sequence */
    if (gl_GlobalInvocationID.x == 0u)
    { 
        ssbo_edge_split_counters_[0].num_split_edges_pass_1 = 0u;
        ssbo_edge_split_counters_[0].num_split_edges = 0u; 
    }
}
#endif




#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_COMPACT)
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 

    uint wedge_id = gl_GlobalInvocationID.x; 
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
    /* reset temporary flags */
    if (pcs_split_iter_ == 0u)
    {
        ef.temp_new_by_split_this_round = false;
        if (valid_thread)
            store_edge_flags(wedge_id, ef); 
    }

    /* Select edges based on split conditions */
    bool is_split_ok = valid_thread 
        && (!ef.dupli)
        && (!ef.border) 
        && (!ef.del_by_split)
        && (!ef.del_by_collapse)
        && (!ef.temp_new_by_split_this_round) /* skip edges generated from current remesh iter */
        && edge_len > get_split_edge_len_min();

    uint split_edge_id = compact_split_select_long_edges(is_split_ok, groupId);
    if (is_split_ok)
    {
        EdgeSplitInfo psei; 
        psei.is_split_ok = true;
        psei.edge_len = edge_len; 
        psei.id = wedge_id;

        store_per_split_edge_info(split_edge_id, psei); 
    }
 
    if (valid_thread)
    {
        EdgeSplitInfo pesi; 
        pesi.is_split_ok = is_split_ok; 
        pesi.edge_len = edge_len; 
        pesi.id = is_split_ok ? split_edge_id : NULL_EDGE; 

        store_per_edge_split_info(wedge_id, pesi); 
    }
}
#endif


#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT)

#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXECUTE)
void Store4_EEAdj(uint w, uvec4 adj_wedges, uvec4 adj_iface_adjs)
{
    uvec4 data;
    data[0] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[0], adj_iface_adjs[0])); 
    data[1] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[1], adj_iface_adjs[1]));
    data[2] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[2], adj_iface_adjs[2]));
    data[3] = encode_adj_wedge_info(AdjWedgeInfo(adj_wedges[3], adj_iface_adjs[3]));

    Store4(ssbo_edge_to_edges_, w, uvec4, data); 
}
#endif

void main()
{
    const uint groupId = gl_LocalInvocationID.x; 

    uint split_edge_id = gl_GlobalInvocationID.x;
    uint num_split_edges = ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges_pass_1; 
    bool valid_thread = split_edge_id < num_split_edges;


    EdgeSplitInfo psei_curr = load_per_split_edge_info(split_edge_id);

    AdjWedgeInfo w[4]; 
    w[0] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 0u]);
    w[1] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 1u]);
    w[2] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 2u]);
    w[3] = decode_adj_wedge_info(ssbo_edge_to_edges_[psei_curr.id*4u + 3u]);



#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_RESOLVE_CONFLICT)
    EdgeSplitInfo pesi[4];
    pesi[0] = load_per_edge_split_info(w[0].wedge_id);
    pesi[1] = load_per_edge_split_info(w[1].wedge_id);
    pesi[2] = load_per_edge_split_info(w[2].wedge_id);
    pesi[3] = load_per_edge_split_info(w[3].wedge_id);

    /* Resolve collision */
    if (psei_curr.is_split_ok && valid_thread)
    {
        EdgeSplitPriorityContext ctx_curr = EdgeSplitPriorityContext(psei_curr, psei_curr.id); 
        psei_curr.is_split_ok = all(
            bvec4(
                split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[0], w[0].wedge_id)),
                split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[1], w[1].wedge_id)),
                split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[2], w[2].wedge_id)),
                split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[3], w[3].wedge_id))
            )
        ); 
    }

    /* Compaction */
    bool split = psei_curr.is_split_ok && valid_thread; 
    uint split_edge_alloc_id = compact_split_resolve_conflict(split, groupId);

    /* Dont need score from now on, so reuse it to store alloc id */
    if (valid_thread)
    {
        psei_curr.edge_len = uintBitsToFloat(split ? split_edge_alloc_id : NULL_EDGE);
        store_per_split_edge_info(split_edge_id, psei_curr); /* no race cond' since only accessed by this thread */
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXECUTE)
    /* Apply edge split */
    if (!psei_curr.is_split_ok || !valid_thread)
        return;

    uint num_existing_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_in_.num_edges; 
    uint num_existing_verts = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts; 

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
    uint v0, v1, v2, v3; 
    v0 = ssbo_edge_to_vert_[psei_curr.id*4u + 0u];
    v1 = ssbo_edge_to_vert_[psei_curr.id*4u + 1u];
    v2 = ssbo_edge_to_vert_[psei_curr.id*4u + 2u];
    v3 = ssbo_edge_to_vert_[psei_curr.id*4u + 3u];

    uint split_edge_alloc_id = floatBitsToUint(psei_curr.edge_len);
    /* alloc x4 new edges */
    uvec4 e = (num_existing_edges + split_edge_alloc_id * 4u).xxxx + uvec4(0, 1, 2, 3); 
    /* alloc x1 new vertex */
    uint v4 = (num_existing_verts + split_edge_alloc_id); 

#define w0 ((w[0].wedge_id))
#define w1 ((w[1].wedge_id))
#define w2 ((w[2].wedge_id))
#define w3 ((w[3].wedge_id))

#define e0 e[0]
#define e1 e[1]
#define e2 e[2]
#define e3 e[3]
    /* Build ee & ev adj. for NEW edges
    * --------------------------------------------------------------                              
    *     Template               E0                    E2
    *        v3                  v0            v1 =E1= v4 =E3= v3                
    *      / | \                /||\             \ <-- || <-- /   
    *    W3 .|  W2             / || \             \ f1 || f0 /    
    *    / /|||  \            W  ||  3             W   E2   2     
    *  v0 f1|W|f0 v2         0   E0   W             1  ||  W      
    *    \  |||/ /          / f1 || f0 \             \ || /       
    *    W0 ||. W1         / --> || --> \             \||/        
    *      \ | /         v1 =E1= v4 =E3= v3            v2         
    *        v1                  
    *    { v0~3 }        { v1, v4, v3, v0 }      { v1, v2, v3, v4 }                 
    *    { w0~3 }        { e1, e3, w3, w0 }      { w1, w2, e3, e1 } 
    */
    Store4(ssbo_edge_to_vert_,  e0, uvec4, uvec4(v1, v4, v3, v0)); // update e0
    Store4_EEAdj(e0, uvec4(e1, e3, w3, w0), uvec4(0, 0, w[3].iface_adj, w[0].iface_adj)); 
    Store4(ssbo_edge_to_vert_,  e2, uvec4, uvec4(v1, v2, v3, v4)); // update e2
    Store4_EEAdj(e2, uvec4(w1, w2, e3, e1), uvec4(w[1].iface_adj, w[2].iface_adj, 1, 1)); 
    /*  ---------------------------------------------------------------                               
    *     Template               E1                    E3
    *                            v0                    v0        
    *                           /||                    ||\       
    *       v0                 / ||                    || \      
    *      /  \               W  ||                    ||  3     
    *    W0 f1 W3            0   E0                    E0   W    
    *    / ---> \           / f1 ||                    || f1 \   
    *  v1-- W  --v3        / --> ||                    || --> \  
    *    \ <--- /        v1 =E1= v4                    v4 =E3= v3
    *    W1 f0  W2         \ <-- ||                    || <-- /   
    *      \  /             \ f0 ||                    || f0 /   
    *       v2               W   E2                    E2   2    
    *                         1  ||                    ||  W     
    *                          \ ||                    || /      
    *                           \||                    ||/       
    *                            v2                    v2    
    *    { v0~3 }        { v0, v1, v2, v4 }      { v0, v4, v2, v3 }
    *    { w0~3 }        { w0, w1, e2, e0 }      { e0, e2, w2, w3 }
    */
    Store4(ssbo_edge_to_vert_,  e1, uvec4, uvec4(v0, v1, v2, v4)); // update e1
    Store4_EEAdj(e1, uvec4(w0, w1, e2, e0), uvec4(w[0].iface_adj, w[1].iface_adj, 0, 0)); 
    Store4(ssbo_edge_to_vert_,  e3, uvec4, uvec4(v0, v4, v2, v3)); // update e3
    Store4_EEAdj(e3, uvec4(e0, e2, w2, w3), uvec4(1, 1, w[2].iface_adj, w[3].iface_adj)); 

    /* Adjust ee & ev adj. for OLD edges 
    * 
    *         v0         
    *        /||\        
    *       / || \       
    *      W  ||  3      
    *     0 f1E0f0 W     
    *    /    ||    \    
    *   / f1->||->f1 \   
    * v1 =E1= v4 =E3= v3 
    *   \ f0<-||<-f0 /   
    *    \    ||    /    
    *     W f1E2f0 2     w0: { {e1,f0}, {e0,f0} } // AdjWedgeInfo data
    *      1  ||  W      w1: { {e2,f0}, {e1,f1} }
    *       \ || /       w2: { {e3,f1}, {e2,f1} }
    *        \||/        w3: { {e0,f1}, {e3,f0} }
    *         v2         
    */
    AdjWedgeInfo awi_updates[8] = {
        AdjWedgeInfo(e1, 0), AdjWedgeInfo(e0, 0), // w0
        AdjWedgeInfo(e2, 0), AdjWedgeInfo(e1, 1), // w1
        AdjWedgeInfo(e3, 1), AdjWedgeInfo(e2, 1), // w2
        AdjWedgeInfo(e0, 1), AdjWedgeInfo(e3, 0)  // w3
    }; 
    
    for (uint i = 0; i < 4; ++i)
    {
        uint iface_overlap = w[i].iface_adj == 0 ? 1 : 0; 
        uvec2 iwedge_update = {
            mark__he_to_wedge(iface_overlap, 1),  
            mark__he_to_wedge(iface_overlap, 2)
        }; 
        ssbo_edge_to_edges_[w[i].wedge_id*4 + iwedge_update[0]] = encode_adj_wedge_info(awi_updates[i*2+0]); 
        ssbo_edge_to_edges_[w[i].wedge_id*4 + iwedge_update[1]] = encode_adj_wedge_info(awi_updates[i*2+1]); 
        
        uint ivert_update = mark__center_wedge_to_oppo_vert__at_face(iface_overlap); 
        ssbo_edge_to_vert_[w[i].wedge_id*4 + ivert_update] = v4; 
    }


    /* Store flags for 4 new edges */
    EdgeFlags ef = init_edge_flags__new_split_edge();
    store_edge_flags(e0, ef);
    store_edge_flags(e1, ef);
    store_edge_flags(e2, ef);
    store_edge_flags(e3, ef);

    /* Mark old edge as deleted */
    EdgeFlags ef_del = load_edge_flags(psei_curr.id);
    ef_del.del_by_split = true; 
    store_edge_flags(psei_curr.id, ef_del);


    /* Build ve adj. for New vert */
    VertWedgeListHeader vwlh; 
    vwlh.wedge_id = e0; 
    vwlh.ivert = 1u; /* check the graph for e0 above */
    ssbo_vert_to_edge_list_header_[v4] = encode_vert_wedge_list_header(vwlh); 


    /* Store new vert pos */ 
    vec3 split_pos = (ld_vpos(v1) + ld_vpos(v3)) / 2.0f; 
    if (valid_thread)
        st_vpos(v4, split_pos); 

    /* Store new vert states */
    VertFlags vf = init_vert_flags__new_split_edge(); 
    store_vert_flags(v4, vf); 


    /* Update dynamic mesh counters */
    if (split_edge_alloc_id == 0u) // !!! don't do this: split_edge_id == 0u because 0th thread might early exit 
    {
        /* accumulate global mesh counters */
        ssbo_dyn_mesh_counters_out_.num_edges = 
            ssbo_dyn_mesh_counters_in_.num_edges + 4u * ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges;
        ssbo_dyn_mesh_counters_out_.num_verts = 
            ssbo_dyn_mesh_counters_in_.num_verts + ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges; 
        
        /* prep local counters for next split iteration */
        ssbo_edge_split_counters_[pcs_split_iter_ + 1].num_split_edges_pass_1 = 0u;
        ssbo_edge_split_counters_[pcs_split_iter_ + 1].num_split_edges = 0u; 
    }

#undef e0
#undef e1
#undef e2
#undef e3

#undef w0
#undef w1
#undef w2
#undef w3

#endif

 }
#endif


