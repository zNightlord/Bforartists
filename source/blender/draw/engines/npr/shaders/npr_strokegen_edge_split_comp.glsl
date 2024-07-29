
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)

#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_loop_subdiv_edge_tree_lib.glsl)

#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)


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

float get_split_edge_len_min(uint v1, uint v3)
{
    float vlen_1 = ld_vtx_remesh_len(v1); 
    float vlen_3 = ld_vtx_remesh_len(v3);
    float targ_len = min(vlen_1, vlen_3); 
    return (4.0f/3.0f) * targ_len; 
    // return pcs_remesh_edge_len_; // use this to stress-test the split kernel
}

float estimate_split_vert_remesh_edge_len(uint v1, uint v3)
{
    float vlen_1 = ld_vtx_remesh_len(v1); 
    float vcurv_1 = average_adaptive_remesh_len_to_vcurv(vlen_1); 

    float vlen_3 = ld_vtx_remesh_len(v3);
    float vcurv_3 = average_adaptive_remesh_len_to_vcurv(vlen_3); 

    /* It is better to calc form curvatures */
    return get_adaptive_remesh_len(.5f * (vcurv_1 + vcurv_3)); 
}


struct EdgeSplitInfo
{
    bool is_split_ok; 
    float edge_len; 
    uint id; 
}; 
uint encode_edge_split_info__id__is_collapse_ok(uint edge_id, bool ok)
{
    return ((edge_id << 1u) | uint(ok));
}
uvec2 encode_edge_split_info(EdgeSplitInfo esi)
{
    uvec2 esi_enc = uvec2(0, 0); 
    esi_enc.x = floatBitsToUint(esi.edge_len); 
    esi_enc.y = encode_edge_split_info__id__is_collapse_ok(esi.id, esi.is_split_ok);

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
void store_per_split_edge_info__id__is_split_ok(uint split_edge_id, uint wedge_id, bool is_split_ok)
{
    uint enc = encode_edge_split_info__id__is_collapse_ok(wedge_id, is_split_ok); 
    ssbo_per_split_edge_info_[split_edge_id*2u+1u] = enc;
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
void store_per_edge_split_info__id__is_split_ok(uint wedge_id, uint split_edge_id, bool is_split_ok)
{
    uint enc = encode_edge_split_info__id__is_collapse_ok(split_edge_id, is_split_ok); 
    ssbo_per_edge_split_info_[wedge_id*2u+1u] = enc;
}
EdgeSplitInfo load_per_edge_split_info(uint wedge_id)
{
    uvec2 pesi_enc;
    Load2(ssbo_per_edge_split_info_, wedge_id, pesi_enc); 
    return decode_edge_split_info(pesi_enc); 
}


uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}





#define EDGE_SPLIT_LONG_EDGES 0u
#define EDGE_SPLIT_CONTOUR_EDGES 1u 
#define EDGE_SPLIT_LOOP_SUBDIV 2u 

bool is_contour_split_pass() { return (pcs_split_mode_ == EDGE_SPLIT_CONTOUR_EDGES); }
bool is_loop_subdiv_pass() { return (pcs_split_mode_ == EDGE_SPLIT_LOOP_SUBDIV); } 

bool is_contour_edge_havent_split(vec3 vpos_0, vec3 vnor_0, VertFlags vf_0, vec3 vpos_1, vec3 vnor_1, VertFlags vf_1, EdgeFlags ef)
{
    if (vf_0.contour || vf_1.contour)
        return false;

    return is_interp_contour_edge__before_tessellation(vf_0, vf_1, ef); 
}

bool is_loop_subdiv_havent_split(EdgeFlags ef, VertFlags vf_0, VertFlags vf_1)
{
    return ((!ef.new_by_split));
}


vec3 calc_split_vpos(uint v1, uint v3, uint wedge_id)
{
    vec3 edge_vpos[2] = { ld_vpos(v1), ld_vpos(v3) };
    if (is_contour_split_pass())
    { /* interpolated contour */
        vec3 vnor[2] = { ld_vnor(v1), ld_vnor(v3) }; 
        mat4 view_to_world = ubo_view_matrices_.viewinv;
        vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
        return calc_interp_contour_vert_pos(vnor, edge_vpos, cam_pos_ws); 
        // float contour_interpo_factor = .5f; 
        // {        
        //     mat4 view_to_world = ubo_view_matrices_.viewinv;
        //     vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

        //     vec2 ndv = vec2(
        //         dot(vnor[0], cam_pos_ws - edge_vpos[0]),  
        //         dot(vnor[1], cam_pos_ws - edge_vpos[1]) 
        //     ); 
        //     contour_interpo_factor = ndv[0] / (ndv[0] - ndv[1]); /* split pos = vpos_0 + interpo * (vpos_1 - vpos_0); */
        // }

        // return edge_vpos[0] + contour_interpo_factor * (edge_vpos[1] - edge_vpos[0]); 
    }
    if (is_loop_subdiv_pass())
    { /* loop subdiv for new verts */
        uvec3 edge_point_pos_enc;
        Load3(ssbo_epos_subd_, wedge_id, edge_point_pos_enc); 
        return uintBitsToFloat(edge_point_pos_enc); 
    }

    // Return mid point by default
    return (edge_vpos[0] + edge_vpos[1]) / 2.0f;
}



/* true/false if esi_0 hase higher/lower priority than esi_1  */
struct EdgeSplitPriorityContext
{
    EdgeSplitInfo pesi; 
    uint wedge_id; 
    bool selected; /* when !selected, pesi is invalid */
}; 
bool split_priority_higher(EdgeSplitPriorityContext ctx_0, EdgeSplitPriorityContext ctx_1)
{
    if ((!ctx_0.pesi.is_split_ok) || (!ctx_0.selected))
        return false;
    if ((!ctx_1.pesi.is_split_ok) || (!ctx_1.selected))
        return true;

    /* -------- Simple Split (no special priority) --------- */
    if (pcs_split_mode_ == EDGE_SPLIT_LOOP_SUBDIV
        || pcs_split_mode_ == EDGE_SPLIT_CONTOUR_EDGES)
    {
        uint key_w0 = pcg(ctx_0.wedge_id * pcs_split_iter_); 
        uint key_w1 = pcg(ctx_1.wedge_id * pcs_split_iter_);  
        if (key_w0 == key_w1)
            return (ctx_0.wedge_id) < (ctx_1.wedge_id); /* use lower id */
        return key_w0 < key_w1; /* use lower id */
    }

    /* -------- Split based on edge length --------- */
    /* Sort Key: [wid | pcg(wid) | edge_len] */
    if (ctx_0.pesi.edge_len == ctx_1.pesi.edge_len)
    { /* use pcg to break tie */
        uint pcg_w0 = pcg(ctx_0.wedge_id + pcs_split_iter_);
        uint pcg_w1 = pcg(ctx_1.wedge_id + pcs_split_iter_);
        if (pcg_w0 == pcg_w1)
            return (ctx_0.wedge_id) < (ctx_1.wedge_id); /* use lower id */
        return pcg_w0 < pcg_w1; /* use lower id */
    }

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
    
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 
    
    if (gl_GlobalInvocationID.x == 0u)
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

   
   
    /* Select edges based on split conditions */
    EdgeFlags ef = load_edge_flags(wedge_id);
    if (valid_thread && pcs_split_iter_ == 0u)
    {
        update_edge_flags__reset_new_split_edge(wedge_id, ef); 
        ef.new_by_split = false; 
    }
    
    bool is_contour_edge = false; 
    if (is_contour_split_pass())
        is_contour_edge = is_contour_edge_havent_split(
            vpos_cwedge[0], ld_vnor(verts_cwedge.x), load_vert_flags(verts_cwedge.x), 
            vpos_cwedge[1], ld_vnor(verts_cwedge.y), load_vert_flags(verts_cwedge.y), 
            ef
        );

    bool is_loop_subd_edge = false;
    if (is_loop_subdiv_pass()) 
        is_loop_subd_edge = is_loop_subdiv_havent_split(
            ef, 
            load_vert_flags(verts_cwedge.x), 
            load_vert_flags(verts_cwedge.y)
        ); 

    bool is_split_ok = valid_thread 
        && (ef.selected)
        && (!ef.sel_border)
        && (!ef.dupli)
        && (!ef.border) 
        && (!ef.del_by_split)
        && (!ef.del_by_collapse)
        && (!ef.new_by_split); 
    if (is_contour_split_pass())
        is_split_ok = is_split_ok && is_contour_edge;
    else if (is_loop_subdiv_pass())
        is_split_ok = is_split_ok && is_loop_subd_edge; 
    else
        is_split_ok = is_split_ok && edge_len > get_split_edge_len_min(verts_cwedge.x, verts_cwedge.y);
    

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

    Store4(ssbo_edge_to_edges_, w, data); 
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



#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXCLUDE_BORDER)
    bool update = false; 
    if (psei_curr.is_split_ok && valid_thread)
    {  
        EdgeFlags ef = load_edge_flags(psei_curr.id); 
        EdgeFlags ef_w0 = load_edge_flags(w[0].wedge_id); 
        EdgeFlags ef_w1 = load_edge_flags(w[1].wedge_id); 
        EdgeFlags ef_w2 = load_edge_flags(w[2].wedge_id); 
        EdgeFlags ef_w3 = load_edge_flags(w[3].wedge_id); 

        /* except border edges */
        if (any(bvec4(ef_w0.border, ef_w1.border, ef_w2.border, ef_w3.border))) 
            psei_curr.is_split_ok = false; 

        /* except edges at selection border */
        if (!all(bvec4(ef_w0.selected, ef_w1.selected, ef_w2.selected, ef_w3.selected)))
            psei_curr.is_split_ok = false; 

        update = !psei_curr.is_split_ok; 
    }

    if (update && valid_thread)
    {
        store_per_edge_split_info__id__is_split_ok(psei_curr.id, split_edge_id, psei_curr.is_split_ok);  
        store_per_split_edge_info__id__is_split_ok(split_edge_id, psei_curr.id, psei_curr.is_split_ok); 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__EDGE_SPLIT_RESOLVE_CONFLICT)
    if (psei_curr.is_split_ok && valid_thread)
    {  
        EdgeSplitInfo pesi[4];
        pesi[0] = load_per_edge_split_info(w[0].wedge_id);
        pesi[1] = load_per_edge_split_info(w[1].wedge_id);
        pesi[2] = load_per_edge_split_info(w[2].wedge_id);
        pesi[3] = load_per_edge_split_info(w[3].wedge_id);

        EdgeFlags ef = load_edge_flags(psei_curr.id); 
        EdgeFlags ef_w0 = load_edge_flags(w[0].wedge_id); 
        EdgeFlags ef_w1 = load_edge_flags(w[1].wedge_id); 
        EdgeFlags ef_w2 = load_edge_flags(w[2].wedge_id); 
        EdgeFlags ef_w3 = load_edge_flags(w[3].wedge_id); 

        // /* except border edges */
        // if (any(bvec4(ef_w0.border, ef_w1.border, ef_w2.border, ef_w3.border))) 
        //     psei_curr.is_split_ok = false; 
        // /* except edges at selection border */
        // if (!all(bvec4(ef_w0.selected, ef_w1.selected, ef_w2.selected, ef_w3.selected)))
        //     psei_curr.is_split_ok = false; 

        /* Resolve collision */
        EdgeSplitPriorityContext ctx_curr = EdgeSplitPriorityContext(psei_curr, psei_curr.id, ef.selected); 
        bool winner = 
            split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[0], w[0].wedge_id, ef_w0.selected)) && 
            split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[1], w[1].wedge_id, ef_w1.selected)) && 
            split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[2], w[2].wedge_id, ef_w2.selected)) && 
            split_priority_higher(ctx_curr, EdgeSplitPriorityContext(pesi[3], w[3].wedge_id, ef_w3.selected)); 
        if (!winner) psei_curr.is_split_ok = false; 
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
    /* replace ssbo binding since we could not squeeze any new slots */
    // if (is_loop_subdiv_pass())
    #define ssbo_subd_edge_tree_node_up_ ssbo_selected_edge_to_edge_ 
    #define ssbo_subd_edge_vert_to_old_edge_ ssbo_vnor_ // swapped when split for loop subdiv
    // if (is_contour_split_pass())
    #define ssbo_contour_vert_to_old_edge_ ssbo_selected_edge_to_edge_ // swapped when split for interpolated contour

    /* Update dynamic mesh counters */
    if (gl_GlobalInvocationID.x == 0u) 
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

    /* Early Exits */
    if (!valid_thread)
        return;
    if (!psei_curr.is_split_ok)
        return; 




    /* Apply edge split */
    uint num_existing_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_in_.num_edges; 
    uint num_existing_verts = pcs_vert_count_ + ssbo_dyn_mesh_counters_in_.num_verts; 
    
    /* Inherit old edge tags */
    EdgeFlags ef_old = load_edge_flags(psei_curr.id); 


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
    Store4(ssbo_edge_to_vert_,  e0, uvec4(v1, v4, v3, v0)); // update e0
    Store4_EEAdj(e0, uvec4(e1, e3, w3, w0), uvec4(0, 0, w[3].iface_adj, w[0].iface_adj)); 
    Store4(ssbo_edge_to_vert_,  e2, uvec4(v1, v2, v3, v4)); // update e2
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
    Store4(ssbo_edge_to_vert_,  e1, uvec4(v0, v1, v2, v4)); // update e1
    Store4_EEAdj(e1, uvec4(w0, w1, e2, e0), uvec4(w[0].iface_adj, w[1].iface_adj, 0, 0)); 
    Store4(ssbo_edge_to_vert_,  e3, uvec4(v0, v4, v2, v3)); // update e3
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
    EdgeFlags ef = init_edge_flags__new_split_edge(false);
    ef.new_by_split_on_old_edge = false; 
    store_edge_flags(e0, ef);
    store_edge_flags(e2, ef);
    ef.new_by_split_on_old_edge = true; 
    ef.crease_level = ef_old.crease_level; 
    store_edge_flags(e1, ef);
    store_edge_flags(e3, ef);

    /* Mark old edge as deleted */
    update_edge_flags__del_by_split(psei_curr.id);

    /* Update ve adj. for Old verts, if thery were connected to the split edge */
    try_update_ve_link(v1, psei_curr.id, VertWedgeListHeader(w0, mark__cwedge_to_beg_vert(w[0].iface_adj))); 
    try_update_ve_link(v3, psei_curr.id, VertWedgeListHeader(w2, mark__cwedge_to_beg_vert(w[2].iface_adj)));
    /* Build ve adj. for New vert */
    VertWedgeListHeader vwlh; 
    vwlh.wedge_id = e0; 
    vwlh.ivert = 1u; /* check the graph for e0 above */
    ssbo_vert_to_edge_list_header_[v4] = encode_vert_wedge_list_header(vwlh); 


    /* Store new vert pos */ 
    vec3 split_pos = calc_split_vpos(v1, v3, psei_curr.id); 
    st_vpos(v4, split_pos); 

    /* Store new vert states */
    VertFlags vf = init_vert_flags__new_split_edge(is_contour_split_pass(), 0u < ef_old.crease_level); 
    store_vert_flags(v4, vf); 

    /* Store adaptive remesh length */
    float edge_len = estimate_split_vert_remesh_edge_len(v1, v3);
    st_vtx_remesh_len(v4, edge_len);

    /* Build subd edge tree */
    if (is_loop_subdiv_pass())
    { 
        LoopSubdEdgeTreeUpNode par_node = decode_loop_subd_tree_node(ssbo_subd_edge_tree_node_up_[psei_curr.id]); 

        /* Calc tree nodes for E1, E3 */
        LoopSubdEdgeTreeUpNode node_e1 = setup_loop_subd_tree_leaf__split_edge(psei_curr.id, par_node, 1u);
        ssbo_subd_edge_tree_node_up_[e1] = encode_loop_subd_tree_node(node_e1); 
        LoopSubdEdgeTreeUpNode node_e3 = setup_loop_subd_tree_leaf__split_edge(psei_curr.id, par_node, 3u);
        ssbo_subd_edge_tree_node_up_[e3] = encode_loop_subd_tree_node(node_e3); 

        /* Init tree nodes for E0, E2, the actual calculation is carried later */
        LoopSubdEdgeTreeUpNode node = init_loop_subd_tree_leaf__face_edge(); // point to nothing 
        ssbo_subd_edge_tree_node_up_[e0] = encode_loop_subd_tree_node(node);
        ssbo_subd_edge_tree_node_up_[e2] = encode_loop_subd_tree_node(node);

        /* In order to eval node for new edge,
         * we cache link from each split vert to its old edge */
        ssbo_subd_edge_vert_to_old_edge_[v4] = psei_curr.id; 
    }
    if (is_contour_split_pass())
    {   /* In order to eval node for contour edge,
         * we cache link from each contour vert to its old edge */
        ssbo_contour_vert_to_old_edge_[v4] = psei_curr.id; 
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


