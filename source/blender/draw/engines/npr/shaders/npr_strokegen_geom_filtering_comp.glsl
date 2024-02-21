#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)


/* Refs: 
 * "Bilateral Normal Filtering for Mesh Denoising"
*/
float gaussian(float dist, float tau)
{
    return exp(-(dist * dist) / (2.0f * tau * tau)); 
}
float vnor_bilateral_filter_weight(vec3 vnor_i, vec3 vpos_i, vec3 vnor_j, vec3 vpos_j)
{
    const float tau_vnor = .2f; 
    float w_vnor = gaussian(distance(vnor_i, vnor_j), tau_vnor); 

    // note: positional weight needs to be added
    // const float tau_vpos = ?;
    float w_vpos = 1.0f; 

    return w_vnor * w_vpos; 
}
float vpos_bilateral_filter_weight(vec3 vnor_i, vec3 vpos_i, vec3 vnor_j, vec3 vpos_j)
{
    const float tau_vnor = 0.1f; 
    float w_vnor = gaussian(distance(vnor_i, vnor_j), tau_vnor); 

    // note: positional weight needs to be added
    // const float tau_vpos = ?;
    float w_vpos = 1.0f; 

    return w_vnor * w_vpos; 
}

/* Inputs: 
ssbo_vbo_full_[]
ssbo_vert_to_edge_list_header_[] 
pcs_vert_count_
ssbo_dyn_mesh_counters_out_
ssbo_bnpr_mesh_pool_counters_
*/

uint get_vert_count()
{
    return pcs_vert_count_ + ssbo_dyn_mesh_counters_out_.num_verts; 
}

void st_vnor_filtered(uint vert_id, vec3 vnor_filtered)
{
    uvec3 vnor_filtered_enc = floatBitsToUint(vnor_filtered); 
    Store3(ssbo_vnor_temp_out_, vert_id, vnor_filtered_enc);
}
vec3 ld_vnor_filtered(uint vert_id)
{
    uvec3 vnor_filtered_enc; 
    Load3(ssbo_vnor_temp_in_, vert_id, vnor_filtered_enc);
    return uintBitsToFloat(vnor_filtered_enc); 
}

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VNOR_FILTERING)
    struct VtxNormalFilteringContext
    {
        vec3 vpos; 
        vec3 vnor;
        vec3 vnor_filtered; 
        float weight_sum; 
        
        bool border; 
        bool sel_border; 
    }; 

    VtxNormalFilteringContext init_vnor_filter_context(
        vec3 vpos, vec3 vnor
    ){
        VtxNormalFilteringContext ctx; 
        
        ctx.vpos = vpos;
        ctx.vnor = vnor; 
        ctx.vnor_filtered = vec3(.0f); 
        ctx.weight_sum = .0f; 
        
        ctx.border = false;
        ctx.sel_border = false; 

        return ctx; 
    }

    /* Example: Rotate fwd(CW) around V3(marked as vc)        
    *        v0  ...  vp             
    *       /  \     /  \            
    *      / f1 \   wp   \           
    *     / ---> \ /      \          
    *   v1 ====== vc--wi--vi         
    *     \ <--- / \<-----/ \        
    *      \ f0 /  wn fi wo  \       
    *       \  /     \  /     \      
    *        v2  ...  vn ----- v_oppo
    */
    bool ve_circulator__vnor_filter(
        CirculatorIterData iter, 
        inout VtxNormalFilteringContext ctx
    ){
        uint wi = iter.awi.wedge_id; 
        uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
        uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 

        vec3 vpos_i = ld_vpos(vi);
        vec3 vnor_i; 
        if (pcs_vnor_filtering_iter_ == 0)
            vnor_i = ld_vnor(vi); 
        else
            vnor_i = ld_vnor_filtered(vi);

        float weight = vnor_bilateral_filter_weight(vnor_i, vpos_i, ctx.vnor, ctx.vpos); 
        ctx.vnor_filtered += vnor_i * weight; 
        ctx.weight_sum    += weight; 

        EdgeFlags ef = load_edge_flags(wi); 
        if (ef.border) ctx.border = true;
        if (!ef.selected) ctx.sel_border = true; 

        return true; 
    }
#endif





#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING)
    struct VtxPositionFilteringContext
    {
        vec3 vpos;
        vec3 vpos_filtered; 
        float weight_sum; 
        vec3 vnor; 
        
        bool border; 
        bool sel_border; 
        float ave_edge_len; 
        int num_edges; 
    }; 

    VtxPositionFilteringContext init_vpos_filter_context(
        vec3 vpos, vec3 vnor
    ){
        VtxPositionFilteringContext ctx; 
        
        ctx.vpos = vpos;
        ctx.vpos_filtered = vec3(.0f); 
        ctx.weight_sum = .0f; 
        ctx.vnor = vnor;
        
        ctx.border = false;
        ctx.sel_border = false; 
        ctx.ave_edge_len = .0f; 
        ctx.num_edges = 0; 

        return ctx; 
    }
    /* Example: Rotate fwd(CW) around V3(marked as vc)        
    *        v0  ...  vp             
    *       /  \     /  \            
    *      / f1 \   wp   \           
    *     / ---> \ /      \          
    *   v1 ====== vc--wi--vi         
    *     \ <--- / \<-----/ \        
    *      \ f0 /  wn fi wo  \       
    *       \  /     \  /     \      
    *        v2  ...  vn ----- v_oppo
    */
    bool ve_circulator__vpos_filter(
        CirculatorIterData iter, 
        inout VtxPositionFilteringContext ctx
    ){
        uint wi = iter.awi.wedge_id; 
        uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
        uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 

        vec3 vpos_i = ld_vpos(vi); 
        vec3 vnor_i = ld_vnor_filtered(vi); 

        float w = vpos_bilateral_filter_weight(vnor_i, vpos_i, ctx.vnor, ctx.vpos); 
        ctx.vpos_filtered += w * vnor_i * dot(vnor_i, vpos_i - ctx.vpos); 
        ctx.weight_sum += 1.0f * w; 

        EdgeFlags ef = load_edge_flags(wi); 
        if (ef.border) ctx.border = true;
        if (!ef.selected) ctx.sel_border = true;

        float edge_len = length(vpos_i - ctx.vpos); 
        ctx.ave_edge_len += edge_len; 
        ctx.num_edges += 1; 

        return true; 
    }
#endif


void main()
{
    uint sel_vert_id = gl_GlobalInvocationID.x; 
    uint vert_id; bool valid_thread; 
    get_vert_id_from_selected_vert(sel_vert_id, /*out*/vert_id, valid_thread);

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VNOR_FILTERING)
    if(!valid_thread) return;

    vec3 vnor = ld_vnor(vert_id); 
    vec3 vpos = ld_vpos(vert_id);

    VtxNormalFilteringContext ctx = init_vnor_filter_context(vpos, vnor);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vnor_filter, ctx, rot_fwd);

    ctx.vnor_filtered /= ctx.weight_sum;

    vec3 vnor_filtered = normalize(ctx.vnor_filtered); 
    if (ctx.border || ctx.sel_border) 
        vnor_filtered = vnor;

    if (valid_thread)
        st_vnor_filtered(vert_id, vnor_filtered);
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING)
    if(!valid_thread) return;
    
    vec3 vpos = ld_vpos(vert_id);
    vec3 vnor = ld_vnor(vert_id); 
    
    VtxPositionFilteringContext ctx = init_vpos_filter_context(vpos, vnor);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true; 
    VE_CIRCULATOR(vwlh, ve_circulator__vpos_filter, ctx, rot_fwd); 

    ctx.vpos_filtered /= ctx.weight_sum;
    ctx.ave_edge_len /= float(ctx.num_edges);

    vec3 vpos_new = vpos + ctx.vpos_filtered; 
    if (ctx.border || ctx.sel_border) /* keep feature verts */
        vpos_new = vpos; 
    // if (length(vpos_new - vpos) > ctx.ave_edge_len * 0.5f)
    //     vpos_new = vpos;
    
    if (valid_thread)
        Store3(ssbo_vpos_temp_, sel_vert_id, floatBitsToUint(vpos_new)); 
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING_FINISH)
    if(!valid_thread) return;
    
    uvec3 vpos_new_enc;
    Load3(ssbo_vpos_temp_, sel_vert_id, vpos_new_enc);
    vec3 vpos_new = uintBitsToFloat(vpos_new_enc);

    if (valid_thread)
        st_vpos(vert_id, vpos_new);
#endif
}