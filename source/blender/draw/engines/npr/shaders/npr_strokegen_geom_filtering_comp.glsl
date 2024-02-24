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
    const float tau_vnor = 0.2f; 
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
    #define TEST_SQRT_3_INTERP 1
    float sqrt_3_subdiv_param(float num_adj_verts)
    {
        float alpha_n = (4.0f - 2.0f * cos((2.0f * 3.14159265359f) / num_adj_verts)) / 9.0f; 
        return alpha_n;
    }

    struct VtxPositionFilteringContext
    {
        vec3 vpos;
        vec3 vpos_p; 
        vec3 vpos_filtered; 
        float weight_sum; 
        vec3 vnor; 
        vec3 vnor_p; 
        
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
        ctx.vpos_p = vec3(.0f); // init in circulator
        ctx.vpos_filtered = vec3(.0f); 
        ctx.weight_sum = .0f; 
        ctx.vnor = vnor;
        ctx.vnor_p = vec3(.0f); // init in circulator 
        
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
        
        // Detect locked verts
        EdgeFlags ef = load_edge_flags(wi); 
        if (ef.border) ctx.border = true;
        if (!ef.selected) ctx.sel_border = true;

        if (iter.rotate_step == 0u)
        { /* init position of vp */
            uint ivert_vp = mark__ve_circ_fwd__get_vp(iter); 
            uint vp = ssbo_edge_to_vert_[wi*4u + ivert_vp]; 
            ctx.vpos_p = ld_vpos(vp); 
            ctx.vnor_p = ld_vnor_filtered(vp); 
        }

        // Filter position
        vec3 filter_pos = vpos_i; 
        vec3 filter_nor = vnor_i; 
        
        float area = length(cross(vpos_i - ctx.vpos, ctx.vpos_p - ctx.vpos)) * 0.5f; 
        float w = /* area */1.0f; // vpos_bilateral_filter_weight(filter_nor, ci, ctx.vnor, ctx.vpos); 
#if defined(TEST_SQRT_3_INTERP)
        ctx.vpos_filtered += vpos_i; 
#else
        ctx.vpos_filtered += w * (ctx.vpos + filter_nor * dot(filter_nor, filter_pos - ctx.vpos)); 
#endif
        ctx.weight_sum += w; 


        // Accumulate edge length 
        float edge_len = length(vpos_i - ctx.vpos); 
        ctx.ave_edge_len += edge_len; 
        ctx.num_edges += 1; 

        ctx.vpos_p = vpos_i;
        ctx.vnor_p = vnor_i; 
        return true; 
    }

    struct VtxPositionValidationContext
    {
        vec3 vpos;
        vec3 vpos_p; 
        vec3 vpos_filtered; 
        vec3 vnor; 
        vec3 vnor_p; 

        float min_dist_proj; 
        vec3 vpos_filtered_proj; // proj to 1-ring neighbors

        bool is_move_ok; 
    }; 
    VtxPositionValidationContext init_filtered_vpos_validate_context(
        vec3 vpos, vec3 vpos_filtered, vec3 vnor
    ){
        VtxPositionValidationContext ctx; 
        
        ctx.vpos = vpos;
        ctx.vpos_p = vec3(.0f); // init in circulator
        ctx.vpos_filtered = ctx.vpos_filtered_proj = vpos_filtered; 
        ctx.vnor = vnor;
        ctx.vnor_p = vec3(.0f); // init in circulator 

        ctx.min_dist_proj = 1000000.0f; 
        ctx.is_move_ok = true; 

        return ctx; 
    }
    bool ve_circulator__vpos_validate(
        CirculatorIterData iter, 
        inout VtxPositionValidationContext ctx
    ){
        uint wi = iter.awi.wedge_id; 
        uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
        uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 
        
        vec3 vpos_i = ld_vpos(vi); 
        vec3 vnor_i = ld_vnor_filtered(vi); 
        
        if (iter.rotate_step == 0u)
        { /* init position of vp */
            uint ivert_vp = mark__ve_circ_fwd__get_vp(iter); 
            uint vp = ssbo_edge_to_vert_[wi*4u + ivert_vp]; 
            ctx.vpos_p = ld_vpos(vp); 
            ctx.vnor_p = ld_vnor_filtered(vp); 
        }

        /* Validate filtered position */
        /* Foldover avoidance, see Ch.3.7 in https://www.cs.cmu.edu/~garland/thesis/thesis-onscreen.pdf */
        vec3 face_normal = normalize(cross(vpos_i-ctx.vpos, ctx.vpos_p-ctx.vpos)); // cross dot in CCW order
        vec3 clip_plane_normal = normalize(cross(face_normal, ctx.vpos_p - vpos_i)); // cross dot in CCW order
        if (dot(ctx.vpos_filtered - ctx.vpos_p, clip_plane_normal) < .0f)
            ctx.is_move_ok = false;

        float proj_dist = abs(dot(ctx.vpos_filtered - ctx.vpos, face_normal)); 
        if (proj_dist < ctx.min_dist_proj)
        {
            ctx.min_dist_proj = proj_dist;
            ctx.vpos_filtered_proj = ctx.vpos_filtered - dot(ctx.vpos_filtered - ctx.vpos, face_normal) * face_normal;
        }

        ctx.vpos_p = vpos_i;
        ctx.vnor_p = vnor_i; 
        return ctx.is_move_ok; 
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
    VertWedgeListHeader vwlh_orig = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    
    // Filtering in 1-ring
    VtxPositionFilteringContext ctx = init_vpos_filter_context(vpos, vnor);
    VertWedgeListHeader vwlh = vwlh_orig; 
    bool rot_fwd = true; 
    VE_CIRCULATOR(vwlh, ve_circulator__vpos_filter, ctx, rot_fwd); 
    ctx.vpos_filtered /= ctx.weight_sum;
    ctx.ave_edge_len /= float(ctx.num_edges);
    #if defined(TEST_SQRT_3_INTERP)
        float alpha_n = sqrt_3_subdiv_param(float(ctx.num_edges)); 
        ctx.vpos_filtered = (1 - alpha_n) * ctx.vpos + alpha_n * ctx.vpos_filtered;
    #endif

    VtxPositionValidationContext ctx_validate = init_filtered_vpos_validate_context(vpos, ctx.vpos_filtered, vnor);
    vwlh = vwlh_orig; 
    rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vpos_validate, ctx_validate, rot_fwd); 
    vec3 vpos_filtered = ctx_validate.vpos_filtered_proj; 

    if (ctx.border || ctx.sel_border) /* keep feature verts */
        ctx_validate.is_move_ok = false; 
    // if (length(vpos_filtered - vpos) > ctx.ave_edge_len * 0.5f)
    //     ctx_validate.is_move_ok = false; 
    
    vec3 vpos_new = vpos;  
    if (ctx_validate.is_move_ok)
        vpos_new = vpos_filtered;

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