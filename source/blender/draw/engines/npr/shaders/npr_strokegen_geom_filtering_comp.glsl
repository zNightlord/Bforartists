#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)


/* Refs: 
 * "Bilateral Normal Filtering for Mesh Denoising"
*/
float vnor_bilateral_filter_weight(vec3 vnor_i, vec3 vpos_i, vec3 vnor_j, vec3 vpos_j)
{
    const float tau_vnor = .2f; 
    float w_vnor = gaussian(distance(vnor_i, vnor_j), tau_vnor); 

    // note: positional weight needs to be added
    // const float tau_vpos = ?;
    float w_vpos = 1.0f; 

    return w_vnor * w_vpos; 
}

float sqrt3_vpos_smooth_weight(float n)
{
    return (4.0f - 2.0f * cos((2.0f * 3.14159265359f) / n)) / 9.0f;
}

float sqrt3_vpos_limit_weight(float n, int step, bool step_inf = false)
{
    float alpha_n = sqrt3_vpos_smooth_weight(n); 
    alpha_n *= 3.0f; 
    float beta_n_inf = alpha_n / (1.0f + alpha_n); 

    if (step_inf) return beta_n_inf; 

    float gamma_n_step = pow(2.0f/3.0f - alpha_n, step); 
    return gamma_n_step; 
}

vec3 sqrt3_limit_vpos(vec3 vpos, vec3 vpos_sum, float n, int step, bool step_inf)
{
    float w = sqrt3_vpos_limit_weight(n, step, step_inf); 
    vec3 lim_pos_inf = (1.0f - w) * vpos + w * (vpos_sum / n); 
    if (step_inf) return lim_pos_inf; 
    return w * vpos + (1.0f - w) * lim_pos_inf;
}

float loop_vpos_smooth_weight(float n)
{ /* mask for interior verts in loop subdiv */
    /* see "Subdivision methods for geometric design - Ch.7.3 Smooth Subdivision for Triangle Meshes" */
    
    if (n == 6.0f) return 3.0f / 8.0f; // fast path, actually same as the formula below

    float alpha_n = (3.0f/8.0f) + (cos((2.0f * 3.14159265359f) / n) / 4.0f);
    alpha_n = (5.0f/8.0f) - alpha_n * alpha_n;
    return alpha_n; 
}


/* Inputs: 
ssbo_vbo_full_[]
ssbo_vert_to_edge_list_header_[] 
pcs_vert_count_
ssbo_dyn_mesh_counters_out_
ssbo_bnpr_mesh_pool_counters_
*/

// TODO: Store data for selected verts in a compacted manner. 

uint get_vert_count()
{
    return pcs_vert_count_ + ssbo_dyn_mesh_counters_out_.num_verts; 
}

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VNOR_FILTERING)
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

    struct VtxNormalFilteringContext
    {
        vec3 vpos; 
        vec3 vnor;
        vec3 vnor_filtered; 
        float weight_sum; 
        
        bool border; 
        bool close_to_sel_border; 
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
        ctx.close_to_sel_border = false; 

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
        if ((!ef.selected) || ef.sel_border) ctx.close_to_sel_border = true; 

        return true; 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING)
    vec3 ld_vpos_filtered(uint vert_id)
    {
        if (pcs_vpos_filtering_iter_ % 2u == 0u)
            return ld_vpos(vert_id); 
        
        uvec3 vpos_filtered_enc;
        Load3(ssbo_vpos_temp_, vert_id, vpos_filtered_enc);
        return uintBitsToFloat(vpos_filtered_enc); 
    }
    void st_vpos_filtered(uint vert_id, vec3 vpos_filtered)
    {
        if (pcs_vpos_filtering_iter_ % 2u == 0u)
        {
            uvec3 vpos_filtered_enc = floatBitsToUint(vpos_filtered); 
            Store3(ssbo_vpos_temp_, vert_id, vpos_filtered_enc);
        }

        st_vpos(vert_id, vpos_filtered); 
    }

    struct VtxPositionFilteringContext
    {
        vec3 vpos; 
        float remesh_len;  

        vec3 vpos_p; 
        float remesh_len_p; 

        vec3 vpos_filtered;
        float weight_sum; 
        
        bool border; 
        bool close_to_sel_border; 
    };

    VtxPositionFilteringContext init_vpos_filter_context(vec3 vpos, float remesh_len)
    {
        VtxPositionFilteringContext ctx; 
        
        ctx.vpos = vpos;
        ctx.remesh_len = remesh_len; 
        ctx.vpos_p = vec3(.0f); /* init in circulator */
        ctx.remesh_len_p = .0f; /* init in circulator */
        ctx.vpos_filtered = vec3(.0f); 
        ctx.weight_sum = .0f; 
        
        ctx.border = false;
        ctx.close_to_sel_border = false; 

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
        vec3 vpos_i = ld_vpos_filtered(vi);
        float remesh_len_i = ld_vtx_remesh_len(vi); 

        if (iter.rotate_step == 0u)
        { // init position of vp
            uint ivert_vp = mark__ve_circ_fwd__get_vp(iter); 
            uint vp = ssbo_edge_to_vert_[wi*4u + ivert_vp]; 
            ctx.vpos_p = ld_vpos_filtered(vp); 
            ctx.remesh_len_p = ld_vtx_remesh_len(vp);
        }

        // TODO: cache this in VertFlags
        EdgeFlags ef = load_edge_flags(wi);
        if (ef.border) ctx.border = true;
        if ((!ef.selected) || ef.sel_border) ctx.close_to_sel_border = true; 

        // float tri_area = abs(length(cross(vpos_i - ctx.vpos, ctx.vpos_p - ctx.vpos)) * 0.5f);
        // float tri_remesh_len = ((remesh_len_i + ctx.remesh_len_p + ctx.remesh_len) / 3.0f);
        // float w = tri_area * tri_remesh_len; 
        // ctx.vpos_filtered += w * ((vpos_i + ctx.vpos_p + ctx.vpos) / 3.0f); 

        float w = 1.0f; 
        ctx.vpos_filtered += w * vpos_i; 


        ctx.weight_sum += w; 

        ctx.vpos_p = vpos_i; 
        ctx.remesh_len_p = remesh_len_i; 
        return true; 
    }
#endif


#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VQUADRICS_DIFFUSION) || defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING)
    bool valid_quadric(Quadric q)
    {
        return (q.area > 0.0f); // invalid quadrics are marked as -1.0f
    }

    void store_vert_quadric(uint vtx_id, Quadric q)
    { 
        vec4 data_0, data_1; 
        vec3 data_2; 
        encode_quadric(q, /*out*/data_0, data_1, data_2); 

        uint st_addr = vtx_id * 11u; 
        for (uint i = 0; i < 4; ++i)
            ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_0[i]);
        st_addr += 4u; 
        for (uint i = 0; i < 4; ++i)
            ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_1[i]);
        st_addr += 4u;
        for (uint i = 0; i < 3; ++i)
            ssbo_vert_quadric_data_out_[st_addr + i] = floatBitsToUint(data_2[i]);
    }

    Quadric load_vert_quadric(uint vtx_id)
    { // TODO: memory perf issue might be here, should opti for coalesced loads
        Quadric q; 

        uint ld_addr = vtx_id * 11u; 
        vec4 data_0, data_1;
        vec3 data_2;
        for (uint i = 0; i < 4; ++i)
            data_0[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);
        ld_addr += 4u;
        for (uint i = 0; i < 4; ++i)
            data_1[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);
        ld_addr += 4u;
        for (uint i = 0; i < 3; ++i)
            data_2[i] = uintBitsToFloat(ssbo_vert_quadric_data_in_[ld_addr + i]);

        q = decode_quadric(data_0, data_1, data_2); 

        return q; 
    }

    struct VtxQuadricsFilteringContext
    {
        vec3 vpos; 
        vec3 vpos_p; 
        Quadric q;  
        Quadric q_filtered;  
        float weight_sum; 
        
        bool border; 
        bool close_to_sel_border; 
    }; 

    VtxQuadricsFilteringContext init_vq_filter_context(
        vec3 vpos, Quadric q 
        // 0th iter: q = 0, 
        // after:    q = load_quadrics(vtx_id), 
    ){
        VtxQuadricsFilteringContext ctx; 
        
        ctx.vpos = vpos;
        ctx.vpos_p = vec3(.0f); 
        ctx.q = q;
        ctx.q_filtered = q;  
        ctx.weight_sum = (pcs_vq_filtering_iter_ == 0) ? .0f : 1.0f;
        
        ctx.border = false;
        ctx.close_to_sel_border = false; 

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
    bool ve_circulator__vq_filter(
        CirculatorIterData iter, 
        inout VtxQuadricsFilteringContext ctx
    ){
        uint wi = iter.awi.wedge_id; 
        uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
        uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 
        vec3 vpos_i = ld_vpos(vi);

        if (iter.rotate_step == 0u)
        { // init position of vp
            uint ivert_vp = mark__ve_circ_fwd__get_vp(iter); 
            uint vp = ssbo_edge_to_vert_[wi*4u + ivert_vp]; 
            ctx.vpos_p = ld_vpos(vp); 
        }

        // TODO: cache this in VertFlags
        EdgeFlags ef = load_edge_flags(wi);
        if (ef.border) ctx.border = true;
        if ((!ef.selected) || ef.sel_border) ctx.close_to_sel_border = true; 

        Quadric vq_i; 
        if (pcs_vq_filtering_iter_ == 0)
        {
            vec3 vci = vpos_i     - ctx.vpos; 
            vec3 vcp = ctx.vpos_p - ctx.vpos; 
            vec3 tri_cross = cross(vci, vcp); 
            vec3 face_normal = normalize(tri_cross); 
            vec3 face_barycenter = (ctx.vpos + vpos_i + ctx.vpos_p) / 3.0f; 
            
            mat4 view_to_world = ubo_view_matrices_.viewinv;
	        vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
            
            vq_i.quadric = compute_filter_quadric(face_normal, face_barycenter, cam_pos_ws, .05f, pcs_position_regularization_scale_); 
            vq_i.area = length(tri_cross) * 0.5f;

            ctx.q_filtered.quadric += vq_i.quadric;
            ctx.q_filtered.area += vq_i.area;  
            ctx.weight_sum += 1.0f; 
        }else
        {
            vq_i = load_vert_quadric(vi);
            if (valid_quadric(vq_i))
            {
                float dev_q = pcs_quadric_deviation_; 
                float dev_g = pcs_geodist_deviation_; 
                float w = compute_vert_quadric_weight(ctx.vpos, ctx.q, vpos_i, vq_i, dev_q, dev_g); 
                
                ctx.q_filtered.quadric += (w * vq_i.quadric);
                ctx.weight_sum += w; 
            }
        }

        ctx.vpos_p = vpos_i; 
        return true; 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VQUADRICS_DIFFUSION) || defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING) || defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING)
    struct VtxPositionValidationContext
    {
        vec3 vpos;
        vec3 vpos_p; 
        vec3 vpos_filtered; 

        float min_dist_proj; 
        vec3 vpos_filtered_proj; // proj to 1-ring neighbors

        bool border; 
        bool close_to_sel_border; 

        bool is_move_ok; 
    }; 
    VtxPositionValidationContext init_filtered_vpos_validate_context(
        vec3 vpos, vec3 vpos_filtered
    ){
        VtxPositionValidationContext ctx; 
        
        ctx.vpos = vpos;
        ctx.vpos_p = vec3(.0f); // init in circulator
        ctx.vpos_filtered = ctx.vpos_filtered_proj = vpos_filtered; 

        ctx.border = false;
        ctx.close_to_sel_border = false; 

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
        
        // Detect locked verts
        EdgeFlags ef = load_edge_flags(wi); 
        if (ef.border) ctx.border = true;
        if ((!ef.selected) || ef.sel_border) ctx.close_to_sel_border = true;

        vec3 vpos_i = ld_vpos(vi); 
        
        if (iter.rotate_step == 0u)
        { /* init position of vp */
            uint ivert_vp = mark__ve_circ_fwd__get_vp(iter); 
            uint vp = ssbo_edge_to_vert_[wi*4u + ivert_vp]; 
            ctx.vpos_p = ld_vpos(vp); 
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
        return ctx.is_move_ok; 
    }
#endif



#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING) || defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING_FINISH)
#define SUBDIV_VPOS_FILTER_TYPE__SQRT3 0u
#define SUBDIV_VPOS_FILTER_TYPE__LOOP 1u 
struct VtxPositionSmoothForSubdivContext
{
    vec3 vpos_filtered; 
    float num_adj_edges; 

    bool border; 
    bool close_to_sel_border; 
}; 
VtxPositionSmoothForSubdivContext init_vpos_subd_smooth_context()
{
    VtxPositionSmoothForSubdivContext ctx; 
    ctx.vpos_filtered = vec3(.0f); 
    ctx.num_adj_edges = .0f; 
    ctx.border = false;
    ctx.close_to_sel_border = false;

    return ctx; 
}
bool ve_circulator__vpos_smooth_for_subdiv(
    CirculatorIterData iter, 
    inout VtxPositionSmoothForSubdivContext ctx
){
    uint wi = iter.awi.wedge_id; 

    // TODO: cache this in VertFlags
    EdgeFlags ef = load_edge_flags(wi);
    if (ef.border) ctx.border = true;
    if ((!ef.selected) || ef.sel_border) ctx.close_to_sel_border = true; 

    uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
    uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 
    vec3 vpos_i = ld_vpos(vi); 
    
    ctx.vpos_filtered += vpos_i; 
    ctx.num_adj_edges += 1.0f; 
    
    return true; 
}
bool is_old_vert_sqrt3(VertFlags vf)
{ 
    bool is_old_vtx = 
        (!vf.new_by_face_split)
        && (!vf.del_by_collapse)
        && (!vf.dupli); 
    return is_old_vtx; 
}
bool is_old_vert_loop(VertFlags vf)
{ 
    bool is_old_vtx = 
        (!vf.new_by_split)
        && (!vf.del_by_collapse)
        && (!vf.dupli); 
    return is_old_vtx; 
}
#endif







#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING)
float ld_vcurv_max_smoothed(uint vert_id)
{
    if (pcs_vcurv_smooth_iter_ % 2u == 0u)
        return ld_vcurv_max(vert_id); 
    
    return uintBitsToFloat(ssbo_vcurv_max_temp_[vert_id]); 
}
void st_vcurv_max_smoothed(uint vert_id, float vcurv_smoothed)
{
    if (pcs_vcurv_smooth_iter_ % 2u == 0u)
        ssbo_vcurv_max_temp_[vert_id] = floatBitsToUint(vcurv_smoothed);

    st_vcurv_max(vert_id, vcurv_smoothed); 
}

struct VtxCurvatureSmoothingContext
{
    float vcurv_summed; 
    float vcurv_weight_sum; 

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    vec3 vpos; 
    float ave_edge_len; 
#endif

    float num_adj_edges; 

}; 
VtxCurvatureSmoothingContext init_vcurv_smooth_context(
#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    vec3 vpos
#endif
)
{
    VtxCurvatureSmoothingContext ctx; 
    ctx.vcurv_summed = .0f; 
    ctx.vcurv_weight_sum = .0f; 

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    ctx.vpos = vpos; 
    ctx.ave_edge_len = .0f; 
#endif

    ctx.num_adj_edges = .0f; 

    return ctx; 
}
bool ve_circulator__vcurv_smooth(
    CirculatorIterData iter, 
    inout VtxCurvatureSmoothingContext ctx
){
    uint wi = iter.awi.wedge_id; 

    uint ivert_vi = mark__ve_circ_fwd__get_vi(iter); 
    uint vi = ssbo_edge_to_vert_[wi*4u + ivert_vi]; 
    float vcurv_i = ld_vcurv_max_smoothed(vi); 
    
    if (valid_vcurv_max(vcurv_i))
    {
        ctx.vcurv_summed += vcurv_i; 
        ctx.vcurv_weight_sum += 1.0f; 
    }

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    vec3 vpos_i = ld_vpos(vi); 
    float edge_len = length(vpos_i - ctx.vpos); 
    ctx.ave_edge_len += edge_len; 
#endif


    ctx.num_adj_edges += 1.0f; 
    return true; 
}
#endif


#if !defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__EDGES)
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
    if (ctx.border || ctx.close_to_sel_border) 
        vnor_filtered = vnor;

    if (valid_thread)
        st_vnor_filtered(vert_id, vnor_filtered);
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING)
    if(!valid_thread) return;

    VertFlags vf = load_vert_flags(vert_id);
    if (pcs_vpos_filtering_iter_ == pcs_num_vpos_filtering_iters_ - 1)
        update_vert_flags__reset_edge_split_new_vert__mov_by_collapse(vert_id, vf);

    vec3 vnor = ld_vnor(vert_id); 
    vec3 vpos_old = ld_vpos(vert_id);
    float remesh_len = ld_vtx_remesh_len(vert_id); 
    vec3 vpos_new = vpos_old; 

    if (pcs_vpos_filtering_iter_ < pcs_num_vpos_filtering_iters_)
    {
        /* Tangential Smoothing, see "Adaptive Remeshing for Real-Time Mesh Deformation" */
        VtxPositionFilteringContext ctx = init_vpos_filter_context(vpos_old, remesh_len);
        VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
        bool rot_fwd = true;
        VE_CIRCULATOR(vwlh, ve_circulator__vpos_filter, ctx, rot_fwd);

        if ((vf.new_by_split || vf.mov_by_collapse))
        {
            /* Putting verts at infinite subdiv level */
            vec3 vpos_subd = sqrt3_limit_vpos(
                vpos_old, ctx.vpos_filtered, ctx.weight_sum, 1, true
            );
            vpos_new = mix(vpos_old, vpos_subd, .1f); 
        }else{
            float alpha_n = sqrt3_vpos_smooth_weight(ctx.weight_sum); 
            vpos_new = (1.0f - alpha_n) * vpos_old + alpha_n * ctx.vpos_filtered / ctx.weight_sum;
            vpos_new = mix(vpos_old, vpos_new, .1f);
        }

        /* Tangential smoothing, weighted by remesh len & tri area */
        // ctx.vpos_filtered /= ctx.weight_sum;
        // vec3 vpos_new = ctx.vpos_filtered - vnor * dot(vnor, ctx.vpos_filtered - vpos_old); 

        /* Avoid fold-over, project to 1-ring fan */
        // VtxPositionValidationContext ctx_validate = init_filtered_vpos_validate_context(vpos_old, vpos_new);
        // rot_fwd = true;
        // VE_CIRCULATOR(vwlh, ve_circulator__vpos_validate, ctx_validate, rot_fwd); 
        // vpos_new = ctx_validate.vpos_filtered_proj; 

        if (ctx.border || ctx.close_to_sel_border) 
            vpos_new = vpos_old;
    }


    if (valid_thread)
        st_vpos_filtered(vert_id, vpos_new);
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VQUADRICS_DIFFUSION)
    if(!valid_thread) return;

    // Load vert quadric 
    Quadric vq;
    if (pcs_vq_filtering_iter_ == 0)
    {   
        vq.quadric = mat4(.0f);
        vq.area = .0f;
    }else
    {
        vq = load_vert_quadric(vert_id); 
    } 
    vec3 vpos = ld_vpos(vert_id);

    // Diffusion
    VtxQuadricsFilteringContext ctx = init_vq_filter_context(vpos, vq);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vq_filter, ctx, rot_fwd);
    if (ctx.weight_sum > .0f)
        ctx.q_filtered.quadric = ctx.q_filtered.quadric / ctx.weight_sum;

    // Fallback on border verts
    if ((ctx.border || ctx.close_to_sel_border) && (pcs_vq_filtering_iter_ > 0))  
    {
        ctx.q_filtered = vq; 
        ctx.q_filtered.area = -1.0f;
    }
    
    if (valid_thread)
        store_vert_quadric(vert_id, ctx.q_filtered); 
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING)
    if(!valid_thread) return;

    /* vert pos from last filtering iter */
    vec3 vpos = ld_vpos(vert_id);
    vec3 vpos_orig = vpos; 
    
    VertFlags vf = load_vert_flags(vert_id);

    /* Project vert pos to quadric surface */
    Quadric q_v = load_vert_quadric(vert_id); 
    bool is_qv_valid = valid_quadric(q_v); 
    vec3 vnor = ld_vnor(vert_id);  

    #if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING__CONSTRAINED_SOLVE)
        mat3 A = mat3(q_v.quadric); /* slice upper 3x3 */
        vec3 b = vec3(q_v.quadric[3][0], q_v.quadric[3][1], q_v.quadric[3][2]); 
        float lambda = -(dot(vpos, A * vnor) + dot(b, vnor)) / dot(vnor, A * vnor); 
        vpos += vnor * lambda; 
    #else /* Directly solve by minimizing the quadrics */
        mat4 q_d = mat4(q_v.quadric); 
 
        mat4 q_reg = compute_plane_quadric(vnor, vpos_orig); 
        q_d = q_d + pcs_position_regularization_scale_ * q_reg; /*damping*/

        dmat3 A = dmat3(q_d); 
        dvec3 b = dvec3(q_d[3][0], q_d[3][1], q_d[3][2]); 
        vpos = vec3(-inverse(A)*b); 
    #endif

    VtxPositionValidationContext ctx_validate = init_filtered_vpos_validate_context(vpos_orig, vpos);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vpos_validate, ctx_validate, rot_fwd); 
    vpos = ctx_validate.vpos_filtered_proj; 

    if ((!is_qv_valid) || any(isnan(vpos)) || any(isinf(vpos)) || 
            ctx_validate.close_to_sel_border || ctx_validate.border ||  (!(vf.new_by_split || vf.mov_by_collapse))) 
        vpos = vpos_orig; 

    if (valid_thread)
        st_vpos(vert_id, vpos);
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING)
    if(!valid_thread) return;

    vec3 vpos = ld_vpos(vert_id);
    vec3 vpos_orig = vpos;

    VtxPositionSmoothForSubdivContext ctx = init_vpos_subd_smooth_context();
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vpos_smooth_for_subdiv, ctx, rot_fwd); 
    ctx.vpos_filtered = ctx.vpos_filtered / ctx.num_adj_edges; 
    
    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__SQRT3)
    {
        float alpha_n = sqrt3_vpos_smooth_weight(ctx.num_adj_edges); 
        vpos = (1.0f - alpha_n) * vpos_orig + alpha_n * ctx.vpos_filtered;
    }
    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__LOOP)
    {
        float alpha_n = loop_vpos_smooth_weight(ctx.num_adj_edges); 
        vpos = (1.0f - alpha_n) * vpos_orig + alpha_n * ctx.vpos_filtered; 
    }

    if (ctx.close_to_sel_border || ctx.border) 
        vpos = vpos_orig; 

    uvec3 vpos_enc = floatBitsToUint(vpos); 
    Store3(ssbo_vpos_subd_, sel_vert_id, vpos_enc); 
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING_FINISH)
    if (!valid_thread) return; 
    
    VertFlags vf = load_vert_flags(vert_id);
    
    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__SQRT3)
        update_vert_flags__reset_face_split_new_vert(vert_id, vf);
    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__SQRT3 && !is_old_vert_sqrt3(vf)) return; 

    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__LOOP)
        update_vert_flags__reset_edge_split_tags(vert_id, vf); 
    if (pcs_subdiv_type_ == SUBDIV_VPOS_FILTER_TYPE__LOOP && !is_old_vert_loop(vf)) return; 

    
    uvec3 vpos_enc;
    Load3(ssbo_vpos_subd_, sel_vert_id, vpos_enc); 
    vec3 vpos = uintBitsToFloat(vpos_enc); 
    
    st_vpos(vert_id, vpos);
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING)
    vert_id = gl_GlobalInvocationID.x; 
    uint num_verts = get_vert_count(); 
    valid_thread = vert_id < num_verts; 

    if(!valid_thread) return;

    float vcurv = ld_vcurv_max(vert_id); 
#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    vec3 vpos = ld_vpos(vert_id); 
#endif

    VtxCurvatureSmoothingContext ctx = init_vcurv_smooth_context(
#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
        vpos
#endif
    );
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 
    bool rot_fwd = true;
    VE_CIRCULATOR(vwlh, ve_circulator__vcurv_smooth, ctx, rot_fwd); 
    
    float vurv_new = vcurv; 
    { /* Calculate new curvature */ 
        bool valid_vcurv = valid_vcurv_max(vcurv); 
        bool valid_vcurv_sum = ctx.vcurv_weight_sum > .0f;  
        if (valid_vcurv_sum)
        { 
            ctx.vcurv_summed = ctx.vcurv_summed / ctx.vcurv_weight_sum;
            if (valid_vcurv)
                vurv_new = mix(vcurv, ctx.vcurv_summed, .1f); 
            else
                vurv_new = ctx.vcurv_summed; 
        }
    }

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
    ctx.ave_edge_len = ctx.ave_edge_len / ctx.num_adj_edges; 
#endif


    if (valid_thread)
    {
        st_vcurv_max_smoothed(vert_id, vurv_new);  
#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN)
        float edge_len = get_adaptive_remesh_len(vurv_new, ctx.ave_edge_len); 
        if (!valid_vcurv_max(vurv_new)) edge_len = ctx.ave_edge_len;
        st_vtx_remesh_len(vert_id, edge_len);
#endif
    }
#endif
}
#endif

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__EDGES)
/* Filter for edges */
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
    uint sel_edge_id = gl_GlobalInvocationID.x; 
    uint wedge_id; bool valid_thread; 
    get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge_id, /*out*/valid_thread); 

#if defined(_KERNEL_MULTICOMPILE__SURF_FILTERING__LOOP_SUBD_EDGE_POINTS)
    if (!valid_thread) return; 

    // EdgeFlags ef = load_edge_flags(wedge_id); 
    // update_edge_flags__reset_new_split_edge(wedge_id, ef); 
    // not needed, will be done at first split operator

    /* Compute subdiv position for new points on each edge */
    uint v[4] = {
        ssbo_edge_to_vert_[wedge_id*4u + 0u], 
        ssbo_edge_to_vert_[wedge_id*4u + 1u], 
        ssbo_edge_to_vert_[wedge_id*4u + 2u], 
        ssbo_edge_to_vert_[wedge_id*4u + 3u]
    };

    vec3 vpos[4] = {
        ld_vpos(v[0]), 
        ld_vpos(v[1]), 
        ld_vpos(v[2]), 
        ld_vpos(v[3])
    };

    vec3 vpos_new = .125f * vpos[0] + .125f * vpos[2] + .375f * vpos[1] + .375f * vpos[3]; 
    uvec3 vpos_new_enc = floatBitsToUint(vpos_new); 
    Store3(ssbo_epos_subd_, wedge_id, vpos_new_enc);
#endif

}
#endif