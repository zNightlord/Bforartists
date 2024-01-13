#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)



/*
.define("WINGED_EDGE_TOPO_INCLUDE", "1")
.define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
.define("VE_CIRCULATOR_INCLUDE", "1")
.define("DYNAMESH_SELECTION_INDEXING_COMMON", "1")

float ssbo_vbo_full_[]
float ssbo_vnor_[]
float ssbo_varea_[]
int pcs_vert_count_
SSBOData_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_out_
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

vec3 ld_vnor(uint vtx_id)
{
	vec3 vnor; 
    Load3(ssbo_vnor_, vtx_id, vnor);
    return vnor; 
}
void st_vnor(uint vtx_id, vec3 vnor)
{
    Store3(ssbo_vnor_, vtx_id, vnor);
}

float ld_varea(uint vtx_id)
{
    float varea; 
    return ssbo_varea_[vtx_id]; 
}
void st_varea(uint vtx_id, float varea)
{
    ssbo_varea_[vtx_id] = varea; 
}

uint get_vert_count()
{
    return pcs_vert_count_ + ssbo_dyn_mesh_counters_out_.num_verts; 
}


float calc_voronoi_area(
    vec3 edge_len_sqr, float cot_x, 
    float ang_c, bool is_obtuse, float face_area
) 
{
    if (is_obtuse) { /* non-obtuse triangle */
        if (ang_c > PI_HALF) 
            return .5f * face_area; 
        else 
            return .25f * face_area; 
    }
    /* obtuse, voronoi safe */
    return .125f * edge_len_sqr * cot_x;
}




struct CalcVertAttrContext_Order0
{
    uint vi; 
    vec3 vpos; 

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    vec3 sum_normal;
    vec3 sum_weight;  
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    float voronoi_area; 
#endif
}; 

CalcVertAttrContext_Order0 init_vert_attr_context_order_0(
    uint vi, vec3 vpos
){
    CalcVertAttrContext_Order0 ctx; 
    ctx.vi = vi; 
    ctx.vpos = vpos; 

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    ctx.sum_normal = .0f; 
    ctx.sum_weight = .0f;
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    ctx.voronoi_area = .0f;
#endif
}


#define PI_HALF 1.57079632679f
/*    Rotate fwd around V3(marked as vc)        
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
bool calc_vert_normal_by_angle(
    CirculatorIterData iter, 
    inout CalcVertAttrContext_Order0 ctx
){
    uint wi = iter.awi.wedge_id; 
    uint vi = ctx.vi; 
    uint ivert_vn = mark__ve_circ_fwd__get_vn(iter); 
    uint vn = ssbo_edge_to_vert_[wi*4u + ivert_vn]; 
    
    vec3 vpos_i = ld_vpos(vi); 
    vec3 vpos_n = ld_vpos(vn); 

    vec3 vci = vpos_i - ctx.vpos; 
    float len_ci = length(vci);
    vec3 vci_dir = vci / len_ci; 
#define vic ((-vci))
#define vic_dir ((-vci_dir))

    vec3 vcn = vpos_n - ctx.vpos;
    float len_cn = length(vcn);
    vec3 vcn_dir = vcn / len_cn; 
#define vnc ((-vcn))
#define vnc_dir ((-vcn_dir))

    vec3 vin = vpos_n - vpos_i; 
    float len_in = length(vin);
    vec3 vin_dir = vin / len_in;
#define vni ((-vin))
#define vni_dir ((-vin_dir))

    float dot_ci_cn = dot(vci_dir, vcn_dir);
    float ang_c = acos(dot_ci_cn);

    float dot_ic_in = dot(vic_dir, vin_dir);
    float ang_i = acos(dot_ic_in);

    float dot_nc_ni = dot(vnc_dir, vni_dir);
    float ang_n = acos(dot_nc_ni);

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    if (len_ci != .0f && len_cn != .0f) 
    { // compute angle-weighted vertex normal
        float dp = dot(vci_dir, vcn_dir); 
        dp = clamp(dp, -1.0f, 1.0f); 

        float angle = acos(dp);
        ctx.sum_weight += angle; 
        ctx.sum_normal += angle * cross(vci_dir, vcn_dir);
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    { // compute voronoi area of the current verts
      // See Section 3.4 in "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds" by Mark Meyer et al. 
      // Here we split the sum of each edge into two parts. 
        bool is_obtuse = PI_HALF < max(max(ang_c, ang_i), ang_n);
        float face_area = .5f * length(cross(vci, vcn)); 
        
        float vci_len_sqr = dot(vci, vci); 
        float cot_n = cot(ang_n); 
        ctx.voronoi_area += calc_voronoi_area(
            vci_len_sqr, cot_n, 
            ang_c,  is_obtuse, face_area
        );  
        
        float vcn_len_sqr = dot(vcn, vcn); 
        float cot_i = cot(ang_i); 
        ctx.voronoi_area += calc_voronoi_area(
            vcn_len_sqr, cot_i, 
            ang_c,  is_obtuse, face_area
        );
    }
#endif

    ctx.vi = vn; 
    return true; 

#undef vic
#undef vic_dir
#undef vnc
#undef vnc_dir
#undef vni
#undef vni_dir
}


#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0)
void main()
{
    uint vert_id = gl_GlobalInvocationID.x; 
    uint num_verts = get_vert_count(); 
    bool valid_thread = vert_id < num_verts; 

    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(vert_wedge_list_headers[vert_id]);

    CalcVertAttrContext_Order0 ctx; 
    ctx.vi = ssbo_edge_to_vert_[vwlh.wedge_id*4u + (vwlh.ivert == 1u) ? 3u : 1u]; 
    ctx.vpos = ld_vpos(vert_id); 
    ctx.sum_weight = ctx.sum_normal = .0f;
    if (valid_thread)
    {
        bool rotate_fwd = true; 
        VE_CIRCULATOR(vwlh, calc_vert_normal_by_angle, ctx, rotate_fwd);
    }
    
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    vec3 vnormal = normalize(ctx.sum_normal / ctx.sum_weight); 
    if (valid_thread)
        st_vnor(vert_id, vnormal);
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    if (valid_thread)
        st_varea(vert_id, ctx.voronoi_area); 
#endif


}
#endif





/* 
 * "Estimating Curvatures and Their Derivatives on Triangle Meshes"
 * Impl based on 
 * https://github.com/yunjay/Apparent-Ridges/tree/master
 * which is a GPU adaptation of Rusinkiewicz curvature estimator from 
 * http://graphics.zcu.cz/curvature.html */
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_CURVATURE)
void main()
{
}
#endif