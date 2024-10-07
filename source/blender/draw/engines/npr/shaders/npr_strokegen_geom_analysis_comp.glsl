
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_contour_geom_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_debug_view_lib.glsl)



/*
.define("WINGED_EDGE_TOPO_INCLUDE", "1")
.define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
.define("VE_CIRCULATOR_INCLUDE", "1")
.define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")
.define("INCLUDE_VERTEX_GEOM", "1")

int pcs_rsc_handle
ObjectMatrices drw_matrix_buf
float ssbo_vbo_full_[]
uint ssbo_vnor_[] // uint so that buffer can be arbitrary
uint ssbo_varea_[] // uint so that buffer can be arbitrary
int pcs_vert_count_
SSBOData_StrokeGenDynamicMeshCounters ssbo_dyn_mesh_counters_out_
*/
#define PI_HALF 1.57079632679f



uint get_vert_count()
{
    return pcs_vert_count_ + ssbo_dyn_mesh_counters_out_.num_verts; 
}






float calc_voronoi_area_corner(
    float edge_len_sqr, float cot_x, 
    float ang_c, bool is_obtuse, float face_area
) 
{
    if (is_obtuse) { /* non-obtuse triangle */
        if (ang_c > (PI_HALF * .99f)) 
            return .5f * face_area; 
        else 
            return .25f * face_area; 
    }
    /* obtuse, voronoi safe */
    return .125f * edge_len_sqr * cot_x;
}


/* 
* Calculate voronoi area of V at face (V, v0, v1)
* V---e0--v0  ang_v:=angle<v0-V-v1>
*  \ <--- /   ang_0:=angle<V-v0-v1>
*  e1    /    ang_1:=angle<V-v1-v0>
*    \  /   
*     v1  
*/ 
float calc_voronoi_area_face(vec3 e0, vec3 e1, float ang_v, float ang_0, float ang_1) 
{
    bool is_obtuse = PI_HALF < max(max(ang_v, ang_0), ang_1);
    float face_area = .5f * length(cross(e0, e1)); 
    
    float e0_len_sqr = dot(e0, e0); 
    float cot_1 = 1.0f/tan(ang_1); 
    float voronoi_area = calc_voronoi_area_corner(e0_len_sqr, cot_1, ang_v, is_obtuse, face_area);
    
    float e1_len_sqr = dot(e1, e1); 
    float cot_0 = 1.0f/tan(ang_0); 
    voronoi_area += calc_voronoi_area_corner(e1_len_sqr, cot_0, ang_v, is_obtuse, face_area);
    
    return voronoi_area;
}

/* 
* Calculate voronoi area of P1 at face (P1, p2, p3)
* P1 ------ p3  
*   \ <--- /   
*    \    /    
*     \  /   
*      p2  
*/ 
float calc_voronoi_area_face(vec3 p1, vec3 p2, vec3 p3)
{
    vec3 e31 = p1 - p3; 
    vec3 e12 = p2 - p1; 
    vec3 e23 = p3 - p2;
    float ang_312 = acos(clamp(dot(normalize(e31), -normalize(e12)),-1.0f, 1.0f));
    float ang_123 = acos(clamp(dot(normalize(e12), -normalize(e23)),-1.0f, 1.0f));
    float ang_231 = acos(clamp(dot(normalize(e23), -normalize(e31)),-1.0f, 1.0f));
    return calc_voronoi_area_face(-e31, e12, ang_312, ang_231, ang_123); 
}

struct CalcVertAttrContext_Order0
{
    uint vi; 
    vec3 vpos; 

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    vec3 sum_normal;
    float sum_weight;  
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    float voronoi_area; 
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER)
    bool is_border; 
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE)
    uint num_adj_crease_edges; 
#endif

}; 

CalcVertAttrContext_Order0 init_vert_attr_context_order_0(
    uint vi, vec3 vpos
){
    CalcVertAttrContext_Order0 ctx; 
    ctx.vi = vi; 
    ctx.vpos = vpos; 

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    ctx.sum_normal = vec3(.0f); 
    ctx.sum_weight = .0f;
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    ctx.voronoi_area = .0f;
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER)
    ctx.is_border = false; 
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE)
    ctx.num_adj_crease_edges = 0u; 
#endif

    return ctx; 
}


/*    Example: Rotate fwd(CW) around V3(marked as vc)        
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
bool calc_vert_attr_order_0(
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
        ctx.sum_normal += angle * cross(vcn_dir, vci_dir/*we are rotating arount CW, but surface oriented CCW*/);
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    { // compute voronoi area of the current verts
      // See Section 3.4 in "Discrete Differential-Geometry Operators for Triangulated 2-Manifolds" by Mark Meyer et al. 
      // Here we split the sum of each edge into two parts. 
        ctx.voronoi_area += calc_voronoi_area_face(vci, vcn, ang_c, ang_i, ang_n); 
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER)
    EdgeFlags ef = load_edge_flags(wi);
    ctx.is_border = ctx.is_border || ef.border;
#endif
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE)
    ctx.num_adj_crease_edges += (0u < ef.crease_level ? 1 : 0); 
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
    uint groupIdx = gl_LocalInvocationID.x; 
    uint num_verts = get_vert_count(); 
    
    uint vert_id; bool valid_thread; 
    if (0 < pcs_only_selected_verts_)
    {
        uint sel_vert_id = gl_GlobalInvocationID.x; 
        get_vert_id_from_selected_vert(sel_vert_id, /*out*/vert_id, valid_thread);
    }
    else
    {
        vert_id = gl_GlobalInvocationID.x; 
        valid_thread = vert_id < get_vert_count(); 
    }
    VertFlags vf = load_vert_flags(vert_id); 


    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]);

    vec3 vpos = ld_vpos(vert_id); 
    CalcVertAttrContext_Order0 ctx = init_vert_attr_context_order_0(
        ssbo_edge_to_vert_[vwlh.wedge_id*4u + ((vwlh.ivert == 1u) ? 3u : 1u)], 
        vpos
    ); 
    if (valid_thread)
    {
        bool rotate_fwd = true; 
        VE_CIRCULATOR(vwlh, calc_vert_attr_order_0, ctx, rotate_fwd);
    }
    
    VertFlags vf_new = vf; 
    bool update_vf = false; 
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER)
    update_vert_flags__border(ctx.is_border, /*inout*/vf_new);
    update_vf = true; 
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE)
    update_vert_flags__crease(2u == ctx.num_adj_crease_edges, /*inout*/vf_new);
    update_vert_flags__corner(2u < ctx.num_adj_crease_edges,  /*inout*/vf_new);
    update_vf = true;
#endif
    
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    vec3 vnormal = normalize(ctx.sum_normal / ctx.sum_weight); 
    if (valid_thread)
    {
        st_vnor(vert_id, vnormal);
	    
        if (0 < pcs_output_vertex_facing_flag_)
        {
            mat4 view_to_world = ubo_view_matrices_.viewinv;
            vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
            vec3 view_dir = vpos - cam_pos_ws; 
            float dot_vn = dot(view_dir, vnormal);
            update_vert_flags__facing_direction(dot_vn >= .0f, dot_vn < .0f, /*inout*/vf_new); 
            update_vf = true; 
        }
    }

    if (0 < pcs_output_dbg_geom_)
    {
        bool dbg_vtx_nor = (!vf.dupli) && valid_thread; 
        uint dbg_line_id = compact_normal_line(dbg_vtx_nor, groupIdx);
        dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 
        if (dbg_vtx_nor)
        {
            float dbg_line_len = pcs_dbg_geom_scale_; 

            vec4 vpos_ws_0 = vec4(ctx.vpos, 1.0f);
            vec4 vpos_ws_1 = vec4(ctx.vpos + vnormal * dbg_line_len * .05f, 1.0f);

            DebugVertData dvd_0 = DebugVertData(vpos_ws_0.xyz, vec3(1.0f), uvec4(0u)); 
            DebugVertData dvd_1 = DebugVertData(vpos_ws_1.xyz, vec3(1.0f), uvec4(0u)); 
            store_debug_line_data(dbg_line_id, dvd_0, dvd_1); 
            dbg_line_id++; 
        }
    }
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA)
    if (valid_thread)
        st_varea(vert_id, ctx.voronoi_area); 
#endif

    if (update_vf)
        store_vert_flags(vert_id, vf_new);
}
#endif





/* 
 * "Estimating Curvatures and Their Derivatives on Triangle Meshes"
 * Impl based on 
 * https://github.com/yunjay/Apparent-Ridges/tree/master
 * which is a GPU adaptation of Rusinkiewicz curvature estimator from 
 * http://graphics.zcu.cz/curvature.html */
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1)

void load_local_vertex_frame(VertWedgeListHeader vwlh, vec3 vpos, vec3 vnor, out vec3 pd1, out vec3 pd2)
{
	mat4 view_to_world = ubo_view_matrices_.viewinv;
	vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */

    // uint ivert_another_on_edge = vwlh.ivert == 1u ? 3u : 1u; 
    // uint another_vid = ssbo_edge_to_vert_[vwlh.wedge_id*4u + ivert_another_on_edge];
    // vec3 vpos_another = ld_vpos(another_vid);
 
    pd1 = cam_pos_ws/* vpos_another */ - vpos;
    pd1 = normalize(cross(pd1, vnor));

    pd2 = cross(vnor, pd1); 
}

void load_local_vertex_frame(uint vid, vec3 vpos, vec3 vnor, out vec3 pd1, out vec3 pd2)
{
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vid]);
    load_local_vertex_frame(vwlh, vpos, vnor, /*out*/pd1, pd2);  
}


/* Ported from https://github.com/Forceflow/trimesh2/blob/main/libsrc/TriMesh_curvature.cc */

// Rotate a coordinate system to be perpendicular to the given normal
void rot_coord_sys(vec3 old_u, vec3 old_v, vec3 new_norm,
                   out vec3 new_u, out vec3 new_v
){
	new_u = old_u;
	new_v = old_v;
	vec3 old_norm = cross(old_u, old_v);
	float ndot = dot(old_norm, new_norm);
	if (/* unlikely */(ndot <= -1.0f)) {
		new_u = -new_u;
		new_v = -new_v;
		return;
	}

	// Perpendicular to old_norm and in the plane of old_norm and new_norm
	vec3 perp_old = new_norm - ndot * old_norm;

	// Perpendicular to new_norm and in the plane of old_norm and new_norm
	// vec3 perp_new = ndot * new_norm - old_norm;

	// perp_old - perp_new, with normalization constants folded in
	vec3 dperp = (1.0f / (1 + ndot)) * (old_norm + new_norm);

	// Subtracts component along perp_old, and adds the same amount along
	// perp_new.  Leaves unchanged the component perpendicular to the
	// plane containing old_norm and new_norm.
	new_u -= dperp * dot(new_u, perp_old);
	new_v -= dperp * dot(new_v, perp_old);
}


// Reproject a curvature tensor from the basis spanned by old_u and old_v
// (which are assumed to be unit-length and perpendicular) to the
// new_u, new_v basis.
void proj_curv(vec3 old_u, vec3 old_v,
               float old_ku, float old_kuv, float old_kv,
               vec3 new_u, vec3 new_v,
               out float new_ku, out float new_kuv, out float new_kv)
{
	vec3 r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, cross(old_u, old_v), r_new_u, r_new_v);

	float u1 = dot(r_new_u, old_u);
	float v1 = dot(r_new_u, old_v);
	float u2 = dot(r_new_v, old_u);
	float v2 = dot(r_new_v, old_v);
	new_ku  = old_ku * u1*u1 + old_kuv * (2.0f  * u1*v1) + old_kv * v1*v1;
	new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
	new_kv  = old_ku * u2*u2 + old_kuv * (2.0f  * u2*v2) + old_kv * v2*v2;
}

// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
void diagonalize_curv(vec3 old_u, vec3 old_v,
                      float ku, float kuv, float kv,
                      vec3 new_norm, 
                      out vec3 pdir1, out vec3 pdir2, out float k1, out float k2)
{
	vec3 r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	float c = 1, s = 0, tt = 0;
    /* jwzw: if kuv == 0, I found it's mostly a planar vertex, 
     * just skip this branch and we will get min/max curvature of 0s */
	if ((kuv != 0.0f)) 
    { // Jacobi rotation to diagonalize
		float h = 0.5f * (kv - ku) / kuv;
		tt = (h < 0.0f) ?
			1.0f / (h - sqrt(1.0f + h*h)) :
			1.0f / (h + sqrt(1.0f + h*h));
		c = 1.0f / sqrt(1.0f + tt*tt);
		s = tt * c;
	}

	k1 = ku - tt * kuv;
	k2 = kv + tt * kuv;

	if (abs(k1) >= abs(k2)) {
		pdir1 = c*r_old_u - s*r_old_v;
	} else {
		// swap(k1, k2);
        float temp = k1; 
        k1 = k2;
        k2 = temp;
		pdir1 = s*r_old_u + c*r_old_v;
	}
	pdir2 = cross(new_norm, pdir1);
}

// LDL^T decomposition of a symmetric positive definite matrix (and some
// other symmetric matrices, but fragile since we don't do pivoting).
// Like Cholesky, but no square roots, which is important for small N.
// Reads diagonal and upper triangle of matrix A.
// On output, lower triangle of A holds LD, while rdiag holds D^-1.
// Algorithm from Golub and van Loan.
#define T float
#define N 3 // only works for 3x3 matrix
T sqr(T x)
{
	return x * x;
}
bool ldltdc(inout T A[N][N], inout T rdiag[N])
{
	// Special case for small N
	{
		T d0 = A[0][0];
		rdiag[0] = 1 / d0;
		
        A[1][0] = A[0][1];
		T l10 = rdiag[0] * A[1][0];
		
        T d1 = A[1][1] - l10 * A[1][0];
		rdiag[1] = 1 / d1;
		
        T d2 = A[2][2] - rdiag[0] * sqr(A[2][0]) - rdiag[1] * sqr(A[2][1]);
		rdiag[2] = 1 / d2;
		
        A[2][0] = A[0][2];
		A[2][1] = A[1][2] - l10 * A[2][0];
		return (d0 != 0 && d1 != 0 && d2 != 0);
	}

	T v[N-1];
	for (uint i = 0; i < N; i++) {
		for (uint k = 0; k < i; k++)
			v[k] = A[i][k] * rdiag[k];
		for (uint j = i; j < N; j++) {
			T sum = A[i][j];
			for (uint k = 0; k < i; k++)
				sum -= v[k] * A[j][k];
			if (i == j) {
				// if (unlikely(sum == 0))
				// 	return false;
				rdiag[i] = 1 / sum;
			} else {
				A[j][i] = sum;
			}
		}
	}

	return true;
}

// Solve Ax=b after ldltdc.  x is allowed to be the same as b.
void ldltsl(T A[N][N],
            T rdiag[N],
            T b[N],
            out T x[N]
){
	for (int i = 0; i < N; i++) {
		T sum = b[i];
		for (int k = 0; k < i; k++)
			sum -= A[i][k] * x[k];
		x[i] = sum * rdiag[i];
	}
	for (int i = N - 1; i >= 0; i--) {
		T sum = 0;
		for (int k = i + 1; k < N; k++)
			sum += A[k][i] * x[k];
		x[i] -= sum * rdiag[i];
	}
}
// Convenience function when output x overwrites input b.
void ldltsl(T A[N][N], T rdiag[N], inout T b[N])
{
	ldltsl(A, rdiag, b, b);
}
#undef N
#undef T



#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURV_PER_FACE)
uvec3 get_iverts_face(uint iface)
{
    if (iface == 0u) return uvec3(2u, 3u, 1u); /* face0 */
    return uvec3(0u, 1u, 3u);  /* face1 */
}

/* Code adopted from "TriMesh::need_curvatures()" in TriMesh_curvature.cc 
 * ----------------------------------------------------------- */
void calc_vert_curv_tensor_at_face( // v012 CCW
    uvec3 vid, vec3 vpos[3], vec3 vnor[3], 
    uint calc_for_vtx_i/*which vertex to solve, we dont need all 3 tensors*/, 
    // Note: we do weighting in the next per-vertex pass
    // float cornerAreasOnFace[3], float pointAreasOnFace[3]
    // Outputs
    out float curv1, out float curv12, out float curv2
)
{
    /* Edges are arranged corresponding to the paper */
    vec3 e[3] = { vpos[2] - vpos[1], vpos[0] - vpos[2], vpos[1] - vpos[0] };
    
    /* (pd1, pd2, normal) coord system per vertex, persistent among connected faces */
    vec3 pd1[3]; vec3 pd2[3];
    for (uint i = 0u; i < 3u; ++i)
        load_local_vertex_frame(vid[i], vpos[i], vnor[i], /*out*/pd1[i], pd2[i]); 

    /* N-T-B coordinate system per face */
    vec3 faceTangent = normalize(e[0]);
    vec3 faceNormal = normalize(cross(e[2], -e[1]));
    if(dot(faceNormal, vnor[0]) < 0.0) { faceNormal = -faceNormal; }
    vec3 faceBitangent = normalize(cross(faceNormal, faceTangent));

    /* ----------------------------------------------------------
    * Estimate curvature on face over normals' finite difference
    * Solving linear least squares problem 
    * ATAx = ATb 
    * by LU decomposition
    * m : Atb, w : ATA 
    * let 
    * ui      := dot(ei, u)
    * vi      := dot(ei, v)
    * dni     := (vnor[(i+2)%3] - vnor[(i+1)%3])
    * dnui    := dot(dni, u)
    * dnvi    := dot(dni, v)
    * => 
    * m = SymmetricMatrix_3x3 [
    *   sum(ui^2),  sum((ui*vi)   ,  0 
    *       ~       sum(ui^2+vi^2),  sum(dot(ui,vi))        
    *       ~             ~          sum(vi^2)        
    * ]
    * w = [sum(dnui*ui), sum(dnui*vi+dnvi*ui), sum(dnvi*vi)]
    * ---------------------------------------- */
    float m[3] = { 0, 0, 0 };
    float w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
    for (uint i = 0; i < 3; i++) 
    { // using the tangent - bitangent as uv coords
        float ui = dot(e[i], faceTangent);
        float vi = dot(e[i], faceBitangent);
        w[0][0] += ui*ui;
        w[0][1] += ui*vi;
        w[2][2] += vi*vi;
        // The below are computed once at the end of the loop
        // w[1][1] += v*v + u*u;
        // w[1][2] += u*v;

        uint prev = (i+2u) % 3u;
        uint next = (i+1u) % 3u;
        /* finite difference for vert normal along ei */
        vec3 dni = vnor[prev] - vnor[next];
        float dnu = dot(dni, faceTangent);
        float dnv = dot(dni, faceBitangent);

        m[0] += dnu*ui;
        m[1] += dnu*vi + dnv*ui;
        m[2] += dnv*vi;
    }
    w[1][1] = w[0][0] + w[2][2];
    w[1][2] = w[0][1];

#define USE_SFU_INVERSE_MATRIX 1
#if defined(USE_SFU_INVERSE_MATRIX)
    dmat3 w_cpy;
    for (uint col = 0u; col < 3u; ++col)
        w_cpy[col] = dvec3(double(w[0][col]), double(w[1][col]), double(w[2][col]));
    w_cpy[0][1] = w_cpy[1][0];
    w_cpy[0][2] = w_cpy[2][0];
    w_cpy[1][2] = w_cpy[2][1]; 
    dvec3 m_cpy = dvec3(m[0], m[1], m[2]);  

    dvec3 m_sol = (inverse(w_cpy) * m_cpy); 
    m[0] = float(m_sol[0]);
    m[1] = float(m_sol[1]);
    m[2] = float(m_sol[2]);
#else
    /* 3x3 LU Solve
    * See "ldltdc" and "ldltsl" in "RusinkiewiczEstimator.cs", from http://graphics.zcu.cz/curvature.html */
    float diag[3] = {0,0,0};
    float old_m[3] = { m[0], m[1], m[2] }; 
    ldltdc(w, diag); /* LU Decomposition */
    ldltsl(w, diag, old_m, /*out*/m);
#endif


    /* Curvature tensor for selected vertex */
    float c1, c12, c2;
    c1 = c12 = c2 = .0f; 
    {
        uint i = calc_for_vtx_i; 

        proj_curv(
            faceTangent, faceBitangent, m[0], m[1], m[2],  
            pd1[i],      pd2[i], /*out*/c1,   c12,  c2 
        ); 
        
        // weight = corner area / point area 
        // Voronoi area weighting. 
        float wt = 1.0f; // jwzw: we do this later;
    }
    curv1  = c1;
    curv12 = c12;
    curv2  = c2;
}

void main()
{
    uint groupIdx = gl_LocalInvocationID.x; 
    
    uint wedge; bool valid_thread; 
    if (0 < pcs_only_selected_verts_/* note: verts on selection border will fuck this up */)
    {
        uint sel_edge_id = gl_GlobalInvocationID.x; 
        get_wedge_id_from_selected_edge(sel_edge_id, /*out*/wedge, /*out*/valid_thread); 
    }else
    {
        wedge = gl_GlobalInvocationID.x; 
        uint num_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_edges; 
        valid_thread = wedge < num_edges; 
    }


    uvec4 v;
    Load4(ssbo_edge_to_vert_, wedge, v);

    /*        v0        
    *        /  \       
    *       /    \      
    *      W      3     
    *     0   f1   W    
    *    /          \   
    *   /  ------->  \  
    * v1 ---- W4 ---- v3
    *   \  <-------  /  
    *    \          /   
    *     W   f0   2    
    *      1      W     
    *       \    /      
    *        \  /       
    *         v2        
    */
    vec3 vpos[4] = {    
        ld_vpos(v.x), 
        ld_vpos(v.y), 
        ld_vpos(v.z), 
        ld_vpos(v.w)
    }; 
    vec3 vnor[4] = {
        ld_vnor(v.x), 
        ld_vnor(v.y), 
        ld_vnor(v.z), 
        ld_vnor(v.w)
    };

    /* calc  2 tensors: v1 at f0, v3 at f1 */
#define sl_iface_for_v1 0u
    const uint v1_id_at_solved_face = 2u; 
#define sl_iface_for_v3 1u
    const uint v3_id_at_solved_face = 2u; 

    float curv1_at_iface[2]; 
    float curv12_at_iface[2]; 
    float curv2_at_iface[2]; 
    for (uint iface = 0u; iface < 2u; ++iface)
    {
        uvec3 iverts = get_iverts_face(iface); 
        uvec3 vids_at_face = uvec3(v[iverts[0]], v[iverts[1]], v[iverts[2]]);
        vec3 vpos_at_face[3] = { vpos[iverts[0]], vpos[iverts[1]], vpos[iverts[2]] };
        vec3 vnor_at_face[3] = { vnor[iverts[0]], vnor[iverts[1]], vnor[iverts[2]] };
        calc_vert_curv_tensor_at_face(
            vids_at_face, vpos_at_face, vnor_at_face, 
            (iface == sl_iface_for_v1) ? v1_id_at_solved_face: v3_id_at_solved_face, 
            /* out */
            curv1_at_iface[iface], curv12_at_iface[iface], curv2_at_iface[iface]
        ); 
    }

    /* Actually, We only calc & store 2 edge verts' tensors */
    float corner_weight_v1; 
    float varea_1 = ld_varea(v[1]); 
    corner_weight_v1 = calc_voronoi_area_face(vpos[1], vpos[2], vpos[3]) / varea_1; 

    float corner_weight_v3; 
    float varea_3 = ld_varea(v[3]); /* TODO: move div-by-varea to next pass */
    corner_weight_v3 = calc_voronoi_area_face(vpos[3], vpos[0], vpos[1]) / varea_3; 

    vec3 v1_tensor_at_f0 = vec3(.0f, .0f, .0f); 
    {    
        float w1 = corner_weight_v1; 
        v1_tensor_at_f0.x = w1 * curv1_at_iface[sl_iface_for_v1]; 
        v1_tensor_at_f0.y = w1 * curv12_at_iface[sl_iface_for_v1];
        v1_tensor_at_f0.z = w1 * curv2_at_iface[sl_iface_for_v1]; 
    }
    
    vec3 v3_tensor_at_f1 = vec3(.0f, .0f, .0f); 
    {
        float w3 = corner_weight_v3;
        v3_tensor_at_f1.x = w3 * curv1_at_iface[sl_iface_for_v3];
        v3_tensor_at_f1.y = w3 * curv12_at_iface[sl_iface_for_v3];
        v3_tensor_at_f1.z = w3 * curv2_at_iface[sl_iface_for_v3];
    }

    uvec3 tensor_enc;
    if (valid_thread)
    { /* note: border edges not handled */
        tensor_enc = floatBitsToUint(v1_tensor_at_f0); 
        Store3(ssbo_edge_vtensors_, (wedge*2u), tensor_enc); 
        tensor_enc = floatBitsToUint(v3_tensor_at_f1); 
        Store3(ssbo_edge_vtensors_, (wedge*2u+1u), tensor_enc); 
    }
}
#endif    



#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__MAIN)

#if defined(INCLUDE_VERTEX_RADIAL_NORMAL)
    #define LOAD_VNOR(vid, vpos) ld_vnor_radial(vid, vpos)
#else
    #define LOAD_VNOR(vid, vpos) ld_vnor(vid)
#endif

struct CalcVertAttrContext_Order1
{
    vec3 vpos_c; 

    /* Rusinkiewicz's Curvature Estimator from 
     * https://github.com/Forceflow/trimesh2/blob/main/libsrc/TriMesh_curvature.cc 
    */
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        vec3 curv_tensor;  
    #endif

    /* JacquesOlivierLachaud's Curvature Estimator from 
     * https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h
     * https://dgtal-team.github.io/doc-nightly/moduleCurvatureMeasures.html 
     * https://github.com/CGAL/cgal/issues/7063
    */
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        float sum_area; 
        vec3 vnor_c; 
        vec3 vpos_i; 
        vec3 vnor_i;  
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
        mat3 curv_tensor; 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
        float mu1;
        float mu2; 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING)
        AdjWedgeInfo wo_prev; 
    #endif

    /* Gradient field of scalar field dot(view_dir, vnor) */
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        vec3 grad_vdotn; 
        float sum_weight_gradvdn; 
    #endif

    // TODO: precalc this
    bool border; 

    // adaptive len for debug lines ---
    float ave_edge_len; 
    float num_adj_edges; 
    // --------------------------------
}; 

CalcVertAttrContext_Order1 init_vert_attr_context_order_1(
    vec3 vpos
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
    , vec3 vnor, vec3 vpos_i, vec3 vnor_i
#endif
){
    CalcVertAttrContext_Order1 ctx; 
    ctx.vpos_c = vpos;  

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        ctx.sum_area = .0f; 
        ctx.vnor_c = vnor; 
        ctx.vpos_i = vpos_i; 
        ctx.vnor_i = vnor_i;
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        ctx.curv_tensor = vec3(.0f); 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
        ctx.curv_tensor = mat3(.0f); 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
        ctx.mu1 = .0f; 
        ctx.mu2 = .0f; 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING)
        ctx.wo_prev = AdjWedgeInfo(0u, 0u);  /* init at beginning iter */
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        ctx.grad_vdotn = vec3(.0f); 
        ctx.sum_weight_gradvdn = .0f;
    #endif

    ctx.border = false; 

    ctx.ave_edge_len = .0f;
    ctx.num_adj_edges = .0f; 

    return ctx; 
}


struct CalcVertAttrContext_Order1_2Ring
{
    float sum_area; /*:=mu0*/
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
    float mu1;
    float mu2;
#endif
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
    mat3 curv_tensor; /*:=muXY*/
#endif
    vec3 vpos_c; 
    vec3 vnor_c; 

    vec3 vpos_i; 
    vec3 vnor_i;  

    AdjWedgeInfo w_end; /* end traversal when wn==w_end */
}; 
bool calc_vert_attr_order_1_2ring(
    CirculatorIterData iter, 
    inout CalcVertAttrContext_Order1_2Ring ctx
)
{
    uint wi = iter.awi.wedge_id; 
    bool at_end = (iter.awi_next.wedge_id == ctx.w_end.wedge_id); 
    if (at_end)
        return false;

    EdgeFlags ef = load_edge_flags(wi); 
    if (ef.border) // note: its an empty face here
        return false; 

    uint ivert_n = mark__ve_circ_bck__get_vn(iter); 
    uint vn = ssbo_edge_to_vert_[wi*4u + ivert_n];

    vec3 vpos_n = ld_vpos(vn);
    vec3 vnor_n = LOAD_VNOR(vn, vpos_n);

    float area_tri = mu0InterpolatedU(ctx.vpos_c, ctx.vpos_i, vpos_n, ctx.vnor_c, ctx.vnor_i, vnor_n); 
    ctx.sum_area += area_tri; 
    
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
    float mu1_tri = mu1InterpolatedU(ctx.vpos_c, ctx.vpos_i, vpos_n, ctx.vnor_c, ctx.vnor_i, vnor_n); 
    ctx.mu1 += mu1_tri;

    float mu2_tri = mu2InterpolatedU(ctx.vpos_c, ctx.vpos_i, vpos_n, ctx.vnor_c, ctx.vnor_i, vnor_n);
    ctx.mu2 += mu2_tri;
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
    mat3 tensor_tri = muXYInterpolatedU(ctx.vpos_c, ctx.vpos_i, vpos_n, ctx.vnor_c, ctx.vnor_i, vnor_n); 
    ctx.curv_tensor += tensor_tri; 
#endif

    ctx.vpos_i = vpos_n; 
    ctx.vnor_i = vnor_n; 

    return true; 
}


/*    Example: Rotate fwd(CW) around V3(marked as vc)        
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
bool calc_vert_attr_order_1(
    CirculatorIterData iter, 
    inout CalcVertAttrContext_Order1 ctx
){
    uint wi = iter.awi.wedge_id; 
    uint ivert_i = mark__ve_circ_fwd__get_vi(iter);
    uint ivert_n = mark__ve_circ_fwd__get_vn(iter); 

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING)
        uint iwedge_wo = mark__vecirc_fwd_get_wo(iter); 
        AdjWedgeInfo wo = decode_adj_wedge_info(ssbo_edge_to_edges_[wi*4u + iwedge_wo]); 
        if (ctx.num_adj_edges < 1e-10f)
        { /* Initialize wo_prev */
            uint iwedge_woprev = mark__vecirc_fwd_get_wop(iter); 
            ctx.wo_prev = decode_adj_wedge_info(ssbo_edge_to_edges_[wi*4u + iwedge_woprev]); 
        }
    #endif

    /* Rusinkiewicz's Curvature Formula */
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        uint iface = iter.awi.iface_adj; 

        uvec3 tensor_enc; 
        Load3(ssbo_edge_vtensors_, (wi*2u+iface), tensor_enc); 
        ctx.curv_tensor += uintBitsToFloat(tensor_enc); 
    #endif

    /* "Interpolated corrected curvature measures for polygonal surfaces" */
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        uint vn = ssbo_edge_to_vert_[wi*4u + ivert_n];
        vec3 vpos_n = ld_vpos(vn);
        vec3 vnor_n = LOAD_VNOR(vn, vpos_n);
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
    { /* visit 2-ring */
        VertWedgeListHeader vwlh_2ring = VertWedgeListHeader(wi, ivert_i); 
        
        CalcVertAttrContext_Order1_2Ring ctx_2ring;
        ctx_2ring.mu1 = (.0f);
        ctx_2ring.mu2 = (.0f);
        ctx_2ring.sum_area = .0f;
        ctx_2ring.vpos_c = ctx.vpos_i; 
        ctx_2ring.vnor_c = ctx.vnor_i;
        ctx_2ring.vpos_i = ctx.vpos_c;
        ctx_2ring.vnor_i = ctx.vnor_c;
        ctx_2ring.w_end = ctx.wo_prev;

        bool rotate_fwd_2ring = false; 
        VE_CIRCULATOR(vwlh_2ring, calc_vert_attr_order_1_2ring, ctx_2ring, rotate_fwd_2ring); 

        ctx.mu1 += ctx_2ring.mu1; 
        ctx.mu2 += ctx_2ring.mu2; 
        ctx.sum_area += ctx_2ring.sum_area; 
    }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
        #if !defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING)
        {            
            float area_tri = mu0InterpolatedU(ctx.vpos_i, ctx.vpos_c, vpos_n, ctx.vnor_i, ctx.vnor_c, vnor_n); 
            mat3 tensor_tri = /* area_tri *  */muXYInterpolatedU(ctx.vpos_i, ctx.vpos_c, vpos_n, ctx.vnor_i, ctx.vnor_c, vnor_n); 
            ctx.sum_area += area_tri;
            ctx.curv_tensor += tensor_tri;
        }
        #else
        { /* visit 2-ring */
            VertWedgeListHeader vwlh_2ring = VertWedgeListHeader(wi, ivert_i); 
            
            CalcVertAttrContext_Order1_2Ring ctx_2ring;
            ctx_2ring.curv_tensor = mat3(.0f);
            ctx_2ring.sum_area = .0f;
            ctx_2ring.vpos_c = ctx.vpos_i; 
            ctx_2ring.vnor_c = ctx.vnor_i;
            ctx_2ring.vpos_i = ctx.vpos_c;
            ctx_2ring.vnor_i = ctx.vnor_c;
            ctx_2ring.w_end = ctx.wo_prev;

            bool rotate_fwd_2ring = false; 
            VE_CIRCULATOR(vwlh_2ring, calc_vert_attr_order_1_2ring, ctx_2ring, rotate_fwd_2ring); 

            ctx.curv_tensor += ctx_2ring.curv_tensor; 
            ctx.sum_area += ctx_2ring.sum_area; 
        }
        #endif
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING)
        ctx.wo_prev = wo; 
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
    {
        vec3 cam_pos_ws = ubo_view_matrices_.viewinv[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
        if (0 < pcs_dbg_ndv_grad_mode_) cam_pos_ws = ubo_view_matrices_last_frame_.viewinv[3].xyz; 
        vec3 vdir_c = normalize(ctx.vpos_c - cam_pos_ws);
        float abs_ndv_c = abs(dot(vdir_c, ctx.vnor_c)); 

        vec3 vdir_i = normalize(ctx.vpos_i - cam_pos_ws);
        float abs_ndv_i = abs(dot(vdir_i, ctx.vnor_i));

        vec3 vdir_n = normalize(vpos_n - cam_pos_ws);
        float abs_ndv_n = abs(dot(vdir_n, vnor_n)); 

        vec3 eic = ctx.vpos_c - ctx.vpos_i; float eic_len = length(eic); 
        vec3 ecn = vpos_n - ctx.vpos_c;     float ecn_len = length(ecn);
        vec3 eni = ctx.vpos_i - vpos_n;     float eni_len = length(eni); 
        float fi_area = length(cross(-eni, eic)) * .5f; 
        vec3 fi_nor = normalize(cross(-eni, eic));

        vec3 grad_fi = cross(fi_nor, (abs_ndv_c * eni + abs_ndv_i * ecn + abs_ndv_n * eic)) / (-2.0f*fi_area);

        float weight_gradvdn = acos(dot(-eni/eni_len, eic/eic_len)); 
        ctx.grad_vdotn          += weight_gradvdn * grad_fi;
        ctx.sum_weight_gradvdn  += weight_gradvdn;    
    }    
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        ctx.vpos_i = vpos_n; 
        ctx.vnor_i = vnor_n; 
    #endif

    EdgeFlags ef = load_edge_flags(wi); 
    if (ef.border) ctx.border = true;

    uint vi = ssbo_edge_to_vert_[wi*4u + ivert_i]; 
    ctx.ave_edge_len += length(ld_vpos(vi) - ctx.vpos_c);
    ctx.num_adj_edges += 1.0f; 
        
    return true; 
}

void main()
{
    uint groupIdx = gl_LocalInvocationID.x; 

    // Selective evaluation of vertices
    // --------------------------------
    uint vert_id; 
    bool valid_thread; 
    // remeshed verts
    if (0 < pcs_only_selected_verts_) {
        uint sel_vert_id = gl_GlobalInvocationID.x; 
        get_vert_id_from_selected_vert(sel_vert_id, /*out*/vert_id, valid_thread);
    }
    else {
        vert_id = gl_GlobalInvocationID.x; 
        uint num_verts = get_vert_count(); 
        valid_thread = vert_id < num_verts; 
    }
    // verts at the interpolated contour
    if (0 < pcs_order_1_eval_only_contour_verts_) {
        VertFlags vf = load_vert_flags(vert_id); 
        if (!vf.contour) valid_thread = false; 
        // unfortunately, we cannot early return since there are compactions to be done
    }
    // verts at the edges that has a interpolated contour vertex waiting to insert
    if (0 < pcs_order_1_eval_only_interpo_contour_adj_verts_) {
        VertFlags vf = load_vert_flags(vert_id); 
        if (!vf.adj_to_contour_vtx) valid_thread = false; 
    }
    
    
    
    // Actual computation starts here
    // ------------------------------
    vec3 vpos = ld_vpos(vert_id); 
    vec3 vnor = LOAD_VNOR(vert_id, vpos);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 

    VertWedgeListHeader vwlh_rot = vwlh;
    CalcVertAttrContext_Order1 ctx; 
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        // Accumulate partial tensors from 1-ring faces
        ctx = init_vert_attr_context_order_1(vpos); 
    #endif
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR) || defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        uint ivert_vi = vwlh_rot.ivert == 1u ? 3u : 1u;
        uint vi = ssbo_edge_to_vert_[vwlh_rot.wedge_id*4u + ivert_vi]; 
        vec3 vpos_vi = ld_vpos(vi);
        vec3 vnor_vi = LOAD_VNOR(vi, vpos_vi);
        ctx = init_vert_attr_context_order_1(
            vpos, vnor, vpos_vi, vnor_vi
        );
    #endif

    if (valid_thread)
    {
        bool rotate_fwd = true; 
        VE_CIRCULATOR(vwlh_rot, calc_vert_attr_order_1, ctx, rotate_fwd);
    }
    ctx.ave_edge_len /= float(ctx.num_adj_edges); 


    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        // Compute principal direction & curvature
        float curv_1 = ctx.curv_tensor.x;
        float curv_12 = ctx.curv_tensor.y;
        float curv_2 = ctx.curv_tensor.z;
        vec3 uaxis, vaxis; 
        load_local_vertex_frame(vwlh, vpos, vnor, /*out*/uaxis, vaxis); 

        vec3 pdir1, pdir2; 
        float curv_1_fin, curv_2_fin; 
        diagonalize_curv(
            uaxis, vaxis, curv_1, curv_12, curv_2, vnor, 
            /*out*/pdir1, pdir2, curv_1_fin, curv_2_fin
        ); 


        float max_curv = max(abs(curv_1_fin), abs(curv_2_fin)); 
        bool valid_curv = !(isnan(max_curv) || isinf(max_curv) || ctx.border);
        if (!valid_curv) max_curv = -1.0f; 

        if (valid_thread)
        {
            if (0 < pcs_output_curv_tensors_)
                st_vcurv_pdirs_k1k2(vert_id, pdir1, curv_1_fin, pdir2, curv_2_fin);
            st_vcurv_max(vert_id, max_curv); 
        }

        // debug lines
        if (0 < pcs_output_dbg_geom_)
        {
            VertFlags vf = decode_vert_flags(ssbo_vert_flags_[vert_id]); 
        
            bool dbg_vtx_curv = (!vf.dupli) && (!vf.del_by_collapse) && valid_thread;
            uint dbg_line_id = compact_general_dbg_lines(dbg_vtx_curv, groupIdx, 3u);
            dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 

            if (dbg_vtx_curv)
            {
                float dbg_line_len = .1f * pcs_dbg_geom_scale_; // min(ctx.ave_edge_len, pcs_dbg_geom_scale_ * .12f);

                vec4 vpos_ws_10 = vec4(vpos /* - curv_1_fin * normalize(pdir1) * dbg_line_len */, 1.0f);
                vec4 vpos_ws_11 = vec4(vpos /* + curv_1_fin * normalize(pdir1) * dbg_line_len */, 1.0f);
                vec3 dvd_col = vec3(1.0, .0, .0); 
                DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dvd_col, uvec4(0u)); 
                DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dvd_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_10, dvd_11); 
                dbg_line_id++; 

                vec4 vpos_ws_20 = vec4(vpos /* - curv_2_fin * normalize(pdir2) * dbg_line_len */, 1.0f);
                vec4 vpos_ws_21 = vec4(vpos /* + curv_2_fin * normalize(pdir2) * dbg_line_len */, 1.0f);
                dvd_col = vec3(.0, 1.0, .0); 
                DebugVertData dvd_20 = DebugVertData(vpos_ws_20.xyz, dvd_col, uvec4(0u)); 
                DebugVertData dvd_21 = DebugVertData(vpos_ws_21.xyz, dvd_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_20, dvd_21); 
                dbg_line_id++; 
                
                #define ssbo_vtx_remesh_len_ ssbo_vcurv_pdirs_k1k2_
                float remesh_edge_len = uintBitsToFloat(ssbo_vtx_remesh_len_[vert_id]); 
                #undef ssbo_vtx_remesh_len_
                dbg_line_len = valid_curv ? max_curv * pcs_dbg_geom_scale_ : .0f;
                // dbg_line_len = valid_curv ? remesh_edge_len * pcs_dbg_geom_scale_ : .0f;

                vec4 vpos_ws_30 = vec4(vpos, 1.0f);
                vec4 vpos_ws_31 = vec4(vpos + vnor * dbg_line_len, 1.0f);
                dvd_col = vec3(.0, 1.0, 1.0); 
                DebugVertData dvd_30 = DebugVertData(vpos_ws_30.xyz, dvd_col, uvec4(0u)); 
                DebugVertData dvd_31 = DebugVertData(vpos_ws_31.xyz, dvd_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_30, dvd_31); 
                dbg_line_id++; 
            }
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE)
        ctx.mu1 /= ctx.sum_area; // Mean
        ctx.mu2 /= ctx.sum_area; // Gaussian

        vec2 curv_fin = vec2(ctx.mu1, ctx.mu2); 
        bool valid_curv = !((any(isnan(curv_fin)) || any(isinf(curv_fin))) || ctx.border);
        if (!valid_curv) ctx.mu1 = ctx.mu2 = .0f;  

        float max_curv = ctx.mu1 + sqrt(max(.0f, ctx.mu1 * ctx.mu1 - abs(ctx.mu2)));
        if (!valid_curv) max_curv = -1.0f; 

        if (valid_thread)
        {
            if (0 < pcs_output_curv_tensors_)
                st_vcurv_pdirs_k1k2(vert_id, pdir1, curv_1_fin, pdir2, curv_2_fin);
            st_vcurv_max(vert_id, max_curv); 
        }

        // debug lines
        if (0 < pcs_output_dbg_geom_)
        {
            VertFlags vf = decode_vert_flags(ssbo_vert_flags_[vert_id]); 
            bool dbg_vtx_curv = (!vf.dupli) && (!vf.del_by_collapse) && valid_thread; 
            uint dbg_line_id = compact_general_dbg_lines(dbg_vtx_curv, groupIdx, 3u);
            dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 

            if (dbg_vtx_curv)
            {
                float dbg_line_len = min(ctx.ave_edge_len * .4f, pcs_dbg_geom_scale_ * .12f);

                vec4 vpos_ws_00 = vec4(vpos, 1.0f);
                vec4 vpos_ws_01 = vec4(vpos/*  + vnor * ctx.mu1 * dbg_line_len */, 1.0f);
                vec3 dbg_col = vec3(1.0, .0, .0); 
                DebugVertData dvd_00 = DebugVertData(vpos_ws_00.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_01 = DebugVertData(vpos_ws_01.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_00, dvd_01); 
                dbg_line_id++; 

                vec4 vpos_ws_10 = vec4(vpos, 1.0f);
                vec4 vpos_ws_11 = vec4(vpos/*  + vnor * ctx.mu2 * dbg_line_len */, 1.0f);
                dbg_col = vec3(.0, 1.0, .0); 
                DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_10, dvd_11); 
                dbg_line_id++; 

                dbg_line_len = pcs_dbg_geom_scale_ * .1f;
                vec4 vpos_ws_20 = vec4(vpos, 1.0f);
                vec4 vpos_ws_21 = vec4(vpos + vnor * max_curv * dbg_line_len, 1.0f);
                dbg_col = vec3(1.0, .0, 1.0); 
                DebugVertData dvd_20 = DebugVertData(vpos_ws_20.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_21 = DebugVertData(vpos_ws_21.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_20, dvd_21); 
                dbg_line_id++; 
            }
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR)
        mat3 vtx_curv_measure = ctx.curv_tensor;
        float vtx_area_measure = ctx.sum_area; 
        vec3 vnor_gt = ld_vnor(vert_id); // NOTE: This must be the geometrically correct normal
        mat3 pdirs; 
        vec3 evals; 
        if (valid_thread) // expensive, so only do it for valid threads
            curvDirFromTensor(vtx_curv_measure, ctx.sum_area, vnor_gt, /*out*/pdirs, evals);

        vec3 pdir0 = vec3(pdirs[0][0], pdirs[1][0], pdirs[2][0]); 
        vec3 pdir1 = vec3(pdirs[0][1], pdirs[1][1], pdirs[2][1]); 
        vec3 pdir2 = vec3(pdirs[0][2], pdirs[1][2], pdirs[2][2]); 
        evals /= vtx_area_measure; 

        // force tangency, there could be tiny deviations from the local tangent plane
		pdir0 = pdir0 - vnor * dot(pdir0, vnor); 
		pdir1 = pdir1 - vnor * dot(pdir0, vnor); 

        
        float max_curv = max(abs(evals[0]), abs(evals[1])); 
        bool valid_curv = !(isnan(max_curv) || isinf(max_curv));
        if (!valid_curv)
            max_curv = evals[0] = evals[1] = -1.0f; 


        float cusp_func = .0f; 
        float cusp_func_alternative = .0f; 
        bool near_contour = false; 
        if (0 < pcs_output_maxcurv_with_cusp_function_ && valid_thread)
        { // Cusp detection from "Illustrating smooth surface" by Hertzmann et al.
            mat4 view_to_world = ubo_view_matrices_.viewinv;
            bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
            vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
            cusp_func = calc_cusp_func(pdir0, pdir1, evals.x, evals.y, vpos, vnor, cam_pos_ws); 
            
            vec3 binormal = calc_cusp_binormal(pdir0, pdir1, evals.x, evals.y, vpos, vnor, cam_pos_ws); 
            float bdotv = dot(binormal, normalize(vpos - cam_pos_ws)); 
            cusp_func_alternative = bdotv; 

            float ndv = dot(normalize(cam_pos_ws - vpos), vnor);
            near_contour = abs(ndv) < .2f; 
        }
        

        if (valid_thread)
        {
            if (0 < pcs_output_curv_tensors_)
                st_vcurv_pdirs_k1k2(vert_id, pdir0, evals[0], pdir1, evals[1]);
            if (0 == pcs_output_maxcurv_with_cusp_function_)
                st_vcurv_max(vert_id, max_curv); 
            else
                st_vcurv_max_with_cusp(vert_id, max_curv, cusp_func_alternative); // cusp_func); 
        }

        
        if (0 < pcs_output_dbg_geom_)
        {
            // debug lines
            VertFlags vf = decode_vert_flags(ssbo_vert_flags_[vert_id]); 
            bool dbg_vtx_curv = (!vf.dupli) && (!vf.del_by_collapse) && valid_thread;

            uint dbg_line_id = compact_general_dbg_lines(dbg_vtx_curv, groupIdx, 3u); /* must run for every thread */
            dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 

            if (dbg_vtx_curv)
            {
                float dbg_line_len = min(ctx.ave_edge_len * .4f, pcs_dbg_geom_scale_ * .12f);
                // dbg_line_len = .0f; 

                vec4 vpos_ws_00 = vec4(vpos - normalize(pdir0) * dbg_line_len, 1.0f);
                vec4 vpos_ws_01 = vec4(vpos + normalize(pdir0) * dbg_line_len, 1.0f);
                vec3 dbg_col = vec3(1.0, .0, .0); 
                DebugVertData dvd_00 = DebugVertData(vpos_ws_00.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_01 = DebugVertData(vpos_ws_01.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_00, dvd_01); 
                dbg_line_id++; 


                // if (!vf.crease) dbg_line_len = .0f; 
                // dbg_line_len = max_curv > pcs_dbg_geom_scale_ ? .05f : .0f; 
                // dbg_line_len = (near_contour && cusp_func.x < .0f) ? pcs_dbg_geom_scale_/*  * cusp_func.x */ : .0f; 
                // dbg_line_len = (vf.front_facing) ? pcs_dbg_geom_scale_/*  * cusp_func.x */ : .0f; 
                // vec4 vpos_ws_10 = vec4(vpos, 1.0f);
                // vec4 vpos_ws_11 = vec4(vpos + normalize(vnor) * dbg_line_len, 1.0f);
                vec4 vpos_ws_10 = vec4(vpos - normalize(pdir1) * dbg_line_len, 1.0f);
                vec4 vpos_ws_11 = vec4(vpos + normalize(pdir1) * dbg_line_len, 1.0f);
                dbg_col = vec3(.0, 1.0, .0); 
                DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_10, dvd_11); 
                dbg_line_id++; 

                // dbg_line_len = pcs_dbg_geom_scale_; 
                // if (!vf.corner) dbg_line_len = .0f; 
                // dbg_line_len = .0f; 
                #define ssbo_vtx_remesh_len_ ssbo_vcurv_pdirs_k1k2_
                float remesh_edge_len = uintBitsToFloat(ssbo_vtx_remesh_len_[vert_id]); 
                // dbg_line_len = pcs_dbg_geom_scale_ * remesh_edge_len;
                #undef ssbo_vtx_remesh_len_
                // dbg_line_len = valid_curv ? .0f : pcs_dbg_geom_scale_;                
                // dbg_line_len = (near_contour && cusp_func.x >= .0f) ? pcs_dbg_geom_scale_/*  * cusp_func.x */ : .0f; 
                // dbg_line_len = (vf.back_facing) ? pcs_dbg_geom_scale_/*  * cusp_func.x */ : .0f; 
                vec4 vpos_ws_20 = vec4(vpos, 1.0f);
                vec4 vpos_ws_21 = vec4(vpos + normalize(vnor) * dbg_line_len, 1.0f);
                dbg_col = vec3(1.0, .0, 1.0); 
                DebugVertData dvd_20 = DebugVertData(vpos_ws_20.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_21 = DebugVertData(vpos_ws_21.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_20, dvd_21); 
                dbg_line_id++; 
            }
        }
    #endif

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN)
        ctx.grad_vdotn /= ctx.sum_weight_gradvdn; 
        
        float grad_vdotn_len = length(ctx.grad_vdotn);
        vec3 bigrad = cross(ctx.grad_vdotn / grad_vdotn_len, vnor); 
        
        if (valid_thread)
        {
            uvec3 grad_enc = floatBitsToUint(ctx.grad_vdotn); 
            Store3(ssbo_vgrad_contour_, vert_id, grad_enc); 
        }

        if (0 < pcs_output_dbg_geom_)
        {
            VertFlags vf = decode_vert_flags(ssbo_vert_flags_[vert_id]); 
            bool dbg_vtx_grad = (!vf.dupli) && (!vf.del_by_collapse) && valid_thread; 
            
            uint dbg_line_id = compact_general_dbg_lines(dbg_vtx_grad, groupIdx, 3u);
            dbg_line_id += get_debug_line_offset(DBG_LINE_TYPE__GENERAL); 
        
            float dbg_line_len = pcs_dbg_geom_scale_ * min(ctx.ave_edge_len * .4f, grad_vdotn_len);

            if (dbg_vtx_grad)
            {
                vec3 grad_vdotn_dir = grad_vdotn_len < 1e-10f ? vec3(.0f) : ctx.grad_vdotn / grad_vdotn_len;

                vec4 vpos_ws_00 = vec4(vpos, 1.0f);
                vec4 vpos_ws_01 = vec4(vpos + (ctx.grad_vdotn / grad_vdotn_len) * dbg_line_len, 1.0f);
                vec3 dbg_col = vec3(1.0, .0, .0); 
                DebugVertData dvd_00 = DebugVertData(vpos_ws_00.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_01 = DebugVertData(vpos_ws_01.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_00, dvd_01); 
                dbg_line_id++; 
                
                vec4 vpos_ws_10 = vec4(vpos, 1.0f);
                vec4 vpos_ws_11 = vec4(vpos + (ctx.grad_vdotn / grad_vdotn_len) * dbg_line_len, 1.0f);
                dbg_col = vec3(.0, 1.0, .0); 
                DebugVertData dvd_10 = DebugVertData(vpos_ws_10.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_11 = DebugVertData(vpos_ws_11.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_10, dvd_11); 
                dbg_line_id++; 

                vec4 vpos_ws_20 = vec4(vpos, 1.0f); // dummy
                vec4 vpos_ws_21 = vec4(vpos, 1.0f);
                dbg_col = vec3(1.0, .0, 1.0); 
                DebugVertData dvd_20 = DebugVertData(vpos_ws_20.xyz, dbg_col, uvec4(0u)); 
                DebugVertData dvd_21 = DebugVertData(vpos_ws_21.xyz, dbg_col, uvec4(0u)); 
                store_debug_line_data(dbg_line_id, dvd_20, dvd_21); 
                dbg_line_id++; 
            }
        }
    #endif
    
}
#endif 
#endif






#if defined(_KERNEL_MULTICOMPILE__CALC_FEATURE_EDGES)

uint get_edge_count()
{
    return pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_edges; 
}

void main()
{
    uint groupIdx = gl_LocalInvocationID.x; 

    uint edge_id; 
    bool valid_thread; 
    if (0 < pcs_only_selected_edges_)
    {
        uint sel_edge_id = gl_GlobalInvocationID.x; 
        get_wedge_id_from_selected_edge(sel_edge_id, /*out*/edge_id, valid_thread);
    }
    else
    {
        edge_id = gl_GlobalInvocationID.x; 
        uint num_edges = get_edge_count(); 
        valid_thread = edge_id < num_edges; 
    }

    EdgeFlags ef = load_edge_flags(edge_id); 
    uvec4 v;
    for (uint i = 0u; i < 4u; i++)
        v[i] = ssbo_edge_to_vert_[edge_id*4u + i]; 

    bool update_edge_flags = false; 

#if defined(_KERNEL_MULTICOMPILE__CALC_FEATURE_EDGES__CREASE_DETECTION)
    if (!valid_thread) return; 

    vec3 vpos[4]; 
    for (uint i = 0u; i < 4u; i++)
        vpos[i] = ld_vpos(v[i]);
    float dihedral = calc_dihedral_angle(vpos[0], vpos[1], vpos[2], vpos[3]); 

    uint crease_level = ((M_PI / 3.0f) < abs(dihedral - M_PI)) ? 3 : 0; 
    if (ef.border) crease_level = 3; 

    update_edge_flags__detect_crease(crease_level, ef);
    update_edge_flags = true; 
#endif

#if defined(_KERNEL_MULTICOMPILE__CALC_FEATURE_EDGES__CONTOUR_SPLIT_DETECTION)
    if (!valid_thread) return; 

    uvec2 iverts_cwedge = mark__cwedge_to_verts(0u);
    uvec2 vids_cwedge = uvec2(v[iverts_cwedge[0]], v[iverts_cwedge[1]]); 
    VertFlags vfs_cwedge[2] = { 
        load_vert_flags(vids_cwedge[0]), 
        load_vert_flags(vids_cwedge[1]) 
    };  

    // Check if the edge is a contour edge - TODO: make this a flag bit, compute only once
    bool is_interp_contour_edge = 
        is_interp_contour_edge__before_tessellation(vfs_cwedge[0], vfs_cwedge[1], ef); 

    // Update edge flags
    update_edge_flags__detect_contour_split(is_interp_contour_edge, /*inout*/ef);
    update_edge_flags = true; 

    // Update vert flags
    // only write when true to avoid racing threads writhing opposite booleans
    // !!! but this requires a clean up pass to set the flags to false !!!
    if (is_interp_contour_edge) 
        for (uint iivert = 0; iivert < 2; ++iivert) { 
            update_vert_flags__adj_to_contour_vtx(true, vfs_cwedge[iivert]); 
            store_vert_flags(vids_cwedge[iivert], vfs_cwedge[iivert]); 
        }
#endif

    if (update_edge_flags)
        store_edge_flags(edge_id, ef); 
}

#endif











#if defined(_KERNEL_MULTICOMPILE__INTERPOLATE_CONTOUR_VERT_ATTRS)

struct InterpVtxAttrs
{
    uint vid; 
    vec3 vpos; 
    vec3 vnor; 
    float max_curv; 
    float cusp_func; 
}; 

vec3 get_cam_pos_ws()
{
    mat4 view_to_world = ubo_view_matrices_.viewinv;
    vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
    return cam_pos_ws; 
}

void main()
{
    // Selective evaluation of vertices
    // --------------------------------
    uint vert_id; 
    bool valid_thread; 
    if (0 < pcs_only_selected_verts_) {
        uint sel_vert_id = gl_GlobalInvocationID.x; 
        get_vert_id_from_selected_vert(sel_vert_id, /*out*/vert_id, valid_thread);
    }
    else {
        vert_id = gl_GlobalInvocationID.x; 
        uint num_verts = get_vert_count(); 
        valid_thread = vert_id < num_verts; 
    }
    // Only for interpolated contour vertices    
    VertFlags vf = load_vert_flags(vert_id); 
    if (!vf.contour) return; 


    // Actual computation starts here
    // ------------------------------ 
    uint old_edge = ssbo_contour_vert_to_old_edge_[vert_id]; // edge that was split to create this contour vertex

    InterpVtxAttrs old_verts[2]; {
        uvec2 iverts_cwedge = mark__cwedge_to_verts(0u);  
        old_verts[0].vid = ssbo_edge_to_vert_[old_edge*4u + iverts_cwedge[0]];
        old_verts[1].vid = ssbo_edge_to_vert_[old_edge*4u + iverts_cwedge[1]];

        for (uint i = 0; i < 2; i++)
        {
            old_verts[i].vpos = ld_vpos(old_verts[i].vid); 
            old_verts[i].vnor = ld_vnor(old_verts[i].vid); 
            ld_vcurv_max_with_cusp(
                old_verts[i].vid, 
                /*out*/old_verts[i].max_curv, old_verts[i].cusp_func
            );
        }
    }

    vec3 edge_vnor[2] = { old_verts[0].vnor, old_verts[1].vnor };
    vec3 edge_vpos[2] = { old_verts[0].vpos, old_verts[1].vpos }; 
    float interp_factor = calc_interp_contour_edge_factor(edge_vnor, edge_vpos, get_cam_pos_ws()); 

    InterpVtxAttrs new_vert;
    new_vert.vid = vert_id;
    new_vert.vpos = mix(old_verts[0].vpos, old_verts[1].vpos, interp_factor);
    new_vert.vnor = mix(old_verts[0].vnor, old_verts[1].vnor, interp_factor);
    new_vert.max_curv = mix(old_verts[0].max_curv, old_verts[1].max_curv, interp_factor);
    new_vert.cusp_func = mix(old_verts[0].cusp_func, old_verts[1].cusp_func, interp_factor);

    if (valid_thread)
    {
        // st_vpos(vert_id, new_vert.vpos); // already done by edge splitting
        st_vnor(vert_id, new_vert.vnor); 
        st_vcurv_max_with_cusp(vert_id, new_vert.max_curv, new_vert.cusp_func); 
    }
}

#endif



