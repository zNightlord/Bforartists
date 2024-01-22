
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_geom_lib.glsl)



/*
.define("WINGED_EDGE_TOPO_INCLUDE", "1")
.define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
.define("VE_CIRCULATOR_INCLUDE", "1")
.define("DYNAMESH_SELECTION_INDEXING_COMMON", "1")
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
        if (ang_c > PI_HALF) 
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
    float ang_312 = acos(dot(normalize(e31), -normalize(e12)));
    float ang_123 = acos(dot(normalize(e12), -normalize(e23)));
    float ang_231 = acos(dot(normalize(e23), -normalize(e31)));
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
    uint vert_id = gl_GlobalInvocationID.x; 
    uint num_verts = get_vert_count(); 
    bool valid_thread = vert_id < num_verts; 

    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]);

    CalcVertAttrContext_Order0 ctx = init_vert_attr_context_order_0(
        ssbo_edge_to_vert_[vwlh.wedge_id*4u + ((vwlh.ivert == 1u) ? 3u : 1u)], 
        ld_vpos(vert_id)
    ); 
    if (valid_thread)
    {
        bool rotate_fwd = true; 
        VE_CIRCULATOR(vwlh, calc_vert_attr_order_0, ctx, rotate_fwd);
    }
    
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL)
    vec3 vnormal = normalize(ctx.sum_normal / ctx.sum_weight); 
    if (valid_thread)
        st_vnor(vert_id, vnormal);

    if (0 < pcs_output_dbg_geom_)
    {
        VertFlags vf = decode_vert_flags(ssbo_vert_flags_[vert_id]); 
       
        bool dbg_vtx_nor = (!vf.dupli) && valid_thread; 
        uint dbg_prim_id = compact_normal_line(dbg_vtx_nor, groupIdx);
        
        if (dbg_vtx_nor)
        {
            vec4 vpos_ws_0 = vec4(ctx.vpos, 1.0f);
            vec4 vpos_ws_1 = vec4(ctx.vpos + vnormal * pcs_dbg_geom_scale_ * .05f, 1.0f);

            uvec3 vpos_enc = floatBitsToUint(vpos_ws_0.xyz); 
            Store3(ssbo_dbg_lines_, dbg_prim_id*6u, vpos_enc);
            vpos_enc = floatBitsToUint(vpos_ws_1.xyz); 
            Store3(ssbo_dbg_lines_, dbg_prim_id*6u + 3u, vpos_enc); 
        }
    }
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
#if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1)
/* Inputs
 * uint ssbo_edge_vtensors_[] // temp buffer
 * uint ssbo_vcurv_tensor_[] // final output 
 * uint ssbo_vcurv_pdirs_k1k2_[]


*/

void load_local_vertex_frame(VertWedgeListHeader vwlh, vec3 vpos, vec3 vnor, out vec3 pd1, out vec3 pd2)
{
    uint ivert_another_on_edge = vwlh.ivert == 1u ? 3u : 1u; 
    uint another_vid = ssbo_edge_to_vert_[vwlh.wedge_id*4u + ivert_another_on_edge];
    vec3 vpos_another = ld_vpos(another_vid);

    pd1 = vpos_another - vpos;
    pd1 = normalize(cross(vnor, pd1));

    pd2 = cross(pd1, vnor); 
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
	vec3 dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);

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
	// if (likely(kuv != 0.0f)) 
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
    *   sum(ui^2), sum(dot(ui,vi)), 0  
    *       ~      sum(ui^2+vi^2),  sum(dot(ui,vi))        
    *       ~            ~          sum(vi^2)        
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

    /* 3x3 LU Solve
    * See "ldltdc" and "ldltsl" in "RusinkiewiczEstimator.cs", from http://graphics.zcu.cz/curvature.html */
    float diag[3] = {0,0,0};
    float old_m[3] = { m[0], m[1], m[2] }; 
    ldltdc(w, diag); /* LU Decomposition */
    ldltsl(w, diag, old_m, /*out*/m);

    // Curvature tensor for each vertex of the face
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
    uint wedge = gl_GlobalInvocationID.x; 
    uint groupIdx = gl_LocalInvocationID.x; 
    uint num_edges = pcs_edge_count_ + ssbo_dyn_mesh_counters_out_.num_edges; 
    bool valid_thread = wedge < num_edges; 

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
    /* TODO: replace this mapping with a function */
    const uvec3 iverts_fk[2] = {
        uvec3(2u, 3u, 1u), /* face0 */
        uvec3(0u, 1u, 3u)  /* face1 */
    }; 
    const uvec2 v1_id_at_iface = uvec2(2u, 1u); 
    const uvec2 v3_id_at_iface = uvec2(1u, 2u); 

    /* calc  2 tensors: v1 at f0, v3 at f1 */
    #define sl_iface_for_v1 0u
    #define sl_iface_for_v3 1u
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
            iface == 0 ? v1_id_at_iface[iface]: v3_id_at_iface[iface], 
            /* out */
            curv1_at_iface[iface], curv12_at_iface[iface], curv2_at_iface[iface]
        ); 
    }

    /* Actually, We only calc & store 2 edge verts' tensors */
    float corner_weight_v1; 
    float varea_1 = ld_varea(v[1]); 
    corner_weight_v1 = calc_voronoi_area_face(vpos[1], vpos[2], vpos[3]) / varea_1; 

    float corner_weight_v3; 
    float varea_3 = ld_varea(v[3]); 
    corner_weight_v3 = calc_voronoi_area_face(vpos[3], vpos[0], vpos[1]) / varea_3; 

    vec3 v1_tensor_at_f0 = vec3(.0f, .0f, .0f); 
    {    
        float w1 = corner_weight_v1; 
        uint v1_index = v1_id_at_iface[sl_iface_for_v1]; 
        v1_tensor_at_f0.x += w1 * curv1_at_iface[sl_iface_for_v1]; 
        v1_tensor_at_f0.y += w1 * curv12_at_iface[sl_iface_for_v1];
        v1_tensor_at_f0.z += w1 * curv2_at_iface[sl_iface_for_v1]; 
    }
    
    vec3 v3_tensor_at_f1 = vec3(.0f, .0f, .0f); 
    {
        float w3 = corner_weight_v3;
        uint v3_index = v3_id_at_iface[sl_iface_for_v3];
        v3_tensor_at_f1.x += w3 * curv1_at_iface[sl_iface_for_v3];
        v3_tensor_at_f1.y += w3 * curv12_at_iface[sl_iface_for_v3];
        v3_tensor_at_f1.z += w3 * curv2_at_iface[sl_iface_for_v3];
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
struct CalcVertAttrContext_Order1
{
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        vec3 curv_tensor;  
    #endif
}; 

CalcVertAttrContext_Order1 init_vert_attr_context_order_1(

){
    CalcVertAttrContext_Order1 ctx; 

    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        ctx.curv_tensor = vec3(.0f); 
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
bool calc_vert_attr_order_1(
    CirculatorIterData iter, 
    inout CalcVertAttrContext_Order1 ctx
){
    uint wi = iter.awi.wedge_id; 
    
    #if defined(_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR)
        uint iface = iter.awi.iface_adj; 

        uvec3 tensor_enc; 
        Load3(ssbo_edge_vtensors_, (wi*2u+iface), tensor_enc); 
        ctx.curv_tensor += uintBitsToFloat(tensor_enc); 
    #endif

    return true; 
}

void main()
{
    uint groupIdx = gl_LocalInvocationID.x; 
    uint vert_id = gl_GlobalInvocationID.x; 
    uint num_verts = get_vert_count(); 
    bool valid_thread = vert_id < num_verts; 

    vec3 vpos = ld_vpos(vert_id); 
    vec3 vnor = ld_vnor(vert_id);
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vert_id]); 

    // Accumulate partial tensors from 1-ring faces
    VertWedgeListHeader vwlh_rot = vwlh;
    CalcVertAttrContext_Order1 ctx = init_vert_attr_context_order_1(); 
    if (valid_thread)
    {
        bool rotate_fwd = true; 
        VE_CIRCULATOR(vwlh_rot, calc_vert_attr_order_1, ctx, rotate_fwd);
    }

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

        if (valid_thread)
        {
            // st_vcurv_tensor(vert_id, ctx.curv_tensor); 
            st_vcurv_pdirs_k1k2(vert_id, pdir1, curv_1_fin, pdir2, curv_2_fin);
        }
    #endif
}
#endif 



#endif