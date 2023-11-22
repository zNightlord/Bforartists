#ifndef BNPR_MESHING_TOPO__INCLUDED
#define BNPR_MESHING_TOPO__INCLUDED


#if !defined(DECODE_IBO_EXCLUDE) && defined(DECODE_IBO_INCLUDE) 
/* line adjacency("edge detection") IBO layout */
/*    v3
 *   /  \
 *  /    \                
 * v1----v2             
 *  \    /                
 *   \  /                     
 *    v0    winding 012, 321                
*/
/** Input
 * #define IBO_BUF XXX   <= index buffer, either 16 bits or 32 bits per index
*/
void load_and_decode_ibo__edge_adj(uint prim_id, uint vertsPerPrim, out uvec4 vid)
{
#if defined(_KERNEL_MULTICOMPILE__INDEX_BUFFER_16BIT)
    uint ibo_items_16 = vertsPerPrim / 2; /* 16 bits per vtx id */
    for (uint i = 0; i < 2; ++i) 
    {
        uint ibo_addr = ibo_items_16 * prim_id + i;
        uint ibo_data = IBO_BUF[ibo_addr];
        /* decode 16 bit index */
        uint ibo_data_16h = (ibo_data >> 16u);
        uint ibo_data_16l = (ibo_data & 0xFFFFu);
        /* fetch vertex pos */
        vid[i * 2] = (ibo_data_16l);
        vid[i * 2 + 1] = (ibo_data_16h);
    }
#else
    uint ibo_items_32 = vertsPerPrim; /* 32 bits per vtx id */
    for (uint i = 0; i < 4; ++i) 
    {
        uint ibo_addr = ibo_items_32 * prim_id + i;
        uint ibo_data = IBO_BUF[ibo_addr];
        /* fetch vertex pos */
        vid[i] = (ibo_data);
    }
#endif
}
#endif





#if defined(WINGED_EDGE_TOPO_INCLUDE)

/* Winged edge structure, 
 * derived from line adjacency("edge detection") IBO layout */

uvec4 line_adj_to_wing_verts(uvec4 line_adj_verts)
{
    uvec4 wing_verts = line_adj_verts.xywz; /* 0123 -> 0132 */
    return wing_verts;
}
uvec4 wing_verts_to_line_adj(uvec4 wing_verts)
{
    uvec4 line_adj_verts = wing_verts.xywz; /* 0123 -> 0132 */
    return line_adj_verts;
}

/* Vertices: x4 vert ids per wing */
/*    v2
 *   /  \
 *  /    \                
 * v1>-->v3 for wedge e, vi => ssbo_edge_to_vert_[e*4+i] , i=0,1,2,3            
 *  \    /  
 *   \  /   Face winding               
 *    v0    013, 123                                  
*/
/* Same to line_adj, we treat the border wedge as two overlapping tris with oppo winding */
bool line_adj_is_border_edge(uvec4 line_adj_verts)
{ /* see "extract_lines_adjacency_finish" */
    bool is_border = line_adj_verts.x == line_adj_verts.w;
    return is_border; 
}
bool wing_verts_is_border_edge(uvec4 wing_verts)
{
    return line_adj_is_border_edge(wing_verts_to_line_adj(wing_verts)); 
}
bool is_border_edge_front_facing(uint ifacefront)
{
    return (ifacefront == 1u); /* := the "actual" face (iface==1 in wedge, v0v1v2 in line_adj) */
}

/* wedges: at most x4 adj. wedges */
/*          v2           
 *        //  \\         v  = vertex     
 *       // /\ \\        he = half-edge     
 *      /1 1  h W\       F  = face     
 *     /W e    e 2\      W  = wedge     
 *    // h  F0  2 \\            
 *   // <--0-eh--\ \\    short cuts 
 *    v1=== W4 ===v3     ----------
 *   \\ \--he-0--> //    iwedge, ivert, ihe, iface: marks of the current wedge, vertex, half-edge, face       
 *    \\ 2  F1  h //     cwedge: center wedge W4
 *     \0 e    e W/      bwedge: border wedges W0~3
 *      \W h  1 3/       
 *       \\ \/ //      
 *        \\  //        
 *          v0         
*/
/*          v0           
 *        //  \\         v  = vertex     
 *       // /\ \\        he = half-edge     
 *      /W h  1 3\       F  = face     
 *     /0 e    e W\      W  = wedge     
 *    // 2  F1  h \\            
 *   // /--eh-0--> \\    short cuts 
 *    v1=== W4 ===v3     ----------
 *   \\ <--he-0--/ //    iwedge, ivert, ihe, iface: marks of the current wedge, vertex, half-edge, face       
 *    \\ h  F0  2 //     cwedge: center wedge W4
 *     \W e    e 2/      bwedge: border wedges W0~3
 *      \1 1  h W/       
 *       \\ \/ //      
 *        \\  //        
 *          v2         
*/
#define NULL_EDGE 0xffffffffu /* non-existing wedges */
uint mark__border_wedge_to_oppo_vert(uint wedge)
{
    return (wedge <= 1u) ? 3u : 1u; /* w0,1->3, w2,3->1 */
}
uint mark__center_wedge_to_oppo_vert__at_face(uint face)
{
    return (face == 0u) ? 2u : 0u; /* f0,w4->2, f1,w4->0 */ 
}
uvec2 mark__wedge_to_verts(uint wedge)
{
    if (wedge == 4u) /* center edge */
        return uvec2(1u, 3u); 
    /* border wedges */
    return uvec2(wedge, (wedge + 1u) % 4u);
}
uint mark__cwedge_to_beg_vert(uint face)
{
    if (face == 0u) return 3u;
    /* if (face == 1u) */ return 1u; 
}
uint mark__cwedge_to_end_vert(uint face)
{
    return (mark__cwedge_to_beg_vert(face) == 1u) ? 3u : 1u; 
}
uint mark__he_to_prev_he(uint he)
{
    return ((he + 2u) % 3u); /* 0->2, 1->0, 2->1 */
}
uint mark__he_to_next_he(uint he)
{
    return ((he + 1u) % 3u); /* 0->1, 1->2, 2->0 */
}
uint mark__he_to_wedge(uint face, uint he)
{
    if (face == 0u)
        return he == 0u ? 4u : he; /* 0->4, 1->1, 2->2 */
    /* if (face == 1u) */
    return he == 0u ? 4u : (he+2u)%4u; /* 0->4, 1->3, 2->0 */
}
uint mark__wedge_to_he(uint face, uint wedge)
{
    if (face == 0u)
        return wedge == 4u ? 0u : wedge; /* 4->0, 1->1, 2->2 */
    /* if (face == 1u) */
    return wedge == 4u ? 0u : (wedge+2u)%4u; /* 4->0, 3->1, 0->2 */
}
uint mark__bwedge_to_face(uint wedge)
{
    return (wedge == 1u || wedge == 2u) ? 0u : 1u; /* 1,2->0, 0,3->1 */
}
uint mark__vert_to_next_vert(uint face, uint vert)
{
    if (face == 0u)
        return (vert == 3u) ? 1u : vert + 1u; /* 1->2, 2->3, 3->1 */
    /* if (face == 1u) */
    return (vert == 1u) ? 3u : ((vert + 1u) % 4u); /* 1->3, 3->0, 0->1 */ 
}
uvec3 mark__face_to_winded_wedges(uint face)
{
    uvec3 res; 
    for (uint he = 0u; he < 3u; ++he)
        res[he] = mark__he_to_wedge(face, he); 
    return res; 
}

/* 
 * Each wedge data also encodes 
 * the mark of non-overlaping adjacent face - F0 or F1
*/
struct AdjWedgeInfo
{
    uint wedge_id;
    uint iface_adj;
};
uint encode_adj_wedge_info(AdjWedgeInfo awi)
{
    uint awi_enc = ((awi.wedge_id << 1u) | (awi.iface_adj & 1u));
    return awi_enc; 
}
AdjWedgeInfo decode_adj_wedge_info(uint awi_enc)
{
    AdjWedgeInfo awi; 
    awi.wedge_id = (awi_enc >> 1u);
    awi.iface_adj = ((awi_enc & 1u));
    return awi; 
}

uint mark__cwedge_rotate_back(uint face)
{
    uint he_pre = mark__he_to_prev_he(0u);
    uint wedge_prev = mark__he_to_wedge(face, he_pre); 
    return wedge_prev; 
}
uint mark__cwedge_rotate_next(uint face)
{
    uint he_next = mark__he_to_next_he(0u); 
    uint wedge_next = mark__he_to_wedge(face, he_next); 
    return wedge_next; 
}


/* Contour-specific functions ------------------- */
struct PerContourWedgeInfo
{
    bool is_border; 
    uint wedge_id;
    uint ifrontface;
};
uint encode_per_contour_wedge_info(PerContourWedgeInfo pcwi)
{
    uint pcwi_enc = ((pcwi.wedge_id << 2u) | ((pcwi.ifrontface & 1u) << 1u) | (uint(pcwi.is_border))); 
    return pcwi_enc; 
} 
PerContourWedgeInfo decode_per_contour_wedge_info(uint pcwi_enc)
{
    PerContourWedgeInfo pcwi; 
    pcwi.wedge_id = (pcwi_enc >> 2u);
    pcwi.ifrontface = ((pcwi_enc >> 1u) & 1u);
    pcwi.is_border = (0u != (pcwi_enc & 1u));
    return pcwi; 
}

struct PerWedgeContourInfo
{
    uint contour_id;
    uint ifrontface;
    bool is_contour;
    bool is_border; 
};
/* Encode for all fields in PerWedgeContourInfo */
uint encode_per_wedge_contour_info(PerWedgeContourInfo peci)
{
    uint peci_enc = ((peci.contour_id << 3u) | ((peci.ifrontface & 1u) << 2u) | (uint(peci.is_contour) << 1u) | (uint(peci.is_border))); 
    return peci_enc; 
}
/* Decode for all fields in PerWedgeContourInfo */
PerWedgeContourInfo decode_per_wedge_contour_info(uint peci_enc)
{
    PerWedgeContourInfo peci; 
    peci.contour_id = (peci_enc >> 3u);
    peci.ifrontface = ((peci_enc >> 2u) & 1u);
    peci.is_contour = (0u != ((peci_enc >> 1u) & 1u));
    peci.is_border = (0u != (peci_enc & 1u));
    return peci; 
}



/* Detect unstable wedges --------------------------------------- */
struct WedgeQuality
{
    float area; 
    float dihedral; 
    bool unstable; 
}; 

WedgeQuality compute_wedge_quality(vec3 p0, vec3 p1, vec3 p2, vec3 p3, vec3 cam_pos_ws) /* world pos of p0~3 */
{ /* impl based on overlay_outline_prepass_vert_no_geom.glsl */
    WedgeQuality res; 
    
	vec3 v10 = p0 - p1;
	vec3 v13 = p3 - p1;
   	vec3 v12 = p2 - p1;

	vec3 n0 = normalize(cross(v13, v10));
	vec3 n2 = normalize(cross(v12, v13));
	
    mat3 m = mat3(v10, v13, vec3(1, 1, 1));
    res.area = abs(determinant(m)); 

    /* 
     * TODO: add dihedral angle, only mark unstable edge with near planar angle 
     * https://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleMeshDerivativesCheatSheet.pdf */
    res.dihedral = atan(dot(v13, cross(n0, n2)), dot(n0, n2));

    vec3 view_dir = normalize(cam_pos_ws - p1); 
	float face_orient_012 = dot(view_dir, n0);
  	float face_orient_321 = dot(view_dir, n2);
    res.unstable = (abs(face_orient_012) < 0.08f) || (abs(face_orient_321) < 0.08f);

	return res; 
}

struct WedgeFloodingPointer
{
    uint next_wedge_id; 
    bool is_border; /* if current wedge is border */
    bool is_seed; /* if wedge pointed by next_wedge_id is a seed */
}; 
uint encode_wedge_flooding_pointer(WedgeFloodingPointer wfp)
{
    uint wfp_enc = ((wfp.next_wedge_id << 2u) | ((uint(wfp.is_border) << 1u) | (uint(wfp.is_seed)))); 
    return wfp_enc; 
}
WedgeFloodingPointer decode_wedge_flooding_pointer(uint wfp_enc)
{
    WedgeFloodingPointer wfp; 
    wfp.next_wedge_id = (wfp_enc >> 2u);
    wfp.is_border = (0u != ((wfp_enc >> 1u) & 1u));
    wfp.is_seed = (0u != (wfp_enc & 1u));
    return wfp; 
}

#endif








#if defined(QUADRICS_FILTERING_INCLUDE)

/* Quardratic Filtering */
struct Quadric
{
    mat4 quadric; /* symmetric */
    float area; 
}; 
/* Pack lower diagonal of 4x4 symmetric matrix "quadric" into 10 floats */
void encode_quadric(Quadric q, out vec4 packed_0, out vec4 packed_1, out vec3 packed_2)
{
    packed_0 = vec4(q.quadric[0][0], q.quadric[1][1], q.quadric[2][2], q.quadric[3][3]);
    packed_1 = vec4(q.quadric[1][0], q.quadric[2][1], q.quadric[3][2], q.quadric[2][0]);
    packed_2 = vec3(q.quadric[3][1], q.quadric[3][0], q.area);
} 
/* Unpack lower diagonal of 4x4 symmetric matrix "quadric" from 10 floats */
Quadric decode_quadric(vec4 packed_0, vec4 packed_1, vec3 packed_2)
{
    Quadric q; 
    q.quadric[0][0] = packed_0.x;
    q.quadric[1][1] = packed_0.y;
    q.quadric[2][2] = packed_0.z;
    q.quadric[3][3] = packed_0.w;
    q.quadric[1][0] = q.quadric[0][1] = packed_1.x;
    q.quadric[2][1] = q.quadric[1][2] = packed_1.y;
    q.quadric[3][2] = q.quadric[2][3] = packed_1.z;
    q.quadric[2][0] = q.quadric[0][2] = packed_1.w;
    q.quadric[3][1] = q.quadric[1][3] = packed_2.x;
    q.quadric[3][0] = q.quadric[0][3] = packed_2.y;
    q.area = packed_2.z; 

    return q; 
}

/* minimizing [nT(x-p)]^2 */
mat4 compute_plane_quadric(vec3 n, vec3 p)
{
    vec4 nq = vec4(n.xyz, 1.0f);
    nq.w = -dot(nq.xyz, p); 
    return mat4(nq.x * nq, nq.y * nq, nq.z * nq, nq.w * nq);
}

mat4 compute_tri_quadric(vec3 p0, vec3 p1, vec3 p2)
{
    mat3 det = mat3(p0, p1, p2); 
    vec4 nq = vec4(cross(p0, p1)+cross(p1, p2)+cross(p2,p0), -determinant(det)); 
    return mat4(nq.x * nq, nq.y * nq, nq.z * nq, nq.w * nq); 
}

/* minimizing ||x-p||^2 */
mat4 compute_point_quadric(vec3 p)
{ 
    mat4 q = mat4(1.0);
    q[3] = vec4(-p.xyz, dot(p, p));
    q[0][3] = -p.x; 
    q[1][3] = -p.y; 
    q[2][3] = -p.z; 

    return q; 
}

mat4 compute_view_binormal_plane_quadric(vec3 n, vec3 p, vec3 cam_pos_ws)
{
    vec3 v = normalize(cam_pos_ws - p); 
    vec3 bn = normalize(cross(n, v));

    return compute_plane_quadric(bn, p);  
}

mat4 compute_view_tagent_plane_quadric(vec3 n, vec3 p, vec3 cam_pos_ws)
{
    vec3 v = normalize(cam_pos_ws - p); 
    vec3 bn = normalize(cross(n, v));
    vec3 t = normalize(cross(bn, n)); 

    return compute_plane_quadric(t, p);  
}

mat4 compute_view_plane_quadric(vec3 n, vec3 p, vec3 cam_pos_ws)
{
    vec3 v = normalize(cam_pos_ws - p); 
    vec3 bn = normalize(cross(n, v));
    vec3 t = normalize(cross(bn, n)); 
    vec3 n_view_proj = normalize(cross(t, bn));

    return compute_plane_quadric(n_view_proj, p);  
}

float compute_view_plane_dist(vec3 n0, vec3 p0, vec3 cam_pos_ws, vec3 p1)
{
    vec3 v = normalize(cam_pos_ws - p0); 
    vec3 bn = normalize(cross(n0, v));
    vec3 t = normalize(cross(bn, n0)); 

    float d = dot(p0, t); 
    return abs(dot(p1, t) - d); 
}

/* Match with "GPUMeshQuadricFilter"
 */ 
#define GEOM_NORMAL_PLANE 0u
#define VIEW_NORMAL_PLANE 1u
#define VIEW_TANGENT_PLANE 2u
#define VIEW_BINORMAL_PLANE 3u
mat4 compute_filter_quadric(vec3 n, vec3 p, vec3 cam_pos_ws, float silouetteness, float pos_reg_weight)
{
    mat4 res; 

    if (pcs_filtered_quadric_type_ == GEOM_NORMAL_PLANE)
        res = compute_plane_quadric(n, p); 
    else if (pcs_filtered_quadric_type_ == VIEW_NORMAL_PLANE)
        res = (0.2f + silouetteness) * compute_view_plane_quadric(n, p, cam_pos_ws); 
    else if (pcs_filtered_quadric_type_ == VIEW_TANGENT_PLANE)
        res = (0.2f + silouetteness) * compute_view_tagent_plane_quadric(n, p, cam_pos_ws);
    else if (pcs_filtered_quadric_type_ == VIEW_BINORMAL_PLANE)
        res = (0.2f + silouetteness) * compute_view_binormal_plane_quadric(n, p, cam_pos_ws);
        
    res += pos_reg_weight * compute_point_quadric(p); /* regularization */

    return res;
}


float compute_wedge_area(vec3 p0, vec3 p1, vec3 p2, vec3 p3)
{
    vec3 v10 = p0 - p1;
	vec3 v13 = p3 - p1;
   	vec3 v12 = p2 - p1;
    
    vec3 v13_cross_v10 = cross(v13, v10); 
    vec3 v12_cross_v13 = cross(v12, v13); 

    /* Averaged area */
    float area = ((length(v13_cross_v10)) + length(v12_cross_v13)) * .25f; 

    return area; 
}

vec3 compute_wedge_normal(vec3 p0, vec3 p1, vec3 p2, vec3 p3)
{
    vec3 v10 = p0 - p1;
	vec3 v13 = p3 - p1;
   	vec3 v12 = p2 - p1;
    
    vec3 v13_cross_v10 = cross(v13, v10); 
    vec3 v12_cross_v13 = cross(v12, v13); 

	vec3 n0 = normalize(v13_cross_v10);
	vec3 n2 = normalize(v12_cross_v13);

    return normalize((n0.xyz + n2.xyz) * .5f); 
}

Quadric compute_wedge_quadric(
    vec3 p0, vec3 p1, vec3 p2, vec3 p3, vec3 cam_pos_ws, float position_regularize_weight, 
    vec3 wedge_normal
) /* world pos of p0~3 */
{ 
    Quadric q; 
    
	vec3 v10 = p0 - p1;
	vec3 v13 = p3 - p1;
   	vec3 v12 = p2 - p1;
    
    vec3 v13_cross_v10 = cross(v13, v10); 
    vec3 v12_cross_v13 = cross(v12, v13); 

    /* Averaged area */
    q.area = ((length(v13_cross_v10)) + length(v12_cross_v13)) * .25f; 

    /* Quadric */
	vec3 n0 = normalize(v13_cross_v10);
    vec3 c0 = (p0 + p1 + p3) / 3.0f;     
    float ndv0 = dot(normalize(cam_pos_ws - c0), n0); 

	vec3 n2 = normalize(v12_cross_v13);
    vec3 c2 = (p1 + p2 + p3) / 3.0f; 
    float ndv2 = dot(normalize(cam_pos_ws - c2), n2);

    float silouetteness = .0f; 
    if (sign(ndv0) != sign(ndv2))
    {
        silouetteness = abs(ndv0 - ndv2) * .5f; 
    }

    mat4 q0 = compute_filter_quadric(n0, c0, cam_pos_ws, silouetteness, position_regularize_weight); 
    mat4 q2 = compute_filter_quadric(n2, c2, cam_pos_ws, silouetteness, position_regularize_weight); 
    q.quadric = (q0 + q2) * .5f;

    wedge_normal = normalize((n0.xyz + n2.xyz) * .5f); 

	return q; 
}

float gaussian(float dist, float tau)
{
    return exp(-(dist * dist) / (2.0f * tau * tau)); 
}

void bilateral_filter_wedge_normal(
    vec3 p0, vec3 p1, vec3 p2, vec3 p3, vec3 np, 
    vec3 q0, vec3 q1, vec3 q2, vec3 q3, vec3 nq, 
    vec3 cam_pos_ws, 
    out vec3 n_filtered, out float weight
){
    float q_area = compute_wedge_area(q0, q1, q2, q3); 

    float dist = distance(p0, p1); 
    float geometry_weight = q_area * gaussian(dist, 1.0f);
    
    vec3 pp = (p1 + p3) / 3.0f;
    vec3 qq = (q1 + q3) / 3.0f; 
    float plane_dist = max(abs(dot(pp - qq, np)), abs(dot(pp - qq, nq)));
    float plane_weight = gaussian(plane_dist, 0.1f);

    // Fast and Effective Feature-Preserving Mesh Denoising
    float normal_weight = dot(np, nq);
    if (normal_weight < 0.8f) normal_weight = 0.0f; 
    normal_weight *= normal_weight; 

    weight = normal_weight/*  * plane_weight *//*  * geometry_weight */; 
    
    
    vec3 vp = normalize(cam_pos_ws - p0); 
    float cos_theta = min(1.0f, max(-1.0f, dot(np, vp))); 
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta); 
    
    vec3 bn = normalize(cross(np, vp));
    vec3 t = normalize(cross(bn, np)); 
    vec3 n_view_proj = normalize(cross(t, bn));

    n_filtered = nq; // normalize(sin_theta * n_view_proj + cos_theta * vp); 
}


float compute_edge_quadric_weight(vec3 p_v, vec3 p_e, Quadric q_e)
{
    float dist = distance(p_v, p_e); 
    float geometry_weight = q_e.area * gaussian(dist, 1.0f); 
    return geometry_weight; 
}

float compute_vert_quadric_weight(vec3 p_v, Quadric q_v, vec3 p_x, Quadric q_x, float dev_q, float dev_g)
{ /* bilateral filtering */
    float dist = distance(p_v, p_x); 
    
    vec4 v = vec4(p_v, 1.0f); 
    float quadric_dist_v2qx = dot(v, q_x.quadric * v); /* vT Q v */
    float quadric_weight = gaussian(sqrt(abs(quadric_dist_v2qx)), dev_q); 

    float geometry_weight = /* q_v.area *  */gaussian(dist, dev_g); 
    
    return geometry_weight * quadric_weight; 
}

#endif








#if defined(VERT_WEDGE_LIST_TOPO_INCLUDE)
/* Non-existing wedges, which should not happen for a 2-manifold mesh with boundaries */
#define NULL_VERT_ADJ_WEDGE 0xffffffffu 
struct VertWedgeListHeader
{
    uint wedge_id; 
    /* the face that contains this vertex that begines at the center wedge */
    uint ivert; /* 1 or 3 */
}; 
uint encode_vert_wedge_list_header(VertWedgeListHeader vwlh)
{
    uint vwlh_enc = ((vwlh.wedge_id << 2u) | (vwlh.ivert & 3u));
    return vwlh_enc; 
}
VertWedgeListHeader decode_vert_wedge_list_header(uint vwlh_enc)
{
    VertWedgeListHeader vwlh; 
    vwlh.wedge_id = (vwlh_enc >> 2u);
    vwlh.ivert = ((vwlh_enc & 3u));
    return vwlh; 
}

#define NULL_FILTERED_VERT 0xffffffffu/* used in vert-to-filtered_vert buffer */
#endif








#endif




