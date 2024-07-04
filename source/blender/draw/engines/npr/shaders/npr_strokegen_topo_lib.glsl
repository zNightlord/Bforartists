#ifndef BNPR_MESHING_TOPO__INCLUDED
#define BNPR_MESHING_TOPO__INCLUDED


/* Code block dependencies */
#if defined(VE_CIRCULATOR_INCLUDE)
    #define WINGED_EDGE_TOPO_INCLUDE 1 
    #define VERT_WEDGE_LIST_TOPO_INCLUDE 1
#endif



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
uint mark__border_iface_mainfold()
{ /* in line_adj triangle 012 has the valid vertex winding */
    return 1u; 
}
uvec2 mark__border_bwedges_mainfold()
{
    return uvec2(0, 3); 
}


/* Vertices: x4 vert ids per wing */
/*    v2
 *   /  \
 *  / f0 \                
 * v1>-->v3 for wedge e, vi => ssbo_edge_to_vert_[e*4+i] , i=0,1,2,3            
 *  \ f1 /  
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
    return (ifacefront == mark__border_iface_mainfold()); /* := the "actual" face (iface==1 in wedge, v0v1v2 in line_adj) */
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
uvec2 mark__cwedge_to_verts(uint face)
{
    if (face == 1u) 
        return uvec2(1u, 3u); 
    return uvec2(3u, 1u); 
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
uint mark__wedge_to_prev_wedge(uint face, uint wedge)
{
    uint he = mark__wedge_to_he(face, wedge);
    uint he_prev = mark__he_to_prev_he(he);
    return mark__he_to_wedge(face, he_prev);
}
uint mark__wedge_to_next_wedge(uint face, uint wedge)
{
    uint he = mark__wedge_to_he(face, wedge);
    uint he_next = mark__he_to_next_he(he);
    return mark__he_to_wedge(face, he_next);
}
uint mark__bwedge_to_face(uint wedge)
{
    return (wedge == 1u || wedge == 2u) ? 0u : 1u; /* 1,2->0, 0,3->1 */
}
uvec2 mark__face_to_bwedges(uint face)
{
    return (face == 0u) ? uvec2(1u, 2u) : uvec2(3u, 0u); 
}  
uint mark__bwedge_to_prev_bwedge(uint wedge)
{ /* Walk along the quad border */
    return ((wedge + 3u) % 4u); 
}
uint mark__bwedge_to_next_bwedge(uint wedge)
{ /* Walk along the quad border */
    return ((wedge + 1u) % 4u); 
}
uint mark__vert_to_next_vert(uint face, uint vert)
{
    if (face == 0u)
        return (vert == 3u) ? 1u : vert + 1u; /* 1->2, 2->3, 3->1 */
    /* if (face == 1u) */
    return (vert == 1u) ? 3u : ((vert + 1u) % 4u); /* 1->3, 3->0, 0->1 */ 
}
uvec3 mark__face_to_winded_verts(uint face)
{
    if (face == 0u) return uvec3(3u, 1u, 2u);
    return uvec3(1u, 3u, 0u); // face == 1u
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
    bool unstable_silouette;
    bool silouette; 
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
    res.unstable_silouette = (abs(face_orient_012) < 0.2f) && (abs(face_orient_321) < 0.2f);
    res.silouette = sign(face_orient_012) != sign(face_orient_321); 

	return res; 
}



/* Edge Selection --------------------------------------------------- 
 * provide function to diffuse from a initial selection
 * construct mapping between selected-edges and edges 
*/
struct WedgeFloodingPointer
{
    uint next_wedge_id; 
    bool is_border; /* if current wedge is border */
    bool is_seed; /* if wedge pointed by next_wedge_id is a seed */
    bool is_unstable; /* if wedge is an unstable sillouette */
}; 
uint encode_wedge_flooding_pointer(WedgeFloodingPointer wfp)
{
    uint wfp_enc = (
        (wfp.next_wedge_id << 3u) 
        | (uint(wfp.is_border) << 2u) 
        | (uint(wfp.is_seed) << 1u)
        | (uint(wfp.is_unstable))
    ); 
    return wfp_enc; 
}
WedgeFloodingPointer decode_wedge_flooding_pointer(uint wfp_enc)
{
    WedgeFloodingPointer wfp; 
    wfp.next_wedge_id = (wfp_enc >> 3u);
    wfp.is_border = (0u != ((wfp_enc >> 2u) & 1u));
    wfp.is_seed = (0u != ((wfp_enc >> 1u) & 1u));
    wfp.is_unstable = (0u != (wfp_enc & 1u)); 
    return wfp; 
}

struct EdgeSelectionInfo
{
    bool is_null_edge; 
    uint edge_id; // mapping between selected-edges and edges
    bool is_unstable; 
    bool is_border; 
}; 
uint encode_edge_selection_info(EdgeSelectionInfo fei)
{
    uint data = 0u; 
    data = (fei.edge_id << 2u); 
    data |= ((fei.is_unstable ? 1u : 0u) << 1u); 
    data |= (fei.is_border   ? 1u : 0u); 

    if (fei.is_null_edge)
        data = NULL_EDGE; 
        
    return data; 
}
EdgeSelectionInfo decode_edge_selection_info(uint data)
{
    EdgeSelectionInfo fei; 
    fei.is_null_edge = (data == NULL_EDGE); 
    fei.edge_id = data >> 2u; 
    fei.is_unstable = (data & 0x2u) != 0u; 
    fei.is_border   = (data & 0x1u) != 0u; 
    return fei; 
}

/* Vertex Selection --------------------------------------------------- 
 * Select vertex by various ways
*/
struct VertSelectionInfo
{
    uint vert_id; 
};
uint encode_vert_selection_info(VertSelectionInfo fvi)
{
    uint data = 0u; 
    data = fvi.vert_id; 

    return data; 
}
VertSelectionInfo decode_vert_selection_info(uint data)
{
    VertSelectionInfo fvi; 
    fvi.vert_id = data; 

    return fvi; 
}
#endif


/* Indexing for selection for dynaremesh */
#if defined(USE_DYNAMESH_EDGE_SELECTION_INDEXING)
    /* Note: DO NOT use wedge_id==0 to do any that must be done by one thread! like initialization or counter-zero-out */
    void get_wedge_id_from_selected_edge(uint sel_edge_id, out uint wedge_id, out bool valid_thread)
    {
        uint num_dyn_edges = ssbo_dyn_mesh_counters_out_.num_edges; 
        uint num_static_edges = pcs_edge_count_;
        
        uint num_presel_edges  = ssbo_bnpr_mesh_pool_counters_.num_filtered_edges;
        uint num_all_sel_edges = num_presel_edges + num_dyn_edges; /* all dyn edges are selected */
        valid_thread = (sel_edge_id < num_all_sel_edges); 

        bool is_dyn_edge = (num_presel_edges <= sel_edge_id) && valid_thread; 
        if (is_dyn_edge)
        {
            wedge_id = num_static_edges + (sel_edge_id - num_presel_edges); 
        }else{
            EdgeSelectionInfo eseli = decode_edge_selection_info(
                ssbo_selected_edge_to_edge_[sel_edge_id]
            ); 
            wedge_id = eseli.edge_id;
        }
    }
#endif


#if defined(USE_DYNAMESH_VERT_SELECTION_INDEXING)
void get_vert_id_from_selected_vert(uint sel_vert_id, out uint vert_id, out bool valid_thread)
{
    uint num_dyn_verts = ssbo_dyn_mesh_counters_out_.num_verts; 
    uint num_static_verts = pcs_vert_count_;

    uint num_presel_verts  = ssbo_bnpr_mesh_pool_counters_.num_filtered_verts;
    uint num_all_sel_verts = num_presel_verts + num_dyn_verts; /* all dyn verts are selected */
    valid_thread = (sel_vert_id < num_all_sel_verts); 

    bool is_dyn_vert = (num_presel_verts <= sel_vert_id) && valid_thread; 
    if (is_dyn_vert)
        vert_id = num_static_verts + (sel_vert_id - num_presel_verts); 
    else{
        VertSelectionInfo eseli = decode_vert_selection_info(
            ssbo_selected_vert_to_vert_[sel_vert_id]
        ); 
        vert_id = eseli.vert_id;
    }
}
#endif



#if defined(VERT_WEDGE_LIST_TOPO_INCLUDE)
/* Non-existing wedges, which should not happen for a 2-manifold mesh with boundaries */
#define NULL_VERT_ADJ_WEDGE 0xffffffffu 
struct VertWedgeListHeader
{
    uint wedge_id; 
    uint ivert; /* 1 or 3, ivert of this vert on #wedge_id */
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
/* Adjust pointer to another edge */
void try_update_ve_link(uint vtx, uint edge_old, VertWedgeListHeader vwlh_new)
{
    VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vtx]); 
    if (vwlh.wedge_id == edge_old)
        ssbo_vert_to_edge_list_header_[vtx] = encode_vert_wedge_list_header(vwlh_new); 
}

#define NULL_FILTERED_VERT 0xffffffffu/* used in vert-to-filtered_vert buffer */
#endif



/* Mesh Element Flags ------------------------------------------------------------------------- */
struct VertFlags
{
    bool dupli; // duplicated vertex, should be ignored, TODO: consider removing in hashing pass
    bool new_by_split; // new vertex created by edge split
    bool new_by_face_split; // new vertex created by face split
    bool mov_by_collapse; // moved vertex by edge collapse
    bool del_by_collapse; // deleted vertex by edge collapse
    /* we provide 4 slots for selection */
    bvec4 selected; // selected for a certain operation
    bool contour; // contour vertex
    bool front_facing; // front facing vertex (dot(n, v) > .0f)
    bool back_facing;  // back facing vertex (dot(n, v) < .0f)
    bool border_eval; // border vertex, only valid after evaluated
    bool crease; // crease vertex for subdivision
    bool corner; // corner vertex for subdivision 
}; 
VertFlags init_vert_flags(bool dupli)
{
    VertFlags vf; 
    vf.dupli = dupli; 
    vf.new_by_split = false; 
    vf.new_by_face_split = false; 
    vf.mov_by_collapse = false; 
    vf.del_by_collapse = false; 
    vf.selected = bvec4(false); 
    vf.contour = false;
    // mark facing as undefined
    vf.front_facing = false;
    vf.back_facing = false;
    vf.border_eval = false; 
    vf.crease = false;
    vf.corner = false;
    
    return vf; 
}
VertFlags init_vert_flags__new_split_edge(bool is_split_for_contour, bool is_crease_edge)
{
    VertFlags vf; 
    vf.dupli = false; 
    vf.new_by_split = true; 
    vf.new_by_face_split = false; 
    vf.mov_by_collapse = false; 
    vf.del_by_collapse = false; 
    vf.selected = bvec4(true); /* always selected for remeshing */
    vf.contour = is_split_for_contour;
    // mark facing as undefined
    vf.front_facing = false;
    vf.back_facing = false;
    vf.border_eval = false;
    vf.crease = is_crease_edge;
    vf.corner = false;
    
    return vf; 
}
VertFlags init_vert_flags__new_split_face()
{
    VertFlags vf; 
    vf.dupli = false; 
    vf.new_by_split = false; 
    vf.new_by_face_split = true; 
    vf.mov_by_collapse = false; 
    vf.del_by_collapse = false; 
    vf.selected = bvec4(true); /* always selected for remeshing */
    vf.contour = false;
    // mark facing as undefined
    vf.front_facing = false;
    vf.back_facing = false;
    vf.border_eval = false;
    vf.crease = false;
    vf.corner = false;
    
    return vf; 
}


uint encode_vert_flags(VertFlags vf)
{
    uint vf_enc = 0u;
    vf_enc = uint(vf.dupli); 
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.new_by_split); 
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.new_by_face_split);
    vf_enc <<= 1u;
    vf_enc |= uint(vf.mov_by_collapse);
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.del_by_collapse); 
    vf_enc <<= 1u; 

    vf_enc |= uint(vf.selected[0]);
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.selected[1]);
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.selected[2]);
    vf_enc <<= 1u; 
    vf_enc |= uint(vf.selected[3]);

    vf_enc <<= 1u;
    vf_enc |= uint(vf.contour);
    vf_enc <<= 1u;
    vf_enc |= uint(vf.front_facing);
    vf_enc <<= 1u;
    vf_enc |= uint(vf.back_facing);
    
    vf_enc <<= 1u;
    vf_enc |= uint(vf.border_eval);
    vf_enc <<= 1u;
    vf_enc |= uint(vf.crease);
    vf_enc <<= 1u;
    vf_enc |= uint(vf.corner);

    return vf_enc; 
}
VertFlags decode_vert_flags(uint vf_enc)
{
    VertFlags vf; 

    vf.corner = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.crease = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.border_eval = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;

    vf.back_facing = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.front_facing = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.contour = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;

    vf.selected[3] = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.selected[2] = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.selected[1] = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.selected[0] = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;

    vf.del_by_collapse = (1u == (vf_enc & 1u)); 
    vf_enc >>= 1u; 
    vf.mov_by_collapse = (1u == (vf_enc & 1u));
    vf_enc >>= 1u; 
    vf.new_by_face_split = (1u == (vf_enc & 1u));
    vf_enc >>= 1u;
    vf.new_by_split = (1u == (vf_enc & 1u)); 
    vf_enc >>= 1u; 
    vf.dupli = (1u == (vf_enc & 1u)); 

    return vf; 
}
#if defined(VERT_FLAGS_INCLUDED)
    VertFlags load_vert_flags(uint vid)
    {
        return decode_vert_flags(ssbo_vert_flags_[vid]); 
    }
    void store_vert_flags(uint vid, VertFlags vf)
    {
        ssbo_vert_flags_[vid] = encode_vert_flags(vf); 
    }

void update_vert_flags_selected(uint vid, uint selection_slot)
{
    VertFlags vf = load_vert_flags(vid); 
    vf.selected[selection_slot] = true; 
    store_vert_flags(vid, vf); 
}
void update_vert_flags__mov_by_collapse(uint vid)
{
    VertFlags vf = load_vert_flags(vid); 
    vf.mov_by_collapse = true;
    store_vert_flags(vid, vf);
}
void update_vert_flags__del_by_collapse(uint vid)
{
    VertFlags vf = load_vert_flags(vid);
    vf.del_by_collapse = true;
    vf.mov_by_collapse = false; 
    vf.selected = bvec4(false); 
    store_vert_flags(vid, vf); 
}
void update_vert_flags__reset_edge_split_new_vert__mov_by_collapse(uint vid, VertFlags vf_old)
{
    vf_old.new_by_split = false; 
    vf_old.mov_by_collapse = false;
    store_vert_flags(vid, vf_old);
}
void update_vert_flags__reset_edge_split_tags(uint vid, VertFlags vf_old)
{
    vf_old.new_by_split = false; 
    store_vert_flags(vid, vf_old);
}
void update_vert_flags__reset_face_split_new_vert(uint vid, VertFlags vf_old)
{
    vf_old.new_by_face_split = false; 
    store_vert_flags(vid, vf_old);
}
void update_vert_flags__facing_direction(bool front_facing, bool back_facing, inout VertFlags vf_update)
{
    vf_update.front_facing = front_facing; 
    vf_update.back_facing = back_facing; 
}
void update_vert_flags__border(bool border, inout VertFlags vf_update)
{
    vf_update.border_eval = border; 
}
void update_vert_flags__crease(bool crease, inout VertFlags vf_update)
{
    vf_update.crease = crease; 
}
void update_vert_flags__corner(bool corner, inout VertFlags vf_update)
{
    vf_update.corner = corner; 
}
#endif



struct EdgeFlags
{
    // Persistent Flags ----------------------------------
    bool dupli; // duplicated edge, should be ignored
    bool border; // border edge
    bool selected; // selected for a certain operation
    bool sel_border; // at the selection border
    bool new_by_split; // new edge created by edge split
    bool new_by_split_on_old_edge; // new edge created by edge split, and it's on the old edge
    bool del_by_split; // deleted edge by edge split
    bool new_by_face_split; // new edge created by edge split
    bool del_by_collapse; // deleted edge by edge collapse
    uint crease_level; // crease level, 2 bits 0~3

    // Temp Flags, can share bits across different passes -
    bool temp_face_split_new_edge; // flipped edge
    bool temp_dbg_draw_edge; // whatever you want to debug draw
}; 

uint encode_edge_flags(EdgeFlags ef)
{
    uint ef_enc = 0u; 
    ef_enc |= uint(ef.dupli); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.border); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.selected); 
    ef_enc <<= 1u;
    ef_enc |= uint(ef.sel_border); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.new_by_split); 
    ef_enc <<= 1u;
    ef_enc |= uint(ef.new_by_split_on_old_edge); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.del_by_split); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.new_by_face_split); 
    ef_enc <<= 1u;
    ef_enc |= uint(ef.del_by_collapse); 
    ef_enc <<= 2u;
    ef_enc |= (ef.crease_level & 0x3u); 
    ef_enc <<= 1u; 
    ef_enc |= uint(ef.temp_face_split_new_edge); 
    ef_enc <<= 1u;
    ef_enc |= uint(ef.temp_dbg_draw_edge); 

    return ef_enc; 
}
EdgeFlags decode_edge_flags(uint ef_enc)
{
    EdgeFlags ef; 
    ef.temp_dbg_draw_edge = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.temp_face_split_new_edge = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.crease_level = (ef_enc & 0x3u);
    ef_enc >>= 2u; 
    ef.del_by_collapse = (1u == (ef_enc & 1u)); 
    ef_enc >>= 1u; 
    ef.new_by_face_split = (1u == (ef_enc & 1u)); 
    ef_enc >>= 1u; 
    ef.del_by_split = (1u == (ef_enc & 1u)); 
    ef_enc >>= 1u; 
    ef.new_by_split_on_old_edge = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.new_by_split = (1u == (ef_enc & 1u)); 
    ef_enc >>= 1u; 
    ef.sel_border = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.selected = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.border = (1u == (ef_enc & 1u));
    ef_enc >>= 1u; 
    ef.dupli = (1u == (ef_enc & 1u)); 

    return ef; 
}

#if defined(EDGE_FLAGS_INCLUDED)
EdgeFlags load_edge_flags(uint wedge_id)
{
    return decode_edge_flags(ssbo_edge_flags_[wedge_id]); 
}
void store_edge_flags(uint wedge_id, EdgeFlags ef)
{
    ssbo_edge_flags_[wedge_id] = encode_edge_flags(ef); 
}

EdgeFlags init_edge_flags(bool dupli, bool border)
{
    EdgeFlags ef; 
    ef.dupli = dupli; 
    ef.border = border; 
    ef.selected = false; 
    ef.sel_border = false; 
    ef.new_by_split = false; 
    ef.new_by_split_on_old_edge = false; 
    ef.del_by_split = false;
    ef.new_by_face_split = false; 
    ef.del_by_collapse = false; 
    ef.crease_level = 0u; 

    ef.temp_face_split_new_edge = false; 
    ef.temp_dbg_draw_edge = false; 

    return ef; 
}

void update_edge_flags__detect_sel_border(uint edge_id, EdgeFlags ef_old)
{
    EdgeFlags ef = ef_old; 
    ef.sel_border = true; 
    store_edge_flags(edge_id, ef); 
}

void update_edge_flags__detect_crease(uint edge_id, EdgeFlags ef_old, uint crease_level)
{
    EdgeFlags ef = ef_old; 
    ef.crease_level = crease_level; 
    store_edge_flags(edge_id, ef); 
}

void update_edge_flags__reset_face_split_new_edge(uint edge_id, EdgeFlags ef_old)
{
    EdgeFlags ef = ef_old; 
    ef.temp_face_split_new_edge = false; 
    store_edge_flags(edge_id, ef);
}

void update_edge_flags__reset_new_split_edge(uint edge_id, EdgeFlags ef_old)
{
    EdgeFlags ef = ef_old; 
    ef.new_by_split = false; 
    ef.new_by_split_on_old_edge = false; 
    store_edge_flags(edge_id, ef);
}


EdgeFlags init_edge_flags__new_split_edge(bool is_on_old_edge)
{
    EdgeFlags ef; 
    ef.dupli = false; 
    ef.border = false; 
    ef.selected = true;
    ef.sel_border = false; 
    ef.new_by_split = true; 
    ef.new_by_split_on_old_edge = is_on_old_edge;
    ef.del_by_split = false; 
    ef.new_by_face_split = false; 
    ef.del_by_collapse = false; 
    ef.crease_level = 0u; 

    ef.temp_face_split_new_edge = false; 
    ef.temp_dbg_draw_edge = false;

    return ef; 
}
void update_edge_flags__del_by_split(uint edge_id)
{
    EdgeFlags ef_del = load_edge_flags(edge_id);
    ef_del.del_by_split = true; 
    ef_del.selected = false; 
    ef_del.sel_border = false; 
    store_edge_flags(edge_id, ef_del);
}

void update_edge_flags__del_by_collapse(uint edge_id)
{
    EdgeFlags ef = load_edge_flags(edge_id); 
    ef.del_by_collapse = true; 
    ef.selected = false; 
    ef.sel_border = false; 
    store_edge_flags(edge_id, ef); 
}

EdgeFlags init_edge_flags__new_face_split_edge()
{
    EdgeFlags ef; 
    ef.dupli = false; 
    ef.border = false; 
    ef.selected = true;
    ef.sel_border = false; 
    ef.new_by_split = false; 
    ef.new_by_split_on_old_edge = false;
    ef.del_by_split = false; 
    ef.new_by_face_split = true; 
    ef.del_by_collapse = false; 
    ef.crease_level = 0u; 

    ef.temp_face_split_new_edge = true; 
    ef.temp_dbg_draw_edge = false;

    return ef; 
}
#endif



#if defined(VE_CIRCULATOR_INCLUDE)

/*      
 * Rotation winding (CW/CCW) is always reverse to global face winding
 * Rotate foward/backward: hedge w always end/begin at v 
 * 
 *    Rotate fwd around V3           Rotate fwd around V1          fi:awi.iface_adj  
 *        v0  ...  vp             v_oppo --- vn  ...  v0           wi:awi.wedge_id  
 *       /  \     /  \                \     /  \     /  \          
 *      / f1 \   wp   \                \  wo fi wn  /    \         When rotating fowards:            
 *     / ---> \ /      \                \ / ---> \ / ---> \        wo = ssbo_edge_to_edges[mark__cwedge_rotate_back(fi)]
 *   v1 ====== v3--wi--vi               vi --wi-- v1 ===== v3      wn = ssbo_edge_to_edges[mark__cwedge_rotate_next(fi)]      
 *     \ <--- / \<-----/ \                \      / \ <--- /        wp can be cached
 *      \ f0 /  wn fi wo  \                \   wp   \    /         
 *       \  /     \  /     \                \  /     \  /          
 *        v2  ...  vn ----- v_oppo           vp  ...  v2           
*/  
struct CirculatorIterData
{
    AdjWedgeInfo awi; // wi
    AdjWedgeInfo awi_next; // wn
    uint rotate_step; 
}; 
#define MAX_WEDGE_ROTATES 36u


/* circulation loop invariant: 
 *     v1 ----- v2  
 *    /  \     /  \     wi:=awi.wedge_id
 *   /    \   /    \    fi:=awi.iface_adj
 *  /      \ /      \  
 * v0 ----- v --wi-- v3 iwedge derivation:       
 *   \__     \<-----/   wo = mark__cwedge_rotate_back(fi)
 *      \__  wn fi wo   wn = mark__cwedge_rotate_next(fi) 
 *         \__ \  /    
 *            \_v4             
*/
/* Parameters
 VertWedgeListHeader vwlh
 bool FUNC(CirculatorIterData, VISIT_DATA)
 CustomData VISIT_DATA
 bool ROTATE_FWD // true then wi always points at v at fi (see graph above), false when wi points from v at fi
*/
#define VE_CIRCULATOR(vwlh, FUNC, VISIT_DATA, ROTATE_FWD) \
{                                                                                                  \
    uint iface_rotate_fwd = (vwlh.ivert == 1u) ? 0 : 1;                                            \
    uint iface_rotate_bcw = (vwlh.ivert == 1u) ? 1 : 0;                                            \
    uint iface_curr = (ROTATE_FWD) ? iface_rotate_fwd : iface_rotate_bcw;                          \
                                                                                                   \
    AdjWedgeInfo awi;                                                                              \
    awi.wedge_id = vwlh.wedge_id; /* current wedge we are traversing */                            \
    awi.iface_adj = iface_curr; /* current iface we are working at */                              \
    uint rotate_step = 0u;                                                                                                 \
    do {                                                                                                                   \
        /* find next wedge */                                                                                              \
        uint iwedge_next = ROTATE_FWD ? mark__cwedge_rotate_next(awi.iface_adj) : mark__cwedge_rotate_back(awi.iface_adj); \
        AdjWedgeInfo awi_next = decode_adj_wedge_info(ssbo_edge_to_edges_[awi.wedge_id*4u + iwedge_next]);                 \
        bool awi_next_border = awi.wedge_id == awi_next.wedge_id;                                                          \
                                                                                                                           \
        /* do something here */                                                                                            \
        if (false == FUNC(CirculatorIterData(awi, awi_next, rotate_step), VISIT_DATA))                                     \
            break;                                                                                                         \
                                                                                                                           \
        if (awi.wedge_id == awi_next.wedge_id)                                                                             \
            break; /* hit the border edge, exit */                                                                         \
        awi = awi_next;                                                                                                    \
        rotate_step++;                                                                                                     \
    }  while (                                                                                                             \
        rotate_step < MAX_WEDGE_ROTATES                                                                                    \
        && awi.wedge_id != vwlh.wedge_id                                                                                   \
    );                                                                                                                     \
}\

/* circulation loop invariant: (rotating foward)
 *     v1 ----- vp  
 *    /  \     /  \     wi:=awi.wedge_id
 *   /    \  wp   wop   fi:=awi.iface_adj
 *  /      \ /      \  
 * v0 ----- v --wi-- vi iwedge derivation:       
 *   \__     \<-----/   wo = mark__cwedge_rotate_back(fi)
 *      \__  wn fi wo   wn = mark__cwedge_rotate_next(fi) 
 *         \__ \  /    
 *            \_vn             
*/
/* get ivert in the context of wi 
 * for example, vi = ssbo_edge_to_vert_[wi*4u + ivert_vi] */
uint mark__ve_circ_fwd__get_vi(CirculatorIterData iter) { return mark__cwedge_to_beg_vert(iter.awi.iface_adj); }
uint mark__ve_circ_fwd__get_vn(CirculatorIterData iter) { return mark__center_wedge_to_oppo_vert__at_face(iter.awi.iface_adj); }
uint mark__ve_circ_fwd__get_vp(CirculatorIterData iter) { return mark__center_wedge_to_oppo_vert__at_face(iter.awi.iface_adj == 1u ? 0u : 1u); }

uint mark__ve_circ_bck__get_vn(CirculatorIterData iter) { return mark__ve_circ_fwd__get_vn(iter); }

uint mark__vecirc_fwd_get_wo(CirculatorIterData iter) { return mark__cwedge_rotate_back(iter.awi.iface_adj); }
uint mark__vecirc_fwd_get_wop(CirculatorIterData iter) { 
    uint iwedge_wo = mark__vecirc_fwd_get_wo(iter); 
    return mark__bwedge_to_next_bwedge(iwedge_wo); 
}
#endif





#if defined(TOPO_DIAGONOSIS_INCLUDE)

/* Diagnose mesh topology */
void validate_wedge_topo(uint wedge_id, out bool valid_ee, out bool valid_ev, out bool valid_ve)
{
    valid_ee = valid_ev = valid_ve = true; 

    AdjWedgeInfo w[4] = {
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 0u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 1u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 2u]),
        decode_adj_wedge_info(ssbo_edge_to_edges_[wedge_id*4u + 3u])
    }; 

    uvec4 v = uvec4(
        ssbo_edge_to_vert_[wedge_id*4u + 0u], 
        ssbo_edge_to_vert_[wedge_id*4u + 1u],
        ssbo_edge_to_vert_[wedge_id*4u + 2u],
        ssbo_edge_to_vert_[wedge_id*4u + 3u]
    ); 
    
    bool border_edge = false; 
    uvec2 iwedges_border = mark__border_bwedges_mainfold(); 

#define w0 ((w[0].wedge_id))
#define w1 ((w[1].wedge_id))
#define w2 ((w[2].wedge_id))
#define w3 ((w[3].wedge_id))
    
    if (w0 == w2 || w1 == w3 || w0 == w3 || w1 == w2) 
        valid_ee = false; 

    for (uint i = 0u; i < 4u; ++i)
    {
        if (border_edge && all(i.xx != iwedges_border.xy))
            continue; 

        uint wi = w[i].wedge_id; 
        uint iface_overlap_at_wi = w[i].iface_adj == 1u ? 0u : 1u; 

        /* Check current wedge's position in neighbor wedge's ee links */
        bool is_cwedge_next = ((i == 0u) || (i == 2u)); /* cwedge == wi.next */
        uint iwedge_at_wi = is_cwedge_next ? 
            mark__cwedge_rotate_next(iface_overlap_at_wi) : 
            mark__cwedge_rotate_back(iface_overlap_at_wi); 

        AdjWedgeInfo reflex = decode_adj_wedge_info(ssbo_edge_to_edges_[wi*4u + iwedge_at_wi]);
        if (reflex.wedge_id != wedge_id)
            valid_ee = false; 

        /* Check ev link coherence */
        uvec2 iverts_wi = mark__wedge_to_verts(i); 
        uvec2 verts_wi = uvec2(v[iverts_wi.x], v[iverts_wi.y]); 
        uvec2 iverts_wi_at_wi = mark__cwedge_to_verts(iface_overlap_at_wi); 
        uvec2 verts_wi_at_wi = uvec2(
            ssbo_edge_to_vert_[wi*4u + iverts_wi_at_wi[0]], 
            ssbo_edge_to_vert_[wi*4u + iverts_wi_at_wi[1]]
        ); 
        if (any(bvec2(verts_wi != verts_wi_at_wi)))
            valid_ev = false;

        uint v_oppo = v[mark__border_wedge_to_oppo_vert(i)]; 
        uint v_oppo_at_wi = ssbo_edge_to_vert_[wi*4u + mark__center_wedge_to_oppo_vert__at_face(iface_overlap_at_wi)];
        if (v_oppo != v_oppo_at_wi)
            valid_ev = false;

        EdgeFlags ef = load_edge_flags(wi); 
        if (ef.del_by_collapse || ef.del_by_split)
            valid_ee = false;
    }


    /* check ve linkage */
    uvec2 iverts = mark__cwedge_to_verts(0u);
    for (uint iivert = 0u; iivert < 2u; ++iivert)
    {
        uint ivert = iverts[iivert]; 
        uint vid = v[ivert]; 

        VertWedgeListHeader vwlh = decode_vert_wedge_list_header(ssbo_vert_to_edge_list_header_[vid]); 
        EdgeFlags ef = load_edge_flags(vwlh.wedge_id); 
        if (ef.del_by_collapse || ef.del_by_split)
            valid_ve = false; 

        uint vid_from_ev_adj = ssbo_edge_to_vert_[vwlh.wedge_id*4u + vwlh.ivert];
        if (vid_from_ev_adj != vid)
            valid_ve = false;
    } 

#undef w0
#undef w1
#undef w2
#undef w3
}

#endif





float gaussian(float dist, float tau)
{
    return exp(-(dist * dist) / (2.0f * tau * tau)); 
}

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


    /* Bilateral weight from "Fast and Effective Feature-Preserving Mesh Denoising" */
    float normal_weight = dot(np, nq);
    if (normal_weight < 0.9f) normal_weight = 0.0f; 
    normal_weight *= normal_weight; 
    weight = normal_weight * plane_weight/*  * geometry_weight */; 
    

    /* Test for view-aligned filtering --------------- */
    vec3 vq = normalize(cam_pos_ws - qq); 
    float cos_theta = dot(nq, vq); 
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta); 
    
    vec3 vp = normalize(cam_pos_ws - pp); 
    vec3 np_view_proj = normalize(np - dot(np, vp) * vp);
    vec3 np_view_filtered = normalize(sin_theta * np_view_proj + cos_theta * vp); 


    n_filtered = nq; 
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
        silouetteness = abs(ndv0 - ndv2) * .5f; 

    vec3 pe = (p1 + p3) * .5f; 
    q.quadric = compute_filter_quadric(
        wedge_normal, pe, cam_pos_ws, silouetteness, position_regularize_weight
    );

	return q; 
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

    float geometry_weight = 1.0f; // q_v.area * gaussian(dist, dev_g); 
    
    return geometry_weight * quadric_weight; 
}

#endif



/* Remeshing ----------------------------------------------------------------------------------- */
/* New vert pos for split/collapse: 
 * interp using Butterfly Scheme, see "Robust Topological Operations for Dynamic Explicit Surfaces"
 *                                  
 *  q0 ----E00----- v0 -----E31---- q3         q0 ----E00----- v0 ----E31---- q3
 *    \  <-------  /  \  <-------  /              \_           |           _/   
 *     \          /    \          0                 \_         |         _/     
 *      E        W      3        3                    \_       |       _/       
 *       0      0   f1   W      E                       \_     |     _/         
 *        1    /          \    /                          \_   |   _/           
 *         \  /  ------->  \  /                             \_ | _/            
 *          v1 ---- W4 ---- v3                                \v1 = (v1+v3)/2 + (v0+v2)/8 - (q0+q1+q3+q4)/16               
 *         /  \  <-------  /  \                             _/ | \_  
 *        E    \          /    \                          _/   |   \_                       
 *       1      W   f0   2      1                       _/     |     \_                     
 *      0        1      W        2                    _/       |       \_                   
 *     /          \    /          E                 _/         |         \_                 
 *    /  ------->  \  /  ------->  \              _/           |           \_               
 *  q1 ----E11----- v2 ----E20----- q2         q1 ----E11----- v2 ----E20---- q2
*/
vec3 interp_mid_point(vec3 v0, vec3 v1, vec3 v2, vec3 v3, vec3 q0, vec3 q1, vec3 q2, vec3 q3)
{
    return (v1+v3)/2.0f + (v0+v2)/8.0f - (q0+q1+q2+q3)/16.0f;  
}

/* Remeshing threshold, 
 * see Ch.6.5.3 "Incremental Remeshing" in
 * http://staff.ustc.edu.cn/~lgliu/Courses/DGP_2014_autumn-winter/References/Book_Polygon%20Mesh%20Processing.pdf
*/
float calc_remesh_edge_len_min(float targ_edge_len) // edge shorter than this will be collapsed
{
    return (4.0f/5.0f) * targ_edge_len; 
}
float calc_remesh_edge_len_max(float targ_edge_len) // edge longer than this will be split
{
    return (4.0f/3.0f) * targ_edge_len; 
}
/* Adaptive Remeshing Length 
 * see "Adaptive Remeshing for Real-Time Mesh Deformation" 
*/
#define EPSI_ADAPTIVE_REMESH_LEN 1e-3f
float get_adaptive_remesh_len(float k/*max_curvature*/, float ref_edge_len = -1.0f)
{
    const float epsi = EPSI_ADAPTIVE_REMESH_LEN; 

    if (ref_edge_len > .0f)
    { // >.0f: use ref_edge_len as the target edge length 
        float targ_edge_len_min = ref_edge_len * .125f; // scales according the how many split/collapse we want 
        float targ_edge_len_max = ref_edge_len * 128.0f;
        if (k < 1e-10f) return targ_edge_len_max; 
    }

    float l = epsi * ((6.0f / k) - (3.0f * epsi)); 
    l = sqrt(max(.0f, l)); 
    // if (ref_edge_len > .0f)
    //     l = clamp(l, targ_edge_len_min, targ_edge_len_max);

    return l; 
}

float average_adaptive_remesh_len_to_vcurv(float l)
{
    const float epsi = EPSI_ADAPTIVE_REMESH_LEN; 
    float k = epsi * (6.0f / (l*l + 3.0f*epsi*epsi));  
    return k;  
}

#endif




