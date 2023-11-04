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
    return line_adj_verts.xywz; /* 0123 -> 0132 */
}
/* Vertices: x4 vert ids per wing */
/*    v2
 *   /  \
 *  /    \                
 * v1>-->v3 for wedge e, vi => ssbo_edge_to_verts_[e*4+i] , i=0,1,2,3            
 *  \    /  
 *   \  /   Face winding               
 *    v0    013, 123                                  
*/

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
#define NULL_EDGE 0xffffffffu /* non-existing wedges, for example at the mesh border */

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

/* contour-specific functions */
struct PerContourWedgeInfo
{
    uint wedge_id;
    uint ifrontface;
};
uint encode_per_contour_wedge_info(PerContourWedgeInfo pcwi)
{
    uint pcwi_enc = ((pcwi.wedge_id << 1u) | (pcwi.ifrontface & 1u));
    return pcwi_enc; 
}
PerContourWedgeInfo decode_per_contour_wedge_info(uint pcwi_enc)
{
    PerContourWedgeInfo pcwi; 
    pcwi.wedge_id = (pcwi_enc >> 1u);
    pcwi.ifrontface = ((pcwi_enc & 1u));
    return pcwi; 
}

struct PerWedgeContourInfo
{
    bool is_contour;
    uint contour_id;
    uint ifrontface;
};
uint encode_per_wedge_contour_info(PerWedgeContourInfo peci)
{
    uint peci_enc = (peci.contour_id << 2u) | ((peci.ifrontface & 1u) << 1u) | (uint(peci.is_contour)); 
    return peci_enc; 
}
PerWedgeContourInfo decode_per_wedge_contour_info(uint peci_enc)
{
    PerWedgeContourInfo peci; 
    peci.contour_id = (peci_enc >> 2u);
    peci.ifrontface = ((peci_enc >> 1u) & 1u);
    peci.is_contour = (0u != (peci_enc & 1u));
    return peci; 
}



#endif













#endif




