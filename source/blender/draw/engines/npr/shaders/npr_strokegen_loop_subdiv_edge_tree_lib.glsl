
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)


#ifndef STROKEGEN_LOOP_SUBDIV_EDGE_TREE_LIB__INCLUDED
#define STROKEGEN_LOOP_SUBDIV_EDGE_TREE_LIB__INCLUDED

/* Loop subdivision tree -------------------------------------------------------------- 
 * Encode sparse signals on a densely subdivided surface. 
 * 
 *                  v0                    Inspired from the "Edge-Tree" encoding in   
 *                 /  \                   "Progressive Geometry Compression".
 *                /    \                  
 *               /      \                 Take edge (v1, v3) as example, 
 *              /        \                After 1 level of loop subd, 
 *             /          \               There is a tree structure between the old edge and its 4 sub-edges:
 *            /            \                                       ... higher nodes
 *          v01= FaceEdge =v30                                     |
 *         /  \            /  \                        .___ old edge (v1, v3) ___.           
 *        /    \          /    \                      /        |            |     \          
 *       /      \        /      \                 (v01,v30),(v12,v23)    (v1,v13),(v13,v3)   
 *      /        \      /        \                      Face Edges           Split Edges     
 *     /          \    /          \   
 *    /            \  /            \      Temporal Consistency -----------------------------------------------------------      
 *  v1= SplitEdge = v13= SplitEdge =v3    The subdivision is non-deterministic.
 *    \            /  \            /      Subdivided edge ids changes in different frames,       
 *     \          /    \          /       We cannot afford an explicit storage for all of them.       
 *      \        /      \        /        
 *       \      /        \      /         We rely on the base edge id and tree code to infer the sub-edge ids.
 *        \    /          \    /                
 *         \  /            \  /           note: base edges are consistent, except for subdivided ones.      
 *          v12= FaceEdge = v23                  
 *            \            /              "subdive face id" is used to determine the code for the two new edges.      
 *             \          /          sf0 <- The face with highest edge id among four adj. edges around (v1,v3).       
 *              \        /           sf1 <- Another face, if there were any.       
 *               \      /                       
 *                \    /                  Use the tag and vertex winding we can number the 4 sub-edges      
 *                 \  /                   The tag is also inherited to the sub-edges for them to number their sub-edges.      
 *                  v2
*/

/*  An example of tree encoding for edge AB after 1 level of loop subd:
 *   
 *           C                     F            
 *          / \                   / \        the quadrant is split into 4 sub-quadrants AGCD,DHBE,DCEF,GIHC
 *         /   \                 /sf0\       The "subdivision face id"(sf) is inherited from parent to child quadrants
 *        /     \               D =2= E                                                           
 *       / sf0   \             / \sf1/ \     The highest bit is 1/0 for face-edge/split-edge, respectively
 *      / ------> \           /sf0\ /sf0\    Split-edges AC and CB has code 00 and 01 (low bit takes winding order from parent quad's sf0)
 *     A  ======== B         A =0= C =1= B   Face-edges  DE and GH has code 10 and 11 (low bit takes which sf this edge is in)
 *      \ <-----  /           \sf1/ \sf1/                                                         
 *       \ sf1   /             \ /sf0\ /                                                          
 *        \     /               G =3= H                                                           
 *         \   /                 \sf1/                                                            
 *          \ /                   \ /                                                             
 *           D                     I                                                              
 *                                                                 
*/
struct LoopSubdEdgeTreeUpNode
{
    uint parent_edge_id; // 29 bits
#define LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID 0x1fffffffu
    uint code; // 2 bits
    uint iface_subd_f0; // 1bit, iface of face tagged as f0 (see explanation above) 
};
uint encode_loop_subd_tree_node(LoopSubdEdgeTreeUpNode edge_tree_ptr)
{
    return ((edge_tree_ptr.parent_edge_id << 3u) | (edge_tree_ptr.code << 1u) | (edge_tree_ptr.iface_subd_f0 & 1u));  
}
LoopSubdEdgeTreeUpNode decode_loop_subd_tree_node(uint enc)
{
    LoopSubdEdgeTreeUpNode ptr; 
    ptr.iface_subd_f0 = enc & 0x01u; 
    enc >>= 1u; 
    ptr.code = enc & 0x03u;
    enc >>= 2u;  
    ptr.parent_edge_id = enc; 
    return ptr; 
}

#if defined(WINGED_EDGE_TOPO_INCLUDE)
LoopSubdEdgeTreeUpNode init_loop_subd_tree_root(uint wedge_id, uvec4 bwedges)
{
    LoopSubdEdgeTreeUpNode node; 
    node.parent_edge_id = wedge_id;
    node.code = 0u;
    uvec2 ibwedges_at_iface[2];
    ibwedges_at_iface[0] = mark__face_to_bwedges(0u); 
    ibwedges_at_iface[1] = mark__face_to_bwedges(1u); 
    
    uvec2 max_bwedge_id_at_iface;
    for (uint iface = 0; iface < 2u; ++iface)
        max_bwedge_id_at_iface[iface] = max(
            bwedges[ibwedges_at_iface[iface].x], 
            bwedges[ibwedges_at_iface[iface].y]
        ); 
    bool is_iface_0_max = max_bwedge_id_at_iface[0] > max_bwedge_id_at_iface[1]; 

    node.iface_subd_f0 = is_iface_0_max ? 0u : 1u; 
    return node; 
}
#endif

LoopSubdEdgeTreeUpNode init_loop_subd_tree_leaf__face_edge()
{
    LoopSubdEdgeTreeUpNode node; 
    node.parent_edge_id = LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID;
    node.code = 0u;
    node.iface_subd_f0 = 0u; 
    return node; 
}
uint get_loop_subd_tree_leaf_code__face_edge(uint subd_face_id)
{ // See the diagram above. The code for two new edges are determined by the iface they reside
    return (2u | (subd_face_id & 1u)); 
}


LoopSubdEdgeTreeUpNode setup_loop_subd_tree_leaf__split_edge(
    uint par_edge_id, LoopSubdEdgeTreeUpNode par_node, 
    uint split_edge_mark/* "E1 or E3", defined in edge-split diagram*/
){
    LoopSubdEdgeTreeUpNode node; 
    node.parent_edge_id = par_edge_id; 
    node.iface_subd_f0 = par_node.iface_subd_f0;  
    if (par_node.iface_subd_f0 == 0u)
        node.code = split_edge_mark == 3u ? 0u : 1u;
    else // par_node.iface_subd_f0 == 1u
        node.code = split_edge_mark == 1u ? 0u : 1u;

    return node; 
}


struct LoopSubdEdgeTreeDwNode
{
    uint wedge_id; // 32 bits
};
uint encode_loop_subd_tree_node_dw(LoopSubdEdgeTreeDwNode edge_tree_ptr)
{
    return edge_tree_ptr.wedge_id; 
}
LoopSubdEdgeTreeDwNode decode_loop_subd_tree_node_dw(uint enc)
{
    LoopSubdEdgeTreeDwNode ptr; 
    ptr.wedge_id = enc; 
    return ptr; 
}


#endif


