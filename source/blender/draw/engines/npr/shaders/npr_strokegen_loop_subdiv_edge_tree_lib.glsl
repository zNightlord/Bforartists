
#ifndef STROKEGEN_LOOP_SUBDIV_EDGE_TREE_LIB__INCLUDED
#define STROKEGEN_LOOP_SUBDIV_EDGE_TREE_LIB__INCLUDED

/* Loop subdivision tree -------------------------------------------------------------- 
 *                  v0                    Inspired from the "Edge-Tree" encoding in   
 *                 /  \                   "Progressive Geometry Compression".
 *                /    \                  
 *               /      \                 Take edge (v1, v3) as example, 
 *              /        \                After 1 level of loop subd, 
 *             /          \               new edges (v12, v23) (v1, v13) (v13, v3) (v01, v30) are parallel to (v1, v3)
 *            /            \              We link them to (v1, v3) with code 0, 1, 2, 3. 
 *          v01==== E0 ====v30        
 *         /  \            /  \       f1  This makes a tree structure, where every old edge is pointed by 4 new edges.
 *        /    \          /    \      
 *       /      \        /      \     
 *      /        \      /        \    
 *     /          \    /          \   
 *    / ---------> \  / ---------> \  
 *  v1 ==== E1 ==== v13 === E3 ==== v3
 *    \ <--------- /  \ <--------- /
 *     \          /    \          /
 *      \        /      \        /
 *       \      /        \      /     f0
 *        \    /          \    /
 *         \  /            \  /
 *          v12==== E2 ==== v23
 *            \            /
 *             \          /
 *              \        /
 *               \      /
 *                \    /
 *                 \  /
 *                  v2
*/
struct LoopSubdEdgeTreeNode
{
    uint parent_edge_id; // 30 bits
    uint code; // 2 bits 
};
uint encode_loop_subd_tree_ptr(LoopSubdEdgeTreeNode edge_tree_ptr)
{
    return ((edge_tree_ptr.parent_edge_id << 2u) | edge_tree_ptr.code);  
}
LoopSubdEdgeTreeNode decode_loop_subd_tree_ptr_edge_id(uint enc)
{
    LoopSubdEdgeTreeNode ptr; 
    ptr.parent_edge_id = enc >> 2u; 
    ptr.code = enc & 0x03u; 
    return ptr; 
}
LoopSubdEdgeTreeNode init_loop_subd_tree_root(uint edge_id)
{
    LoopSubdEdgeTreeNode node; 
    node.parent_edge_id = edge_id;
    node.code = 0u;
    return node; 
}
LoopSubdEdgeTreeNode init_loop_subd_tree_leaf()
{
#define LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID 0x3fffffffu
    LoopSubdEdgeTreeNode node; 
    node.parent_edge_id = LOOP_SUBD_TREE_INVALID_PARENT_EDGE_ID;
    node.code = 0u;
    return node; 
}
uint get_loop_subd_tree_leaf_code(uint iface)
{ // See the diagram above. The code for two new edges are determined by the iface they reside
    return (iface == 0u) ? 2u : 0u; 
}
LoopSubdEdgeTreeNode setup_loop_subd_tree_node__edge_split_at_parent(uint edge_id, uint code)
{
    LoopSubdEdgeTreeNode node; 
    node.parent_edge_id = edge_id; 
    node.code = code; 
    return node; 
}

#endif


