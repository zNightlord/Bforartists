#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)

uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING)
/*
 * uint pcs_edge_id_offset_; 
 * uint pcs_edge_count_; 
 * uint ssbo_wedge_flooding_pointers_in_[]; 
 * uint ssbo_wedge_flooding_pointers_out_[];
 * uint ssbo_edge_to_edges_[]; 
*/
void main()
{
    const uint groupId = gl_LocalInvocationID.x; 
	const uint idx = gl_GlobalInvocationID.x; 
    /* Do not use EdgeID here since it's offseted with current mesh batch */
    bool valid_thread = (idx.x < pcs_edge_count_); 
    
    const uint EdgeID = idx.x + pcs_edge_id_offset_; 


#if defined(_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER)
    WedgeFloodingPointer wfptr = decode_wedge_flooding_pointer(ssbo_wedge_flooding_pointers_in_[EdgeID]); 
    if (wfptr.is_seed) return; /* already a seed, quit */ 

    /* Now try to find a seed */
    AdjWedgeInfo awi[4]; 
    WedgeFloodingPointer wfptr_bwedge[4]; 
    for (int ibwedge = 0; ibwedge < 4; ++ibwedge)
    {
        awi[ibwedge]          = decode_adj_wedge_info(
            ssbo_edge_to_edges_[EdgeID*4 + ibwedge]
        ); 
        wfptr_bwedge[ibwedge] = decode_wedge_flooding_pointer(
            ssbo_wedge_flooding_pointers_in_[awi[ibwedge].wedge_id]
        );

        if (wfptr_bwedge[ibwedge].is_seed)
        {
            wfptr.next_wedge_id = wfptr_bwedge[ibwedge].next_wedge_id; 
            wfptr.is_seed = true; 
            break; /* found a seed, quit */
        }
    }

    if (false == wfptr.is_seed)
    { /* no seed is found in adjacent wedges, jump with the long range ptr  */
        WedgeFloodingPointer wfptr_next = decode_wedge_flooding_pointer(
            ssbo_wedge_flooding_pointers_in_[wfptr.next_wedge_id]
        );
        if (wfptr_next.next_wedge_id != EdgeID)
        {
            wfptr.next_wedge_id = wfptr_next.next_wedge_id;
            wfptr.is_seed = wfptr_next.is_seed; 
        }else
        { /* long range ptr is pointing to this wedge */
            uint random_select = pcg(EdgeID) % 4u; 
            for (uint roll = 0; roll < 4u; ++ roll)
            {
                uint ibwedge = (random_select + roll) % 4u; 
                if (wfptr_bwedge[ibwedge].next_wedge_id != EdgeID)
                {
                    wfptr.next_wedge_id = wfptr_bwedge[ibwedge].next_wedge_id;
                    wfptr.is_seed = wfptr_bwedge[ibwedge].is_seed;
                    break; 
                }
            }
        }
    }

    if (valid_thread)
        ssbo_wedge_flooding_pointers_out_[EdgeID] = encode_wedge_flooding_pointer(wfptr); 
#endif




}
#endif