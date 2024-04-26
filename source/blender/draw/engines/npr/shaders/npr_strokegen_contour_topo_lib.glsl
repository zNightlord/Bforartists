#ifndef BNPR_CONTOUR_TOPO__INCLUDED
#define BNPR_CONTOUR_TOPO__INCLUDED

struct ContourFlags
{
    bool seg_head; 
    bool seg_tail; 
    bool looped_curve; 
}; 

uint encode_contour_flags(ContourFlags cf)
{
    uint cf_enc = 0u; 
    cf_enc |= uint(cf.seg_head); 
    cf_enc <<= 1u; 
    cf_enc |= uint(cf.seg_tail); 
    cf_enc <<= 1u;
    cf_enc |= uint(cf.looped_curve);

    return cf_enc; 
}

ContourFlags decode_contour_flags(uint cf_enc)
{
    ContourFlags cf; 
    cf.looped_curve = (1u == (cf_enc & 1u));
    cf_enc >>= 1u;
    cf.seg_tail = (1u == (cf_enc & 1u)); 
    cf_enc >>= 1u; 
    cf.seg_head = (1u == (cf_enc & 1u)); 

    return cf; 
}

ContourFlags init_contour_flags(bool seg_head)
{
    ContourFlags cf;
    cf.seg_head = seg_head;
    cf.seg_tail = false;
    cf.looped_curve = true; 
    return cf; 
}

#if defined(INCLUDE_CONTOUR_FLAGS_LOAD_STORE)
void store_contour_flags(uint contour_id, ContourFlags cf)
{
    ssbo_contour_edge_flags_[contour_id] = encode_contour_flags(cf); 
}
ContourFlags load_contour_flags(uint contour_id)
{
    return decode_contour_flags(ssbo_contour_edge_flags_[contour_id]); 
}
#endif



// Covers full knowledge of edge-loop that current edge lies in
struct ContourCurveTopo
{
    bool looped_curve; 
	uint head_contour_id;
	uint tail_contour_id;
	uint len;
};

uint move_elem_along_loop(uint curr_elem_id, int offset, uint loop_head_elem_id, uint loop_len)
{
	bool move_left = offset < 0;

	uint d = uint(abs(offset));
	d = d % loop_len;

	curr_elem_id -= loop_head_elem_id;
	curr_elem_id += (move_left ? (loop_len - d) : d);
	curr_elem_id = (curr_elem_id % loop_len);
	curr_elem_id += loop_head_elem_id;

	return curr_elem_id;
}
uint move_contour_id_along_loop(ContourCurveTopo curve_topo, uint start_contour_id, float offset)
{
	return move_elem_along_loop(
		start_contour_id, (int(offset + 1e-10f)), 
		curve_topo.head_contour_id, curve_topo.len
	);
}
#if defined(INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD)
ContourCurveTopo load_contour_curve_topo(uint contour_id, ContourFlags cf)
{
    ContourCurveTopo cct; 
    cct.looped_curve    = cf.looped_curve; 
    cct.head_contour_id = ssbo_contour_snake_list_head_[contour_id]; 
    cct.len             = ssbo_contour_snake_list_len_[contour_id]; 
    cct.tail_contour_id = cct.head_contour_id + cct.len - 1u;  
    return cct; 
}
#endif























#endif