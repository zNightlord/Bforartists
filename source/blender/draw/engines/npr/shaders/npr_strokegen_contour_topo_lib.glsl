#ifndef BNPR_CONTOUR_TOPO__INCLUDED
#define BNPR_CONTOUR_TOPO__INCLUDED



#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE__EXTRACT_MESH_CONTOUR_DATA)
/* Rambling between global and local id when extracting contours from curren mesh */
uint calc_local_contour_edge_id(uint global_contour_edge_id)
{
	uint num_contour_edges_prev = ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges;
	return global_contour_edge_id - num_contour_edges_prev; 
}

uint calc_global_contour_edge_id(uint local_contour_edge_id)
{
	uint num_contour_edges_prev = ssbo_bnpr_mesh_pool_counters_prev_.num_contour_edges;
	return local_contour_edge_id + num_contour_edges_prev; 
}
#endif



struct ContourFlags
{
    bool seg_head; 
    bool seg_tail; 
    bool looped_curve; 
    bool cusp_func_pstv; /*has positive cusp function*/
	bool occluded; /* valid when visibility test activated */
}; 

uint encode_contour_flags(ContourFlags cf)
{
    uint cf_enc = 0u; 
    cf_enc |= uint(cf.seg_head); 
    cf_enc <<= 1u; 
    cf_enc |= uint(cf.seg_tail); 
    cf_enc <<= 1u;
    cf_enc |= uint(cf.looped_curve);
    cf_enc <<= 1u;
    cf_enc |= uint(cf.cusp_func_pstv);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.occluded);

    return cf_enc; 
}

ContourFlags decode_contour_flags(uint cf_enc)
{
    ContourFlags cf; 
	cf.occluded = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
    cf.cusp_func_pstv = (1u == (cf_enc & 1u));
    cf_enc >>= 1u;
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
    cf.seg_tail = false;       // needs further setup
    cf.looped_curve = true;    // needs further setup
    cf.cusp_func_pstv = false; // needs further setup
	cf.occluded = false;       // needs further setup
    return cf; 
}

void set_contour_seg_head(bool seg_head, inout ContourFlags cf)
{
    cf.seg_head = seg_head; 
}

void set_contour_seg_tail(bool seg_tail, inout ContourFlags cf)
{
    cf.seg_tail = seg_tail; 
}

void init_contour_looped_curve(bool looped, inout ContourFlags cf)
{
    cf.looped_curve = looped; 
}

void init_contour_cusp_flags(bool cusp_func_pstv, inout ContourFlags cf)
{
    cf.cusp_func_pstv = cusp_func_pstv; 
}

void set_contour_flags_occluded(inout ContourFlags cf)
{
	cf.occluded = true; 
}

#if defined(INCLUDE_CONTOUR_FLAGS_LOAD_STORE)
void store_contour_flags(uint contour_id, ContourFlags cf)
{
    ssbo_contour_snake_flags_[contour_id] = encode_contour_flags(cf); 
}
ContourFlags load_contour_flags(uint contour_id)
{
    return decode_contour_flags(ssbo_contour_snake_flags_[contour_id]); 
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






#if defined(USE_CONTOUR_TRANSFER_DATA_BUFFER)
// Packed per-contour-edge data.
// Used as minimum cache as we iterate through meshes. 
// used for processing all extracted contours after all meshes are processed
struct ContourEdgeTransferData
{
	vec3 vpos_ws[2]; 
	ContourFlags cf; 
	vec2 cusp_funcs; // cusp function at 2 verts
};

void store_contour_edge_transfer_data_(uint contour_edge_id, ContourEdgeTransferData cetd)
{
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+0] = floatBitsToUint(cetd.vpos_ws[0].x);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+1] = floatBitsToUint(cetd.vpos_ws[0].y);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+2] = floatBitsToUint(cetd.vpos_ws[0].z);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+3] = floatBitsToUint(cetd.vpos_ws[1].x);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+4] = floatBitsToUint(cetd.vpos_ws[1].y);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+5] = floatBitsToUint(cetd.vpos_ws[1].z);

	ssbo_contour_edge_transfer_data_[contour_edge_id*8+6] = encode_contour_flags(cetd.cf);
	ssbo_contour_edge_transfer_data_[contour_edge_id*8+7] = packHalf2x16(
		vec2(
			cetd.cusp_funcs[0] > .0f ? 1.0f : -1.0f, 
			cetd.cusp_funcs[1] > .0f ? 1.0f : -1.0f
		)
	);
}

ContourEdgeTransferData load_contour_edge_transfer_data(uint contour_edge_id)
{
	ContourEdgeTransferData cetd; 
	uvec4 enc_data[2];
	Load4(ssbo_contour_edge_transfer_data_, contour_edge_id*2u,    enc_data[0]); 
	Load4(ssbo_contour_edge_transfer_data_, contour_edge_id*2u+1u, enc_data[1]);

	uvec3 vpos_0_enc, vpos_1_enc; 
	vpos_0_enc = enc_data[0].xyz;
	vpos_1_enc = uvec3(enc_data[0].w, enc_data[1].xy);
	cetd.vpos_ws[0] = uintBitsToFloat(vpos_0_enc); 
	cetd.vpos_ws[1] = uintBitsToFloat(vpos_1_enc); 

	cetd.cf = decode_contour_flags(enc_data[1].z); 
	cetd.cusp_funcs = unpackHalf2x16(enc_data[1].w);

    return cetd;  
}
#endif
















#endif


