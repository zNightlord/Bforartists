#ifndef BNPR_CONTOUR_TOPO__INCLUDED
#define BNPR_CONTOUR_TOPO__INCLUDED



#if defined(_KERNEL_MULTICOMPILE__GEOM_EXTRACT) || defined(_KERNEL_MULTICOMPILE__EXTRACT_MESH_CONTOUR_DATA)
/* Rambling between global and local id when extracting contours from current mesh */
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
	/* the start/end snake(i.e. vertex) of the contour curve? */ 
	bool curve_head; 
	bool curve_tail; 

    bool seg_head; /* dynamically set for segmentation */
    bool seg_tail; /* dynamically set for segmentation */
    
	bool looped_curve; 
    bool cusp_func_pstv; /*has positive cusp function*/
	bool occluded; /* valid when visibility test activated */

	/* flags used in 2d-resampled curves -------------------------- */
	bool curve_clipped;  // TODO: may be buggy, needs further testing
	bool is_corner;
	bool seg_head_contour; // seg-head based on contour segmentation 
	bool seg_head_clipped; // seg-head based on clipping
	bool seg_tail_clipped; // seg-tail based on clipping 
	bool is_curv_minima; // local curvature minima during corner detection
	bool occluded_filtered; 
}; 

uint encode_contour_flags(ContourFlags cf)
{
    uint cf_enc = 0u;
	cf_enc |= uint(cf.curve_head);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.curve_tail);
	cf_enc <<= 1u; 
    cf_enc |= uint(cf.seg_head); 
    cf_enc <<= 1u; 
    cf_enc |= uint(cf.seg_tail); 
    cf_enc <<= 1u;
    cf_enc |= uint(cf.looped_curve);
    cf_enc <<= 1u;
    cf_enc |= uint(cf.cusp_func_pstv);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.occluded);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.curve_clipped);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.is_corner);
	cf_enc <<= 1u; 
	cf_enc |= uint(cf.seg_head_contour); 
	cf_enc <<= 1u; 
	cf_enc |= uint(cf.seg_head_clipped); 
	cf_enc <<= 1u;
	cf_enc |= uint(cf.seg_tail_clipped);
	cf_enc <<= 1u;
	cf_enc |= uint(cf.is_curv_minima); 
	cf_enc <<= 1u;
	cf_enc |= uint(cf.occluded_filtered);

    return cf_enc; 
}

ContourFlags decode_contour_flags(uint cf_enc)
{
    ContourFlags cf; 
	cf.occluded_filtered = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.is_curv_minima = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.seg_tail_clipped = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.seg_head_clipped = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.seg_head_contour = (1u == (cf_enc & 1u)); 
	cf_enc >>= 1u; 
	cf.is_corner = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.curve_clipped = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.occluded = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
    cf.cusp_func_pstv = (1u == (cf_enc & 1u));
    cf_enc >>= 1u;
    cf.looped_curve = (1u == (cf_enc & 1u));
    cf_enc >>= 1u;
    cf.seg_tail = (1u == (cf_enc & 1u)); 
    cf_enc >>= 1u; 
    cf.seg_head = (1u == (cf_enc & 1u)); 
	cf_enc >>= 1u;
	cf.curve_tail = (1u == (cf_enc & 1u));
	cf_enc >>= 1u;
	cf.curve_head = (1u == (cf_enc & 1u));

    return cf; 
}

ContourFlags init_contour_flags(bool seg_head)
{
    ContourFlags cf;
	cf.curve_head = true; // needs further setup
	cf.curve_tail = true; // needs further setup
    cf.seg_head = seg_head;
    cf.seg_tail = false;       // needs further setup
    cf.looped_curve = true;    // needs further setup
    cf.cusp_func_pstv = false; // needs further setup
	cf.occluded = false;       // needs further setup
	cf.curve_clipped = false;  // needs further setup 
	cf.is_corner = false;      // needs further setup 
	cf.seg_head_contour = false;  // needs further setup 
	cf.seg_head_clipped = false;  // needs further setup
	cf.seg_tail_clipped = false;  // needs further setup 
	cf.is_curv_minima = false;  // needs further setup
	cf.occluded_filtered = false;  // needs further setup 
    return cf; 
}

void init_contour_curve_head(bool curve_head, inout ContourFlags cf)
{ /* should happen only once for all contour snakes */
	cf.curve_head = curve_head; 
}

void init_contour_curve_tail(bool curve_tail, inout ContourFlags cf)
{ /* should happen only once for all contour snakes */
	cf.curve_tail = curve_tail; 
}

void init_seg_head_contour(bool sample_2d_is_seg_head, inout ContourFlags cf)
{ /* should happen only once for all 2d samples */
	cf.seg_head_contour = sample_2d_is_seg_head; 
}

void init_seg_head_clipped(bool sample_2d_is_seg_head, inout ContourFlags cf)
{ /* should happen only once for all 2d samples */
	cf.seg_head_clipped = sample_2d_is_seg_head; 
}

void set_seg_tail_clipped(bool sample_2d_is_seg_tail, inout ContourFlags cf)
{ 
	cf.seg_tail_clipped = sample_2d_is_seg_tail; 
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
{ /* should happen only once for all contour snakes */
    cf.looped_curve = looped; 
}

void set_contour_cusp_flags(bool cusp_func_pstv, inout ContourFlags cf)
{
    cf.cusp_func_pstv = cusp_func_pstv; 
}

void set_contour_flags_occluded(inout ContourFlags cf)
{
	cf.occluded = true; 
}

void set_contour_flags_curve_clipped(bool clipped, inout ContourFlags cf)
{
	cf.curve_clipped = clipped; 
}

void set_contour_flags_is_corner(inout ContourFlags cf)
{
	cf.is_corner = true; 
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
	uint head_id;
	uint tail_id;
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
uint move_contour_id_along_loop(ContourCurveTopo cct, uint start_contour_id, float offset)
{
	return move_elem_along_loop(
		start_contour_id, (int(offset + 1e-10f)), 
		cct.head_id, cct.len
	);
}
#if defined(INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD)
ContourCurveTopo load_contour_curve_topo(uint contour_id, ContourFlags cf)
{
    ContourCurveTopo cct; 
    cct.looped_curve    = cf.looped_curve; 
    cct.head_id = ssbo_contour_snake_list_head_[contour_id]; 
    cct.len             = ssbo_contour_snake_list_len_[contour_id]; 
    cct.tail_id = cct.head_id + cct.len - 1u;  
    return cct; 
}
#endif



void FixLoopedJumps(
	inout uint seg_rank, inout uint seg_len, 
    ContourCurveTopo cct,
    uint scanseg_head, uint scanseg_tail, 
    ContourFlags cf_curv_head, ContourFlags cf_curv_tail,
    SSBOData_SegScanType_uint scan_res_step_0_tailseg, // segscan_uint_hf_decode(scan_output_buf_0_[cct.tail_id]);
    SSBOData_SegScanType_uint scan_res_step_1_headseg  // segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(cct.head_id)]);       
)
{
  // Fix seg rank & seg len for looped curve
  if (cct.looped_curve) {
    bool first_seg_in_loop = scanseg_head == cct.head_id;
    bool last_seg_in_loop = scanseg_tail == cct.tail_id;
    bool is_self_loop = seg_len == cct.len;

    if (first_seg_in_loop && !is_self_loop) {
      bool downflow = !(cf_curv_head.seg_head);

      if (downflow) {
        // prev half segment, at the end of loop
        // SSBOData_SegScanType_uint scan_res_step_1_tailseg =
        // 	segscan_uint_hf_decode(scan_output_buf_1_[REVERSE_ID(cct.tail_id)]); // should
        // equal to 0
        uint tail_seg_len = scan_res_step_0_tailseg.val + 1u;
        seg_rank += tail_seg_len; // fix seg rank
        seg_len += tail_seg_len;  // fix seg len
      }
    }
    if (last_seg_in_loop && !is_self_loop) {
      bool overflow = !(cf_curv_tail.seg_tail);

      if (overflow) {
        // next half segment, at the head of loop
        // SSBOData_SegScanType_uint scan_res_step_0_headseg =
        // 	segscan_uint_hf_decode(scan_output_buf_0_[cct.head_id]); // should equal to 0
        uint head_seg_len = scan_res_step_1_headseg.val + 1u;
        seg_len += head_seg_len; // fix seg len
      }
    }
  }
}





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



#if defined(USE_CONTOUR_2D_SAMPLE_TOPOLOGY_BUFFER)
void store_ssbo_contour_2d_sample_topology__flags(uint sample_id, ContourFlags flags)
{
	ssbo_contour_2d_sample_topology_[sample_id] = encode_contour_flags(flags);  
}
ContourFlags load_ssbo_contour_2d_sample_topology__flags(uint sample_id)
{
	return decode_contour_flags(ssbo_contour_2d_sample_topology_[sample_id]); 
}
void store_ssbo_contour_2d_sample_topology__curve_rank(uint sample_id, uint rank, uint num_samples)
{
	uint subbuff_offset = num_samples; 
	ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id] = rank; 
}
uint load_ssbo_contour_2d_sample_topology__curve_rank(uint sample_id, uint num_samples)
{
	uint subbuff_offset = num_samples; 
	return ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id]; 
}
void store_ssbo_contour_2d_sample_topology__curve_len(uint sample_id, uint len, uint num_samples)
{
	uint subbuff_offset = num_samples * 2u; 
	ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id] = len; 
}
uint load_ssbo_contour_2d_sample_topology__curve_len(uint sample_id, uint num_samples)
{
	uint subbuff_offset = num_samples * 2u; 
	return ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id]; 
}
void store_ssbo_contour_2d_sample_topology__seg_rank(uint sample_id, uint rank, uint num_samples)
{
	uint subbuff_offset = num_samples * 3u; 
	ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id] = rank; 
}
uint load_ssbo_contour_2d_sample_topology__seg_rank(uint sample_id, uint num_samples)
{
	uint subbuff_offset = num_samples * 3u; 
	return ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id]; 
}
void store_ssbo_contour_2d_sample_topology__seg_len(uint sample_id, uint len, uint num_samples)
{
	uint subbuff_offset = num_samples * 4u; 
	ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id] = len; 
}
uint load_ssbo_contour_2d_sample_topology__seg_len(uint sample_id, uint num_samples)
{
	uint subbuff_offset = num_samples * 4u; 
	return ssbo_contour_2d_sample_topology_[subbuff_offset + sample_id]; 
}
ContourCurveTopo load_contour_2d_sample_curve_topo(uint contour_id, ContourFlags cf, uint num_samples)
{
	uint curve_rank = load_ssbo_contour_2d_sample_topology__curve_rank(contour_id, num_samples); 

    ContourCurveTopo cct; 
    cct.looped_curve    = cf.looped_curve; 
    cct.head_id = contour_id - curve_rank; 
    cct.len             = load_ssbo_contour_2d_sample_topology__curve_len(contour_id, num_samples); 
    cct.tail_id = cct.head_id + cct.len - 1u;  
    return cct; 
}
bool is_2d_sample_curve_looped(bool contour_looped, bool contour_crve_clipped, bool single_sub_seg)
{
	bool is_sub_seg_loop = contour_looped && (!contour_crve_clipped) && single_sub_seg; 
	return is_sub_seg_loop; 
}
#endif






#endif


