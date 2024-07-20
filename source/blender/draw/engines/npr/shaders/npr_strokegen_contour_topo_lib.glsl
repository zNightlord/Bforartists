

#pragma BLENDER_REQUIRE(npr_strokegen_encode_lib.glsl)


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
	uint temporal_rec_id;  
};

void store_contour_edge_transfer_data_(uint contour_edge_id, ContourEdgeTransferData cetd)
{
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+0] = floatBitsToUint(cetd.vpos_ws[0].x);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+1] = floatBitsToUint(cetd.vpos_ws[0].y);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+2] = floatBitsToUint(cetd.vpos_ws[0].z);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+3] = floatBitsToUint(cetd.vpos_ws[1].x);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+4] = floatBitsToUint(cetd.vpos_ws[1].y);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+5] = floatBitsToUint(cetd.vpos_ws[1].z);

	ssbo_contour_edge_transfer_data_[contour_edge_id*9+6] = encode_contour_flags(cetd.cf);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+7] = packHalf2x16(
		vec2(
			cetd.cusp_funcs[0] > .0f ? 1.0f : -1.0f, 
			cetd.cusp_funcs[1] > .0f ? 1.0f : -1.0f
		)
	);
	ssbo_contour_edge_transfer_data_[contour_edge_id*9+8] = cetd.temporal_rec_id; 
}

ContourEdgeTransferData load_contour_edge_transfer_data(uint contour_edge_id)
{
	ContourEdgeTransferData cetd; 
	uvec4 enc_data[2];
	cetd.vpos_ws[0].x = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+0u]);
	cetd.vpos_ws[0].y = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+1u]);
	cetd.vpos_ws[0].z = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+2u]);
	cetd.vpos_ws[1].x = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+3u]);
	cetd.vpos_ws[1].y = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+4u]);
	cetd.vpos_ws[1].z = uintBitsToFloat(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+5u]);

	cetd.cf = decode_contour_flags(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+6u]);
	cetd.cusp_funcs = unpackHalf2x16(ssbo_contour_edge_transfer_data_[contour_edge_id*9u+7u]);
	cetd.temporal_rec_id = ssbo_contour_edge_transfer_data_[contour_edge_id*9u+8u];

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



#if defined(_KERNEL_MULTICOMPILE__CALC_TEMPORAL_CONTOUR_RECORDS__MAIN)
 // reuse ssbo slots
#define ssbo_contour_temporal_records_old_ ssbo_selected_edge_to_edge_
#endif

#define MAX_TEMPORAL_FRAMES 2u 			 // Match to MAX_TEMPORAL_FRAMES 		  in bnpr_defines.hh
#define MAX_TEMPORAL_TRACKED_OBJECTS 64u // Match to MAX_TEMPORAL_TRACKED_OBJECTS in bnpr_defines.hh
// note: in the furture, we could try to match against multiple history frames 
#define temporal_record_counter(obj_id, frame_id) ssbo_temporal_record_counters_[(obj_id * MAX_TEMPORAL_FRAMES) + (frame_id % MAX_TEMPORAL_FRAMES)] 


/* Temporal contour record, generated per contour edge. */
struct TemporalRecordFlags
{
	uint subd_tree_code_chain; // 9 bits = 3 bits x 3 subd levels 
	bool valid_contour_data; 
}; 
uint encode_temporal_record_flags(TemporalRecordFlags trf)
{
	uint enc = 0u; 
	enc <<= 9u; 
	enc |= (trf.subd_tree_code_chain & 0x01ffu); 
	enc <<= 1u; 
	enc |= uint(trf.valid_contour_data); 

	return enc; 
}
TemporalRecordFlags decode_temporal_record_flags(uint enc)
{
	TemporalRecordFlags trf; 
	
	trf.valid_contour_data = (1u == (enc & 1u)); 
	enc >>= 1u; 
	trf.subd_tree_code_chain = enc & 0x01ffu; 
	enc >>= 9u; 

	return trf; 
}
TemporalRecordFlags init_temporal_record_flags(uint subd_code_chain)
{
	TemporalRecordFlags trf; 
	trf.subd_tree_code_chain = subd_code_chain; 
	trf.valid_contour_data = false; 
	return trf; 
}

void append_subd_tree_node_to_code(uint tree_code, inout uint subd_tree_code_chain)
{ // Traverse tree upwards and concatenate node code, 3 bits each level
	subd_tree_code_chain <<= 3u; 
	subd_tree_code_chain |= ((tree_code << 1u) | 1u/*valid code*/); 
}
void pop_subd_tree_code_chain(inout uint subd_tree_code_chain, out uint tree_code, out bool null_node)
{
	tree_code = (subd_tree_code_chain >> 1u) & 0x03u;
	null_node = (0u == (subd_tree_code_chain & 0x01u)); 
	subd_tree_code_chain >>= 3u; 
}


// Temporal Record - Flags | offset stride 0 
// ------------------------------------------------
#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_NEW)
void store_ssbo_contour_temporal_records_new__flags(uint enc_id, TemporalRecordFlags flags)
{
	ssbo_contour_temporal_records_new_[enc_id] = encode_temporal_record_flags(flags);  
}
TemporalRecordFlags load_ssbo_contour_temporal_records_new__flags(uint enc_id)
{
	return decode_temporal_record_flags(ssbo_contour_temporal_records_new_[enc_id]); 
}
#endif
#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_OLD)
TemporalRecordFlags load_ssbo_contour_temporal_records_old__flags(uint enc_id)
{
	return decode_temporal_record_flags(ssbo_contour_temporal_records_old_[enc_id]); 
}
#endif

// Temporal Record - Base Edge Id | offset stride 1
// ------------------------------------------------
#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_NEW)
void store_ssbo_contour_temporal_records_new__subd_root_edge_id(uint enc_id, uint root_edge_id, uint num_recs)
{
	uint subbuff_offset = num_recs * 1u; 
	ssbo_contour_temporal_records_new_[subbuff_offset + enc_id] = root_edge_id; 
}
uint load_ssbo_contour_temporal_records_new__subd_root_edge_id(uint enc_id, uint num_recs)
{
	uint subbuff_offset = num_recs * 1u; 
	return ssbo_contour_temporal_records_new_[subbuff_offset + enc_id]; 
}
#endif
#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_OLD)
uint load_ssbo_contour_temporal_records_old__subd_root_edge_id(uint enc_id, uint num_recs)
{
	uint subbuff_offset = num_recs * 1u; 
	return ssbo_contour_temporal_records_old_[subbuff_offset + enc_id]; 
}
#endif


// Temporal Record - Contour Data | offset stride 2
struct TemporalRecordContourData
{
	uint curve_key; 	 // 24 bits
	uint seg_key; 		 // 24 bits 
	uint curve_rank; 	 // 24 bits
	uint seg_rank; 		 // 24 bits

	ContourFlags cf; 	 // 32 bits
}; 
TemporalRecordContourData init_temporal_record_contour_data()
{
	TemporalRecordContourData trcd; 
	trcd.curve_key = trcd.seg_key = trcd.curve_rank = trcd.seg_rank = 0u; 
	trcd.cf = init_contour_flags(false); 
	return trcd; 
}
void encode_temporal_record_contour_data(TemporalRecordContourData trcd, out uvec4 enc_0)
{
	enc_0.xyz = pack_u24_x4(
		uvec4(trcd.curve_key, trcd.seg_key, trcd.curve_rank, trcd.seg_rank)
	); 
	enc_0.w = encode_contour_flags(trcd.cf); 
} 
TemporalRecordContourData decode_temporal_record_contour_data(uvec4 enc_0)
{
	TemporalRecordContourData trcd; 

	uvec4 dec_0_xyz = unpack_u24_x4(enc_0.xyz); 
	trcd.curve_key = dec_0_xyz.x; 
	trcd.seg_key = dec_0_xyz.y; 
	trcd.curve_rank = dec_0_xyz.z; 
	trcd.seg_rank = dec_0_xyz.w; 

	trcd.cf = decode_contour_flags(enc_0.w); 

	return trcd; 
}

#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_NEW)
void store_ssbo_contour_temporal_records_new__contour_data(uint rec_id, TemporalRecordContourData trcd, uint num_recs)
{
	uint subbuff_offset = num_recs * 2u; 
	uvec4 enc;
	encode_temporal_record_contour_data(trcd, /*out*/enc); 
	ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 0u] = enc[0]; 
	ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 1u] = enc[1]; 
	ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 2u] = enc[2]; 
	ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 3u] = enc[3]; 
}
TemporalRecordContourData load_ssbo_contour_temporal_records_new__contour_data(uint rec_id, uint num_recs)
{
	uint subbuff_offset = num_recs * 2u; 
	uvec4 enc; 
	enc[0] = ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 0u]; 
	enc[1] = ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 1u]; 
	enc[2] = ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 2u]; 
	enc[3] = ssbo_contour_temporal_records_new_[subbuff_offset + rec_id * 4u + 3u]; 
	return decode_temporal_record_contour_data(enc); 
}
#endif
#if defined(USE_CONTOUR_TEMPORAL_RECORD_BUFFER_OLD)
TemporalRecordContourData load_ssbo_contour_temporal_records_old__contour_data(uint rec_id, uint num_recs)
{
	uint subbuff_offset = num_recs * 2u; 
	uvec4 enc; 
	enc[0] = ssbo_contour_temporal_records_old_[subbuff_offset + rec_id * 4u + 0u]; 
	enc[1] = ssbo_contour_temporal_records_old_[subbuff_offset + rec_id * 4u + 1u]; 
	enc[2] = ssbo_contour_temporal_records_old_[subbuff_offset + rec_id * 4u + 2u]; 
	enc[3] = ssbo_contour_temporal_records_old_[subbuff_offset + rec_id * 4u + 3u]; 
	return decode_temporal_record_contour_data(enc); 
}
#endif


// TODO: put the actual logic here, for now only for debug purposes
bool match_history_rec(
	TemporalRecordFlags rec_flags, TemporalRecordContourData rec_contour_data
	/* float curr_cusp_func */)
{
	return rec_flags.valid_contour_data; 
}



/* Each edge holds a pointer to its records, 
 * the ptr is null if it wasn't a contour */
#define PER_EDGE_TEMPORAL_REC_ID_NULL 0xffffffffu
#if defined(USE_EDGE_TO_TEMPORAL_RECORD_BUFFER)
void store_ssbo_edge_to_new_temporal_record_(uint wedge_id, uint rec_id)
{ // not using interleaved here for easier debugging
	ssbo_edge_to_temporal_record_[wedge_id * 2u] = rec_id; 
}
uint load_ssbo_edge_to_new_temporal_record_(uint wedge_id)
{
	return ssbo_edge_to_temporal_record_[wedge_id * 2u]; 
}
void store_ssbo_edge_to_old_temporal_record_(uint wedge_id, uint rec_id)
{
	ssbo_edge_to_temporal_record_[wedge_id * 2u + 1u] = rec_id; 
}
uint load_ssbo_edge_to_old_temporal_record_(uint wedge_id)
{
	return ssbo_edge_to_temporal_record_[wedge_id * 2u + 1u]; 
}
#endif

#define PER_CONTOUR_TEMPORAL_REC_ID_NULL 0xffffffffu


#endif


