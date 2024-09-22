
#pragma BLENDER_REQUIRE(npr_strokegen_encode_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)

#ifndef BNPR_CONTOUR_GEOM__INCLUDED
#define BNPR_CONTOUR_GEOM__INCLUDED



// ---------------------------------------------------------------------------------------------------
// Functions for Detecting Contour Edges
// ---------------------------------------------------------------------------------------------------
bool is_back_face(float ndv)
{
	return ndv <= .0f; 
}
bool is_contour_edge_common(EdgeFlags ef)
{
	bool is_contour = true; 
	is_contour = is_contour && (!ef.del_by_split) && (!ef.del_by_collapse) && (!ef.dupli); 
	is_contour = is_contour && (!ef.border); // TODO: support for border edges
	return is_contour; 
}
/*     v0
 *    /  \
 *   / f1 \ 
 *  v1 -- v3 Winding:CCW
 *   \ f0 /  
 *    \  /    
 *     v2                   
*/
bool is_contour_edge(
	vec3 v0, vec3 v1, vec3 v2, vec3 v3, vec3 cam_pos, EdgeFlags ef, 
	out float face_orient_123, out float face_orient_013
)
{ /* impl based on overlay_outline_prepass_vert_no_geom.glsl */
	vec3 v10 = v0 - v1;
   	vec3 v12 = v2 - v1;
	vec3 v13 = v3 - v1;

	vec3 n0 = cross(v12, v13);
	vec3 n3 = cross(v13, v10);

	vec3 view_dir = cam_pos - v1; 

	face_orient_123 = dot(view_dir, n0);
  	face_orient_013 = dot(view_dir, n3);
	bool is_contour = is_back_face(face_orient_123) != is_back_face(face_orient_013); 
	
	is_contour = is_contour && is_contour_edge_common(ef); 

	return is_contour; 
}

bool is_interp_contour_edge__before_tessellation(VertFlags vf_0, VertFlags vf_1, EdgeFlags ef)
{
    bool is_contour = (
        (vf_0.front_facing && vf_1.back_facing)
        || (vf_0.back_facing && vf_1.front_facing)
    ); 
	
	is_contour = is_contour && is_contour_edge_common(ef); 

	return is_contour; 
}

bool is_interp_contour_edge__after_tessellation(
	VertFlags vf_0, VertFlags vf_1, VertFlags vf_2, VertFlags vf_3, EdgeFlags ef, 
	out float face_orient_123, out float face_orient_013
)
{
	bool is_contour = vf_1.contour && vf_3.contour;  
	face_orient_123 = vf_2.front_facing ? 1.0f : -1.0f; 
	face_orient_013 = vf_0.front_facing ? 1.0f : -1.0f; 

	is_contour = is_contour && is_contour_edge_common(ef); 

	return is_contour; 
}

float calc_interp_contour_edge_factor(vec3 edge_vnor[2], vec3 edge_vpos[2], vec3 cam_pos_ws)
{ // edge_vnor, edge_vpos: normal/pos of 2 verts on the edge; 
	float contour_interpo_factor = .5f; 
	vec2 ndv = vec2(
		dot(edge_vnor[0], cam_pos_ws - edge_vpos[0]),  
		dot(edge_vnor[1], cam_pos_ws - edge_vpos[1]) 
	); 
	contour_interpo_factor = ndv[0] / (ndv[0] - ndv[1]); 
	/* split pos = mix(vpos_0, vpos_1, interpo) = vpos_0 + interpo * (vpos_1 - vpos_0) */

	return contour_interpo_factor; 
}
vec3 calc_interp_contour_vert_pos(vec3 edge_vnor[2], vec3 edge_vpos[2], vec3 cam_pos_ws)
{ 
	float contour_interpo_factor = calc_interp_contour_edge_factor(edge_vnor, edge_vpos, cam_pos_ws); 
	return edge_vpos[0] + contour_interpo_factor * (edge_vpos[1] - edge_vpos[0]); 
}





// ---------------------------------------------------------------------------------------------------
// Functions for Rasterization
// ---------------------------------------------------------------------------------------------------
// uniform vec2 pcs_line_raster_resolution_;

// Note: Some drivers do not implement uniform initializers correctly.
// https://www.khronos.org/opengl/wiki/Uniform_(GLSL)
uniform vec4 g_homogenous_clip_planes[6] =
{
	vec4(1, 0, 0, 1), // -w <= x
	vec4(-1, 0, 0, 1), // x <= w
	vec4(0, -1, 0, 1), // y <= w 
	vec4(0, 1, 0, 1), // -w <= y
	vec4(0, 0, -1, 1), // z <= w
	vec4(0, 0, 1, 0) // 0 <= z
};

vec3 hclip_to_ndc(vec4 pos_cs)
{
	pos_cs /= (pos_cs.w);
	pos_cs.xy = pos_cs.xy * 0.5 + 0.5;

	return pos_cs.xyz;
}

/* Line clipping in 3D homogeneous space ----------------------------- */
vec2 frustum_clip_line_segment(
	vec4 vpos_hcs_0, vec4 vpos_hcs_1, // homogenous clip space positions
	out bool reject, out bool inside)
{
	float alpha0 = 0; // v0 = lerp(v0, v1, alpha0)
	float alpha1 = 1; // v1 = lerp(v0, v1, alpha1)
	reject = false; // Is this edge totally out of frustum
	inside = true;  // Is this edge inside view frustum

	for (uint i = 0u; i < 6u; ++i)
	{
		// Compute clip factor
		float d0 = dot(vpos_hcs_0, g_homogenous_clip_planes[i]);
		float d1 = dot(vpos_hcs_1, g_homogenous_clip_planes[i]);
		float alpha = d0 / (d0 - d1);

		// Clip edge
		bool out0 = d0 < 0; // if v0 is outside ith clip plane
		bool out1 = d1 < 0; // if v1 is outside ith clip plane
		inside = (inside && ((!out0) && (!out1)));

		if (out0 && (!out1))
		{
			// v0 outside, v1 inside, clip v0 to new position
			alpha0 = max(alpha, alpha0);
		}
		if ((!out0) && (out1))
		{
			// clip v1 to new position
			alpha1 = min(alpha, alpha1);
		}

		// Reject conditions:
		// 1. both verts of edge are on the outside of same clip plane
		// 2. verts "crosses" clip planes, 
		// -- where alpha1 > alpha0:
		//             .  V1
		//             . /
		//             X <--- alpha0
		//  alpha1   / .
		// .....->.X...+-----------------+
		//       /     |                 |
		//      V0     |  Screen Frustum |
		//             |     (in 2D)     |
		//             |                 |
		//             +-----------------+
		if ((out0 && out1) || (alpha0 > alpha1))
		{
			reject = true;
		}
	}

	return vec2(alpha0, alpha1);
}


vec2 ndc_to_viewport(vec2 coord_ndc, vec2 raster_resolution)
{
	vec2 coordVP =
		max(
			vec2(0.5, 0.5),
			min(
				coord_ndc * raster_resolution.xy,
				raster_resolution.xy - vec2(1.0f)
			)
		);
	return coordVP;
}

vec2 viewport_to_ndc(vec2 coord_ss, vec2 raster_resolution)
{
	vec2 coordNDC = coord_ss / raster_resolution.xy;
	return clamp(coordNDC, vec2(.0f), vec2(1.0f));
}



struct LineRasterResult
{
    uint num_frags; 
	// note that the vertex order can be flipped ----------------
    vec4 begend_uvs; // clipped beg/end uvs. within range of 0-1, 
    vec2 begend_wclips; // homogeneous w
	// ----------------------------------------------------------
    
    bool is_line_clip_rejected; // totally outside frustum
	bool is_line_clipped; 		// partially within frustum 
    bool is_x_major_line; 
    bool beg_from_p0; // raster begins at the first vert on the original line 
}; 

void encode_line_raster_result(LineRasterResult res, out uvec4 d0123, out uint d4)
{
    uint d0 = res.num_frags; // 28 bits for num_frags should be sufficient
    d0 <<= 1u; 
    d0 |= (res.is_line_clip_rejected ? 1u : 0u);
	d0 <<= 1u;
	d0 |= (res.is_line_clipped ? 1u : 0u); 
	d0 <<= 1u; 
	d0 |= (res.is_x_major_line ? 1u : 0u); 
	d0 <<= 1u; 
	d0 |= (res.beg_from_p0 ? 1u : 0u); 

	uvec2 d12 = uvec2(0u); 
	res.begend_uvs = clamp(res.begend_uvs, vec4(.0f), vec4(1.0f)); 
	float fixed_factor = float((1u << 16u) - 1u); 
	uvec4 uv_fixed = uvec4(round(res.begend_uvs * fixed_factor)) & 0xffffu; 
	d12 = uv_fixed.xy | (uv_fixed.zw << 16u); 

	uvec2 d34 = uvec2(0u);
	d34.x = floatBitsToUint(res.begend_wclips.x);
	d34.y = floatBitsToUint(res.begend_wclips.y);

	d0123 = uvec4(d0, d12, d34.x);
	d4 = d34.y;
}

LineRasterResult decode_line_raster_result(uvec4 d0123, uint d4)
{
	LineRasterResult res; 

	uint d0 = d0123.x; 
	res.num_frags = (d0 >> 4u); 
	res.beg_from_p0           = (d0 & 0x1u) != 0u; 
	res.is_x_major_line       = (d0 & 0x2u) != 0u; 
	res.is_line_clipped       = (d0 & 0x4u) != 0u; 
	res.is_line_clip_rejected = (d0 & 0x8u) != 0u; 

	uvec2 d12 = uvec2(d0123.yz);
	uvec4 uv_fixed = uvec4(
		d12 & 0xffffu, 
		d12 >> 16u
	); 
	float fixed_factor = float((1u << 16u) - 1u); 
	res.begend_uvs = clamp(vec4(uv_fixed) / fixed_factor, vec4(.0f), vec4(1.0f)); 

	uvec2 d34 = uvec2(d0123.w, d4); 
	res.begend_wclips.x = uintBitsToFloat(d34.x); 
	res.begend_wclips.y = uintBitsToFloat(d34.y); 

	return res; 
}



#if defined(INCLUDE_LOAD_STORE_CONTOUR_RASTER_DATA)

#define LINE_FRAGS_OFFSET_NONE 0xffffffffu // no fragments are allocated for this line

void store_contour_edge_raster_data(uint contour_edge_id, LineRasterResult line_raster_data, uint frag_alloc_offset)
{
	uvec4 data_enc_0123; 
	uint data_enc_4; 
	encode_line_raster_result(line_raster_data, /*out*/data_enc_0123, data_enc_4); 

	ssbo_contour_raster_data_[contour_edge_id * 6 + 0u] = data_enc_0123[0];
	ssbo_contour_raster_data_[contour_edge_id * 6 + 1u] = data_enc_0123[1];
	ssbo_contour_raster_data_[contour_edge_id * 6 + 2u] = data_enc_0123[2];
	ssbo_contour_raster_data_[contour_edge_id * 6 + 3u] = data_enc_0123[3];
	ssbo_contour_raster_data_[contour_edge_id * 6 + 4u] = data_enc_4; 

	ssbo_contour_raster_data_[contour_edge_id * 6 + 5u] = 
			0 < line_raster_data.num_frags ? frag_alloc_offset : LINE_FRAGS_OFFSET_NONE; 
}

LineRasterResult load_contour_edge_raster_data(uint contour_edge_id)
{
	LineRasterResult line_raster_data = decode_line_raster_result(
		uvec4(
			ssbo_contour_raster_data_[contour_edge_id * 6u + 0u],
			ssbo_contour_raster_data_[contour_edge_id * 6u + 1u],
			ssbo_contour_raster_data_[contour_edge_id * 6u + 2u],
			ssbo_contour_raster_data_[contour_edge_id * 6u + 3u]
		), 
		ssbo_contour_raster_data_[contour_edge_id * 6u + 4u]
	);
	
	return line_raster_data; 
}

LineRasterResult load_contour_edge_raster_data(uint contour_edge_id, out uint head_frag_id)
{
	LineRasterResult line_raster_data = load_contour_edge_raster_data(contour_edge_id);
	
	head_frag_id = ssbo_contour_raster_data_[contour_edge_id * 6u + 5u]; // start frag rastered from this edge
	
	return line_raster_data; 
}

LineRasterResult load_contour_edge_raster_data(uint contour_edge_id, out uint head_frag_id, out uint tail_frag_id)
{
	LineRasterResult line_raster_data = load_contour_edge_raster_data(contour_edge_id, /*out*/head_frag_id);
	
	tail_frag_id = head_frag_id + line_raster_data.num_frags - 1u; // end frag rastered from this edge
	
	return line_raster_data; 
}
#endif

LineRasterResult raster_line_segment(
    vec4 vpos_0, vec4 vpos_1, 
    vec2 raster_resolution, /* mat4 Matrix_M,  */mat4 Matrix_V, mat4 Matrix_P
) {
    // OS --> WS (nothing here since all contours are in world space) 

    // WS --> VS
    // -----------------------------------------------------------------
    vpos_0 = (Matrix_V * vpos_0); // mul(Matrix_M, vpos_0));
    vpos_1 = (Matrix_V * vpos_1); // mul(Matrix_M, vpos_1));

    // VS --> HClip
    // --------------------------------------------------------------------
    vpos_0 = (Matrix_P * vec4(vpos_0.xyz, 1.0f)); // HClip coord  
    vpos_1 = (Matrix_P * vec4(vpos_1.xyz, 1.0f)); // HClip coord

    // HClip --> HClip(Clipped) 
    // -----------------------------------------------------------------------------
    bool reject, inside;
    vec2 interpolate = frustum_clip_line_segment(
        vpos_0, vpos_1, // in
        reject, inside  // out
    );
    vpos_0 = mix(vpos_0, vpos_1, interpolate.x);
    vpos_1 = mix(vpos_0, vpos_1, interpolate.y);
    vec2 whclip_v0v1 = vec2(vpos_0.w, vpos_1.w);

    // HClip(Clipped)->NDC
    // -----------------------------------------------------------------------------
    // .w component is preserved & for later use
    vpos_0.xyz = reject ? vec3(0, 0, 0) : hclip_to_ndc(vpos_0);
    vpos_1.xyz = reject ? vec3(0, 0, 0) : hclip_to_ndc(vpos_1);


    // NDC->Viewport
    // -------------------------------------------------------------------------
    vec4 frags_coord = vec4(
		ndc_to_viewport(vpos_0.xy, raster_resolution).xy, 
		ndc_to_viewport(vpos_1.xy, raster_resolution).xy
	);

    // Compute fragment(== segment, in our context) count
    // Is edge X or Y major?
    // \          /  | A & C: X-Major
    //   \   B  /    | B & D: Y-Major 
    //     \  /      | 
    //  A   O    C   | 
    //    /    \     | 
    //  /   D    \   | 
    ///            \ | 
    vec2 dxdy = abs(frags_coord.zw - frags_coord.xy);
    bool is_x_major_line = (dxdy.y < dxdy.x);

	// raster begins at the lower point of the line
    bvec2 is_p0_lower = lessThan(frags_coord.xy, frags_coord.zw);
    bool beg_from_p0 = is_x_major_line ? is_p0_lower.x : is_p0_lower.y; 

	// Reorder two points to fit rasterization order
    frags_coord = beg_from_p0 ? frags_coord.xyzw : frags_coord.zwxy;
	vec2 frags_whclip = beg_from_p0 ? whclip_v0v1.xy : whclip_v0v1.yx; 

    // Compute fragment count according to X or Y major
    // Only accept contours that passes frustum clipping 
    uint num_frags = reject ? 0 : uint(max(0, ceil(max(dxdy.x, dxdy.y) + 0.0001f)));


	// Assembly the result
	LineRasterResult res; 
	res.num_frags = num_frags; 
    res.begend_uvs = frags_coord.xyzw / raster_resolution.xyxy; // in original shader, I encoded coord instead of 01uv
    res.begend_wclips = frags_whclip; 
    
    res.is_line_clip_rejected = reject;
	res.is_line_clipped = (!reject) && (!inside); 
    res.is_x_major_line = is_x_major_line; 
    res.beg_from_p0 = beg_from_p0; 

	return res; 
}


vec2 calc_frag_screen_pos(
	vec4 line_beg_end_coords, 
	uint head_frag_id, uint curr_frag_id,
	bool is_x_major,
	out float linear_factor,
	out float linear_step
)
{
	// line_beg_end_coords.x = floor(line_beg_end_coords.x);
	// line_beg_end_coords.z = ceil(line_beg_end_coords.z);
	// Flip xy/zw to yx/wz if edge is y-major,
	// so that we always have 'major' axis coord in first slot
	line_beg_end_coords = is_x_major ? line_beg_end_coords.xyzw : line_beg_end_coords.yxwz;
	float majorAxisBegPos = floor(line_beg_end_coords.x);
	float majorAxisOffset = float(curr_frag_id - head_frag_id) + 0.5;

	vec2 targTexel; // .x: Major axis, .y: Dual axis
	targTexel.x = majorAxisBegPos + majorAxisOffset;
	linear_step = 1.0 / (line_beg_end_coords.z - line_beg_end_coords.x);
	linear_factor = clamp((targTexel.x - line_beg_end_coords.x) * linear_step, .0f, 1.0f);
	targTexel.y = mix(line_beg_end_coords.y, line_beg_end_coords.w, linear_factor);

	
	targTexel = is_x_major ? targTexel.xy : targTexel.yx; // Flip to actual coord

	linear_step =
		abs(line_beg_end_coords.w - line_beg_end_coords.y) / 
		(ceil(line_beg_end_coords.z) - floor(line_beg_end_coords.x));
	
	return targTexel;
}

float interpolate_frag_depth(float whclip0, float whclip1, float linear_factor)
{
	return 1.0f / mix(
			1.0f / whclip0, 
			1.0f / whclip1, 
			linear_factor
		);
}

struct FragVisibilityTestResult
{ 
	bool visible; 
	uint head_frag_id; 
}; 

uint encode_frag_visibility_test_result(FragVisibilityTestResult res)
{
	uint d0 = res.head_frag_id; 
	d0 <<= 1u; 
	d0 |= (res.visible ? 1u : 0u); 
	return d0; 
}

FragVisibilityTestResult decode_frag_visibility_test_result(uint d0)
{
	FragVisibilityTestResult res; 
	res.visible = (d0 & 0x1u) != 0u; 
	res.head_frag_id = (d0 >> 1u); 
	return res; 
}

struct ContourVisibilitySplitInfo
{
	// d0
	uint parent_contour_id; 
	bool is_new_contour; 
	bool is_visible; 
	bool no_rastered_frags; 
	bool split_by_occlusion;
	// d123, 24 bits for each
	vec2 begend_ratios; // range [0, 1], the ratio of start/end of this segment at the parent contour
#define FRAG_SEG_ID_MASK 0x00ffffffu
	uint prev_frag_seg_head_id; // 24 bits 
	bool no_rastered_prev_contour; // prev contour is clipped away from view
	uint next_frag_seg_head_id; // 24 bits
	bool no_rastered_next_contour; // next contour is clipped away from view
}; 


#define BIT_TYPE_CONTOUR_VIS_SPLIT__NEW_EDGE 0x1u
#define BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_NOT_RASTERED 0x2u
#define BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_SPLIT 0x4u
#define BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_RASTERED_NO_SPLIT 0x8u
uint calc_contour_visibility_states(ContourVisibilitySplitInfo cvsi)
{ // note: partially clipped edges are considered old_edge_rastered_no_split
	bool new_edge                   = cvsi.is_new_contour;
	bool old_edge_not_rastered      = cvsi.no_rastered_frags;
	bool old_edge_split             = !(new_edge || old_edge_not_rastered) && cvsi.split_by_occlusion; 
	bool old_edge_rastered_no_split = !(new_edge || old_edge_split || old_edge_not_rastered); 

	uint res = 0u; 
	res |= (new_edge                   ? BIT_TYPE_CONTOUR_VIS_SPLIT__NEW_EDGE : 0u);
	res |= (old_edge_not_rastered      ? BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_NOT_RASTERED : 0u);
	res |= (old_edge_split             ? BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_SPLIT : 0u);
	res |= (old_edge_rastered_no_split ? BIT_TYPE_CONTOUR_VIS_SPLIT__OLD_EDGE_RASTERED_NO_SPLIT : 0u);

	return res; 
}



void init_contour_visibility_split_info(uint contour_id, inout ContourVisibilitySplitInfo cvsi)
{
	cvsi.parent_contour_id = contour_id; 
	cvsi.is_new_contour = false; 
	cvsi.is_visible = true; // match the default state in ContourFlags 
	cvsi.no_rastered_frags = true;
	
	cvsi.begend_ratios = vec2(.0f, 1.0f); 
	cvsi.prev_frag_seg_head_id = 0u;
	cvsi.next_frag_seg_head_id = 0u;
	cvsi.no_rastered_prev_contour = true; 
	cvsi.no_rastered_next_contour = true; 
}

uvec4 encode_contour_visibility_split_info(ContourVisibilitySplitInfo cvsi)
{
	uvec4 d0123 = uvec4(0u); 
	d0123.x = cvsi.parent_contour_id; 
	d0123.x <<= 1u; 
	d0123.x |= (cvsi.is_new_contour 	? 1u : 0u);
	d0123.x <<= 1u;
	d0123.x |= (cvsi.is_visible     	? 1u : 0u);
	d0123.x <<= 1u;
	d0123.x |= (cvsi.no_rastered_frags    	? 1u : 0u);
	d0123.x <<= 1u;
	d0123.x |= (cvsi.split_by_occlusion 	? 1u : 0u);

	uvec2 begend_ratios_fixed; 
	{ // Pack into 2x24 bits
		cvsi.begend_ratios = clamp(cvsi.begend_ratios, vec2(.0f), vec2(1.0f));
		float fixed_factor = float((1u << 24u) - 1u);
		begend_ratios_fixed = uvec2(round(cvsi.begend_ratios * fixed_factor)) & FRAG_SEG_ID_MASK; 
	}
	if (cvsi.no_rastered_prev_contour) cvsi.prev_frag_seg_head_id = FRAG_SEG_ID_MASK; 
	if (cvsi.no_rastered_next_contour) cvsi.next_frag_seg_head_id = FRAG_SEG_ID_MASK;
	
	d0123.yzw = pack_u24_x4(
		uvec4(
			begend_ratios_fixed.x, 
			begend_ratios_fixed.y, 
			cvsi.prev_frag_seg_head_id & FRAG_SEG_ID_MASK, 
			cvsi.next_frag_seg_head_id & FRAG_SEG_ID_MASK
		)
	); 

	return d0123; 
}

ContourVisibilitySplitInfo decode_contour_visibility_split_info(uvec4 d0123)
{
	ContourVisibilitySplitInfo cvsi; 
	cvsi.split_by_occlusion = ((d0123.x & 0x1u) != 0u); 
	cvsi.no_rastered_frags      = ((d0123.x & 0x2u) != 0u); 
	cvsi.is_visible        = ((d0123.x & 0x4u) != 0u); 
	cvsi.is_new_contour    = ((d0123.x & 0x8u) != 0u); 
	cvsi.parent_contour_id = (d0123.x >> 4u); 

	uvec4 d123_dec = unpack_u24_x4(d0123.yzw); 
	cvsi.begend_ratios = vec2(d123_dec.xy) / float((1u << 24u) - 1u); 
	cvsi.prev_frag_seg_head_id = d123_dec.z; 
	cvsi.no_rastered_prev_contour = cvsi.prev_frag_seg_head_id == FRAG_SEG_ID_MASK; 
	cvsi.next_frag_seg_head_id = d123_dec.w; 
	cvsi.no_rastered_next_contour = cvsi.next_frag_seg_head_id == FRAG_SEG_ID_MASK; 

	return cvsi; 
}









// ----------------------------------------------------------------------------------------------
// 2D resampling of 3D contour curves 
// ----------------------------------------------------------------------------------------------
struct Contour2DResampleRasterData
{
	bool has_samples; 
	vec4 begend_uvs; // the same as LineRasterResult.begend_uvs
};


uvec3 encode_contour_2d_resample_data(Contour2DResampleRasterData data)
{
	uvec3 d = uvec3(0u); 
	d.x = (data.has_samples ? 1u : 0u); 
	
	data.begend_uvs = clamp(data.begend_uvs, vec4(.0f), vec4(1.0f)); 
	d.yz = pack_2d_uv_x2(data.begend_uvs); 
	
	return d; 
}

Contour2DResampleRasterData decode_contour_2d_resample_data(uvec3 d)
{
	Contour2DResampleRasterData data; 
	data.has_samples = (d.x != 0u); 
	data.begend_uvs = unpack_2d_uv_x2(d.yz); 
	
	return data; 
}

#define CONTOUR_2D_SAMPLES_OFFSET_NONE 0xffffffffu // no samples are allocated for this contour



// Tangent Estimation ------------------------------------------------
// Visit each neighboring sample in the local window 
// and accumulate weighted least squares parameters 
struct Curve2DWLSFitParams
{
	vec3 a_123; 
	vec4 b_x1y1_x2y2; 
}; 

Curve2DWLSFitParams curve_2d_wls_fit_init_params()
{
	Curve2DWLSFitParams wls; 
	wls.a_123 = vec3(.0f, .0f, .0f); 
	wls.b_x1y1_x2y2 = vec4(.0f, .0f, .0f, .0f); 
	return wls; 
}
void curve_2d_wls_fit_acc_params(
	uint windowSize, float l, vec2 dl,
	inout Curve2DWLSFitParams wls_params
){
	float wi = abs(l) / float(windowSize);
	
	float val = wi * l * l;
	wls_params.a_123[0] += val; // a1 <- a1 + wi * li^2
	val *= (0.5 * l);
	wls_params.a_123[1] += val; // a2 <- a2 + wi/2 * li^3
	val *= (0.5 * l);
	wls_params.a_123[2] += val; // a3 <- a3 + wi/4 * li^4
	
	dl = (wi * l) * dl;
	wls_params.b_x1y1_x2y2.xyzw += vec4(dl.xy, 0.5 * l * dl.xy); // bx1, by1, bx2, by2 
}

void curve_2d_wls_fit_iter(
	uint windowSize/* radius*2+1 */, bool left_samples, 
	vec2 coord_center_sample, vec2 coord_neigh_sample, 
	inout float l/*init as 0*/, inout Curve2DWLSFitParams fit_params
)
{
	vec2 dl = coord_neigh_sample - coord_center_sample;
	if (left_samples) l -= length(dl); // negative arc-length
	else l += length(dl);
	
	curve_2d_wls_fit_acc_params(windowSize, l, dl, /*inout*/fit_params);
}

void curve_2d_wls_fit_solve(Curve2DWLSFitParams fit_params, out vec2 tangent, out float curv)
{
#define A1  ((fit_params.a_123[0]))
#define A2  ((fit_params.a_123[1]))
#define A3  ((fit_params.a_123[2]))
#define BX1 ((fit_params.b_x1y1_x2y2.x))
#define BY1 ((fit_params.b_x1y1_x2y2.y))
#define BX2 ((fit_params.b_x1y1_x2y2.z))
#define BY2 ((fit_params.b_x1y1_x2y2.w))
	float dInv = 1.0f / ((A1 * A3) - (A2 * A2));
	vec2 dr = vec2(		// (dr/dx, dr/dy)
		(A3 * BX1) - (A2 * BX2),
		(A3 * BY1) - (A2 * BY2)
	);
	dr *= dInv;
	float dr_len = length(dr); 
	tangent = dr / dr_len;
	// Note: we can also output the curvature. 

	float ddInv = 1.0f / ((A2 * A2) - (A1 * A3)); 
	vec2 ddr = vec2(
		(A2 * BX1) - (A1 * BX2), 
		(A2 * BY1) - (A1 * BY2)
	); 
	ddr *= ddInv;
	curv = (dr.y * ddr.x - dr.x * ddr.y) / (dr_len * dr_len * dr_len);
#undef A1
#undef A2
#undef A3
#undef BX1
#undef BY1
#undef BX2
#undef BY2
}





#if defined(USE_CONTOUR_2D_SAMPLE_GEOMETRY_BUFFER)
void store_ssbo_contour_2d_sample_geometry__position(uint sample_id, vec2 P)
{ // stride : 1
	ssbo_contour_2d_sample_geometry_[sample_id] = pack_2d_uv(P); 
}
vec2 load_ssbo_contour_2d_sample_geometry__position(uint sample_id)
{
	return unpack_2d_uv(ssbo_contour_2d_sample_geometry_[sample_id]); 
}
void store_ssbo_contour_2d_sample_geometry__tangent(uint sample_id, vec2 tangent, uint num_samples)
{
	uint subbuff_offset = num_samples; 
	tangent = tangent * .5f + .5f; 
	ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id] = pack_2d_uv(tangent.xy); 
}
vec2 load_ssbo_contour_2d_sample_geometry__tangent(uint sample_id, uint num_samples)
{
	uint subbuff_offset = num_samples; 
	vec2 tangent = unpack_2d_uv(ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id]); 
	tangent = tangent * 2.0f - 1.0f; 
	return tangent; 
}
// subbuff 2, can be reused ----------------------------------
void store_ssbo_contour_2d_sample_geometry__curv_arclen_param(uint sample_id, float arclen_param, uint num_samples)
{ 
	uint subbuff_offset = 2u * num_samples; 
	ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id] = floatBitsToUint(arclen_param); 
}
float load_ssbo_contour_2d_sample_geometry__curv_arclen_param(uint sample_id, uint num_samples)
{
	uint subbuff_offset = 2u * num_samples; 
	return uintBitsToFloat(ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id]); 
}
// subbuff 3, can be reused ----------------------------------
void store_ssbo_contour_2d_sample_geometry__corner_curvature(uint sample_id, float corner_curvature, uint num_samples)
{ // stride : 1
	uint subbuff_offset = 3u * num_samples; 
	ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id] = floatBitsToUint(corner_curvature); 
}
float load_ssbo_contour_2d_sample_geometry__corner_curvature(uint sample_id, uint num_samples)
{
	uint subbuff_offset = 3u * num_samples; 
	return uintBitsToFloat(ssbo_contour_2d_sample_geometry_[subbuff_offset + sample_id]); 
}
#endif




#endif


