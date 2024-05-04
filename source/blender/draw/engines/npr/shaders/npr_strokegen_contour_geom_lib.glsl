#ifndef BNPR_CONTOUR_GEOM__INCLUDED
#define BNPR_CONTOUR_GEOM__INCLUDED

uvec3 pack_u24_x4(uvec4 u)
{
    uvec3 p = uvec3(0u);
    p.xyz = (u.xyz << 8);

    p.x |= (u.w & 0x000000ff);
    u.w >>= 8;
    p.y |= (u.w & 0x000000ff);
    u.w >>= 8;
    p.z |= (u.w & 0x000000ff);

    return p;
}

uvec4 unpack_u24_x4(uvec3 p)
{
    uvec4 u = uvec4(0u);
    u.xyz = (p.xyz >> 8);
    p.xyz &= 0x000000ff;
    u.w = ((p.z << 16) | (p.y << 8) | p.x);

    return u;
}


// ----------------------------------------------------------
// Functions for Rasterization
// ----------------------------------------------------------
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
// #if UNITY_UV_STARTS_AT_TOP
// 	// Our world space, view space, screen space and NDC space are Y-up.
// 	// Our clip space is flipped upside-down due to poor legacy Unity design.
// 	// The flip is baked into the projection matrix, so we only have to flip
//         // manually when going from CS to NDC and back.
//         pos_cs.y = -pos_cs.y;
// #endif

	pos_cs /= (pos_cs.w);
	pos_cs.xy = pos_cs.xy * 0.5 + 0.5;

	return pos_cs.xyz;
}

vec2 frustum_clip_line_segment(
	vec4 vpos_hcs_0, vec4 vpos_hcs_1, // honogenous clip space positions
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
    vec4 begend_uvs; // note: must be range of 0-1
    vec2 begend_wclips; 
    
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
	bool no_occlusion; 
	bool is_contour_header;
	// d123, 24 bits for each
	vec2 begend_ratios; // range [0, 1], the ratio of start/end of this segment at the parent contour
	uint prev_frag_seg_head_id; 
	uint next_frag_seg_head_id; 
}; 

void init_contour_visibility_split_info(uint contour_id, inout ContourVisibilitySplitInfo cvsi)
{
	cvsi.parent_contour_id = contour_id; 
	cvsi.is_new_contour = false; 
	cvsi.is_visible = true; 
	cvsi.no_occlusion = true;
	cvsi.is_contour_header = true;  
	
	cvsi.begend_ratios = vec2(.0f, 1.0f); 
	cvsi.prev_frag_seg_head_id = 0u;
	cvsi.next_frag_seg_head_id = 0u;
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
	d0123.x |= (cvsi.no_occlusion    	? 1u : 0u);
	d0123.x <<= 1u;
	d0123.x |= (cvsi.is_contour_header 	? 1u : 0u);

	uvec2 begend_ratios_fixed; 
	{ // Pack into 2x24 bits
		cvsi.begend_ratios = clamp(cvsi.begend_ratios, vec2(.0f), vec2(1.0f));
		float fixed_factor = float((1u << 24u) - 1u);
		begend_ratios_fixed = uvec2(round(cvsi.begend_ratios * fixed_factor)) & 0x00ffffffu; 
	}
	d0123.yzw = pack_u24_x4(
		uvec4(
			begend_ratios_fixed.x, 
			begend_ratios_fixed.y, 
			cvsi.prev_frag_seg_head_id & 0x00ffffffu, 
			cvsi.next_frag_seg_head_id & 0x00ffffffu
		)
	); 

	return d0123; 
}

ContourVisibilitySplitInfo decode_contour_visibility_split_info(uvec4 d0123)
{
	ContourVisibilitySplitInfo cvsi; 
	cvsi.is_contour_header = ((d0123.x & 0x1u) != 0u);
	cvsi.no_occlusion      = ((d0123.x & 0x2u) != 0u);
	cvsi.is_visible        = ((d0123.x & 0x4u) != 0u); 
	cvsi.is_new_contour    = ((d0123.x & 0x8u) != 0u); 
	cvsi.parent_contour_id = (d0123.x >> 4u); 

	uvec4 d123_dec = unpack_u24_x4(d0123.yzw);
	cvsi.begend_ratios = vec2(d123_dec.xy) / float((1u << 24u) - 1u);
	cvsi.prev_frag_seg_head_id = d123_dec.z;
	cvsi.next_frag_seg_head_id = d123_dec.w;

	return cvsi; 
}

#endif


