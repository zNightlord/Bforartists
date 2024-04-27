#ifndef BNPR_CONTOUR_GEOM__INCLUDED
#define BNPR_CONTOUR_GEOM__INCLUDED






// ----------------------------------------------------------
// Functions for Rasterization
// ----------------------------------------------------------
uniform vec2 pcs_screen_size_;

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


vec2 ndc_to_viewport(vec2 coord_ndc)
{
	vec2 coordVP =
		max(
			vec2(0.5, 0.5),
			min(
				coord_ndc * pcs_screen_size_.xy,
				pcs_screen_size_.xy - 1
			)
		);
	return coordVP;
}

vec2 viewport_to_ndc(vec2 coord_ss)
{
	vec2 coordNDC = coord_ss / pcs_screen_size_.xy;
	return saturate(coordNDC);
}


struct LineRasterResult
{
    uint frag_count; 
    vec4 begend_uvs; // note: must be range of 0-1
    vec2 begend_wclips; 
    
    bool is_outside_frust_line; 
    bool is_x_major_line; 
    bool beg_from_p0; // raster begins at the first vert on the original line 
}; 
void encode_line_raster_result(LineRasterResult res, out uvec4 d0123, out uint d4)
{
    uint d0 = frag_count; // 28 bits for frag_count should be sufficient
    d0 <<= 1u; 
    d0 |= (is_outside_frust_line ? 1u : 0u);
	d0 <<= 1u; 
	d0 |= (is_x_major_line ? 1u : 0u); 
	d0 <<= 1u; 
	d0 |= (beg_from_p0 ? 1u : 0u)

	uvec2 d12 = uvec2(0u); 
	begend_uvs = clamp(begend_uvs, vec4(.0f), vec4(1.0f)); 
	float fixed_factor = float((1u << 16u) - 1u); 
	uvec4 uv_fixed = uvec4(begend_uvs * fixed_factor) & 0xffffu; 
	d12 = uv_fixed.xy | (uv_fixed.zw << 16u); 

	uvec2 d34 = uvec2(0u);
	d34.x = floatBitsToUint(begend_wclips.x);
	d34.y = floatBitsToUint(begend_wclips.y);

	d123 = uvec4(d0, d12, d34.x);
	d4 = d34.y;
}

void raster_line_segment(
    vec4 vpos_0, vec4 vpos_1, vec4 vnor_0, vec4 vnor_1, 
    /* mat4 Matrix_M,  */mat4 Matrix_V, mat4 Matrix_P, 
    bool is_valid_thread
) {
    // OS --> WS (nothing here since all contours are in world space) 

    // WS --> VS
    // -----------------------------------------------------------------
    vpos_0 = mul(Matrix_V, vpos_0); // mul(Matrix_M, vpos_0));
    vpos_1 = mul(Matrix_V, vpos_1); // mul(Matrix_M, vpos_1));

    // VS --> HClip
    // --------------------------------------------------------------------
    vpos_0 = mul(Matrix_P, vpos_0); // HClip coord  
    vpos_1 = mul(Matrix_P, vpos_1); // HClip coord

    // HClip --> HClip(Clipped) 
    // -----------------------------------------------------------------------------
    bool reject, inside;
    vec2 interpolate = frustum_clip_line_segment(
        vpos_0, vpos_1, // in
        reject, inside  // out
    );
    vpos_0 = lerp(vpos_0, vpos_1, interpolate.x);
    vpos_1 = lerp(vpos_0, vpos_1, interpolate.y);
    vec2 whclip_v0v1 = vec2(vpos_0.w, vpos_1.w);

    // HClip(Clipped)->NDC
    // -----------------------------------------------------------------------------
    // .w component is preserved & for later use
    vpos_0.xyz = reject ? vec3(0, 0, 0) : hclip_to_ndc(vpos_0);
    vpos_1.xyz = reject ? vec3(0, 0, 0) : hclip_to_ndc(vpos_1);


    // NDC->Viewport
    // -------------------------------------------------------------------------
    vec4 frags_pos = // easier to swizzle this way
        vec4(ndc_to_viewport(vpos_0.xy).xy, ndc_to_viewport(vpos_1.xy).xy);

    // Compute fragment(== segment, in our context) count
    // Is edge X or Y major?
    // \          /  | A & C: X-Major
    //   \   B  /    | B & D: Y-Major 
    //     \  /      | 
    //  A   O    C   | 
    //    /    \     | 
    //  /   D    \   | 
    ///            \ | 
    vec2 dxdy = abs(frags_pos.zw - frags_pos.xy);
    bool is_x_major_line = (dxdy.y < dxdy.x);

    bvec2 is_p0_lower = frags_pos.xy < frags_pos.zw;
    bool beg_from_p0 = is_x_major_line ? is_p0_lower.x : is_p0_lower.y;

    frags_pos = beg_from_p0 ? frags_pos.xyzw : frags_pos.zwxy;

    // Compute fragment count according to X or Y major
    // Only accept contours that passes frustum clipping 
    uint fragCount = reject ? 0 : (uint) max(0, ceil(max(dxdy.x, dxdy.y) + 0.0001f));
    fragCount = (is_valid_thread) ? fragCount : 0;
}

    

#endif