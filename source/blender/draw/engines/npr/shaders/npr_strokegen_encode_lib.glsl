#ifndef BNPR_ENCODE_LIB_INCLUDED
#define BNPR_ENCODE_LIB_INCLUDED

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

uint pack_2d_uv(vec2 uv)
{
	float fixed_factor = float((1u << 16u) - 1u); 
	uvec2 uv_fixed = uvec2(round(uv * fixed_factor)); 
	return ((uv_fixed.y << 16u) | uv_fixed.x);
}

vec2 unpack_2d_uv(uint uv_packed)
{
	float fixed_factor = float((1u << 16u) - 1u); 
	uvec2 uv_fixed = uvec2(uv_packed & 0xffffu, uv_packed >> 16u); 
	return clamp(vec2(uv_fixed) / fixed_factor, vec2(.0f), vec2(1.0f));
}

uvec2 pack_2d_uv_x2(vec4 begend_uvs)
{
	begend_uvs = clamp(begend_uvs, vec4(.0f), vec4(1.0f)); 

	float fixed_factor = float((1u << 16u) - 1u); 
	uvec4 uv_fixed = uvec4(round(begend_uvs * fixed_factor)) & 0xffffu; 
	return ((uv_fixed.zw << 16u) | uv_fixed.xy); 
}

vec4 unpack_2d_uv_x2(uvec2 uv_packed)
{
	uvec4 uv_fixed = uvec4(
		uv_packed & 0xffffu, 
		uv_packed >> 16u
	); 
	float fixed_factor = float((1u << 16u) - 1u); 
	return clamp(vec4(uv_fixed) / fixed_factor, vec4(.0f), vec4(1.0f)); 
}

#endif


