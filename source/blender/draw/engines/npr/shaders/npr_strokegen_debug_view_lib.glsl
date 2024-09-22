 

#ifndef BNPR_STROKEGEN_DEBUG_VIEW__INCLUDED
#define BNPR_STROKEGEN_DEBUG_VIEW__INCLUDED

/* https://www.shadertoy.com/view/7tlXR4 */
vec3 hash32(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}
vec3 hsl2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
vec3 rand_col_rgb(uint seed0, uint seed1)
{
    vec3 rnd = hash32(vec2(float(seed0 * 17), float(seed1 * 33)));

    float hue = rnd.x;
    float saturation = 0.6 + rnd.z*0.4;
    float luminosity  = 0.6 + rnd.y*0.4;

    return hsl2rgb(vec3(hue, saturation, luminosity));
}
/* https://www.pcg-random.org/ */
uint pcg(uint v) /* pcg hash */
{
	uint state = v * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}
/* Hashing 3d vert position
 * Converting 1D hashing func into ND via nesting.
 * this is better than linear combination of N hashes, 
 * see Ch 4.2 of https://jcgt.org/published/0009/03/02/paper.pdf */
uint pcg_nested_3d(vec3 v)
{
    uvec3 v3 = floatBitsToUint(v);
    return pcg(v3.z + pcg(v3.y + pcg(v3.x))); 
}



#if defined(INCLUDE_DEBUG_LINE_CONFIG)
    /* Definitions for Debugging Lines */
    #define DBG_LINE_TYPE__GENERAL 2u

    uint get_debug_line_counter(uint line_type)
    {
        // if (line_type == DBG_LINE_TYPE__VNOR)
        //     return ssbo_bnpr_mesh_pool_counters_.num_dbg_vnor_lines;
        if (line_type == DBG_LINE_TYPE__GENERAL)
            return ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines;
        // else if (line_type == DBG_LINE_TYPE__EDGES)
        //     return ssbo_bnpr_mesh_pool_counters_.num_dbg_edge_lines; 
        return 0u; 
    }

    uint get_debug_line_offset(uint line_type)
    {
        /* memory offset line_id based on line type */
        uint line_offset = 0u; 
        // if (line_type == DBG_LINE_TYPE__VNOR)
        //     return line_offset; 
        // line_offset += get_debug_line_counter(DBG_LINE_TYPE__VNOR); 

        if (line_type == DBG_LINE_TYPE__GENERAL)
            return line_offset; 
        line_offset += get_debug_line_counter(DBG_LINE_TYPE__GENERAL); 

        // if (line_type == DBG_LINE_TYPE__EDGES)
        //     return line_offset;
        // line_offset += get_debug_line_counter(DBG_LINE_TYPE__EDGES);
        
        return line_offset; 
    }

    uint pack_r11_g11_b10(vec3 x)
    {
        const vec3 normalization = vec3(
            float(1u << 11u) - 1.0f, 
            float(1u << 11u) - 1.0f, 
            float(1u << 10u) - 1.0f
        ); 
        x = clamp(x, vec3(.0f), vec3(1.0f)); 
        x = round(x * normalization); 
        
        uvec3 x_u = uvec3(x.rgb);
    #define PACK_X_MASK 0x000007ffu
    #define PACK_Y_MASK 0x000007ffu 
    #define PACK_Z_MASK 0x000003ffu 
        uint enc = x_u.r & PACK_X_MASK;
        enc <<= 11u;
        enc |= (x_u.g & PACK_Y_MASK); 
        enc <<= 10u; 
        enc |= (x_u.b & PACK_Z_MASK); 

        return enc; 
    }
    vec3 unpack_r11_g11_b10(uint enc)
    {
        uvec3 x_u; 
        x_u.b = enc & PACK_Z_MASK; 
        enc >>= 10u; 
        x_u.g = enc & PACK_Y_MASK; 
        enc >>= 11u; 
        x_u.r = enc & PACK_X_MASK; 
    #undef PACK_Z_MASK
    #undef PACK_Y_MASK
    #undef PACK_X_MASK

        const vec3 normalization = vec3(
            float(1u << 11u) - 1.0f, 
            float(1u << 11u) - 1.0f, 
            float(1u << 10u) - 1.0f
        ); 
        vec3 x = vec3(x_u) / normalization; 
        return clamp(x, vec3(.0f), vec3(1.0f)); 
    }

    struct DebugVertData
    {
        vec3 pos; /* world space vertex pos */
        vec3 col;
        uvec4 dbg_data; 
// Special flag values for dbg_data ---
// handle overlapped line with different layers of depth bias
#define DBG_LINE_SPEC_DATA_X__BOTTOM_LAYER 0xfefefefeu 
// ------------------------------------
    }; 
    void encode_debug_vert_data(DebugVertData vtx, out uvec4 enc_0, out uvec4 enc_1)
    {
        enc_0 = uvec4(floatBitsToUint(vtx.pos), pack_r11_g11_b10(vtx.col)); 
        enc_1 = vtx.dbg_data; 
    }
    DebugVertData decode_debug_vert_data(uvec4 enc_0, uvec4 enc_1)
    {
        DebugVertData vtx; 
        vtx.pos = uintBitsToFloat(enc_0.xyz); 
        vtx.col = unpack_r11_g11_b10(enc_0.w); 
        vtx.dbg_data = enc_1; 
        return vtx; 
    }

    #if defined(INCLUDE_DEBUG_LINE_CONFIG_LOAD_STORE)
        void store_debug_line_data(uint dbg_line_id, DebugVertData v0, DebugVertData v1)
        {
            uvec4 enc_0, enc_1;
            encode_debug_vert_data(v0, /*out*/enc_0, enc_1); 
            Store4(ssbo_dbg_lines_, dbg_line_id*4u,      enc_0); 
            Store4(ssbo_dbg_lines_, dbg_line_id*4u + 1u, enc_1); 

            encode_debug_vert_data(v1, /*out*/enc_0, enc_1); 
            Store4(ssbo_dbg_lines_, dbg_line_id*4u + 2u, enc_0); 
            Store4(ssbo_dbg_lines_, dbg_line_id*4u + 3u, enc_1); 
        }
        DebugVertData load_debug_vtx_data(uint dbg_line_id, uint drw_vtx_id)
        {
            uint vtx_offset = dbg_line_id * 4u + (drw_vtx_id % 2u) * 2u;
            uvec4 enc_0, enc_1; 
            Load4(ssbo_dbg_lines_, vtx_offset + 0u, enc_0); 
            Load4(ssbo_dbg_lines_, vtx_offset + 1u, enc_1); 
            return decode_debug_vert_data(enc_0, enc_1); 
        }
        void store_debug_line_color(uint dbg_line_id, vec3 col)
        {
            uint enc_0_w = pack_r11_g11_b10(col);

            uint st_addr = dbg_line_id * 16u + 3u; // x16 u32s per line 
            ssbo_dbg_lines_[st_addr] = enc_0_w; // store color in vtx0
            st_addr += 8u; 
            ssbo_dbg_lines_[st_addr] = enc_0_w; // store color in vtx1 
        }
    #endif
#endif


#endif


