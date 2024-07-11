#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_linear_algebra_lib.glsl)


#ifndef BNPR_MESHING_GEOM__INCLUDED
#define BNPR_MESHING_GEOM__INCLUDED


#if defined(INCLUDE_VERTEX_POSITION)
    vec3 ld_vpos(uint vtx_id)
    {
        vec3 vpos; 
        Load3(ssbo_vbo_full_, vtx_id, vpos);

        return vpos; 
    }
    void st_vpos(uint vtx_id, vec3 vpos)
    {
        Store3(ssbo_vbo_full_, vtx_id, vpos);
    }
#endif

#if defined(INCLUDE_VERTEX_NORMAL)
    vec3 ld_vnor(uint vtx_id)
    {
        uvec3 vnor_enc; 
        Load3(ssbo_vnor_, vtx_id, vnor_enc);

        return uintBitsToFloat(vnor_enc); 
    }
    void st_vnor(uint vtx_id, vec3 vnor)
    {
        uvec3 vnor_enc = floatBitsToUint(vnor);
        Store3(ssbo_vnor_, vtx_id, vnor_enc);
    }

    #if defined(INCLUDE_VERTEX_RADIAL_NORMAL)
        vec3 ld_vnor_radial(uint vid, vec3 vpos)
        {
            vec3 n = ld_vnor(vid); 

            mat4 view_to_world = ubo_view_matrices_.viewinv;
            bool is_persp = (ubo_view_matrices_.winmat[3][3] == 0.0);
            vec3 cam_pos_ws = view_to_world[3].xyz; /* see "#define cameraPos ViewMatrixInverse[3].xyz" */
            
            vec3 v = cam_pos_ws - vpos; 
            vec3 v_ = normalize(v); 

            // vec3 n_view_proj = normalize(n - dot(n, v_) * v_); 
            // if (dot(n_view_proj, n) < .0f) n_view_proj = -n_view_proj; 

            // return n_view_proj; 

            vec3 bn = normalize(cross(n, v));
            vec3 t = normalize(cross(bn, n)); 
            vec3 n_view_proj = normalize(cross(t, bn));

            vec3 res = t; 
            float cos_theta = dot(v_, n); 
            res = cos_theta < -.01f ? -res : res; 
            if (abs(cos_theta) <= .01f) res = n_view_proj; 

            return n_view_proj; 
        }
    #endif
#endif



#if defined(INCLUDE_VERTEX_VORONOI_AREA)
    float ld_varea(uint vtx_id)
    {
        uint varea_enc = ssbo_varea_[vtx_id]; 
        return uintBitsToFloat(varea_enc); 
    }
    void st_varea(uint vtx_id, float varea)
    {
        uint varea_enc = floatBitsToUint(varea);
        ssbo_varea_[vtx_id] = varea_enc;
    }
#endif

#if defined(INCLUDE_VERTEX_CURV_TENSOR)

    // vec3 ld_vcurv_tensor(uint vtx_id)
    // {
    //     uvec3 vcurv_tensor_enc;
    //     Load3(ssbo_vcurv_tensor_, vtx_id, vcurv_tensor_enc); 
    //     return uintBitsToFloat(vcurv_tensor_enc); 
    // }
    // void st_vcurv_tensor(uint vtx_id, vec3 vcurv_tensor)
    // {
    //     uvec3 vcurv_tensor_enc = floatBitsToUint(vcurv_tensor);
    //     Store3(ssbo_vcurv_tensor_, vtx_id, vcurv_tensor_enc); 
    // }

    void st_vcurv_pdirs_k1k2(uint vtx_id, vec3 pd1, float k1, vec3 pd2, float k2)
    {
        uvec4 pd1_k1_enc = floatBitsToUint(vec4(pd1.xyz, k1));
        Store4(ssbo_vcurv_pdirs_k1k2_, (vtx_id*2u), pd1_k1_enc); 

        uvec4 pd2_k2_enc = floatBitsToUint(vec4(pd2.xyz, k2));
        Store4(ssbo_vcurv_pdirs_k1k2_, (vtx_id*2u+1u), pd2_k2_enc); 
    }
    void ld_vcurv_pdir1_k1(uint vtx_id, out vec3 pdir, out float k)
    {
        uvec4 pd1_k1_enc;
        Load4(ssbo_vcurv_pdirs_k1k2_, (vtx_id*2u), pd1_k1_enc); 
        
        vec4 pd1_k1 = uintBitsToFloat(pd1_k1_enc); 
        pdir = pd1_k1.xyz;
        k    = pd1_k1.w;
    }
    void ld_vcurv_pdir2_k2(uint vtx_id, out vec3 pdir, out float k)
    {
        uvec4 pd2_k2_enc;
        Load4(ssbo_vcurv_pdirs_k1k2_, (vtx_id*2u+1u), pd2_k2_enc); 
        
        vec4 pd2_k2 = uintBitsToFloat(pd2_k2_enc); 
        pdir = pd2_k2.xyz;
        k    = pd2_k2.w;
    }
#endif

bool valid_vcurv_max(float vcurv_max)
{
    return vcurv_max >= 0.0f; 
}

#if defined(INCLUDE_VERTEX_CURV_MAX)
    float ld_vcurv_max(uint vtx_id)
    {
        uint vcurv_max_enc = ssbo_vcurv_max_[vtx_id]; 
        return uintBitsToFloat(vcurv_max_enc); 
    }
    void st_vcurv_max(uint vtx_id, float vcurv_max)
    {
        uint vcurv_max_enc = floatBitsToUint(vcurv_max);
        ssbo_vcurv_max_[vtx_id] = vcurv_max_enc;
    }
    void ld_vcurv_max_with_cusp(uint vtx_id, out float vcurv_max, out float cusp_func)
    {
        uint vcurv_max_enc = ssbo_vcurv_max_[vtx_id*2u]; 
        vcurv_max = uintBitsToFloat(vcurv_max_enc); 

        uint cusp_enc = ssbo_vcurv_max_[vtx_id*2u+1u];
        cusp_func = uintBitsToFloat(cusp_enc);
    }
    void st_vcurv_max_with_cusp(uint vtx_id, float vcurv_max, float cusp_func)
    {
        ssbo_vcurv_max_[vtx_id*2u] = floatBitsToUint(vcurv_max);
        ssbo_vcurv_max_[vtx_id*2u+1u] = floatBitsToUint(cusp_func);
    }
#endif

#if defined(INCLUDE_VERTEX_REMESH_LEN)
#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT) || defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE)
    #define ssbo_vtx_remesh_len_ ssbo_per_vert_collapse_wedge_id_ // ssbo slots are full (x16)
#endif
    void st_vtx_remesh_len(uint vtx_id, float len)
    {
        uint len_enc = floatBitsToUint(len);
        ssbo_vtx_remesh_len_[vtx_id] = len_enc;
    }
    float ld_vtx_remesh_len(uint vtx_id)
    {
        uint len_enc = ssbo_vtx_remesh_len_[vtx_id];
        return uintBitsToFloat(len_enc);
    }
#endif





#if defined(INCLUDE_DEBUG_LINE_CONFIG)
    /* Definitions for Debugging Lines */
    #define DBG_LINE_TYPE__VNOR  0u
    #define DBG_LINE_TYPE__GENERAL 1u
    #define DBG_LINE_TYPE__EDGES 2u

    uint get_debug_line_counter(uint line_type)
    {
        if (line_type == DBG_LINE_TYPE__VNOR)
            return ssbo_bnpr_mesh_pool_counters_.num_dbg_vnor_lines;
        else if (line_type == DBG_LINE_TYPE__GENERAL)
            return ssbo_bnpr_mesh_pool_counters_.num_dbg_general_lines;
        else if (line_type == DBG_LINE_TYPE__EDGES)
            return ssbo_bnpr_mesh_pool_counters_.num_dbg_edge_lines; 
        return 0u; 
    }

    uint get_debug_line_offset(uint line_type)
    {
        /* memory offset line_id based on line type */
        uint line_offset = 0u; 
        if (line_type == DBG_LINE_TYPE__VNOR)
            return line_offset; 
        line_offset += get_debug_line_counter(DBG_LINE_TYPE__VNOR); 

        if (line_type == DBG_LINE_TYPE__GENERAL)
            return line_offset; 
        line_offset += get_debug_line_counter(DBG_LINE_TYPE__GENERAL); 

        if (line_type == DBG_LINE_TYPE__EDGES)
            return line_offset;
        line_offset += get_debug_line_counter(DBG_LINE_TYPE__EDGES);
        
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
    }; 
    uvec4 encode_debug_vert_data(DebugVertData vtx)
    {
        return uvec4(floatBitsToUint(vtx.pos), pack_r11_g11_b10(vtx.col)); 
    }
    DebugVertData decode_debug_vert_data(uvec4 enc)
    {
        DebugVertData vtx; 
        vtx.pos = uintBitsToFloat(enc.xyz); 
        vtx.col = unpack_r11_g11_b10(enc.w); 
        return vtx; 
    }

    #if defined(INCLUDE_DEBUG_LINE_CONFIG_LOAD_STORE)
        void store_debug_line_data(uint dbg_line_id, DebugVertData v0, DebugVertData v1)
        {
            uvec4 enc_v0 = encode_debug_vert_data(v0); 
            Store4(ssbo_dbg_lines_, dbg_line_id*2u,    enc_v0); 
            uvec4 enc_v1 = encode_debug_vert_data(v1); 
            Store4(ssbo_dbg_lines_, dbg_line_id*2u+1u, enc_v1); 
        }
        DebugVertData load_debug_vtx_data(uint dbg_line_id, uint drw_vtx_id)
        {
            uvec4 enc; 
            Load4(ssbo_dbg_lines_, (dbg_line_id*2u + (drw_vtx_id%2u)), enc); 
            return decode_debug_vert_data(enc); 
        }
    #endif
#endif



// Subdivision Surfaces -------------------------------------------------------
float sqrt3_vpos_smooth_weight(float n)
{
    return (4.0f - 2.0f * cos((2.0f * 3.14159265359f) / n)) / 9.0f;
}

float sqrt3_vpos_limit_weight(float n, int step, bool step_inf = false)
{
    float alpha_n = sqrt3_vpos_smooth_weight(n); 
    alpha_n *= 3.0f; 
    float beta_n_inf = alpha_n / (1.0f + alpha_n); 

    if (step_inf) return beta_n_inf; 

    float gamma_n_step = pow(2.0f/3.0f - alpha_n, step); 
    return gamma_n_step; 
}

vec3 sqrt3_limit_vpos(vec3 vpos, vec3 vpos_sum, float n, int step, bool step_inf)
{
    float w = sqrt3_vpos_limit_weight(n, step, step_inf); 
    vec3 lim_pos_inf = (1.0f - w) * vpos + w * (vpos_sum / n); 
    if (step_inf) return lim_pos_inf; 
    return w * vpos + (1.0f - w) * lim_pos_inf;
}


float loop_vpos_smooth_weight(float n, bool crease = false, bool corner = false)
{ /* mask for interior verts in loop subdiv */
    /* see "Subdivision methods for geometric design - Ch.7.3 Smooth Subdivision for Triangle Meshes" */
    if (corner) return .0f; 
    if (crease) return .25f; 

    if (n == 6.0f) return 3.0f / 8.0f; // fast path, actually same as the formula below

    float alpha_n = (3.0f/8.0f) + (cos((2.0f * 3.14159265359f) / n) / 4.0f);
    alpha_n = (5.0f/8.0f) - alpha_n * alpha_n;
    return alpha_n; 
}




/* Scalable Curvature Computations on Triangle Meshes ----------------------- */
#ifndef M_PI
#define M_PI           (3.1415926535897932384626433832795)
#endif
#ifndef M_PI_2
#define M_PI_2         (1.5707963267948966192313216916398)
#endif
// struct SphericalTriangle ----------------------------------------------------
// https://www.cs.cmu.edu/~kmcrane/Projects/ScalableCurvature/
/// @param[in] a the first vertex of spherical triangle ABC (ccw oriented)
/// @param[in] b the second vertex of spherical triangle ABC (ccw oriented)
/// @param[in] c the third vertex of spherical triangle ABC (ccw oriented)
/// @return 'true' if ABC is made of almost collinear points.
bool isDegenerate(vec3 a, vec3 b, vec3 c)
{
    vec3 d = { length(a - b),length(a - c),length(b - c) };
    // Checks that the spherical triangle is small or thin.
    if ( ( d[0] < 1e-8f ) || ( d[1] < 1e-8f ) || ( d[2] < 1e-8f ) )
        return true;
    // Checks that the spherical triangle is flat.
    uint m = 0;
    if ( d[ 1 ] > d[ m ] ) m = 1;
    if ( d[ 2 ] > d[ m ] ) m = 2;
    return ( abs( d[ m ] - d[ (m+1)%3 ] - d[ (m+2)%3 ] ) < 1e-8f );
}

/// Computes the polar triangle associated with triangle ABC.
/// @param[in] a the first vertex of spherical triangle ABC (ccw oriented)
/// @param[in] b the second vertex of spherical triangle ABC (ccw oriented)
/// @param[in] c the third vertex of spherical triangle ABC (ccw oriented)
/// @param[out] Ap the first vertex of its polar triangle A'B'C'
/// @param[out] Bp the second vertex of its polar triangle A'B'C'
/// @param[out] Cp the third vertex of its polar triangle A'B'C'
void polarTriangle(in vec3 a, in vec3 b, in vec3 c,
                   out vec3 Ap, out vec3 Bp, out vec3 Cp)
{
    Ap = cross(b, c);
    Bp = cross(c, a);
    Cp = cross(a, b);
    // Reorient points.
    if ( dot(Ap, a) < 0.0f ) Ap = -Ap;
    if ( dot(Bp, b) < 0.0f ) Bp = -Bp;
    if ( dot(Cp, c) < 0.0f ) Cp = -Cp;
}

/// Computes the interior angles of the spherical triangle ABC.
/// @param[in] a the first vertex of spherical triangle ABC (ccw oriented)
/// @param[in] b the second vertex of spherical triangle ABC (ccw oriented)
/// @param[in] c the third vertex of spherical triangle ABC (ccw oriented)
/// @param[out] alpha the interior angle at vertex A.
/// @param[out] beta  the interior angle at vertex B.
/// @param[out] gamma the interior angle at vertex C.
void interiorAngles(in vec3 a, in vec3 b, in vec3 c,
                    out float alpha, out float beta, out float gamma )
{
    vec3 Ta, Tb, Tc;
    polarTriangle(a, b, c, Ta, Tb, Tc);
    Ta /= length(Ta);
    Tb /= length(Tb);
    Tc /= length(Tc);
    if ( all(Ta.xyz == .0f) || all(Tb.xyz == .0f) || all(Tc.xyz == .0f) )
        alpha = beta = gamma = 0.0; // TODO: change ""==.0f" to "<1e-10f"
    else
    {
        float ca = max( -1.0, min( 1.0, dot(Tb, Tc) ) );
        float cb = max( -1.0, min( 1.0, dot(Tc, Ta) ) );
        float cc = max( -1.0, min( 1.0, dot(Ta, Tb) ) );
        alpha     = float(acos( (ca) ));
        beta      = float(acos( (cb) ));
        gamma     = float(acos( (cc) ));
    }
}

/// @param[in] a the first vertex of spherical triangle ABC (ccw oriented)
/// @param[in] b the second vertex of spherical triangle ABC (ccw oriented)
/// @param[in] c the third vertex of spherical triangle ABC (ccw oriented)
/// @return the (unsigned) area of the spherical triangle (below 2pi).
float area(in vec3 a, in vec3 b, in vec3 c)
{
    float alpha, beta, gamma;
    if ( isDegenerate(a,b,c) ) return 0.0;
    interiorAngles( a, b, c, alpha, beta, gamma );
    return ( (alpha == 0.0) || (beta == 0.0) || (gamma == 0.0) )
        ? 0.0 : (2.0*M_PI - alpha - beta - gamma);
}

float l1Norm(in vec3 v)
{
    return abs(v.x) + abs(v.y) + abs(v.z);
} 
/// @param[in] a the first vertex of spherical triangle ABC (ccw oriented)
/// @param[in] b the second vertex of spherical triangle ABC (ccw oriented)
/// @param[in] c the third vertex of spherical triangle ABC (ccw oriented)
/// @return the (signed) area of the spherical triangle (below 2pi).
float algebraicArea(in vec3 a, in vec3 b, in vec3 c) 
{
    float S = area(a,b,c);
    vec3 M = a + b + c;
    vec3 X = cross(( b - a ), ( c - a )); 
    // if ( M.lpNorm<1>() <= 1e-8 || X.lpNorm<1>() <= 1e-8 ) return 0.0;
    // TODO: dangerous here, need to check
    if ( l1Norm(M) <= 1e-8f || l1Norm(X) <= 1e-8f ) return 0.0;
    return dot(M, X) < 0.0 ? -S : S;
}
// struct SphericalTriangle----------------------------------------------------





///---------------------- Main functions ----------------

/// Computes mu0 measure (area) of triangle abc given an interpolated
/// corrected normal vector \a ua, \a \ub, \a uc.
/// @param a any point
/// @param b any point
/// @param c any point
/// @param ua the corrected normal vector at point a
/// @param ub the corrected normal vector at point b
/// @param uc the corrected normal vector at point c
/// @param unit_u when 'true' considers that interpolated
/// corrected normals should be made unitary, otherwise
/// interpolated corrected normals may have smaller norms.
/// @return the mu0-measure of triangle abc, i.e. its area.
float mu0InterpolatedU(in vec3 a,
                        in vec3 b,
                        in vec3 c,
                        in vec3 ua,
                        in vec3 ub,
                        in vec3 uc,
                        bool unit_u = false)
{
    // MU0=1/2*det( uM, B-A, C-A )
    //    =  1/2 < ( (u_A + u_B + u_C)/3.0 ) | (AB x AC ) >
    vec3 uM = ( ua+ub+uc ) / 3.0f;
    if ( unit_u )
    {
        float uM_norm = length(uM);
        uM = uM_norm == 0.0 ? uM : uM / uM_norm;
    }
    return 0.5f * dot((cross(( b - a ), ( c - a ))), uM );
}

 /// Computes mu1 measure (mean curvature) of triangle abc given an interpolated
/// corrected normal vector \a ua, \a \ub, \a uc.
/// @param a any point
/// @param b any point
/// @param c any point
/// @param ua the corrected normal vector at point a
/// @param ub the corrected normal vector at point b
/// @param uc the corrected normal vector at point c
/// @param unit_u when 'true' considers that interpolated
/// corrected normals should be made unitary, otherwise
/// interpolated corrected normals may have smaller norms.
/// @return the mu1-measure of triangle abc, i.e. \b twice its mean curvature.
float mu1InterpolatedU(in vec3 a,
                    in vec3 b,
                    in vec3 c,
                    in vec3 ua,
                    in vec3 ub,
                    in vec3 uc,
                    bool unit_u = false)
  {
    // MU1=1/2( | uM u_C-u_B A | + | uM u_A-u_C B | + | uM u_B-u_A C |
    vec3 uM = ( ua+ub+uc ) / 3.0f;
    if ( unit_u ) uM /= length(uM);
    return 0.5f * dot(cross(uM, ( uc - ub )), a)
                   + dot(cross(uM, ( ua - uc )), b)
                   + dot(cross(uM, ( ub - ua )), c);
  }

    
  /// Computes mu2 measure (Gaussian curvature) of triangle abc given an interpolated
  /// corrected normal vector \a ua, \a \ub, \a uc.
  /// @param a any point
  /// @param b any point
  /// @param c any point
  /// @param ua the corrected normal vector at point a
  /// @param ub the corrected normal vector at point b
  /// @param uc the corrected normal vector at point c
  /// @param unit_u when 'true' considers that interpolated
  /// corrected normals should be made unitary, otherwise
  /// interpolated corrected normals may have smaller norms.
  /// @return the mu2-measure of triangle abc, i.e. its Gaussian curvature.
  float mu2InterpolatedU(in vec3 a,
                    in vec3 b,
                    in vec3 c,
                    in vec3 ua,
                    in vec3 ub,
                    in vec3 uc,
                    bool unit_u = false)
  {
    
    // Using non unitary interpolated normals give
    // MU2=1/2*det( uA, uB, uC )
    // When normals are unitary, it is the area of a spherical triangle.
    if ( unit_u )
      return algebraicArea(ua,ub,uc);
    else
      return 0.5f * dot(cross(ua, ub), uc);
  }

/// Computes muXY measure (anisotropic curvature) of triangle abc
/// given an interpolated corrected normal vector \a ua, \a \ub, \a
/// uc.
///
/// @param a any point
/// @param b any point
/// @param c any point
/// @param ua the corrected normal vector at point a
/// @param ub the corrected normal vector at point b
/// @param uc the corrected normal vector at point c
/// @return the muXY-measure of triangle abc, i.e. its anisotropic curvature.
mat3 muXYInterpolatedU( in vec3 a,
                        in vec3 b,
                        in vec3 c,
                        in vec3 ua,
                        in vec3 ub,
                        in vec3 uc,
                        bool unit_u = false)
{
    mat3 T = mat3(.0f, .0f, .0f, .0f, .0f, .0f, .0f, .0f, .0f);
    vec3 uM = ( ua+ub+uc ) / 3.0;
    if ( unit_u ) uM /= length(uM);
    const vec3 uac = uc - ua;
    const vec3 uab = ub - ua;
    const vec3  ab = b - a;
    const vec3  ac = c - a;
    for ( uint i = 0u; i < 3u; ++i )
    {
        vec3 X = vec3(.0f);
        X[i] = 1.0 ;
        for ( uint j = 0u; j < 3u; ++j )
        {
            // Since RealVector Y = RealVector::base( j, 1.0 );
            // < Y | uac > = uac[ j ]
            const float tij =
            0.5 * dot(uM, ( uac[ j ] * cross(X, ab)
                           - uab[ j ] * cross(X, ac) ));
            T[i][j] = tij;
        }
    }
    return T;
}


///---------------------- Helper functions ----------------
///
///
/// Computing principal curvatures k1 and k2 from tensor
/// @param tensor The muXY integrated tensor
/// @param area Area of the face
/// @param N the normal vector
/// @return a pair of principal directions.
void curvDirFromTensor(in mat3 tensor,
                       in float area,
                       in vec3 N, 
                       out mat3 pdirs, 
                       out vec3 evals)
{
    mat3 Mt = transpose(tensor);
    mat3 M = tensor;
    M += Mt;
    M *= 0.5f;
    const float coef_N = 1000.0f * area;
    // Adding regulated term == 1000 area n x n 
    // force the principal direction eigenvectors to be tangential to the surface
    // (see @cite lachaud2020interpolated, section 2) 
    for ( uint j = 0u; j < 3u; j++ )
        for ( uint k = 0u; k < 3u; k++ )
            M[ j ][ k ] += coef_N * N[ j ] * N[ k ];

    pdirs = mat3(.0f); 
    evals = vec3(.0f); 
    dsyevv3(M, pdirs, evals); // TODO: we might need to transpose since dsyevv3 is for hlsl row major matrix


    // Sort eigenvalues and eigenvectors
#define SWAP_EVAL_EVEC(i, j) \
    if (evals[i] > evals[j])      \
    {                             \
        float temp = evals[i];    \
        evals[i] = evals[j];      \
        evals[j] = temp;          \
                                \
        vec3 temp_vec = vec3(pdirs[0][i], pdirs[1][i], pdirs[2][i]); \
        pdirs[0][i] = pdirs[0][j];   \
        pdirs[1][i] = pdirs[1][j];   \
        pdirs[2][i] = pdirs[2][j];   \
        pdirs[0][j] = temp_vec[0];   \
        pdirs[1][j] = temp_vec[1];   \
        pdirs[2][j] = temp_vec[2];   \
    }                                \

    SWAP_EVAL_EVEC(0, 1);
    SWAP_EVAL_EVEC(1, 2); 
    SWAP_EVAL_EVEC(0, 1);
#undef SWAP_EVAL_EVEC

}

/* Compute cusp function *
 * Cusp detection from "Illustrating smooth surface" by Hertzmann et al.
 * Surface point p is a cusp when the tangent to the silhouette at p is parallel to the viewing direction v. 
 * 
 * using min&max pricipal dirs w1, w2 & curvatures k1, k2, 
 * the local vicinity of a contour point P is parameterized as a quadartic surface,
 * and the cusp function C(P) = k1(dot(v, w1))^2 + k2(dot(v, w2))^2, where v = p - camera_pos
*/
float calc_cusp_func(vec3 pdir0, vec3 pdir1, float curv0, float curv1, vec3 vpos, vec3 vnor, vec3 cam_pos_ws)
{
    vec3 v = cam_pos_ws - vpos; 
    vec3 v_ = normalize(v); 
    float ndv = dot(v_, vnor);

    vec2 cusp_func = vec2(dot(v, normalize(pdir0)), dot(v, normalize(pdir1))); 
    cusp_func *= cusp_func;
    cusp_func.x = dot(cusp_func, vec2(curv0, curv1)); 

    return cusp_func.x; 
}


/*
 * Compute the dihedral angle between two triangles sharing an edge v1-v3.
 *     v1 ------- v3                       
 *    /->\ <----_/                                 
 *   / __/\   _/                                 
 *  /_/    \_/  we're using quad edge layout in npr_strokegen_topo_lib.glsl.                              
 * v0      v2     
 *                                        
 * Calculated Dihedral Angle                         
 *  0                                v0    v2        v0v2 
 *  v1        v1    v0 -- v1 -- v2    \    /          ||                   
 *  ||       /..\        '..'         .\  /.        .'||'.                 
 *  ||      /pi/2\        pi          . v1 .        . v1 . 2pi             
 * v0v2    v0    v2                    '''' 3pi/2    ''''      
*/
float calc_dihedral_angle(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{ 
	vec3 v10 = v0 - v1;
	vec3 v13 = v3 - v1;
   	vec3 v12 = v2 - v1;

	vec3 n0 = normalize(cross(v13, v10));
	vec3 n2 = normalize(cross(v12, v13));
	
    /* https://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleMeshDerivativesCheatSheet.pdf */
    float dihedral = acos(dot(n0, n2)); /* [0, pi] */
    float convex = .0f < dot(v13, cross(n0, n2)) ? -1.0 : 1.0f; 
    dihedral *= convex; /* [-pi, pi] */
    dihedral += M_PI; /* [0, 2pi] */

    return dihedral; 
}







/* Simple 3D Math ----------------------- */
// Project D onto the plane of the triangle q012(q0-q1-q2 has CCW winding), 
// parameterized with D_proj = u*Q1 + v*Q2, where Q1 = q1 - q0, Q2 = q2 - q0
vec2 proj_vec_to_triangle_plane(vec3 D, vec3 Q1, vec3 Q2)
{
    vec3 N = cross(Q1, Q2); 
    N = normalize(N); 
    D = D - dot(D, N)*N; 
    D = normalize(D); 
    // Solve for underdetermined system 
    // D = [Q1 Q2] * [u v]^T
    float Q1dotQ2 = dot(Q1, Q2);
    vec2 MT_mul_D = vec2(
        dot(Q1, D), dot(Q2, D)
    );
    mat2 MT_mul_M = mat2(
        dot(Q1, Q1), Q1dotQ2,
        Q1dotQ2, dot(Q2, Q2)
    );
    // [uD vD] = (M^T * M)^-1 * M^T * D
    vec2 uvD = inverse(MT_mul_M) * MT_mul_D;

    return uvD; 
}





#endif


