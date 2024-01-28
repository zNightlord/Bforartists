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









// Scalable Curvature Computations on Triangle Meshes
#ifndef M_PI
#define M_PI           (3.14159265358979323846)
#endif
#ifndef M_PI_2
#define M_PI_2         (1.57079632679489661923)
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
    // Adding 1000 area n x n to anisotropic measure
    // force the principal direction eigenvectors to be tangential to the surface
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

















#endif


