#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)


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

    vec3 ld_vcurv_tensor(uint vtx_id)
    {
        uvec3 vcurv_tensor_enc;
        Load3(ssbo_vcurv_tensor_, vtx_id, vcurv_tensor_enc); 
        return uintBitsToFloat(vcurv_tensor_enc); 
    }
    void st_vcurv_tensor(uint vtx_id, vec3 vcurv_tensor)
    {
        uvec3 vcurv_tensor_enc = floatBitsToUint(vcurv_tensor);
        Store3(ssbo_vcurv_tensor_, vtx_id, vcurv_tensor_enc); 
    }

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



#endif


