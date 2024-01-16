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





#endif


