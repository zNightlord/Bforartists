
#pragma BLENDER_REQUIRE(npr_strokegen_compaction_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_topo_lib.glsl)
#pragma BLENDER_REQUIRE(npr_strokegen_load_store_lib.glsl)



/* Inputs: 
int pcs_collapse_iter_

!!! Remember to zero out these counters !!!
* SSBO_StrokeGenEdgeCollapseCounters    ssbo_edge_collapse_counters_[]
* SSBO_StrokeGenDynamicMeshCounters     ssbo_dyn_mesh_counters_in_

#ifdef _KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT
.define("GLOBAL_COUNTER", "ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapse_edges_pass_1")
.define("CP_TAG", "collapse_inital_selection")
#endif

#ifdef _KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT_PASS_1
.define("GLOBAL_COUNTER", "ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_split_edges")
.define("CP_TAG", "collapse_resolve_conflict")
#endif

 float ssbo_vbo_full_[]
 float pcs_remesh_edge_len_

 uint ssbo_per_collapse_edge_info_[]
 uint ssbo_per_edge_collapse_info_in_[]
 uint ssbo_per_edge_collapse_info_out_[]
 uint ssbo_per_vert_collapse_wedge_id_[] // !!!clean this buffer to NULL_EDGEs
*/

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

float get_collapse_edge_len_max()
{
    return calc_remesh_edge_len_min(pcs_remesh_edge_len_); 
}

struct EdgeCollapseInfo
{
    bool is_collapse_ok; 
    float edge_len; 
    uint id; 
}; 
uvec2 encode_edge_collapse_info(EdgeCollapseInfo eci)
{
    uvec2 eci_enc = uvec2(0, 0); 
    eci_enc.x = floatBitsToUint(eci.edge_len); 
    eci_enc.y = ((eci.id << 1u) | uint(eci.is_collapse_ok));

    return eci_enc;  
}
EdgeCollapseInfo decode_edge_collapse_info(uvec2 eci_enc)
{
    EdgeCollapseInfo eci; 
    eci.edge_len = uintBitsToFloat(eci_enc.x);
    eci.is_collapse_ok = (0u != (eci_enc.y & 1u));
    eci_enc.y >>= 1u; 
    eci.id = eci_enc.y;  

    return eci; 
}
void store_per_collapse_edge_info(uint collapse_edge_id, EdgeCollapseInfo pcei)
{
    uvec2 pcei_enc = encode_edge_collapse_info(pcei); 
    Store2(ssbo_per_collapse_edge_info_, collapse_edge_id, pcei_enc);
}
EdgeCollapseInfo load_per_collapse_edge_info(uint collapse_edge_id)
{
    uvec2 pcei_enc; 
    Load2(ssbo_per_collapse_edge_info_, collapse_edge_id, pcei_enc); 
    return decode_edge_collapse_info(pcei_enc); 
}
void store_per_edge_collapse_info(uint wedge_id, EdgeCollapseInfo peci)
{
    uvec2 peci_enc = encode_edge_collapse_info(peci); 
    Store2(ssbo_per_edge_collapse_info_out_, wedge_id, peci_enc); 
}
EdgeCollapseInfo load_per_edge_collapse_info(uint wedge_id)
{
    uvec2 peci_enc;
    Load2(ssbo_per_edge_collapse_info_in_, wedge_id, peci_enc); 
    return decode_edge_collapse_info(peci_enc); 
}

bool collapse_score_larger(EdgeCollapseInfo eci_0, uint wedge_id_0, EdgeCollapseInfo eci_1, uint wedge_id_1)
{
    if (!eci_0.is_collapse_ok) return false;
    if (!eci_1.is_collapse_ok) return true;
    if (eci_0.edge_len == eci_0.edge_len)
        return wedge_id_0 < wedge_id_1; 
    return eci_0.edge_len > eci_0.edge_len; 
}


void store_per_vert_collapse_wedge_id(uint vid, uint wedge_id)
{
    ssbo_per_vert_collapse_wedge_id_[vid] = wedge_id; 
}
uint load_per_vert_collapse_wedge_id(uint vid)
{
    return ssbo_per_vert_collapse_wedge_id_[vid]; 
}



#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_INIT)
void main()
{ /* Granularity: per vertex */
    uint vid = gl_GlobalInvocationID.x; 

    if (pcs_collapse_iter_ == 0 && gl_GlobalInvocationID.x == 0u)
    { 
        ssbo_edge_collapse_counters_[0].num_split_edges_pass_1 = 0u;
        ssbo_edge_collapse_counters_[0].num_split_edges = 0u; 
    }
    ssbo_per_vert_collapse_wedge_id_[vid] = NULL_EDGE;
}
#endif







#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT)
void main()
{
    uint wedge_id = gl_GlobalInvocationID.x; 

}
#endif






#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT_PASS_0)
void main()
{
    uint collapse_edge_id = gl_GlobalInvocationID.x; 

}
#endif





#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT_PASS_1)
void main()
{
    uint collapse_edge_id = gl_GlobalInvocationID.x; 

}
#endif






#if defined(_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE)
void main()
{
    uint collapse_edge_id = gl_GlobalInvocationID.x; 

}
#endif