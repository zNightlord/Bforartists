
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

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


void main()
{
    uint vid = gl_VertexID;
    uint contour_edge_id = vid / 2u; 
    
    uint base_addr = vid * 6; 
    uint num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 
 
    uint addr_uv = mesh_pool_addr__edgeuv(contour_edge_id); 
    if (vid % 2u == 1u) 
        addr_uv += 2u;
    vec2 uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+1])
    );

    uint addr_zhclip = mesh_pool_addr__zwhclip(contour_edge_id); 
    if (vid % 2u == 1u) 
        addr_zhclip += 2u; 
    float zhclip = uintBitsToFloat(buf_strokegen_mesh_pool[addr_zhclip+0]); 
    float whclip = uintBitsToFloat(buf_strokegen_mesh_pool[addr_zhclip+1]); 


    vec3 pos_view = get_view_space_from_depth(uv, zhclip);

    vec4 pos_hclip; 
    pos_hclip.xy = uv * 2.0f - 1.0f;
    pos_hclip.z = zhclip;
    pos_hclip.w = whclip; 
    pos_hclip.xy *= pos_hclip.w; 

    /* Fetch edge dir from ssbo */
    uint addr_edgedir = mesh_pool_addr__edgedir(contour_edge_id); 
    vec2 edgedir_uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+1])
    );
    vec2 edgenor_uv = vec2(edgedir_uv.y, -edgedir_uv.x); 

    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    gl_Position = pos_hclip;
    /* Apply depth bias to curve breaks, 
     * which can happen due to Z-fighting artifacts. */
    /* TODO: Further imporve this */
    gl_Position.z -= 2.0e-5 * whclip;
    gl_Position.xy += edgenor_uv * whclip * pcs_screen_size_inv_ * 1.0f; 
	vec2 edge_dir_ext = (vid % 2u == 1u) ? edgedir_uv : -edgedir_uv; 
	gl_Position.xy += edge_dir_ext * whclip * pcs_screen_size_inv_ * 1.5f; 


    uint contour_edge_rank = ssbo_contour_edge_rank_[contour_edge_id]; 
    uint contour_edge_list_len = ssbo_contour_edge_list_len_[contour_edge_id]; 
    float contour_edge_param = float(contour_edge_rank) / float(contour_edge_list_len); 
    uint contour_edge_list_head = ssbo_contour_edge_list_head_[contour_edge_id]; 

    uvec2 dbg_contour_prev_next = uvec2(
        ssbo_contour_to_contour_[2*contour_edge_id], 
        ssbo_contour_to_contour_[2*contour_edge_id+1u]
    ); 

    color = /* min(1.0f, contour_edge_param * 1.5f) *  */
        /* vec4(float(contour_edge_id), vec2(dbg_contour_prev_next.xy), 1.0f);   */
        vec4(.5f * rand_col_rgb(contour_edge_list_head, contour_edge_list_len).rgb, 1.0f); 
        /* vec4((float(contour_edge_rank) / float(contour_edge_list_len)).xxx, 1.0f); */
        /* vec4(edgedir_uv.xy * 0.5f + .5f, 0, 1); */ 
}
