
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
    // goes CCW on screen. 
    vec2 edgenor_uv = vec2(edgedir_uv.y, -edgedir_uv.x); 


    // Debug stroke topo
    uint contour_edge_rank = ssbo_contour_edge_rank_[contour_edge_id]; 
    uint contour_edge_list_len = ssbo_contour_edge_list_len_[contour_edge_id]; 
    float contour_edge_param = float(contour_edge_rank) / float(contour_edge_list_len); 
    uint contour_edge_list_head = ssbo_contour_edge_list_head_[contour_edge_id]; 
    uvec2 link = uvec2(
        ssbo_contour_to_contour_[2*contour_edge_id], 
        ssbo_contour_to_contour_[2*contour_edge_id+1u]
    ); 
    uvec2 link_prev = uvec2(
        ssbo_contour_to_contour_[2*link.x], 
        ssbo_contour_to_contour_[2*link.x+1u]
    ); 
    uvec2 link_next = uvec2(
        ssbo_contour_to_contour_[2*link.y], 
        ssbo_contour_to_contour_[2*link.y+1u]
    ); 
    bool invalid_link = (contour_edge_id != link_prev.y) && (contour_edge_id != link_next.x);

    
    uint end_node_id = contour_edge_list_head; 
    uvec2 link_end = uvec2(
        ssbo_contour_to_contour_[2*end_node_id], 
        ssbo_contour_to_contour_[2*end_node_id+1u]
    ); 
	uint end_node_rank = ssbo_contour_edge_rank_[end_node_id]; 


	uint jump_node_id = contour_edge_id; 
	bool found_end_node = false; 
	uint jumps_find_end_node = 0; 
	if (contour_edge_list_len < contour_edge_rank)
	{
		for (uint jump = 1; jump < 10 * min(65536, contour_edge_list_len - 1); ++jump)
		{
			jumps_find_end_node++; 
            if (jump_node_id == end_node_id)
            {
                found_end_node = true; 
                break;
            }
			
            uint jump_node_id_next = ssbo_contour_to_contour_[2*jump_node_id+1u]; 
			if (jump_node_id == jump_node_id_next) break; 
			
			jump_node_id = jump_node_id_next; 
		}
        if (!found_end_node)
        {
            jump_node_id = contour_edge_id;  
            for (uint jump = 1; jump < 10 * min(65536, contour_edge_list_len - 1); ++jump)
            {
				jumps_find_end_node++; 
                if (jump_node_id == end_node_id)
                {
                    found_end_node = true; 
                    break;
                }
                
                uint jump_node_id_prev = ssbo_contour_to_contour_[2*jump_node_id]; 
				if (jump_node_id == jump_node_id_prev) break; 
				
				jump_node_id = jump_node_id_prev; 
            }
        }
	}

    bool found_curr_edge = false; 
    jump_node_id = end_node_id; 
	uint jumps_find_curr_edge = 0; 
    for (uint jump = 1; jump < min(65536, contour_edge_list_len - 2); ++jump)
    {
		jumps_find_curr_edge++; 
        jump_node_id = ssbo_contour_to_contour_[2*jump_node_id + 1]; 
        if (jump_node_id == contour_edge_id)
        {
            found_curr_edge = true; 
            break; 
        }
    }

	uint code = (uint(found_curr_edge) << 1) | uint(found_end_node); 

    color = 
        /* min(1.0f, contour_edge_param * 1.5f) *  */
        // vec4(edgedir_uv.xy * 0.5f + .5f, 0, 1); 

        // vec4(
        //     /* contour_edge_param *  */.5f * rand_col_rgb(contour_edge_list_len, contour_edge_list_len).rgb, 
        //     1.0f/* contour_edge_rank *//* contour_edge_list_head */
        // ); 

        contour_edge_list_len < contour_edge_rank ? 
			vec4(
					code, 
					jumps_find_end_node, 
					jumps_find_curr_edge, 
					contour_edge_list_len
			)
            // vec4(contour_edge_rank, contour_edge_list_len, contour_edge_list_head, end_node_rank + 1) 
            : vec4(.0f); // vec4(contour_edge_rank, contour_edge_list_len, contour_edge_list_head, 0);
        
		// vec4(float(contour_edge_id), vec2(link.xy), 1.0f);  
        /* vec4((float(contour_edge_rank) / float(contour_edge_list_len)).xxx, 1.0f); */
        /* vec4(edgedir_uv.xy * 0.5f + .5f, 0, 1); */ 


    normal = vec3(0, 0, 1);
    tangent.xyz = vec3(0, 0, 1);

    gl_Position = pos_hclip;
    /* Apply depth bias to curve breaks, 
     * which can happen due to Z-fighting artifacts. */
    /* TODO: Further imporve this */
    gl_Position.z -= 5.0e-7 * whclip;
    gl_Position.xy += edgenor_uv * whclip * pcs_screen_size_inv_ * 1.0f; 
	vec2 edge_dir_ext = (vid % 2u == 1u) ? edgedir_uv : -edgedir_uv; 
	gl_Position.xy += edge_dir_ext * whclip * pcs_screen_size_inv_ * 1.5f; 
}
