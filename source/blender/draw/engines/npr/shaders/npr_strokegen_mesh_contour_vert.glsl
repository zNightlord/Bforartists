
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

struct ListValidationContext
{
    uint list_len; 
    uint list_head;

    uint curr_node_id;
    uint jump_node_id; 
    uint num_jumps; 
    uint search_node_id;  

    bool found_search_node; 
    bool consistent_list_len;
    bool consistent_list_head; 
    bool reach_tail; 
    bool reach_head; 

    uint curr_node_rank; 
    uint prev_node_rank; 
    uint num_rank_fails; 
}; 


void main()
{
    uint vid = gl_VertexID;
    uint contour_edge_id = vid / 2u; 
    
    uint base_addr = vid * 6; 
    uint num_contour_edges = ssbo_bnpr_mesh_pool_counters_.num_contour_edges; 




    // Debug stroke topo
    uint contour_edge_rank = ssbo_contour_edge_rank_[contour_edge_id]; 
    uint contour_edge_list_len = ssbo_contour_edge_list_len_[contour_edge_id]; 
    float contour_edge_param = float(contour_edge_rank) / float(contour_edge_list_len); 
    uint contour_edge_list_head = ssbo_contour_edge_list_head_[contour_edge_id]; 

    
    uint end_node_id = contour_edge_list_head; 
    uvec2 link_end = uvec2(
        ssbo_contour_to_contour_[2*end_node_id], 
        ssbo_contour_to_contour_[2*end_node_id+1u]
    ); 
	uint end_node_rank = ssbo_contour_edge_rank_[end_node_id]; 


	
    ListValidationContext ctx; 
    ctx.list_len = contour_edge_list_len; 
    ctx.list_head = contour_edge_list_head;
    ctx.curr_node_id = contour_edge_id;
    ctx.jump_node_id = contour_edge_id; 
    ctx.search_node_id = end_node_id;  
    ctx.found_search_node = false; 
    ctx.consistent_list_len = true;
    ctx.consistent_list_head = true; 
    ctx.reach_tail = false; 
    ctx.reach_head = false; 
    ctx.num_jumps = 0; 
    ctx.prev_node_rank = int(contour_edge_rank); 
    ctx.num_rank_fails = 0; 

    ListValidationContext ctx_another; 
    ctx_another.list_len = ssbo_contour_edge_list_len_[end_node_id]; 
    ctx_another.list_head = ssbo_contour_edge_list_head_[end_node_id];
    ctx_another.curr_node_id = end_node_id;
    ctx_another.jump_node_id = end_node_id; 
    ctx_another.search_node_id = contour_edge_id;  
    ctx_another.found_search_node = false; 
    ctx_another.consistent_list_len = true;
    ctx_another.consistent_list_head = true; 
    ctx_another.reach_tail = false; 
    ctx_another.reach_head = false; 
    ctx_another.num_jumps = 0; 
    ctx_another.prev_node_rank = int(ssbo_contour_edge_rank_[end_node_id]); 
    ctx_another.num_rank_fails = 0; 

    bool bad_contour = contour_edge_list_len < contour_edge_rank; 
    bool bad_list = false; 

    uint dbg_draw_node_id = end_node_id; 
	bool found_dbg_draw_node = false; 
	uint dbg_draw_node_rank = 0; 
	uint max_rank = (ssbo_contour_edge_rank_[end_node_id]); 
	uint first_fail_node_id = 0; 
	if (bad_contour)
	{
        uint head_node_id = contour_edge_id; 
        ctx.jump_node_id = contour_edge_id; 
        ctx.prev_node_rank = contour_edge_rank;
		for (uint jump = 1; jump < 100 * max(800, ctx.list_len); ++jump)
		{
            uint curr_node_id = ssbo_contour_to_contour_[2*ctx.jump_node_id];
            if (curr_node_id == ctx.jump_node_id)
            {
                head_node_id = curr_node_id; 
                break; 
            }

            uint curr_node_rank = ssbo_contour_edge_rank_[curr_node_id]; 
            if (curr_node_rank < ctx.prev_node_rank)
            {
                head_node_id = ctx.jump_node_id;
                break; 
            }

            ctx.prev_node_rank = curr_node_rank; 
            ctx.jump_node_id = curr_node_id; 
		}

        bool loop_list = true; 
        uint list_len = 0; 
        ctx.curr_node_id = head_node_id; // header node
        ctx.jump_node_id = head_node_id;
        ctx.curr_node_rank = ssbo_contour_edge_rank_[head_node_id]; 
        ctx.consistent_list_head = true; 
        for (uint jump = 1; jump < 100 * max(800, ctx.list_len); ++jump)
		{
            list_len++; 

            uint next_node = ssbo_contour_to_contour_[2*ctx.jump_node_id + 1];

            if (next_node == ctx.jump_node_id)
            {
                loop_list = false; 
                break; 
            }
            if (next_node == head_node_id)
                break; 

			uint next_node_rank = ssbo_contour_edge_rank_[next_node];
            if (ctx.curr_node_rank != next_node_rank + jump)
            {
                ctx.num_rank_fails++; 
            }

            uint next_node_list_len = ssbo_contour_edge_list_len_[next_node]; 
            if (ctx.list_len != next_node_list_len) ctx.consistent_list_len = false; 

            uint next_node_list_head = ssbo_contour_edge_list_head_[next_node];
            if (ctx.list_head != next_node_list_head) ctx.consistent_list_head = false;

            uint next_node_prev = ssbo_contour_to_contour_[2*next_node];
            if (next_node_prev != ctx.jump_node_id)
            {
                bad_list = true;  
            }

            ctx.jump_node_id = next_node;
		}
			
		ctx_another.jump_node_id = end_node_id; 
		ctx_another.num_rank_fails = 0; 
		ctx_another.reach_tail = false;
		ctx_another.curr_node_rank = int(ssbo_contour_edge_rank_[end_node_id]); 
        for (uint jump = 1; jump < 100 * max(800, ctx_another.list_len); ++jump)
        {
            ctx_another.num_jumps++; 
            uint next_node = ssbo_contour_to_contour_[2*ctx_another.jump_node_id + 1];
			
			if (next_node == end_node_id || next_node == ctx_another.jump_node_id) 
			{
				ctx_another.reach_tail = true; 
				break; 
			}
            if (ctx_another.jump_node_id == ctx_another.search_node_id)
                ctx_another.found_search_node = true; 

            uint list_len_ = ssbo_contour_edge_list_len_[next_node]; 
            if (ctx_another.list_len != list_len_) ctx_another.consistent_list_len = false; 
			
			uint next_node_rank = ssbo_contour_edge_rank_[next_node]; 
            if (ctx_another.curr_node_rank != next_node_rank + jump)
            {
                ctx_another.num_rank_fails++; 
            }
            
            if (((contour_edge_rank % 1733) == (jump % 1733))
                && (!found_dbg_draw_node) 
                && 0 < ctx_another.num_rank_fails
            ){
                dbg_draw_node_id = next_node; 
                found_dbg_draw_node = true; 
				dbg_draw_node_rank = next_node_rank; 
            }



			max_rank = min(max_rank, next_node_rank); 

			ctx_another.jump_node_id = next_node; 
        }
	}
    
    
    uint draw_contour_edge_id = contour_edge_id;
	bool draw_dbg_edge = false; // found_dbg_draw_node && bad_contour; 
    if (draw_dbg_edge)
        draw_contour_edge_id = dbg_draw_node_id; 

    uint addr_uv = mesh_pool_addr__edgeuv(draw_contour_edge_id); 
    uint wedge_id = buf_strokegen_mesh_pool[addr_uv+4u];

    if (vid % 2u == 1u) 
        addr_uv += 2u;
    vec2 uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv+1])
    );

    uint addr_uv_another = mesh_pool_addr__edgeuv(draw_contour_edge_id);
    if (vid % 2u == 0u) 
        addr_uv_another += 2u;
    vec2 uv_another_vtx = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv_another+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_uv_another+1])
    );
    float edge_len = length((uv - uv_another_vtx) / pcs_screen_size_inv_.xy); 



    uint addr_zhclip = mesh_pool_addr__zwhclip(draw_contour_edge_id); 
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
    uint addr_edgedir = mesh_pool_addr__edgedir(draw_contour_edge_id); 
    vec2 edgedir_uv = vec2(
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+0]), 
        uintBitsToFloat(buf_strokegen_mesh_pool[addr_edgedir+1])
    );
    // goes CCW on screen. 
    vec2 edgenor_uv = vec2(edgedir_uv.y, -edgedir_uv.x); 


    color = 
        vec4(rand_col_rgb(contour_edge_list_len, contour_edge_list_len), wedge_id);

        // edge_len < 100.0f ? vec4(rand_col_rgb(contour_edge_list_len, contour_edge_list_len), 1.0f) : vec4(.0f);

        // contour_edge_list_len < contour_edge_rank ? 
        //     vec4(1.0f)
		// 	// vec4( 
		// 	// 	first_fail_node_id, 
		// 	// 	max_rank, 
		// 	// 	ctx_another.num_rank_fails, 
		// 	// 	ctx_another.num_jumps
		// 	// )
        //     // vec4(contour_edge_rank, contour_edge_list_len, contour_edge_list_head, end_node_rank + 1) 
        //     : vec4(.0f); 
        //     // : vec4(contour_edge_list_len.xxx, .0f); 
        //     // vec4(contour_edge_rank, contour_edge_list_len, contour_edge_list_head, 0);


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
