
#pragma BLENDER_REQUIRE(npr_strokegen_encode_lib.glsl)

#ifndef BNPR_BRUSH_TOOLBOX_INCLUIDED
#define BNPR_BRUSH_TOOLBOX_INCLUIDED

// Skeletal Stroke -------------------------------
// Wing-Quad Anatomy
// V2 --------- V5
// \7         8/10\
//  \	    _/     \	... ...
//   \	  _/		\
//    \6 /9	       11\
// ==> V1 ========= V4 ==> Stroke Dir ==>
//     |1		  2/4|  
//     |	    _/	 |
//     |	  _/	 |
//     |    _/		 |	... ...
//     |0 _/		 |
//     |_/3		    5|
//	  V0 ----------- V3
// each Pixel-Edge forms a Spine: (V1->V4)
// and spans a "wing" (V1-{V0, V2})
// each Wing has 3 vertices
// point list(x12) is generated from vertices shared by multiple quad triangles
#define POINTS_PER_WING_QUAD 12
uniform uint wing_quad_vertex_buffer[POINTS_PER_WING_QUAD] =
{
	0, 1, 4, // v0, v1, v4
	0, 4, 3,
	1, 2, 5,
	1, 5, 4
};

struct SpineTopoInfo
{
	bool looped; // if last spine is connected to the first spine 
	bool is_tail_spine; 
	uint head_sample_id; 
}; 

uint get_wing_quad_spine_index(uint drw_vtx_id)
{
	return drw_vtx_id / POINTS_PER_WING_QUAD;
}

uint get_wing_quad_vertex_index(uint drw_vtx_id, SpineTopoInfo sti)
{
	uint spine_id = drw_vtx_id / POINTS_PER_WING_QUAD;
	uint point_id = drw_vtx_id % POINTS_PER_WING_QUAD;

	uint vert_id = wing_quad_vertex_buffer[point_id];

	if (sti.is_tail_spine)
	{
		if (sti.looped && 3u <= vert_id)
		{ // jump back to header spine
			spine_id = sti.head_sample_id;
			vert_id -= 3u; 
		} 
		else if (3u <= vert_id)
		{ // collapse the last piece so that it draws nothing
			vert_id -= 3u; 
		}
	}
	
	return 3u * spine_id + vert_id;
}

vec2 get_wing_quad_vertex_uv(uint vert_id)
{
//V2- --------- V5
//(0,1)        /(1,1)
//  \	    _/     \	... ...
//   \	  _/		\
//    \  /9	         \
// ==> V1 ========= V4 ==> Stroke Dir ==>
//    (0,.5)       (1,.5)
//     |	    _/	 |
//     |	  _/	 |
//     |    _/		 |	... ...
//     |  _/		 |
//    (0,0)		   (1,0)
//	  V0 ----------- V3
	float offset_u = float(vert_id / 3u);
	float offset_v = .5f * float(vert_id % 3u);
	
	return vec2(offset_u, offset_v);
}

mat3x2 compute_wing_quad_verts(
	vec2 center_coord,
	vec2 normal, float width
)
{
	mat3x2 verts;
	verts[0] = center_coord - width * normal;
	verts[1] = center_coord;
	verts[2] = center_coord + width * normal;

	return verts;
}




#if defined(USE_CONTOUR_2D_STROKE_MESH_BUFFER)
void store_ssbo_stroke_mesh_pool__skeletal_VB(uint curve_sample_id, mat3x2 verts)
{
    ssbo_stroke_mesh_pool_[curve_sample_id * 3u + 0u] = pack_2d_uv(verts[0]); 
    ssbo_stroke_mesh_pool_[curve_sample_id * 3u + 1u] = pack_2d_uv(verts[1]); 
    ssbo_stroke_mesh_pool_[curve_sample_id * 3u + 2u] = pack_2d_uv(verts[2]); 
}
vec2 load_ssbo_stroke_mesh_pool__skeletal_VB(uint drw_vtx_id, SpineTopoInfo sti)
{
    uint wq_vtx_id = get_wing_quad_vertex_index(drw_vtx_id, sti); 
    return unpack_2d_uv(ssbo_stroke_mesh_pool_[wq_vtx_id]); 
}
void store_ssbo_stroke_mesh_pool__skeletal_color(uint curve_sample_id, vec4 col, uint num_samples)
{
    uint offset = num_samples * 3u;
    ssbo_stroke_mesh_pool_[offset + curve_sample_id*2u + 0u] = pack_2d_uv(col.rg); 
    ssbo_stroke_mesh_pool_[offset + curve_sample_id*2u + 1u] = pack_2d_uv(col.ba); 
}
vec4 load_ssbo_stroke_mesh_pool__skeletal_color(uint drw_vtx_id, uint num_samples)
{
    uint offset = num_samples * 3u;
    uint curve_sample_id = get_wing_quad_spine_index(drw_vtx_id); 

    vec4 col; 
    col.rg = unpack_2d_uv(ssbo_stroke_mesh_pool_[offset + curve_sample_id*2u + 0u]); 
    col.ba = unpack_2d_uv(ssbo_stroke_mesh_pool_[offset + curve_sample_id*2u + 1u]); 
    return col; 
}
void store_ssbo_stroke_mesh_pool__skeletal_topo_info(uint curve_sample_id, SpineTopoInfo sti, uint num_samples)
{
	uint sti_enc = sti.head_sample_id;
	sti_enc <<= 1u; 
	sti_enc |= uint(sti.is_tail_spine);
	sti_enc <<= 1u;
	sti_enc |= uint(sti.looped);

    uint offset = num_samples * 5u;
    ssbo_stroke_mesh_pool_[offset + curve_sample_id] = sti_enc; 
}
SpineTopoInfo load_ssbo_stroke_mesh_pool__skeletal_topo_info(uint drw_vtx_id, uint num_samples)
{
	uint offset = num_samples * 5u;
	uint curve_sample_id = get_wing_quad_spine_index(drw_vtx_id); 
	uint sti_enc = ssbo_stroke_mesh_pool_[offset + curve_sample_id]; 

	SpineTopoInfo sti; 
	sti.looped = bool(sti_enc & 1u);
	sti_enc >>= 1u;
	sti.is_tail_spine = bool(sti_enc & 1u);
	sti_enc >>= 1u;
	sti.head_sample_id = sti_enc;

	return sti; 
}

#endif


#endif


