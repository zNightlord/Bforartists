
#pragma BLENDER_REQUIRE(common_view_clipping_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

/*
ubo_view_matrices_
ssbo_face_to_vert_draw_depth_
*/

void main()
{
    mat4 world_to_view = ubo_view_matrices_.viewmat;
	mat4 mat_camera_proj = ubo_view_matrices_.winmat; 

    uint face_id = gl_VertexID / 3u; 
    uint vtx_order = gl_VertexID % 3u; 
    uint data_start = face_id * 4u; 
    uint vert_id = ssbo_face_to_vert_draw_depth_[data_start + vtx_order];
    uint wedge_id = ssbo_face_to_vert_draw_depth_[data_start + 3u]; 

    vec3 vpos_ws; 
    vpos_ws.x = (ssbo_vbo_full_[vert_id * 3u + 0u]); 
    vpos_ws.y = (ssbo_vbo_full_[vert_id * 3u + 1u]);
    vpos_ws.z = (ssbo_vbo_full_[vert_id * 3u + 2u]);

    vec3 vpos_vs = (world_to_view * vec4(vpos_ws, 1.0)).xyz;
    vec4 vpos_cs = mat_camera_proj * vec4(vpos_vs, 1.0);

    gl_Position = vpos_cs;
    color = vec4(float(wedge_id));
}
