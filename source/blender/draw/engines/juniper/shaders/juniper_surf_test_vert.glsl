
#pragma BLENDER_REQUIRE(juniper_nodetree_lib.glsl)
#pragma BLENDER_REQUIRE(juniper_attributes_lib.glsl)

Closure out_surface, out_volume, out_displacement, out_thickness;

void init_interface()
{
    #ifdef GPU_VERTEX_SHADER
    interp.P = vec3(0.0);
    interp.N = vec3(0.0);
    // interp.barycentric_coords = vec2(0.0); TODO barycentric support
    #endif
}

void main()
{
    PASS_RESOURCE_ID;
    init_interface();
    interp.N = normal_object_to_world(nor).xyz;
    interp.P = point_object_to_world(pos);

    attrib_load();

#ifdef MAT_DISPLACEMENT_BUMP
    interp.P += nodetree_displacement();

    gl_Position = point_world_to_ndc(interp.P);
#else
    gl_Position = point_object_to_ndc(pos);
#endif
}
