
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)

#if !defined(SRC_OVERRIDE_FRAG) // User overrides

void main()
{
    out_color = texture(color_buffer, uvcoordsvar.xy);
    gl_FragDepth = texture(depth_buffer, uvcoordsvar.xy).x;
}

#endif
