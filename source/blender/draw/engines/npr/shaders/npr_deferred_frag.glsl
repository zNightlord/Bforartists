
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  ivec2 uv = ivec2(gl_FragCoord.xy);
  uint id = texelFetch(id_tx, uv, 0).r;

  vec4 contour_data = texelFetch(contour_tx, uv, 0).rgba; 
  out_color.rgba = contour_data.rgba;
  out_line_color.rgba = contour_data.rgba;
  out_line_data = 0;

/*   out_color.rgb = id != 0 ? vec3(1) : vec3(0.5);
  out_line_color = vec4(0);
  out_line_data = 0;

  ivec2 offsets[4] = ivec2[4](ivec2(-1, 0), ivec2(1, 0), ivec2(0, 1), ivec2(0, -1));
  for (int i = 0; i < 4; i++) {
    uint offset_id = texelFetch(id_tx, uv + offsets[i], 0).r;
    if (id > offset_id) {
      
      out_line_color.a = 1;
      out_line_data = 1;
      out_color.rgb = vec3(0);
      return;
    }
  } */
}
