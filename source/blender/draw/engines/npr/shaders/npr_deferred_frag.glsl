
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  ivec2 uv = ivec2(gl_FragCoord.xy);
  // uint id = texelFetch(id_tx, uv, 0).r;

  vec4 contour_data = texelFetch(contour_tx, uv, 0).rgba; 
  
  out_color.rgba = contour_data.rgba;
  if (0 < pcs_full_overlay_) out_color.a = 1; 

  out_line_color.rgba = contour_data.rgba;
  out_line_data = 0;
}
