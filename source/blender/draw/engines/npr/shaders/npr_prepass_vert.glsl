
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  id = uint(resource_id);
  // normal = normal_object_to_world(nor);
  // /* TODO: Load pre-computed tangents */
  // tangent.xyz = normal_object_to_world(normalize(cross(vec3(0, 0, 1), nor)));

  vec3 wP = point_object_to_world(pos);
  gl_Position = point_world_to_ndc(wP);
}
