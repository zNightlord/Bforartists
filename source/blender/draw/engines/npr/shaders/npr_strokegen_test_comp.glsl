
/**
* Testing compute shader for bnpr engine.
 */

#pragma BLENDER_REQUIRE(gpu_shader_codegen_lib.glsl)

void main()
{
  buf_test[0] = buf_ibo[0];
}
