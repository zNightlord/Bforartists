#pragma BLENDER_REQUIRE(npr_strokegen_fill_indirect_args_inputs_lib.glsl)

void main()
{
    if (gl_GlobalInvocationID.x == 0)
    {
#if !defined(_KERNEL_MULTICOMPILE__FILL_DRAW_ARGS__REMESHING)
        uvec3 dispatch_args = uvec3(1, 1, 1); 
        GetDispatchArgs(/*out*/ dispatch_args); 
        FillDispatchArgsBuffer(dispatch_args); 
#else
        DrawCommand drw_args; 
        GetDrawArgs(/*out*/ drw_args); 
        FillDrawArgsBuffer(drw_args); 
#endif
    } 
}