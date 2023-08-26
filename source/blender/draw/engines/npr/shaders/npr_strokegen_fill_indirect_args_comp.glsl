#pragma BLENDER_REQUIRE(npr_strokegen_fill_indirect_args_inputs_lib.glsl)


#ifdef _KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS

void main()
{
    if (gl_GlobalInvocationID.x == 0)
    {
        uvec3 dispatch_args = uvec3(1, 1, 1); 
        GetDispatchArgs(/*out*/ dispatch_args); 
        FillDispatchArgsBuffer(dispatch_args); 
    }
}

#endif /* _KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS */