#pragma BLENDER_REQUIRE(npr_strokegen_fill_indirect_args_inputs_lib.glsl)


#ifdef _KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS

void main()
{
    if (gl_GlobalInvocationID.x == 0)
    {
        uvec3 args = GetDispatchArgs(); 
        
        SSBO_ARGS.num_groups_x = args.x; 
        SSBO_ARGS.num_groups_y = args.y;
        SSBO_ARGS.num_groups_z = args.z;    
    }
}

#endif /* _KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS */