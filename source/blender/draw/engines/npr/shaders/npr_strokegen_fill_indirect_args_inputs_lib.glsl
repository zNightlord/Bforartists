#ifndef NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H
#define NPR_STROKEGEN_FILL_INDIRECT_ARGS_INPUTS_LIB_H



#if defined(_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__LIST_RANKING_ANCHORS)

uvec3 GetDispatchArgs()
{
    return uvec3(
        (ssbo_list_ranking_atomics_.counter_anchors + INDIRECT_THREAD_GROUP_SIZE - 1u) / INDIRECT_THREAD_GROUP_SIZE, 
        1, 
        1
    ); 
}
#define GET_DISPATCH_ARGS() GetDispatchArgs()

#endif 



#endif

