#ifndef BNPR_SEGLOOPCONV1D_LIB_INPUTS_INCLUDED
#define BNPR_SEGLOOPCONV1D_LIB_INPUTS_INCLUDED


#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE)
void func_device_store_loopconv1d_patch_id(uint stAddr, uint val)
{
    ssbo_segloopconv1d_patch_table_[stAddr] = val;
}
#endif


#if defined(_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION)
uint func_device_load_loopconv1d_patch_id(uint ldAddr)
{
    return ssbo_segloopconv1d_patch_table_[ldAddr];
}

uint func_device_load_loopconv1d_data(uint elemId)
{
    return wang_hash(elemId * 73u) % 711u;
}
#endif


#endif


