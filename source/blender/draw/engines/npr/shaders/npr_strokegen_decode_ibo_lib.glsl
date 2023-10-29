#ifndef BNPR_DECODE_IBO__INCLUDED
#define BNPR_DECODE_IBO__INCLUDED


#if !defined(DECODE_IBO_EXCLUDE)


/** Input
 * #define buf_ibo XXX   <= index buffer, either 16 bits or 32 bits per index
*/
void load_and_decode_ibo(uvec4 vid, uint prim_id, uint vertsPerPrim)
{
#if defined(_KERNEL_MULTICOMPILE__INDEX_BUFFER_16BIT)
    uint ibo_items_16 = vertsPerPrim / 2; /* 16 bits per vtx id */
    for (uint i = 0; i < 2; ++i) 
    {
        uint ibo_addr = ibo_items_16 * prim_id + i;
        uint ibo_data = buf_ibo[ibo_addr];
        /* decode 16 bit index */
        uint ibo_data_16h = (ibo_data >> 16u);
        uint ibo_data_16l = (ibo_data & 0xFFFFu);
        /* fetch vertex pos */
        vid[i * 2] = (ibo_data_16l);
        vid[i * 2 + 1] = (ibo_data_16h);
    }
#else
    uint ibo_items_32 = vertsPerPrim; /* 32 bits per vtx id */
    for (uint i = 0; i < 4; ++i) 
    {
        uint ibo_addr = ibo_items_32 * prim_id + i;
        uint ibo_data = buf_ibo[ibo_addr];
        /* fetch vertex pos */
        vid[i] = (ibo_data);
    }
#endif
}

#endif


#endif


