
#define PIXELS_PER_THREAD 32u /* compress 1bit per pixel */
#define COMPRESS_RECT_SIZE CONTOUR_PIX_MARK_COMPRESS_RECT_SIZE /* (uvec2(8u, 4u)) */
void calc_coord_tex2d_contour_pix_marks(uvec2 screen_coord, out uvec2 block_coord, out uint bit_mask)
{
    block_coord = screen_coord.xy / COMPRESS_RECT_SIZE.xy; 
    uvec2 block_offset = screen_coord.xy % COMPRESS_RECT_SIZE.xy; 
    uint bit_id = block_offset.x + block_offset.y * COMPRESS_RECT_SIZE.x; 
    bit_mask = 1u << bit_id; 
}

#if defined(_KERNEL_MULTICOMPILE__COMPRESS_CONTOUR_IMG)
void main()
{
    const uvec2 groupId = gl_LocalInvocationID.xy;
	const uvec2 idx = gl_GlobalInvocationID.xy;
    
#if !defined(_KERNEL_MULTICOMPILE__COMPRESS_CONTOUR_IMG_VALIDATE)
    const uvec2 image_res = uvec2(pcs_screen_size_ + 1e-10f); 
    const uvec2 grid_beg = idx.xy * COMPRESS_RECT_SIZE; 

    /* Compressed block for this thread: 
     * Y3 | byte_3  
     * Y2 | byte_2  
     * Y1 | byte_1  
     * Y0 | byte_0
     *    +--------> X
     *   bit_0 ~ bit_7
     *      X0 ~ X7
    */
    uint compressed_32bits = 0u; 
    uint curr_bit_mask = 1u; 
    for (uint off_y = 0u; off_y < COMPRESS_RECT_SIZE.y; off_y++)
        for (uint off_x = 0u; off_x < COMPRESS_RECT_SIZE.x; off_x++)
        {
            const uvec2 grid_pos = grid_beg + uvec2(off_x, off_y);
            if (grid_pos.x >= image_res.x || grid_pos.y >= image_res.y)
                continue;
            const vec4 pixel = imageLoad(tex2d_contour_image_, ivec2(grid_pos)).rgba;
            bool b_contour_pixel = !all(lessThan(abs(pixel.xy), vec2(1e-10f, 1e-10f)));  

            if (b_contour_pixel)
                compressed_32bits |= curr_bit_mask; 
            curr_bit_mask <<= 1u; 
        }

    imageStore(tex2d_contour_pix_marks_, ivec2(idx.xy), compressed_32bits.xxxx);
#else /* debug kernel */
    const uvec2 screen_res = uvec2(pcs_screen_size_ + 1e-10f); 
    const uvec2 screen_coord = idx.xy;

    if (any(greaterThan(screen_coord, screen_res - uvec2(1u, 1u)))) return; 

    uvec2 block_coord; uint bit_mask; 
    calc_coord_tex2d_contour_pix_marks(screen_coord, /*out*/block_coord, bit_mask); 
    const uint compressed_32bits = imageLoad(tex2d_contour_pix_marks_, ivec2(block_coord.xy)).x;
    uint contour_pixel = uint((compressed_32bits & bit_mask) != 0u);
    
    imageStore(tex2d_contour_pix_marks_dbg_, ivec2(screen_coord.xy), uvec4(contour_pixel, contour_pixel, contour_pixel, 1u));
#endif
}
#endif



/* const uvec2 grid_size_xy = gl_WorkGroupSize.xy;
const uvec2 num_grids_xy = gl_NumWorkGroups.xy; */


