
#if defined(_KERNEL_MULTICOMPILE_DRAW_CONTOUR__EDGES)
void main()
{
    out_col = color;  
    out_normal = normalize(normal);
    out_tangent.xyz = normalize(tangent.xyz);
}
#endif

#if defined(_KERNEL_MULTICOMPILE_DRAW_CONTOUR__2D_SAMPLES)
void main()
{
    out_col = color;  
    // out_normal = normalize(normal);
    // out_tangent.xyz = normalize(tangent.xyz);
    if (out_col.a < .1f) discard; 
}
#endif