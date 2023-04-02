
void main()
{
	const uint groupId =  gl_LocalInvocationID.x;
	const uint idx = gl_GlobalInvocationID.x;
	
	const bool PCS_ib_fmt_u16_ = pcs_ib_fmt_u16;
	
	uint EdgeId = idx;
	
	uint ibo_addr = get_prim_base_addr_ibo(PCS_ib_fmt_u16_, EdgeId, 4u);
	uint ibo_data = buf_ibo[ibo_addr];
	uint vbo_addr = get_vbo_addr(PCS_ib_fmt_u16_, ibo_data, EdgeId, 4u);
	
	buf_strokegen_mesh_pool[idx] = ibo_data; // only for testing
}