
/* Inputs: 
float ssbo_vbo_full_[]
uint ssbo_vert_flags_[]
uint ssbo_edge_flags_[]
*/

vec3 ld_vpos(uint vtx_id)
{
	uint base_addr = vtx_id * 3; 
	return vec3(
        ssbo_vbo_full_[base_addr], 
        ssbo_vbo_full_[base_addr+1], 
        ssbo_vbo_full_[base_addr+2]
    ); 
}
void st_vpos(uint vtx_id, vec3 vpos)
{
    uint base_addr = vtx_id * 3; 
    ssbo_vbo_full_[base_addr]   = vpos.x;
    ssbo_vbo_full_[base_addr+1] = vpos.y;
    ssbo_vbo_full_[base_addr+2] = vpos.z;
}







void main()
{

}
