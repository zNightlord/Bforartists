#include "strokegen_mesh_raster_pass.hh"

void StrokegenMeshRasterPass::append_draw_subpass(IndirectDrawGfxResource &gfx_rsc)
{
  draw::PassMain::Sub* subpass = &sub("strokegen raster pass");
  subpass->bind_ssbo(0, *(gfx_rsc.ssbo_mesh));
  subpass->draw_procedural_indirect(GPUPrimType::GPU_PRIM_LINES, *gfx_rsc.ssbo_draw_args); 
}
