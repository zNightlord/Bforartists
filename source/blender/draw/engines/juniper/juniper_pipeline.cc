
#include "juniper_pipeline.hh"
#include "juniper_instance.hh"

namespace blender::juniper {

void NPRPipeline::sync() {
  {
    npr_opaque_ps_.init();
    npr_opaque_ps_.bind_texture(SH_TEX_NORMAL_SLOT, jnpr_.buffers.normal_tx);
    npr_opaque_ps_.bind_texture(SH_TEX_DEPTH_SLOT, jnpr_.buffers.depth_copy_tx);
    jnpr_.lights.bind_resources(&npr_opaque_ps_);

    npr_opaque_sub_ = &npr_opaque_ps_.sub("Test Subpass");
    npr_opaque_sub_->state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_EQUAL);
  }

  {
    npr_prepass_ps_.init();

    npr_prepass_sub_ = &npr_prepass_ps_.sub("Prepass Sub");
    npr_prepass_sub_->state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS);
  }
}

void NPRPipeline::render(draw::View &view) {
  jnpr_.manager->submit(npr_opaque_ps_, view);
}

void NPRPipeline::prepass(draw::View &view) {
  jnpr_.manager->submit(npr_prepass_ps_, view);
}

} // namespace blender::juniper

