#pragma once

#include "DRW_render.h"
#include "draw_shader_shared.h"
#include "draw_pass.hh"
#include "juniper_shaders.hh"
#include "GPU_common_types.h"
#include "GPU_texture.h"
#include "juniper_defines.hh"
#include "DRW_gpu_wrapper.hh"

#include "juniper_shader_shared.hh"

namespace blender::juniper {
class JNPR;

// CPU version of shader light object.
struct Light : public LightData {
public:
    bool used = false;

    void sync(const Object* ob);
};

class LightManager {
private:
    JNPR &jnpr_;
    Map<Object*, Light> scene_lights;
    LightDataBuf light_buf;
    LightMetaBuf meta_buf;
public:
    LightManager(JNPR &jnpr) : jnpr_(jnpr){};

    void begin_sync();
    void sync_light(Object *ob);
    void end_sync();

    template<typename T> void bind_resources(draw::detail::PassBase<T> *pass)
    {
      pass->bind_ssbo(SH_BUF_LIGHTDATA_SLOT, &light_buf);
      pass->bind_ssbo(SH_BUF_LIGHTMETA_SLOT, &meta_buf);
    }
};

}