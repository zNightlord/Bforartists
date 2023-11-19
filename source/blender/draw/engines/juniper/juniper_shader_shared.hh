#ifndef USE_GPU_SHADER_CREATE_INFO
#  pragma once

namespace blender::juniper {

    using namespace draw;

#endif

struct LightMeta {
    int num_lights;
    float _pad[3];
};
BLI_STATIC_ASSERT_ALIGN(LightMeta, 16)

struct LightData {
    float4x4 object_mat;
#ifndef USE_GPU_SHADER_CREATE_INFO
#  define _right object_mat[0]
#  define _up object_mat[1]
#  define _back object_mat[2]
#  define _position object_mat[3]
#else
#  define _right object_mat[0].xyz
#  define _up object_mat[1].xyz
#  define _back object_mat[2].xyz
#  define _position object_mat[3].xyz
#endif
    float3 color;
    float power;

    float radius;
    int type;
    int2 _pad;

    int4 light_groups;
};
BLI_STATIC_ASSERT_ALIGN(LightData, 16)

#if defined(__cplusplus) && !defined(GPU_SHADER)
using LightDataBuf = draw::StorageArrayBuffer<LightData, 1>;
using LightMetaBuf = draw::StorageBuffer<LightMeta>;

}
#endif
