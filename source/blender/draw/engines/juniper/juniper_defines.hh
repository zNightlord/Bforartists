// Shared CPU/GPU defines.
#ifndef GPU_SHADER
#  pragma once
#endif

#define SH_TEX_DEPTH_SLOT 1
#define SH_TEX_COLOR_SLOT 2
#define SH_TEX_NORMAL_SLOT 3

#define SH_BUF_LIGHTDATA_SLOT 0
#define SH_BUF_LIGHTMETA_SLOT 1

#define LIGHT_TYPE_SUN 0
#define LIGHT_TYPE_POINT 1

// 0-3 reserved by draw manager
#define UBO_SLOT_MATERIAL_INFO 4
