#pragma once

/* Forward-declare in blender:: to match what armature_skinning.cc sees */
namespace blender {
struct Object;
struct Mesh;
struct bDeformGroup;
} /* namespace blender */

#define GVB_DEFAULT_RES_X 256
#define GVB_DEFAULT_RES_Y 256
#define GVB_DEFAULT_RES_Z 128

void geodesic_voxel_bone_weighting(blender::Object *ob,
                                   blender::Mesh *mesh,
                                   float (*verts)[3],
                                   int numbones,
                                   blender::bDeformGroup **dgrouplist,
                                   blender::bDeformGroup **dgroupflip,
                                   float (*root)[3],
                                   float (*tip)[3],
                                   const bool *selected,
                                   const char **r_error_str);