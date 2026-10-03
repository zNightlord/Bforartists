/** \file
 * \ingroup edarmature
 *
 * Geodesic Voxel Binding — CPU implementation.
 * Reference: Dionne & de Lasa, "Geodesic Voxel Binding for Production
 * Character Meshes", SCA 2013.
 */

#pragma once


struct Object;
struct Mesh;
struct bDeformGroup;

/**
 * Default voxel grid resolution from the paper (Section 8).
 * 256×256×128 worked well for all 11 production test meshes.
 */
#define GVB_DEFAULT_RES_X 256
#define GVB_DEFAULT_RES_Y 256
#define GVB_DEFAULT_RES_Z 128
#define GVB_DEFAULT_ALPHA 0.7f

void geodesic_voxel_bone_weighting(Object *ob,
                                   Mesh *mesh,
                                   float (*verts)[3],
                                   int numbones,
                                   bDeformGroup **dgrouplist,
                                   bDeformGroup **dgroupflip,
                                   float (*root)[3],
                                   float (*tip)[3],
                                   const bool *selected,
                                   const char **r_error_str);