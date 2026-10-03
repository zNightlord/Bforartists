/* SPDX-FileCopyrightText: 2024 Bforartists. All rights reserved.
 * SPDX-License-Identifier: GPL-2.0-or-later */

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

/**
 * Default alpha bind-smoothness parameter (Section 6, Eq. 8).
 * Paper used alpha = 0.7 across all test characters.
 */
#define GVB_DEFAULT_ALPHA 0.7f

/**
 * Compute skinning weights via geodesic voxel distances.
 *
 * \param ob         Mesh object receiving the weights.
 * \param mesh       Rest-pose mesh data block.
 * \param verts      World-space vertex positions [mesh->verts_num][3].
 * \param numbones   Number of bones.
 * \param dgrouplist Deform group per bone (nullptr = skip that bone).
 * \param dgroupflip Mirror deform group per bone (nullptr = no mirror).
 * \param root       World-space bone root positions [numbones][3].
 * \param tip        World-space bone tip positions  [numbones][3].
 * \param selected   Per-bone selection flags.
 * \param alpha      Bind smoothness [0,1]. Higher = stiffer local bind.
 * \param res        Voxel grid resolution {x, y, z}.
 * \param error_str  Set to a static error string on failure; nullptr = success.
 */
void geodesic_voxel_bone_weighting(Object *ob,
                                   Mesh *mesh,
                                   float (*verts)[3],
                                   int numbones,
                                   bDeformGroup **dgrouplist,
                                   bDeformGroup **dgroupflip,
                                   float (*root)[3],
                                   float (*tip)[3],
                                   const bool *selected,
                                   float alpha,
                                   const int res[3],
                                   const char **error_str);
