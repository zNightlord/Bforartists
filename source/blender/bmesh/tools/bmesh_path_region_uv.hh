/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bmesh
 */

#include "BLI_linklist.h"

#include "bmesh_class.hh"

LinkNode *BM_mesh_calc_path_uv_region_vert(BMesh *bm,
                                           BMElem *ele_src,
                                           BMElem *ele_dst,
                                           int cd_loop_uv_offset,
                                           bool (*filter_fn)(BMLoop *, void *user_data),
                                           void *user_data) ATTR_WARN_UNUSED_RESULT
    ATTR_NONNULL(1, 2, 3);

LinkNode *BM_mesh_calc_path_uv_region_edge(BMesh *bm,
                                           BMElem *ele_src,
                                           BMElem *ele_dst,
                                           int cd_loop_uv_offset,
                                           bool (*filter_fn)(BMLoop *, void *user_data),
                                           void *user_data) ATTR_WARN_UNUSED_RESULT
    ATTR_NONNULL(1, 2, 3);

LinkNode *BM_mesh_calc_path_uv_region_face(BMesh *bm,
                                           BMElem *ele_src,
                                           BMElem *ele_dst,
                                           int cd_loop_uv_offset,
                                           bool (*filter_fn)(BMFace *, void *user_data),
                                           void *user_data) ATTR_WARN_UNUSED_RESULT
    ATTR_NONNULL(1, 2, 3);
