/* SPDX-FileCopyrightText: 2021 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup draw
 */

#include "DNA_meshdata_types.h"

#include "BLI_array_utils.hh"

#include "BKE_deform.hh"

#include "draw_subdivision.hh"
#include "extract_mesh.hh"

namespace blender::draw {

static float3 blended_vgroup_color(const MDeformVert *dvert,
                                   const Span<float3> defgroup_colors,
                                   const Span<bool> validmap,
                                   eV3D_Overlay_WPaint_VGroupColorMode mode,
                                   int active_index)
{
  if (!dvert || dvert->totweight == 0) {
    return float3(0.0f);
  }
  const uint defgroup_len = uint(defgroup_colors.size());

  bool use_validmap = false;
  if (!validmap.is_empty()) {
    for (const bool valid : validmap) {
      if (valid) {
        use_validmap = true;
        break;
      }
    }
  }

  float3 result(0.0f);
  float total_weight = 0.0f;

  for (int i = 0; i < dvert->totweight; i++) {
    const uint def_nr = uint(dvert->dw[i].def_nr);

    if (def_nr >= defgroup_len) {
      continue;
    }
    /* Skip non deform vertex groups. */
    if (use_validmap && !validmap[def_nr]) {
      continue;
    }

    /* Active mode only contribute the current active vertex group. */
    if (mode == V3D_OVERLAY_WPAINT_VGROUP_COLOR_ACTIVE && def_nr != active_index) {
      continue;
    }

    const float w = float(dvert->dw[i].weight);
    result += defgroup_colors[def_nr] * w;
    total_weight += w;
  }

  if (mode == V3D_OVERLAY_WPAINT_VGROUP_COLOR_ALL && total_weight > 1.0f) {
    result /= total_weight;
  }

  return result;
}

static float evaluate_vertex_weight(const MDeformVert *dvert, const DRW_MeshWeightState *wstate)
{
  /* Error state. */
  if ((wstate->defgroup_active < 0) && (wstate->defgroup_len > 0)) {
    return -2.0f;
  }

  float input = 0.0f;
  if (wstate->flags & DRW_MESH_WEIGHT_STATE_MULTIPAINT) {
    /* Multi-Paint feature */
    bool is_normalized = (wstate->flags & (DRW_MESH_WEIGHT_STATE_AUTO_NORMALIZE |
                                           DRW_MESH_WEIGHT_STATE_LOCK_RELATIVE));
    input = BKE_defvert_multipaint_collective_weight(dvert,
                                                     wstate->defgroup_len,
                                                     wstate->defgroup_sel,
                                                     wstate->defgroup_sel_count,
                                                     is_normalized);
    /* make it black if the selected groups have no weight on a vertex */
    if (input == 0.0f) {
      return -1.0f;
    }
  }
  else {
    /* default, non tricky behavior */
    input = BKE_defvert_find_weight(dvert, wstate->defgroup_active);

    if (input == 0.0f) {
      switch (wstate->alert_mode) {
        case OB_DRAW_GROUPUSER_ACTIVE:
          return -1.0f;
          break;
        case OB_DRAW_GROUPUSER_ALL:
          if (BKE_defvert_is_weight_zero(dvert, wstate->defgroup_len)) {
            return -1.0f;
          }
          break;
      }
    }
  }

  /* Lock-Relative: display the fraction of current weight vs total unlocked weight. */
  if (wstate->flags & DRW_MESH_WEIGHT_STATE_LOCK_RELATIVE) {
    input = BKE_defvert_lock_relative_weight(
        input, dvert, wstate->defgroup_len, wstate->defgroup_locked, wstate->defgroup_unlocked);
  }

  CLAMP(input, 0.0f, 1.0f);
  return input;
}

static void extract_weights_mesh(const MeshRenderData &mr,
                                 const DRW_MeshWeightState &weight_state,
                                 MutableSpan<float> vbo_data)
{
  const Mesh &mesh = *mr.mesh;
  const Span<MDeformVert> dverts = mesh.deform_verts();
  if (dverts.is_empty()) {
    vbo_data.fill(weight_state.alert_mode == OB_DRAW_GROUPUSER_NONE ? 0.0f : -1.0f);
    return;
  }

  threading::parallel_for(vbo_data.index_range(), 1024, [&](const IndexRange range) {
    for (const int vert : range) {
      vbo_data[vert] = evaluate_vertex_weight(&dverts[vert], &weight_state);
    }
  });
}

static void extract_weights_bm(const MeshRenderData &mr,
                               const DRW_MeshWeightState &weight_state,
                               MutableSpan<float> vbo_data)
{
  const BMesh &bm = *mr.bm;
  const int offset = CustomData_get_offset(&bm.vdata, CD_MDEFORMVERT);
  if (offset == -1) {
    vbo_data.fill(weight_state.alert_mode == OB_DRAW_GROUPUSER_NONE ? 0.0f : -1.0f);
    return;
  }

  threading::parallel_for(IndexRange(bm.totvert), 2048, [&](const IndexRange range) {
    for (const int vert_index : range) {
      const BMVert &vert = *BM_vert_at_index(&const_cast<BMesh &>(bm), vert_index);
      vbo_data[vert_index] = evaluate_vertex_weight(
          static_cast<const MDeformVert *>(BM_ELEM_CD_GET_VOID_P(&vert, offset)), &weight_state);
    }
  });
}

gpu::VertBufPtr extract_weights(const MeshRenderData &mr,
                                const MeshBatchCache &batch_cache,
                                MeshBufferCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("weight",
                                                              gpu::VertAttrType::SFLOAT_32);

  gpu::VertBufPtr vert_vbo = gpu::VertBufPtr(GPU_vertbuf_create_with_format(format));
  GPU_vertbuf_data_alloc(*vert_vbo, mr.verts_num);
  MutableSpan<float> vbo_data = vert_vbo->data<float>();

  const DRW_MeshWeightState &weight_state = batch_cache.weight_state;
  if (weight_state.defgroup_active == -1) {
    vbo_data.fill(weight_state.alert_mode == OB_DRAW_GROUPUSER_NONE ? 0.0f : -1.0f);
  }
  else {
    if (mr.extract_type == MeshExtractType::Mesh) {
      extract_weights_mesh(mr, weight_state, vbo_data);
    }
    else {
      extract_weights_bm(mr, weight_state, vbo_data);
    }
  }

  gpu::VertBufPtr vbo = gpu::VertBufPtr(GPU_vertbuf_create_on_device(format, mr.corners_num));

  gather_vert_to_corner_gpu(mr, cache, *vert_vbo, *vbo);

  return vbo;
}

gpu::VertBufPtr extract_weights_subdiv(const MeshRenderData &mr,
                                       const DRWSubdivCache &subdiv_cache,
                                       const MeshBatchCache &batch_cache,
                                       MeshBufferCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("weight",
                                                              gpu::VertAttrType::SFLOAT_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(
      GPU_vertbuf_create_on_device(format, subdiv_cache.num_subdiv_loops));

  gpu::VertBufPtr coarse_weights = extract_weights(mr, batch_cache, cache);
  draw_subdiv_interp_custom_data(subdiv_cache, *coarse_weights, *vbo, GPU_COMP_F32, 1, 0);
  return vbo;
}

gpu::VertBufPtr extract_weight_vgroup_blended_color(const MeshRenderData &mr,
                                                    const MeshBatchCache &batch_cache,
                                                    MeshBufferCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("vgroup_color_blended",
                                                              gpu::VertAttrType::SFLOAT_32_32_32);

  gpu::VertBufPtr vert_vbo = gpu::VertBufPtr(GPU_vertbuf_create_with_format(format));
  GPU_vertbuf_data_alloc(*vert_vbo, mr.verts_num);
  MutableSpan<float3> vert_data = vert_vbo->data<float3>();

  const DRW_MeshWeightState &weight_state = batch_cache.weight_state;
  const int active_index = weight_state.defgroup_active;
  const eV3D_Overlay_WPaint_VGroupColorMode mode = weight_state.vgroup_color_mode;
  /* Nothing to compute. */
  if (mode == V3D_OVERLAY_WPAINT_VGROUP_COLOR_OFF) {
    vert_data.fill(float3(0.0f));
    return vert_vbo;
  }

  const Span<float3> defgroup_colors(weight_state.defgroup_colors, weight_state.defgroup_len);
  const Span<bool> validmap(weight_state.defgroup_validmap,
                            weight_state.defgroup_validmap ? weight_state.defgroup_len : 0);

  /* Handles Mesh or BMesh API cases.*/
  if (mr.extract_type == MeshExtractType::Mesh) {
    const Mesh &mesh = *mr.mesh;
    const Span<MDeformVert> dverts = mesh.deform_verts();
    if (dverts.is_empty()) {
      vert_data.fill(float3(0.0f));
      return vert_vbo;
    }

    threading::parallel_for(vert_data.index_range(), 1024, [&](const IndexRange range) {
      for (const int vert : range) {
        vert_data[vert] = blended_vgroup_color(
            &dverts[vert], defgroup_colors, validmap, mode, active_index);
      }
    });
  }
  else {
    const BMesh &bm = *mr.bm;
    const int offset = CustomData_get_offset(&bm.vdata, CD_MDEFORMVERT);
    if (offset == -1) {
      vert_data.fill(float3(0.0f));
      return vert_vbo;
    }
    threading::parallel_for(IndexRange(bm.totvert), 2048, [&](const IndexRange range) {
      for (const int vert_index : range) {
        const BMVert &vert = *BM_vert_at_index(&const_cast<BMesh &>(bm), vert_index);
        vert_data[vert_index] = blended_vgroup_color(
            static_cast<const MDeformVert *>(BM_ELEM_CD_GET_VOID_P(&vert, offset)),
            defgroup_colors,
            validmap,
            mode,
            active_index);
      }
    });
  }

  gpu::VertBufPtr vbo = gpu::VertBufPtr(GPU_vertbuf_create_on_device(format, mr.corners_num));
  gather_vert_to_corner_gpu(mr, cache, *vert_vbo, *vbo);
  return vbo;
}

gpu::VertBufPtr extract_weight_vgroup_blended_color_subdiv(const MeshRenderData &mr,
                                                           const DRWSubdivCache &subdiv_cache,
                                                           const MeshBatchCache &batch_cache,
                                                           MeshBufferCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("vgroup_color_blended",
                                                              gpu::VertAttrType::SFLOAT_32_32_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(
      GPU_vertbuf_create_on_device(format, subdiv_cache.num_subdiv_loops));

  gpu::VertBufPtr coarse = extract_weight_vgroup_blended_color(mr, batch_cache, cache);
  draw_subdiv_interp_custom_data(subdiv_cache, *coarse, *vbo, GPU_COMP_F32, 3, 0);
  return vbo;
}

}  // namespace blender::draw
