/* SPDX-FileCopyrightText: 2021 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup draw
 */

#include "DNA_meshdata_types.h"

#include "BLI_array_utils.hh"
#include "BLI_noise.hh"

#include "BKE_deform.hh"

#include "draw_subdivision.hh"
#include "extract_mesh.hh"

namespace blender::draw {

static float3 hash_group_color(int def_nr, int random_id)
{
  return blender::noise::hash_float_to_float3(float(def_nr + random_id + 1));
}

static float3 blended_vgroup_color(const MDeformVert *dvert,
                                   int random_id,
                                   int effective_mode,
                                   int active_index,
                                   const bool *validmap,
                                   int defgroup_len)
{
  if (!dvert || dvert->totweight == 0) {
    return float3(0.0f);
  }
  float3 result(0.0f);
  // float total_weight = 0.0f;

  for (int i = 0; i < dvert->totweight; i++) {
    const int def_nr = dvert->dw[i].def_nr;

    if (def_nr < 0 || def_nr >= defgroup_len) {
      continue;
    }
    /* Filter by deform bones when armature exists. */
    if (validmap && !validmap[def_nr]) {
      continue;
    }
    /* ACTIVE mode: only contribute the active group.
     * When effective_mode is ALL (including colored_vertex override),
     * this filter is skipped so all groups contribute to vertex dot color. */
    if (effective_mode == 1 && def_nr != active_index) {
      continue;
    }

    const float w = float(dvert->dw[i].weight);
    result += hash_group_color(def_nr, random_id) * w;
    // total_weight += w;
  }

  /* Normalize for ALL mode so colors show proportional group influence.
   * ACTIVE mode stays un-normalized — intensity encodes weight magnitude. */
  // if (effective_mode != VGROUP_COLOR_ACTIVE && total_weight > 0.0f) {
  //   result /= total_weight;
  // }

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

  Array<float> weights(dverts.size());
  threading::parallel_for(weights.index_range(), 1024, [&](const IndexRange range) {
    for (const int vert : range) {
      weights[vert] = evaluate_vertex_weight(&dverts[vert], &weight_state);
    }
  });
  array_utils::gather(weights.as_span(), mr.corner_verts, vbo_data);
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

  threading::parallel_for(IndexRange(bm.totface), 2048, [&](const IndexRange range) {
    for (const int face_index : range) {
      const BMFace &face = *BM_face_at_index(&const_cast<BMesh &>(bm), face_index);
      const BMLoop *loop = BM_FACE_FIRST_LOOP(&face);
      for ([[maybe_unused]] const int i : IndexRange(face.len)) {
        const int index = BM_elem_index_get(loop);
        vbo_data[index] = evaluate_vertex_weight(
            static_cast<const MDeformVert *>(BM_ELEM_CD_GET_VOID_P(loop->v, offset)),
            &weight_state);
        loop = loop->next;
      }
    }
  });
}

gpu::VertBufPtr extract_weights(const MeshRenderData &mr, const MeshBatchCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("weight",
                                                              gpu::VertAttrType::SFLOAT_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(GPU_vertbuf_create_with_format(format));
  GPU_vertbuf_data_alloc(*vbo, mr.corners_num);
  MutableSpan<float> vbo_data = vbo->data<float>();

  const DRW_MeshWeightState &weight_state = cache.weight_state;
  if (weight_state.defgroup_active == -1) {
    vbo_data.fill(weight_state.alert_mode == OB_DRAW_GROUPUSER_NONE ? 0.0f : -1.0f);
    return vbo;
  }

  if (mr.extract_type == MeshExtractType::Mesh) {
    extract_weights_mesh(mr, weight_state, vbo_data);
  }
  else {
    extract_weights_bm(mr, weight_state, vbo_data);
  }
  return vbo;
}

gpu::VertBufPtr extract_weights_subdiv(const MeshRenderData &mr,
                                       const DRWSubdivCache &subdiv_cache,
                                       const MeshBatchCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("weight",
                                                              gpu::VertAttrType::SFLOAT_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(
      GPU_vertbuf_create_on_device(format, subdiv_cache.num_subdiv_loops));

  gpu::VertBufPtr coarse_weights = extract_weights(mr, cache);
  draw_subdiv_interp_custom_data(subdiv_cache, *coarse_weights, *vbo, GPU_COMP_F32, 1, 0);
  return vbo;
}

gpu::VertBufPtr extract_weight_vgroup_blended_color(const MeshRenderData &mr,
                                                    const MeshBatchCache &cache)
{
  static GPUVertFormat format = GPU_vertformat_from_attribute("vgroup_color_blended",
                                                              gpu::VertAttrType::SFLOAT_32_32_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(GPU_vertbuf_create_with_format(format));
  GPU_vertbuf_data_alloc(*vbo, mr.corners_num);
  MutableSpan<float3> vbo_data = vbo->data<float3>();

  const DRW_MeshWeightState &weight_state = cache.weight_state;
  const int random_id = weight_state.vgroup_color_random_id;
  const int active_index = weight_state.defgroup_active;
  const bool *validmap = weight_state.defgroup_validmap;
  const int defgroup_len = weight_state.defgroup_len;
  const bool colored_vertex = weight_state.colored_vertex;

  /* When colored_vertex is on but surface mode is NONE,
   * treat as ALL for VBO computation — point shader uses the color,
   * surface shader ignores it via its own push constant. */
  const int effective_mode = colored_vertex ?
                                 V3D_OVERLAY_WPAINT_VGROUP_COLOR_ALL :
                                 weight_state.vgroup_color_mode;

  /* Nothing to compute */
  if (effective_mode == 0) {
    vbo_data.fill(float3(0.0f));
    return vbo;
  }

  if (mr.extract_type == MeshExtractType::Mesh) {
    const Mesh &mesh = *mr.mesh;
    const Span<MDeformVert> dverts = mesh.deform_verts();
    if (dverts.is_empty()) {
      vbo_data.fill(float3(0.0f));
      return vbo;
    }
    Array<float3> colors(dverts.size());
    threading::parallel_for(colors.index_range(), 1024, [&](const IndexRange range) {
      for (const int vert : range) {
        colors[vert] = blended_vgroup_color(
            &dverts[vert], random_id, effective_mode, active_index, validmap, defgroup_len);
      }
    });
    array_utils::gather(colors.as_span(), mr.corner_verts, vbo_data);
  }
  else {
    const BMesh &bm = *mr.bm;
    const int offset = CustomData_get_offset(&bm.vdata, CD_MDEFORMVERT);
    if (offset == -1) {
      vbo_data.fill(float3(0.0f));
      return vbo;
    }
    threading::parallel_for(IndexRange(bm.totface), 2048, [&](const IndexRange range) {
      for (const int face_index : range) {
        const BMFace &face = *BM_face_at_index(&const_cast<BMesh &>(bm), face_index);
        const BMLoop *loop = BM_FACE_FIRST_LOOP(&face);
        for ([[maybe_unused]] const int i : IndexRange(face.len)) {
          const int index = BM_elem_index_get(loop);
          vbo_data[index] = blended_vgroup_color(
              static_cast<const MDeformVert *>(BM_ELEM_CD_GET_VOID_P(loop->v, offset)),
              random_id,
              effective_mode,
              active_index,
              validmap,
              defgroup_len);
          loop = loop->next;
        }
      }
    });
  }
  return vbo;
}

gpu::VertBufPtr extract_weight_vgroup_blended_color_subdiv(const MeshRenderData &mr,
                                                           const DRWSubdivCache &subdiv_cache,
                                                           const MeshBatchCache &cache)
{
  GPUVertFormat format{};
  GPU_vertformat_attr_add(&format, "vgroup_color_blended", gpu::VertAttrType::SFLOAT_32_32_32);

  gpu::VertBufPtr vbo = gpu::VertBufPtr(
      GPU_vertbuf_create_on_device(format, subdiv_cache.num_subdiv_loops));

  gpu::VertBufPtr coarse_colors = extract_weight_vgroup_blended_color(mr, cache);
  draw_subdiv_interp_custom_data(subdiv_cache, *coarse_colors, *vbo, GPU_COMP_F32, 3, 0);

  return vbo;
}

}  // namespace blender::draw
