/* SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 *
 * Deform coordinates by a armature object (used by modifier).
 */

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "MEM_guardedalloc.h"

#include "BLI_array_utils.hh"
#include "BLI_index_mask.hh"
#include "BLI_listbase.h"
#include "BLI_listbase_wrapper.hh"
#include "BLI_math_matrix.h"
#include "BLI_math_quaternion.hh"
#include "BLI_math_rotation.h"
#include "BLI_math_vector.h"
#include "BLI_task.h"
#include "BLI_task.hh"

#include "DNA_armature_types.h"
#include "DNA_lattice_types.h"
#include "DNA_listBase.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BKE_action.hh"
#include "BKE_armature.hh"
#include "BKE_customdata.hh"
#include "BKE_deform.hh"
#include "BKE_editmesh.hh"
#include "BKE_lattice.hh"
#include "BKE_mesh.hh"

#include "CLG_log.h"

static CLG_LogRef LOG = {"geom.armature_deform"};

/* -------------------------------------------------------------------- */
/** \name Armature Deform Internal Utilities
 * \{ */

static float bone_envelope_falloff(const float distance_squared,
                                   const float closest_radius,
                                   const float falloff_distance)
{
  using namespace blender;

  if (distance_squared < closest_radius * closest_radius) {
    return 1.0f;
  }

  /* Zero influence beyond falloff distance. */
  if (falloff_distance == 0.0f ||
      distance_squared >= math::square(closest_radius + falloff_distance))
  {
    return 0.0f;
  }

  /* Compute influence from envelope over the falloff distance. */
  const float dist_envelope = sqrtf(distance_squared) - closest_radius;
  return 1.0f - (dist_envelope * dist_envelope) / (falloff_distance * falloff_distance);
}

float distfactor_to_bone(const blender::float3 &position,
                         const blender::float3 &head,
                         const blender::float3 &tail,
                         const float radius_head,
                         const float radius_tail,
                         const float falloff_distance)
{
  using namespace blender;

  float bone_length;
  const float3 bone_axis = math::normalize_and_get_length(tail - head, bone_length);
  /* Distance along the bone axis from head. */
  const float height = math::dot(position - head, bone_axis);

  if (height < 0.0f) {
    /* Below the start of the bone use the head radius. */
    const float distance_squared = math::distance_squared(position, head);
    return bone_envelope_falloff(distance_squared, radius_head, falloff_distance);
  }
  else if (height > bone_length) {
    /* After the end of the bone use the tail radius. */
    const float distance_squared = math::distance_squared(tail, position);
    return bone_envelope_falloff(distance_squared, radius_tail, falloff_distance);
  }
  else {
    /* Interpolate radius. */
    const float distance_squared = math::distance_squared(position, head) - height * height;
    const float closest_radius = bone_length != 0.0f ? math::interpolate(radius_head,
                                                                         radius_tail,
                                                                         height / bone_length) :
                                                       radius_head;
    return bone_envelope_falloff(distance_squared, closest_radius, falloff_distance);
  }
}

namespace blender::bke {

/**
 * Utility class for accumulating linear bone deformation.
 * If full_deform is true the deformation matrix is also computed.
 */
template<bool full_deform> struct BoneDeformLinearMixer {

  float3 position_delta = float3(0.0f);
  float3x3 deform = float3x3::zero();

  void accumulate(const bPoseChannel &pchan, const float3 &co, const float weight)
  {
    const float4x4 &pose_mat = float4x4(pchan.chan_mat);

    position_delta += weight * (math::transform_point(pose_mat, co) - co);
    if constexpr (full_deform) {
      deform += weight * pose_mat.view<3, 3>();
    }
  }

  void accumulate_bbone(const bPoseChannel &pchan,
                        const float3 &co,
                        const float weight,
                        const int index)
  {
    const Span<float4x4> pose_mats = Span<Mat4>(pchan.runtime.bbone_deform_mats,
                                                pchan.runtime.bbone_segments + 2)
                                         .cast<float4x4>();
    const float4x4 &pose_mat = pose_mats[index + 1];

    position_delta += weight * (math::transform_point(pose_mat, co) - co);
    if constexpr (full_deform) {
      deform += weight * pose_mat.view<3, 3>();
    }
  }

  void finalize(const float3 & /*co*/,
                float total,
                float armature_weight,
                float3 &r_delta_co,
                float3x3 &r_deform_mat)
  {
    const float scale_factor = armature_weight / total;
    r_delta_co = position_delta * scale_factor;
    r_deform_mat = deform * scale_factor;
  };
};

/**
 * Utility class for accumulating dual quaternion bone deformation.
 * If full_deform is true the deformation matrix is also computed.
 */
template<bool full_deform> struct BoneDeformDualQuaternionMixer {
  DualQuat dq = {};

  void accumulate(const bPoseChannel &pchan, const float3 &co, const float weight)
  {
    const DualQuat &deform_quat = pchan.runtime.deform_dual_quat;

    add_weighted_dq_dq_pivot(&dq, &deform_quat, co, weight, full_deform);
  }

  void accumulate_bbone(const bPoseChannel &pchan,
                        const float3 &co,
                        const float weight,
                        const int index)
  {
    const Span<DualQuat> quats = {pchan.runtime.bbone_dual_quats,
                                  pchan.runtime.bbone_segments + 1};
    const DualQuat &deform_quat = quats[index];

    add_weighted_dq_dq_pivot(&dq, &deform_quat, co, weight, full_deform);
  }

  void finalize(const float3 &co,
                float total,
                float armature_weight,
                float3 &r_delta_co,
                float3x3 &r_deform_mat)
  {
    normalize_dq(&dq, total);
    float3 dco = co;
    float3x3 dmat;
    mul_v3m3_dq(dco, dmat.ptr(), &dq);
    r_delta_co = (dco - co) * armature_weight;
    /* Quaternion already is scale corrected. */
    r_deform_mat = dmat;
  }
};

/* Add interpolated deformation along a b-bone segment of the pose channel. */
template<typename MixerT>
static void b_bone_deform(const bPoseChannel &pchan,
                          const float3 &co,
                          const float weight,
                          MixerT &mixer)
{
  /* Calculate the indices of the 2 affecting b_bone segments. */
  int index;
  float blend;
  BKE_pchan_bbone_deform_segment_index(&pchan, co, &index, &blend);

  mixer.accumulate_bbone(pchan, co, weight * (1.0f - blend), index);
  mixer.accumulate_bbone(pchan, co, weight * blend, index + 1);
}

/* Add bone deformation based on envelope distance. */
template<typename MixerT>
static float dist_bone_deform(const bPoseChannel &pchan, const float3 &co, MixerT &mixer)
{
  const Bone *bone = pchan.bone;
  if (bone == nullptr || bone->weight == 0.0f) {
    return 0.0f;
  }

  const float fac = distfactor_to_bone(co,
                                       float3(bone->arm_head),
                                       float3(bone->arm_tail),
                                       bone->rad_head,
                                       bone->rad_tail,
                                       bone->dist);
  if (fac == 0.0f) {
    return 0.0f;
  }

  const float weight = fac * bone->weight;
  if (bone->segments > 1 && pchan.runtime.bbone_segments == bone->segments) {
    b_bone_deform(pchan, co, weight, mixer);
  }
  else {
    mixer.accumulate(pchan, co, weight);
  }

  return weight;
}

/* Add bone deformation based on vertex group weight. */
template<typename MixerT>
static float pchan_bone_deform(const bPoseChannel &pchan,
                               const float weight,
                               const float3 &co,
                               MixerT &mixer)
{
  const Bone *bone = pchan.bone;

  if (!weight) {
    return 0.0f;
  }

  if (bone->segments > 1 && pchan.runtime.bbone_segments == bone->segments) {
    b_bone_deform(pchan, co, weight, mixer);
  }
  else {
    mixer.accumulate(pchan, co, weight);
  }

  return weight;
}

}  // namespace blender::bke

/** \} */

/* -------------------------------------------------------------------- */
/** \name Armature Deform #BKE_armature_deform_coords API
 *
 * #BKE_armature_deform_coords and related functions.
 * \{ */

namespace blender::bke {

struct ArmatureDeformParams {
  MutableSpan<float3> vert_coords;
  std::optional<MutableSpan<float3x3>> vert_deform_mats;
  std::optional<Span<float3>> vert_coords_prev;

  bool use_envelope = false;
  bool invert_vgroup = false;
  bool use_dverts = false;

  const IndexMask *selection = nullptr;
  std::optional<Span<float>> vert_influence;

  /* List of all pose channels on the target object. */
  ConstListBaseWrapper<bPoseChannel> pose_channels = {{nullptr, nullptr}};
  /* Maps vertex group index (def_nr) to pose channels, if vertex groups are used.
   * Vertex groups used for deform can be different from the target object vertex groups list,
   * the def_nr needs to be mapped to the correct pose channel first. */
  Array<bPoseChannel *> pose_channel_by_vertex_group;

  float4x4 target_to_armature = float4x4::identity();
  float4x4 armature_to_target = float4x4::identity();
};

static ArmatureDeformParams get_armature_deform_params(
    const Object &ob_arm,
    const Object &ob_target,
    const ListBase *defbase,
    MutableSpan<float3> vert_coords,
    std::optional<Span<float3>> vert_coords_prev,
    std::optional<MutableSpan<float3x3>> vert_deform_mats,
    const int deformflag,
    std::optional<Span<float>> vert_influence,
    const bool try_use_dverts)
{
  const bool dverts_supported = BKE_object_supports_vertex_groups(&ob_target);

  ArmatureDeformParams deform_params;
  deform_params.vert_coords = vert_coords;
  deform_params.vert_deform_mats = vert_deform_mats;
  deform_params.vert_coords_prev = vert_coords_prev;
  deform_params.use_envelope = bool(deformflag & ARM_DEF_ENVELOPE);
  deform_params.invert_vgroup = bool(deformflag & ARM_DEF_INVERT_VGROUP);
  deform_params.vert_influence = vert_influence;

  deform_params.pose_channels = {ob_arm.pose->chanbase};
  deform_params.use_dverts = try_use_dverts && dverts_supported && (deformflag & ARM_DEF_VGROUP);
  if (deform_params.use_dverts) {
    const int defbase_len = BLI_listbase_count(defbase);
    deform_params.pose_channel_by_vertex_group.reinitialize(defbase_len);
    /* TODO(sergey): Some considerations here:
     *
     * - Check whether keeping this consistent across frames gives speedup.
     */
    int i;
    LISTBASE_FOREACH_INDEX (bDeformGroup *, dg, defbase, i) {
      bPoseChannel *pchan = BKE_pose_channel_find_name(ob_arm.pose, dg->name);
      /* Exclude non-deforming bones. */
      deform_params.pose_channel_by_vertex_group[i] = (pchan &&
                                                       !(pchan->bone->flag & BONE_NO_DEFORM)) ?
                                                          pchan :
                                                          nullptr;
    }
  }

/* TODO using the existing matrices directly is better, but fails tests because old code was
 * doing a double-inverse of the object matrix, leading to small differences on the order of 10^-5.
 * Test data needs to be updated if the transforms change. */
#if 0
  deform_params.target_to_armature = ob_arm.world_to_object() * ob_target.object_to_world();
  deform_params.armature_to_target = ob_target.world_to_object() * ob_arm.object_to_world();
#else
  deform_params.armature_to_target = ob_target.world_to_object() * ob_arm.object_to_world();
  deform_params.target_to_armature = math::invert(deform_params.armature_to_target);
#endif

  return deform_params;
}

/* Alternative for defining parameters without a target object. */
static ArmatureDeformParams get_armature_deform_params(
    const Object &ob_arm,
    const float4x4 &target_to_world,
    MutableSpan<float3> vert_coords,
    std::optional<Span<float3>> vert_coords_prev,
    std::optional<MutableSpan<float3x3>> vert_deform_mats,
    const bool use_envelope,
    const std::optional<ListBase> vertex_groups,
    const bool invert_vertex_groups,
    const IndexMask &selection,
    std::optional<Span<float>> vert_influence)
{
  ArmatureDeformParams deform_params;
  deform_params.vert_coords = vert_coords;
  deform_params.vert_deform_mats = vert_deform_mats;
  deform_params.vert_coords_prev = vert_coords_prev;
  deform_params.use_envelope = use_envelope;
  deform_params.invert_vgroup = invert_vertex_groups;
  deform_params.selection = &selection;
  deform_params.vert_influence = vert_influence;

  deform_params.pose_channels = {ob_arm.pose->chanbase};
  if (vertex_groups) {
    deform_params.use_dverts = true;

    const ListBase *defbase = &(*vertex_groups);
    const int defbase_len = BLI_listbase_count(defbase);
    deform_params.pose_channel_by_vertex_group.reinitialize(defbase_len);
    /* TODO(sergey): Some considerations here:
     *
     * - Check whether keeping this consistent across frames gives speedup.
     */
    int i;
    LISTBASE_FOREACH_INDEX (bDeformGroup *, dg, defbase, i) {
      bPoseChannel *pchan = BKE_pose_channel_find_name(ob_arm.pose, dg->name);
      /* Exclude non-deforming bones. */
      deform_params.pose_channel_by_vertex_group[i] = (pchan &&
                                                       !(pchan->bone->flag & BONE_NO_DEFORM)) ?
                                                          pchan :
                                                          nullptr;
    }
  }

/* TODO using the existing matrices directly is better, but fails tests because old code was
 * doing a double-inverse of the object matrix, leading to small differences on the order of 10^-5.
 * Test data needs to be updated if the transforms change. */
#if 0
  deform_params.target_to_armature = ob_arm.world_to_object() * target_to_world;
  deform_params.armature_to_target = math::invert(target_to_world) * ob_arm.object_to_world();
#else
  deform_params.armature_to_target = math::invert(target_to_world) * ob_arm.object_to_world();
  deform_params.target_to_armature = math::invert(deform_params.armature_to_target);
#endif

  return deform_params;
}

struct VertexGroupMask {
  int def_nr;
  IndexMask mask;
  /* Compressed weights, matches positions in the index mask. */
  Array<float> weights;
};

static Vector<VertexGroupMask> vertex_groups_to_masks(const Span<MDeformVert> dverts,
                                                      const IndexMask &selection,
                                                      const IndexRange def_nr_range,
                                                      const std::optional<float> weight_threshold,
                                                      const bool discard_empty_groups,
                                                      IndexMaskMemory &memory)
{
  Array<int> def_nr_counts(def_nr_range.size() + 1, 0);
  selection.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (!def_nr_range.contains(dw.def_nr)) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      ++def_nr_counts[def_nr_range[dw.def_nr]];
    }
  });
  const OffsetIndices verts_by_def_nr = offset_indices::accumulate_counts_to_offsets(
      def_nr_counts);

  Array<int> all_indices(verts_by_def_nr.total_size());
  Array<float> all_weights(verts_by_def_nr.total_size());
  Array<int> next_entry = def_nr_counts.as_span().drop_back(1);
  selection.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (!def_nr_range.contains(dw.def_nr)) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      const int def_nr_i = dw.def_nr - def_nr_range.start();
      int &entry = next_entry[def_nr_i];
      BLI_assert(entry < def_nr_counts[def_nr_i + 1]);
      all_indices[entry] = index;
      all_weights[entry] = dw.weight;
      ++entry;
    }
  });

  Vector<VertexGroupMask> group_masks;
  group_masks.reserve(def_nr_range.size());
  for (const int def_nr_i : def_nr_range.index_range()) {
    const IndexRange entries = verts_by_def_nr[def_nr_i];
    if (discard_empty_groups && entries.is_empty()) {
      continue;
    }

    const int def_nr = def_nr_range[def_nr_i];
    IndexMask mask = IndexMask::from_indices(all_indices.as_span().slice(entries), memory);
    Array<float> weights = all_weights.as_span().slice(entries);
    group_masks.append_unchecked({def_nr, std::move(mask), std::move(weights)});
  }

  return group_masks;
}

/* Accumulate bone deformations using the mixer implementation. */
template<typename MixerT>
static void armature_vert_task_with_mixer(const ArmatureDeformParams &params,
                                          const int i,
                                          const MDeformVert *dvert,
                                          MixerT &mixer)
{
  const bool full_deform = params.vert_deform_mats.has_value();

  /* Overall influence, can change by masking with a vertex group. */
  float armature_weight = 1.0f;
  float prevco_weight = 0.0f; /* weight for optional cached vertexcos */
  if (params.vert_influence) {
    const float mask_weight = (*params.vert_influence)[i];
    /* On multi-modifier the mask is used to blend with previous coordinates. */
    if (params.vert_coords_prev) {
      prevco_weight = params.invert_vgroup ? mask_weight : 1.0f - mask_weight;
      if (prevco_weight == 1.0f) {
        return;
      }
    }
    else {
      armature_weight = params.invert_vgroup ? 1.0f - mask_weight : mask_weight;
      if (armature_weight == 0.0f) {
        return;
      }
    }
  }

  /* Input coordinates to start from. */
  float3 co = params.vert_coords_prev ? (*params.vert_coords_prev)[i] : params.vert_coords[i];
  /* Transform to armature space. */
  co = math::transform_point(params.target_to_armature, co);

  float contrib = 0.0f;
  bool deformed = false;
  /* Apply vertex group deformation if enabled. */
  if (params.use_dverts && dvert) {
    /* Range of valid def_nr in MDeformWeight. */
    const IndexRange def_nr_range = params.pose_channel_by_vertex_group.index_range();
    const Span<MDeformWeight> dweights(dvert->dw, dvert->totweight);
    for (const auto &dw : dweights) {
      const bPoseChannel *pchan = def_nr_range.contains(dw.def_nr) ?
                                      params.pose_channel_by_vertex_group[dw.def_nr] :
                                      nullptr;
      if (pchan == nullptr) {
        continue;
      }

      float weight = dw.weight;

      /* Bone option to mix with envelope weight. */
      const Bone *bone = pchan->bone;
      if (bone && bone->flag & BONE_MULT_VG_ENV) {
        weight *= distfactor_to_bone(co,
                                     float3(bone->arm_head),
                                     float3(bone->arm_tail),
                                     bone->rad_head,
                                     bone->rad_tail,
                                     bone->dist);
      }

      contrib += pchan_bone_deform(*pchan, weight, co, mixer);
      deformed = true;
    }
  }
  /* Use envelope if enabled and no bone deformed the vertex yet. */
  if (!deformed && params.use_envelope) {
    for (const bPoseChannel *pchan : params.pose_channels) {
      if (!(pchan->bone->flag & BONE_NO_DEFORM)) {
        contrib += dist_bone_deform(*pchan, co, mixer);
      }
    }
  }

  /* TODO Actually should be EPSILON? Weight values and contrib can be like 10e-39 small. */
  constexpr float contrib_threshold = 0.0001f;
  if (contrib > contrib_threshold) {
    float3 delta_co;
    float3x3 local_deform_mat;
    mixer.finalize(co, contrib, armature_weight, delta_co, local_deform_mat);

    co += delta_co;
    if (full_deform) {
      float3x3 &deform_mat = (*params.vert_deform_mats)[i];
      const float3x3 armature_to_target = params.armature_to_target.view<3, 3>();
      const float3x3 target_to_armature = params.target_to_armature.view<3, 3>();
      deform_mat = armature_to_target * local_deform_mat * target_to_armature * deform_mat;
    }
  }

  /* Transform back to target object space. */
  co = math::transform_point(params.armature_to_target, co);

  /* Multi-modifier: Interpolate with previous modifier position using the vertex group mask. */
  if (params.vert_coords_prev) {
    copy_v3_v3(params.vert_coords[i], math::interpolate(co, params.vert_coords[i], prevco_weight));
  }
  else {
    copy_v3_v3(params.vert_coords[i], co);
  }
}

/* Accumulate bone deformations using the mixer implementation. */
template<typename MixerT>
static void armature_vert_task_with_mixer(const ArmatureDeformParams &params,
                                          const IndexMask &selection,
                                          const std::optional<Span<MDeformVert>> dverts)
{
  constexpr GrainSize grain_size = GrainSize(1024);
  const bool full_deform = params.vert_deform_mats.has_value();

  /* Input coordinates to start from. */
  Array<float3> positions = params.vert_coords_prev ? *params.vert_coords_prev :
                                                      params.vert_coords;
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    float3 &position = positions[index];
    position = math::transform_point(params.target_to_armature, position);
  });

  Array<MixerT> mixers(positions.size(), MixerT{});
  Array<float> contrib(positions.size(), 0.0f);
  Array<bool> deformed(positions.size(), false);
  /* Apply vertex group deformation if enabled. */
  if (params.use_dverts && dverts) {
    /* Range of valid def_nr in MDeformWeight. */
    const IndexRange def_nr_range = params.pose_channel_by_vertex_group.index_range();
    IndexMaskMemory memory;
    /* Build index masks for all vertex groups.
     * TODO It may be worthwile to cache this - unless in the future the sparse attribute storage
     * will provide vertex weights in this order anyway, which would make index masks a lot
     * cheaper to create. */
    const Vector<VertexGroupMask> group_masks = vertex_groups_to_masks(
        *dverts, selection, def_nr_range, std::nullopt, true, memory);

    for (const VertexGroupMask &group_mask : group_masks) {
      BLI_assert(def_nr_range.contains(group_mask.def_nr));
      const bPoseChannel *pchan = params.pose_channel_by_vertex_group[group_mask.def_nr];
      if (pchan == nullptr) {
        continue;
      }
      /* Bone option to mix with envelope weight. */
      const Bone *bone = pchan->bone;
      const bool use_envelope_multiply = (bone && bone->flag & BONE_MULT_VG_ENV);

      group_mask.mask.foreach_index(grain_size, [&](const int index, const int pos) {
        const float3 &position = positions[index];
        /* Weights are stored as compressed array. */
        float weight = group_mask.weights[pos];
        /* Multiply with bone envelope falloff if enabled for this bone. */
        if (use_envelope_multiply) {
          weight *= distfactor_to_bone(position,
                                       float3(bone->arm_head),
                                       float3(bone->arm_tail),
                                       bone->rad_head,
                                       bone->rad_tail,
                                       bone->dist);
        }

        contrib[index] += pchan_bone_deform(*pchan, weight, position, mixers[index]);
        deformed[index] = true;
      });
    }
  }
  /* Use envelope if enabled and no bone deformed the vertex yet. */
  if (params.use_envelope) {
    selection.foreach_index(grain_size, [&](const int index) {
      if (deformed[index]) {
        return;
      }

      const float3 &position = positions[index];
      for (const bPoseChannel *pchan : params.pose_channels) {
        if (!(pchan->bone->flag & BONE_NO_DEFORM)) {
          contrib[index] += dist_bone_deform(*pchan, position, mixers[index]);
        }
      }
    });
  }

  selection.foreach_index(grain_size, [&](const int index) {
    float3 &position = positions[index];

    /* Overall influence, can change by masking with a vertex group. */
    float armature_weight = 1.0f;
    float prevco_weight = 0.0f; /* weight for optional cached vertexcos */
    if (params.vert_influence) {
      const float mask_weight = (*params.vert_influence)[index];
      /* On multi-modifier the mask is used to blend with previous coordinates. */
      if (params.vert_coords_prev) {
        prevco_weight = params.invert_vgroup ? mask_weight : 1.0f - mask_weight;
        if (prevco_weight == 1.0f) {
          return;
        }
      }
      else {
        armature_weight = params.invert_vgroup ? 1.0f - mask_weight : mask_weight;
        if (armature_weight == 0.0f) {
          return;
        }
      }
    }

    /* TODO Actually should be EPSILON? Weight values and contrib can be like 10e-39 small. */
    constexpr float contrib_threshold = 0.0001f;
    if (contrib[index] > contrib_threshold) {
      float3 delta_co;
      float3x3 local_deform_mat;
      mixers[index].finalize(
          position, contrib[index], armature_weight, delta_co, local_deform_mat);

      position += delta_co;
      if (full_deform) {
        float3x3 &deform_mat = (*params.vert_deform_mats)[index];
        const float3x3 armature_to_target = params.armature_to_target.view<3, 3>();
        const float3x3 target_to_armature = params.target_to_armature.view<3, 3>();
        deform_mat = armature_to_target * local_deform_mat * target_to_armature * deform_mat;
      }
    }

    /* Transform back to target object space. */
    position = math::transform_point(params.armature_to_target, position);

    /* Multi-modifier: Interpolate with previous modifier position using the vertex group mask. */
    if (params.vert_coords_prev) {
      copy_v3_v3(params.vert_coords[index],
                 math::interpolate(position, params.vert_coords[index], prevco_weight));
    }
    else {
      copy_v3_v3(params.vert_coords[index], position);
    }
  });
}

/* Accumulate bone deformations for a vertex. */
static void armature_vert_task_with_dvert(const ArmatureDeformParams &deform_params,
                                          const int i,
                                          const MDeformVert *dvert,
                                          const bool use_quaternion)
{
  const bool full_deform = deform_params.vert_deform_mats.has_value();
  if (use_quaternion) {
    if (full_deform) {
      bke::BoneDeformDualQuaternionMixer<true> mixer;
      armature_vert_task_with_mixer(deform_params, i, dvert, mixer);
    }
    else {
      bke::BoneDeformDualQuaternionMixer<false> mixer;
      armature_vert_task_with_mixer(deform_params, i, dvert, mixer);
    }
  }
  else {
    if (full_deform) {
      bke::BoneDeformLinearMixer<true> mixer;
      armature_vert_task_with_mixer(deform_params, i, dvert, mixer);
    }
    else {
      bke::BoneDeformLinearMixer<false> mixer;
      armature_vert_task_with_mixer(deform_params, i, dvert, mixer);
    }
  }
}

/* Accumulate bone deformations. */
static void armature_vert_task_with_dvert(const ArmatureDeformParams &deform_params,
                                          const IndexMask &selection,
                                          const std::optional<Span<MDeformVert>> dverts,
                                          const bool use_quaternion)
{
  const bool full_deform = deform_params.vert_deform_mats.has_value();
  if (use_quaternion) {
    if (full_deform) {
      armature_vert_task_with_mixer<bke::BoneDeformDualQuaternionMixer<true>>(
          deform_params, selection, dverts);
    }
    else {
      armature_vert_task_with_mixer<bke::BoneDeformDualQuaternionMixer<false>>(
          deform_params, selection, dverts);
    }
  }
  else {
    if (full_deform) {
      armature_vert_task_with_mixer<bke::BoneDeformLinearMixer<true>>(
          deform_params, selection, dverts);
    }
    else {
      armature_vert_task_with_mixer<bke::BoneDeformLinearMixer<false>>(
          deform_params, selection, dverts);
    }
  }
}

static void armature_deform_coords(const MutableSpan<float3> vert_coords,
                                   const bool use_quaternion,
                                   const std::optional<Span<MDeformVert>> dverts,
                                   const ArmatureDeformParams &deform_params)
{
  const IndexMask &masked_selection = deform_params.selection ?
                                          *deform_params.selection :
                                          IndexMask(vert_coords.index_range());
  const bool use_dverts = deform_params.use_dverts && dverts;
  /* Split selection to avoid checking each item for valid dvert individually. */
  const IndexMask dvert_selection = use_dverts ?
                                        masked_selection.slice_content(dverts->index_range()) :
                                        IndexMask();
  const IndexMask no_dvert_selection = use_dverts ?
                                           masked_selection.slice_content(
                                               dverts->size(), masked_selection.min_array_size()) :
                                           masked_selection;

  armature_vert_task_with_dvert(deform_params, dvert_selection, dverts, use_quaternion);
  armature_vert_task_with_dvert(deform_params, no_dvert_selection, std::nullopt, use_quaternion);
}

static std::optional<Array<float>> get_vert_influence_from_dverts(const ListBase &defbase,
                                                                  StringRefNull defgrp_name,
                                                                  const Span<MDeformVert> dverts)
{
  const int armature_def_nr = BKE_defgroup_name_index(&defbase, defgrp_name);
  if (armature_def_nr < 0) {
    return std::nullopt;
  }

  Array<float> vert_influence(dverts.size());
  threading::parallel_for(vert_influence.index_range(), 1024, [&](const IndexRange range) {
    for (const int i : range) {
      const float mask_weight = BKE_defvert_find_weight(&dverts[i], armature_def_nr);
      /* On multi-modifier the mask is used to blend with previous coordinates. */
      vert_influence[i] = mask_weight;
    }
  });

  return vert_influence;
}

static std::optional<Array<float>> get_vert_influence_from_bmesh(const ListBase &defbase,
                                                                 StringRefNull defgrp_name,
                                                                 const BMEditMesh &em_target)
{
  const int armature_def_nr = BKE_defgroup_name_index(&defbase, defgrp_name);
  if (armature_def_nr < 0) {
    return std::nullopt;
  }

  const int cd_dvert_offset = CustomData_get_offset(&em_target.bm->vdata, CD_MDEFORMVERT);
  Array<float> vert_influence(em_target.bm->totvert);

  struct UserData {
    int cd_dvert_offset;
    int armature_def_nr;
    MutableSpan<float> vert_influence;
  };
  UserData userdata = {cd_dvert_offset, armature_def_nr, vert_influence};

  /* While this could cause an extra loop over mesh data, in most cases this will
   * have already been properly set. */
  BM_mesh_elem_index_ensure(em_target.bm, BM_VERT);

  TaskParallelSettings settings;
  BLI_parallel_mempool_settings_defaults(&settings);
  BLI_task_parallel_mempool(
      em_target.bm->vpool,
      &userdata,
      [](void *__restrict userdata_v,
         MempoolIterData *iter,
         const TaskParallelTLS *__restrict /*tls*/) {
        const UserData &userdata = *static_cast<UserData *>(userdata_v);

        BMVert *v = reinterpret_cast<BMVert *>(iter);
        const MDeformVert *dvert = static_cast<const MDeformVert *>(
            BM_ELEM_CD_GET_VOID_P(v, userdata.cd_dvert_offset));
        const float mask_weight = BKE_defvert_find_weight(dvert, userdata.armature_def_nr);

        userdata.vert_influence[BM_elem_index_get(v)] = mask_weight;
      },
      &settings);

  return vert_influence;
}

struct ArmatureEditMeshUserdata {
  bool use_quaternion = false;
  int cd_dvert_offset = -1;

  ArmatureDeformParams deform_params;
};

template<bool use_dvert>
static void armature_vert_task_editmesh(void *__restrict userdata,
                                        MempoolIterData *iter,
                                        const TaskParallelTLS *__restrict /*tls*/)
{
  const ArmatureEditMeshUserdata &data = *static_cast<const ArmatureEditMeshUserdata *>(userdata);
  BMVert *v = reinterpret_cast<BMVert *>(iter);
  const MDeformVert *dvert = use_dvert ? static_cast<const MDeformVert *>(
                                             BM_ELEM_CD_GET_VOID_P(v, data.cd_dvert_offset)) :
                                         nullptr;
  armature_vert_task_with_dvert(
      data.deform_params, BM_elem_index_get(v), dvert, data.use_quaternion);
}

static void armature_deform_editmesh(const Object &ob_arm,
                                     const Object &ob_target,
                                     const ListBase *defbase,
                                     const MutableSpan<float3> vert_coords,
                                     const std::optional<MutableSpan<float3x3>> vert_deform_mats,
                                     const int deformflag,
                                     const std::optional<Span<float3>> vert_coords_prev,
                                     blender::StringRefNull defgrp_name,
                                     const BMEditMesh &em_target,
                                     const int cd_dvert_offset)
{
  std::optional<Array<float>> vert_influence;
  /* Note: For legacy reasons vertex masking is only enabled if vertex group deform is enabled.
   * It actually works fine and is enabled for regular meshes even if vertex group deform is
   * disabled. */
  if (defbase && !defgrp_name.is_empty() && cd_dvert_offset >= 0 && (deformflag & ARM_DEF_VGROUP))
  {
    vert_influence = get_vert_influence_from_bmesh(*defbase, defgrp_name, em_target);
  }
  ArmatureDeformParams deform_params = get_armature_deform_params(ob_arm,
                                                                  ob_target,
                                                                  defbase,
                                                                  vert_coords,
                                                                  vert_coords_prev,
                                                                  vert_deform_mats,
                                                                  deformflag,
                                                                  vert_influence,
                                                                  cd_dvert_offset >= 0);

  ArmatureEditMeshUserdata data{};
  data.use_quaternion = bool(deformflag & ARM_DEF_QUATERNION);
  data.cd_dvert_offset = cd_dvert_offset;
  data.deform_params = std::move(deform_params);

  /* While this could cause an extra loop over mesh data, in most cases this will
   * have already been properly set. */
  BM_mesh_elem_index_ensure(em_target.bm, BM_VERT);

  TaskParallelSettings settings;
  BLI_parallel_mempool_settings_defaults(&settings);

  if (deform_params.use_dverts) {
    BLI_task_parallel_mempool(
        em_target.bm->vpool, &data, armature_vert_task_editmesh<true>, &settings);
  }
  else {
    BLI_task_parallel_mempool(
        em_target.bm->vpool, &data, armature_vert_task_editmesh<false>, &settings);
  }
}

static bool verify_armature_deform_valid(const Object &ob_arm)
{
  /* Not supported in armature edit mode or without pose data. */
  const bArmature *arm = static_cast<const bArmature *>(ob_arm.data);
  if (arm->edbo || (ob_arm.pose == nullptr)) {
    return false;
  }
  if ((ob_arm.pose->flag & POSE_RECALC) != 0) {
    CLOG_ERROR(&LOG,
               "Trying to evaluate influence of armature '%s' which needs Pose recalc!",
               ob_arm.id.name);
    BLI_assert_unreachable();
  }
  return true;
}

}  // namespace blender::bke

void BKE_armature_deform_coords_with_curves(
    const Object &ob_arm,
    const Object &ob_target,
    const ListBase *defbase,
    blender::MutableSpan<blender::float3> vert_coords,
    std::optional<blender::Span<blender::float3>> vert_coords_prev,
    std::optional<blender::MutableSpan<blender::float3x3>> vert_deform_mats,
    blender::Span<MDeformVert> dverts,
    int deformflag,
    blender::StringRefNull defgrp_name)
{
  using namespace blender;

  if (!bke::verify_armature_deform_valid(ob_arm)) {
    return;
  }

  /* Vertex groups must be provided explicitly, cannot rely on object vertex groups since this is
   * used for Grease Pencil layers as well. */
  BLI_assert(dverts.size() == vert_coords.size());

  std::optional<Array<float>> vert_influence;
  if (defbase && !defgrp_name.is_empty()) {
    vert_influence = bke::get_vert_influence_from_dverts(*defbase, defgrp_name, dverts);
  }
  bke::ArmatureDeformParams deform_params = bke::get_armature_deform_params(ob_arm,
                                                                            ob_target,
                                                                            defbase,
                                                                            vert_coords,
                                                                            vert_coords_prev,
                                                                            vert_deform_mats,
                                                                            deformflag,
                                                                            vert_influence,
                                                                            true);
  const bool use_quaternion = bool(deformflag & ARM_DEF_QUATERNION);
  bke::armature_deform_coords(vert_coords, use_quaternion, dverts, deform_params);
}

void BKE_armature_deform_coords_with_mesh(
    const Object &ob_arm,
    const Object &ob_target,
    blender::MutableSpan<blender::float3> vert_coords,
    std::optional<blender::Span<blender::float3>> vert_coords_prev,
    std::optional<blender::MutableSpan<blender::float3x3>> vert_deform_mats,
    int deformflag,
    blender::StringRefNull defgrp_name,
    const Mesh *me_target)
{
  using namespace blender;

  if (!bke::verify_armature_deform_valid(ob_arm)) {
    return;
  }

  /* Note armature modifier on legacy curves calls this, so vertex groups are not guaranteed to
   * exist. */
  const ID *id_target = static_cast<const ID *>(ob_target.data);
  const ListBase *defbase = nullptr;
  if (me_target) {
    /* Use the vertex groups from the evaluated mesh that is being deformed. */
    defbase = BKE_id_defgroup_list_get(&me_target->id);
  }
  else if (BKE_id_supports_vertex_groups(id_target)) {
    /* Take the vertex groups from the original object data. */
    defbase = BKE_id_defgroup_list_get(id_target);
  }

  blender::Span<MDeformVert> dverts;
  if (ob_target.type == OB_MESH) {
    if (me_target == nullptr) {
      me_target = static_cast<const Mesh *>(ob_target.data);
    }
    dverts = me_target->deform_verts();
  }
  else if (ob_target.type == OB_LATTICE) {
    const Lattice *lt = static_cast<const Lattice *>(ob_target.data);
    if (lt->dvert != nullptr) {
      dverts = blender::Span<MDeformVert>(lt->dvert, lt->pntsu * lt->pntsv * lt->pntsw);
    }
  }

  std::optional<Span<MDeformVert>> dverts_opt;
  if ((me_target && !me_target->deform_verts().is_empty()) || dverts.size() == vert_coords.size())
  {
    dverts_opt = dverts;
  }
  if (me_target) {
    BLI_assert(vert_coords.size() <= me_target->verts_num);
  }

  std::optional<Array<float>> vert_influence;
  if (defbase && !defgrp_name.is_empty() && dverts_opt) {
    vert_influence = bke::get_vert_influence_from_dverts(*defbase, defgrp_name, *dverts_opt);
  }
  bke::ArmatureDeformParams deform_params = bke::get_armature_deform_params(
      ob_arm,
      ob_target,
      defbase,
      vert_coords,
      vert_coords_prev,
      vert_deform_mats,
      deformflag,
      vert_influence,
      dverts_opt.has_value());
  const bool use_quaternion = bool(deformflag & ARM_DEF_QUATERNION);
  bke::armature_deform_coords(vert_coords, use_quaternion, dverts_opt, deform_params);
}

void BKE_armature_deform_coords_with_editmesh(
    const Object &ob_arm,
    const Object &ob_target,
    blender::MutableSpan<blender::float3> vert_coords,
    std::optional<blender::Span<blender::float3>> vert_coords_prev,
    std::optional<blender::MutableSpan<blender::float3x3>> vert_deform_mats,
    int deformflag,
    blender::StringRefNull defgrp_name,
    const BMEditMesh &em_target)
{
  using namespace blender;

  if (!bke::verify_armature_deform_valid(ob_arm)) {
    return;
  }

  const ListBase *defbase = BKE_id_defgroup_list_get(static_cast<const ID *>(ob_target.data));
  const int cd_dvert_offset = CustomData_get_offset(&em_target.bm->vdata, CD_MDEFORMVERT);
  bke::armature_deform_editmesh(ob_arm,
                                ob_target,
                                defbase,
                                vert_coords,
                                vert_deform_mats,
                                deformflag,
                                vert_coords_prev,
                                defgrp_name,
                                em_target,
                                cd_dvert_offset);
}

/** \} */
