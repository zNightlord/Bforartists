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
#include "BLI_index_ranges_builder.hh"
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

namespace blender::bke {

ArmatureDeformGroup build_deform_group_for_vertex_group(
    const IndexMask &universe,
    const Span<float3> positions,
    const Span<MDeformVert> dverts,
    const int def_nr,
    const Bone &bone,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    IndexMaskMemory &memory)
{
  BLI_assert(positions.size() >= universe.min_array_size());
  BLI_assert(dverts.size() >= universe.min_array_size());

  int def_nr_count = 0;
  universe.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr != def_nr) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      ++def_nr_count;
    }
  });
  if (def_nr_count == 0) {
    return {};
  }

  // XXX TODO THIS IS NOT THREADSAFE! NEEDS THREAD-LOCAL DATA AND REDUCE!
  Vector<int> indices;
  indices.reserve(def_nr_count);
  Array<float> weights(def_nr_count);

  universe.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr != def_nr) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      indices.append_unchecked(index);
      weights[indices.size() - 1] = dw.weight;
    }
  });

  /* Bone option to mix with envelope weight. */
  if (use_envelope_multiply && (bone.flag & BONE_MULT_VG_ENV)) {
    for (const int i : indices.index_range()) {
      const int index = indices[i];
      weights[i] *= distfactor_to_bone(
          positions[index], bone.head, bone.tail, bone.rad_head, bone.rad_tail, bone.dist);
    }
  }

  IndexMask mask = IndexMask::from_indices(indices.as_span(), memory);
  return {std::move(mask), std::move(weights)};
}

ArmatureDeformGroup build_deform_group_for_envelope(const IndexMask &universe,
                                                    const Span<float3> positions,
                                                    const Bone &bone,
                                                    const std::optional<float> weight_threshold,
                                                    IndexMaskMemory &memory)
{
  // TODO
  UNUSED_VARS(universe, positions, bone, weight_threshold, memory);
  return {};
}

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

  void finalize(const float total, float3 &r_position) const
  {
    const float scale_factor = 1.0f / total;
    r_position = r_position + position_delta * scale_factor;
  };

  void finalize(const float total, float3 &r_position, float3x3 &r_deform_mat)
  {
    BLI_assert(full_deform);
    const float scale_factor = 1.0f / total;
    r_position = r_position + position_delta * scale_factor;
    r_deform_mat = (deform * scale_factor) * r_deform_mat;
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

  void finalize(const float total, float3 &r_position)
  {
    normalize_dq(&dq, total);
    float3x3 local_deform_mat;
    mul_v3m3_dq(r_position, local_deform_mat.ptr(), &dq);
  }

  void finalize(const float total, float3 &r_position, float3x3 &r_deform_mat)
  {
    BLI_assert(full_deform);
    normalize_dq(&dq, total);
    float3x3 local_deform_mat;
    mul_v3m3_dq(r_position, local_deform_mat.ptr(), &dq);
    r_deform_mat = local_deform_mat * r_deform_mat;
  }
};

static std::optional<int> find_def_nr_from_pose_channel(const bPoseChannel &pchan,
                                                        const ListBase &vertex_groups)
{
  const int def_nr = BLI_findstringindex(&vertex_groups, pchan.name, offsetof(bDeformGroup, name));
  if (def_nr >= 0) {
    return def_nr;
  }
  return std::nullopt;
}

template<typename MixerT>
static void pchan_bone_deform(const bPoseChannel &pchan,
                              const float3 &position,
                              const float weight,
                              MixerT &mixer)
{
  BLI_assert(pchan.bone != nullptr);
  const Bone &bone = *pchan.bone;
  BLI_assert(!(bone.flag & BONE_NO_DEFORM));

  if (bone.segments > 1 && pchan.runtime.bbone_segments == bone.segments) {
    /* Calculate the indices of the 2 affecting b_bone segments. */
    int index;
    float blend;
    BKE_pchan_bbone_deform_segment_index(&pchan, position, &index, &blend);

    mixer.accumulate_bbone(pchan, position, weight * (1.0f - blend), index);
    mixer.accumulate_bbone(pchan, position, weight * blend, index + 1);
  }
  else {
    mixer.accumulate(pchan, position, weight);
  }
}

static IndexMask filter_index_mask(const IndexMask &universe,
                                   const Span<bool> deny_list,
                                   IndexMaskMemory &memory)
{
  constexpr GrainSize grain_size = GrainSize(4096);
  return IndexMask::from_batch_predicate(
      universe,
      grain_size,
      memory,
      [&](const IndexMaskSegment universe_segment, IndexRangesBuilder<int16_t> &builder) {
        for (const int16_t local_index : universe_segment.base_span()) {
          if (!deny_list[universe_segment.offset() + local_index]) {
            builder.add(local_index);
          }
        }
        return universe_segment.offset();
      });
}

template<typename MixerT>
static void deform_with_mixer(
    const ListBase &pose_channels,
    const blender::IndexMask &selection,
    const std::optional<Span<float>> point_weights,
    const Span<PoseChannelDeformGroup> custom_groups,
    const std::optional<ArmatureDeformVertexGroupParams> &vertex_group_params,
    const bool use_envelope,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    const MutableSpan<float3> positions,
    const std::optional<MutableSpan<float3x3>> deform_mats)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  BLI_assert(positions.size() >= selection.min_array_size());

  IndexMaskMemory memory;

  Array<MixerT> mixers(selection.min_array_size(), MixerT{});
  Array<float> contrib(selection.min_array_size(), 0.0f);
  Array<bool> deformed(selection.min_array_size(), false);

  for (const PoseChannelDeformGroup &pchan_group : custom_groups) {
    if (pchan_group.pose_channel == nullptr) {
      continue;
    }
    pchan_group.deform_group.mask.foreach_index(grain_size, [&](const int index, const int pos) {
      const float3 &position = positions[index];
      /* Weights are stored as compressed array. */
      const float weight = pchan_group.deform_group.weights[pos];
      pchan_bone_deform(*pchan_group.pose_channel, position, weight, mixers[index]);
      contrib[index] += weight;
      deformed[index] = true;
    });
  }

  if (vertex_group_params) {
    const IndexMask undeformed_mask = filter_index_mask(selection, deformed, memory);

    int pchan_i;
    LISTBASE_FOREACH_INDEX (bPoseChannel *, pchan, &pose_channels, pchan_i) {
      const Bone *bone = pchan->bone;
      if (!bone || (bone->flag & BONE_NO_DEFORM)) {
        continue;
      }

      const std::optional<int> def_nr = find_def_nr_from_pose_channel(
          *pchan, vertex_group_params->vertex_groups);
      if (!def_nr) {
        continue;
      }

      ArmatureDeformGroup deform_group = build_deform_group_for_vertex_group(
          undeformed_mask,
          positions,
          vertex_group_params->dverts,
          *def_nr,
          *bone,
          weight_threshold,
          use_envelope_multiply,
          memory);
      deform_group.mask.foreach_index(grain_size, [&](const int index, const int pos) {
        const float3 &position = positions[index];
        /* Weights are stored as compressed array. */
        const float weight = deform_group.weights[pos];
        pchan_bone_deform(*pchan, position, weight, mixers[index]);
        contrib[index] += weight;
        deformed[index] = true;
      });
    }
  }

  if (use_envelope) {
    const IndexMask undeformed_mask = filter_index_mask(selection, deformed, memory);
    int pchan_i;
    LISTBASE_FOREACH_INDEX (bPoseChannel *, pchan, &pose_channels, pchan_i) {
      const Bone *bone = pchan->bone;
      if (!bone || (bone->flag & BONE_NO_DEFORM)) {
        continue;
      }

      ArmatureDeformGroup deform_group = build_deform_group_for_envelope(
          undeformed_mask, positions, *bone, weight_threshold, memory);
      // TODO
    }
  }

  const bool full_deform = deform_mats.has_value();
  selection.foreach_index(grain_size, [&](const int index) {
    /* Overall influence. */
    const float mask_weight = point_weights ? (*point_weights)[index] : 1.0f;
    const float weight = contrib[index];
    const float total_weight = weight * mask_weight;
    if (weight_threshold && total_weight <= weight_threshold) {
      return;
    }

    if (full_deform) {
      float3 &position = positions[index];
      float3x3 &deform_mat = (*deform_mats)[index];
      mixers[index].finalize(weight, position, deform_mat);
    }
    else {
      float3 &position = positions[index];
      mixers[index].finalize(weight, position);
    }
  });
}

/* Statically determine the correct mixer type. */
static void deform_values(
    const ListBase &pose_channels,
    const blender::IndexMask &selection,
    const std::optional<Span<float>> point_weights,
    const Span<PoseChannelDeformGroup> custom_groups,
    const std::optional<ArmatureDeformVertexGroupParams> &vertex_group_params,
    const bool use_envelope,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    const ArmatureDeformSkinningMode skinning_mode,
    const MutableSpan<float3> positions,
    const std::optional<MutableSpan<float3x3>> deform_mats)
{
  if (deform_mats) {
    switch (skinning_mode) {
      case ArmatureDeformSkinningMode::Linear:
        deform_with_mixer<BoneDeformLinearMixer<true>>(pose_channels,
                                                       selection,
                                                       point_weights,
                                                       custom_groups,
                                                       vertex_group_params,
                                                       use_envelope,
                                                       weight_threshold,
                                                       use_envelope_multiply,
                                                       positions,
                                                       deform_mats);
        break;
      case ArmatureDeformSkinningMode::DualQuatenrion:
        deform_with_mixer<BoneDeformDualQuaternionMixer<true>>(pose_channels,
                                                               selection,
                                                               point_weights,
                                                               custom_groups,
                                                               vertex_group_params,
                                                               use_envelope,
                                                               weight_threshold,
                                                               use_envelope_multiply,
                                                               positions,
                                                               deform_mats);
        break;
    }
  }
  else {
    switch (skinning_mode) {
      case ArmatureDeformSkinningMode::Linear:
        deform_with_mixer<BoneDeformLinearMixer<false>>(pose_channels,
                                                        selection,
                                                        point_weights,
                                                        custom_groups,
                                                        vertex_group_params,
                                                        use_envelope,
                                                        weight_threshold,
                                                        use_envelope_multiply,
                                                        positions,
                                                        std::nullopt);
        break;
      case ArmatureDeformSkinningMode::DualQuatenrion:
        deform_with_mixer<BoneDeformDualQuaternionMixer<false>>(pose_channels,
                                                                selection,
                                                                point_weights,
                                                                custom_groups,
                                                                vertex_group_params,
                                                                use_envelope,
                                                                weight_threshold,
                                                                use_envelope_multiply,
                                                                positions,
                                                                std::nullopt);
        break;
    }
  }
}

void armature_deform_positions(
    const Object &ob_arm,
    const blender::float4x4 &target_to_world,
    const blender::IndexMask &selection,
    const std::optional<Span<float>> point_weights,
    const Span<PoseChannelDeformGroup> custom_groups,
    const std::optional<ArmatureDeformVertexGroupParams> vertex_group_params,
    bool use_envelope,
    const ArmatureDeformSkinningMode skinning_mode,
    blender::MutableSpan<blender::float3> positions)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  const float4x4 target_to_armature = ob_arm.world_to_object() * target_to_world;
  const float4x4 armature_to_target = math::invert(target_to_world) * ob_arm.object_to_world();

  /* Input coordinates to start from. */
  Array<float3> armature_positions(positions.size());
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    armature_positions[index] = math::transform_point(target_to_armature, positions[index]);
  });

  deform_values(ob_arm.pose->chanbase,
                selection,
                point_weights,
                custom_groups,
                vertex_group_params,
                use_envelope,
                0.0f,
                true,
                skinning_mode,
                armature_positions,
                std::nullopt);

  /* Transform back to target object space. */
  selection.foreach_index(grain_size, [&](const int index) {
    positions[index] = math::transform_point(armature_to_target, armature_positions[index]);
  });
}

void armature_deform_matrices(
    const Object &ob_arm,
    const blender::float4x4 &target_to_world,
    const blender::IndexMask &selection,
    const std::optional<Span<float>> point_weights,
    const Span<PoseChannelDeformGroup> custom_groups,
    const std::optional<ArmatureDeformVertexGroupParams> vertex_group_params,
    const bool use_envelope,
    const ArmatureDeformSkinningMode skinning_mode,
    blender::MutableSpan<blender::float4x4> matrices)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  const float4x4 target_to_armature = ob_arm.world_to_object() * target_to_world;
  const float4x4 armature_to_target = ob_arm.object_to_world() * math::invert(target_to_world);

  /* Input coordinates to start from. */
  Array<float3> armature_positions(matrices.size());
  Array<float3x3> deform_mats(matrices.size());
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    const float4x4 armature_matrix = target_to_armature * matrices[index];
    armature_positions[index] = armature_matrix.location();
    deform_mats[index] = armature_matrix.view<3, 3>();
  });

  deform_values(ob_arm.pose->chanbase,
                selection,
                point_weights,
                custom_groups,
                vertex_group_params,
                use_envelope,
                0.0f,
                true,
                skinning_mode,
                armature_positions,
                std::nullopt);

  /* Transform back to target object space. */
  selection.foreach_index(grain_size, [&](const int index) {
    const float3x3 &deform_mat = deform_mats[index];
    const float4x4 armature_matrix = float4x4(float4(deform_mat.x, 0),
                                              float4(deform_mat.y, 0),
                                              float4(deform_mat.z, 0),
                                              float4(armature_positions[index], 1));
    matrices[index] = armature_to_target * armature_matrix;
  });
}

}  // namespace blender::bke

/** \} */
