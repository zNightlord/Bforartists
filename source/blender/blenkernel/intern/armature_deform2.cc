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
#include "BLI_enumerable_thread_specific.hh"
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
  /* dverts span is empty if no vertex groups are defined. */
  if (dverts.is_empty()) {
    return {};
  }
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
      break;
    }
  });
  if (def_nr_count == 0) {
    return {};
  }

  struct GroupIndices {
    Vector<int> indices;
    Vector<float> weights;
    GroupIndices(const int max_size)
    {
      indices.reserve(max_size);
      weights.reserve(max_size);
    }
    GroupIndices(const GroupIndices &other)
    {
      indices.reserve(other.indices.capacity());
      weights.reserve(other.weights.capacity());
    }
  };

  threading::EnumerableThreadSpecific<GroupIndices> tls(
      [&]() { return GroupIndices(def_nr_count); });
  universe.foreach_index(GrainSize(1024), [&](const int index) {
    GroupIndices &data = tls.local();
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr != def_nr) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      data.indices.append_unchecked(index);
      data.weights.append_unchecked(dw.weight);
      break;
    }
  });
  GroupIndices result(def_nr_count);
  for (const GroupIndices &data : tls) {
    result.indices.extend_unchecked(data.indices);
    result.weights.extend_unchecked(data.weights);
  }

  /* Bone option to mix with envelope weight. */
  if (use_envelope_multiply && (bone.flag & BONE_MULT_VG_ENV)) {
    for (const int i : result.indices.index_range()) {
      const int index = result.indices[i];
      result.weights[i] *= distfactor_to_bone(
          positions[index], bone.head, bone.tail, bone.rad_head, bone.rad_tail, bone.dist);
    }
  }

  IndexMask mask = IndexMask::from_indices(result.indices.as_span(), memory);
  return {std::move(mask), Array<float>(result.weights.as_span())};
}

static void build_deform_groups_for_all_vertex_groups(
    const IndexMask &universe,
    const Span<float3> positions,
    const Span<MDeformVert> dverts,
    const Span<const Bone *> bone_by_def_nr,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    IndexMaskMemory &memory,
    MutableSpan<ArmatureDeformGroup> r_group_by_def_nr)
{
  /* dverts span is empty if no vertex groups are defined. */
  if (dverts.is_empty()) {
    return;
  }
  BLI_assert(positions.size() >= universe.min_array_size());
  BLI_assert(dverts.size() >= universe.min_array_size());

  const int bones_num = bone_by_def_nr.size();
  Array<int> def_nr_offsets(bones_num + 1, 0);
  universe.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr >= bones_num) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      ++def_nr_offsets[dw.def_nr];
    }
  });
  const OffsetIndices verts_by_def_nr = offset_indices::accumulate_counts_to_offsets(
      def_nr_offsets);

  struct GroupIndices {
    Array<int> all_indices;
    Array<float> all_weights;
    Array<int> next_entry;
    GroupIndices(const OffsetIndices<int> verts_by_def_nr)
    {
      all_indices.reinitialize(verts_by_def_nr.total_size());
      all_weights.reinitialize(verts_by_def_nr.total_size());
      next_entry = verts_by_def_nr.data().drop_back(1);
    }
    GroupIndices(const GroupIndices &other)
    {
      all_indices.reinitialize(other.all_indices.size());
      all_weights.reinitialize(other.all_weights.size());
      next_entry = other.next_entry;
    }
  };

  threading::EnumerableThreadSpecific<GroupIndices> tls(
      [&]() { return GroupIndices(verts_by_def_nr); });

  universe.foreach_index(GrainSize(1024), [&](const int index) {
    GroupIndices &data = tls.local();
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr >= bones_num) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      int &entry = data.next_entry[dw.def_nr];
      BLI_assert(entry < def_nr_offsets[dw.def_nr + 1]);
      data.all_indices[entry] = index;
      data.all_weights[entry] = dw.weight;
      ++entry;
    }
  });
  GroupIndices result(verts_by_def_nr);
  for (const GroupIndices &data : tls) {
    for (const int def_nr : result.next_entry.index_range()) {
      /* `next_entry` in the source data is the end positions of the data range. */
      const IndexRange src_range = IndexRange::from_begin_end(def_nr_offsets[def_nr],
                                                              data.next_entry[def_nr]);
      /* `next_entry` in the result data is current write position of collected data. */
      const IndexRange dst_range = IndexRange::from_begin_size(result.next_entry[def_nr],
                                                               src_range.size());

      result.all_indices.as_mutable_span().slice(dst_range).copy_from(
          data.all_indices.as_span().slice(src_range));
      result.all_weights.as_mutable_span().slice(dst_range).copy_from(
          data.all_weights.as_span().slice(src_range));
      result.next_entry[def_nr] = dst_range.one_after_last();
    }
  }
#ifndef NDEBUG
  /* Verify that the whole range has been filled. */
  for (const int def_nr : result.next_entry.index_range()) {
    BLI_assert(result.next_entry[def_nr] == def_nr_offsets[def_nr + 1]);
  }
#endif

  for (const int def_nr : IndexRange(bones_num)) {
    const IndexRange entries = verts_by_def_nr[def_nr];
    const Span<int> indices = result.all_indices.as_span().slice(entries);
    MutableSpan<float> weights = result.all_weights.as_mutable_span().slice(entries);

    /* Bone option to mix with envelope weight. */
    if (use_envelope_multiply) {
      const Bone *bone = bone_by_def_nr[def_nr];
      if (bone && bool(bone->flag & BONE_MULT_VG_ENV)) {
        for (const int i : indices.index_range()) {
          weights[i] *= distfactor_to_bone(positions[indices[i]],
                                           bone->head,
                                           bone->tail,
                                           bone->rad_head,
                                           bone->rad_tail,
                                           bone->dist);
        }
      }
    }

    IndexMask mask = IndexMask::from_indices(indices, memory);
    r_group_by_def_nr[def_nr] = {std::move(mask), Array<float>(weights.as_span())};
  }
}

static std::optional<int> find_def_nr_from_pose_channel(const bPoseChannel &pchan,
                                                        const ListBase &vertex_groups)
{
  const int def_nr = BLI_findstringindex(&vertex_groups, pchan.name, offsetof(bDeformGroup, name));
  if (def_nr >= 0) {
    return def_nr;
  }
  return std::nullopt;
}

Vector<PoseChannelDeformGroup> pose_channel_groups_from_vertex_groups(
    const ListBase &pose_channels,
    const blender::IndexMask &mask,
    const std::optional<ArmatureDeformVertexGroupParams> &vertex_group_params,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    const Span<float3> positions,
    IndexMaskMemory &memory)
{
  const int vertex_group_num = BLI_listbase_count(&vertex_group_params->vertex_groups);
  Array<const bPoseChannel *> pose_channel_by_def_nr(vertex_group_num, nullptr);
  Array<const Bone *> bone_by_def_nr(vertex_group_num, nullptr);
  LISTBASE_FOREACH (bPoseChannel *, pchan, &pose_channels) {
    const Bone *bone = pchan->bone;
    if (!bone || (bone->flag & BONE_NO_DEFORM)) {
      continue;
    }

    const std::optional<int> def_nr = find_def_nr_from_pose_channel(
        *pchan, vertex_group_params->vertex_groups);
    if (!def_nr) {
      continue;
    }

    pose_channel_by_def_nr[*def_nr] = pchan;
    bone_by_def_nr[*def_nr] = bone;
  }

  Array<ArmatureDeformGroup> deform_group_by_def_nr(vertex_group_num);
  build_deform_groups_for_all_vertex_groups(mask,
                                            positions,
                                            vertex_group_params->dverts,
                                            bone_by_def_nr,
                                            weight_threshold,
                                            use_envelope_multiply,
                                            memory,
                                            deform_group_by_def_nr);

  Vector<PoseChannelDeformGroup> pchan_groups;
  pchan_groups.reserve(vertex_group_num);
  for (const int def_nr : deform_group_by_def_nr.index_range()) {
    const bPoseChannel *pchan = pose_channel_by_def_nr[def_nr];
    ArmatureDeformGroup &deform_group = deform_group_by_def_nr[def_nr];
    if (!pchan || deform_group.mask.is_empty()) {
      continue;
    }
    pchan_groups.append_unchecked({std::move(deform_group), pchan});
  }
  return pchan_groups;
}

static void build_deform_groups_for_envelopes(const IndexMask &universe,
                                              const Span<float3> positions,
                                              const Span<const Bone *> bones,
                                              const std::optional<float> weight_threshold,
                                              IndexMaskMemory &memory,
                                              MutableSpan<ArmatureDeformGroup> r_groups)
{
  constexpr GrainSize grain_size(1024);

  for (const int i : bones.index_range()) {
    const Bone &bone = *bones[i];

    Array<bool> inside_falloff(universe.size());
    Array<float> all_weights(universe.size());
    universe.foreach_index(grain_size, [&](const int index, const int pos) {
      const float weight = distfactor_to_bone(
          positions[index], bone.head, bone.tail, bone.rad_head, bone.rad_tail, bone.dist);

      inside_falloff[pos] = (weight > 0.0f);
      all_weights[pos] = weight;
    });

    IndexMask mask = IndexMask::from_bools(inside_falloff, memory);
    Array<float> weights(mask.size());
    array_utils::gather(all_weights.as_span(), mask, weights.as_mutable_span());
    r_groups[i] = {std::move(mask), std::move(weights)};
  }
}

Vector<PoseChannelDeformGroup> pose_channel_groups_from_envelopes(
    const ListBase &pose_channels,
    const blender::IndexMask &mask,
    const std::optional<float> weight_threshold,
    const Span<float3> positions,
    IndexMaskMemory &memory)
{
  const int pose_channel_num = BLI_listbase_count(&pose_channels);
  Array<const Bone *> bones(pose_channel_num, nullptr);
  int pchan_i;
  LISTBASE_FOREACH_INDEX (bPoseChannel *, pchan, &pose_channels, pchan_i) {
    const Bone *bone = pchan->bone;
    if (!bone || (bone->flag & BONE_NO_DEFORM)) {
      continue;
    }
    bones[pchan_i] = pchan->bone;
  }

  Array<ArmatureDeformGroup> deform_groups(pose_channel_num);
  build_deform_groups_for_envelopes(
      mask, positions, bones, weight_threshold, memory, deform_groups);

  Vector<PoseChannelDeformGroup> pchan_groups;
  pchan_groups.reserve(pose_channel_num);
  LISTBASE_FOREACH_INDEX (bPoseChannel *, pchan, &pose_channels, pchan_i) {
    ArmatureDeformGroup &deform_group = deform_groups[pchan_i];
    if (!pchan || deform_group.mask.is_empty()) {
      continue;
    }
    pchan_groups.append_unchecked({std::move(deform_group), pchan});
  }
  return pchan_groups;
}

/**
 * Utility class for accumulating linear bone deformation.
 * If full_deform is true the deformation matrix is also computed.
 */
template<bool full_deform> struct BoneDeformLinearMixer {
  Array<float3> position_delta;
  Array<float3x3> deform;

  BoneDeformLinearMixer(const int size)
      : position_delta(size, float3(0.0f)), deform(size, float3x3::identity())
  {
  }

  void accumulate(const bPoseChannel &pchan, const int index, const float3 &co, const float weight)
  {
    const float4x4 &pose_mat = float4x4(pchan.chan_mat);

    position_delta[index] += weight * (math::transform_point(pose_mat, co) - co);
    if constexpr (full_deform) {
      deform[index] += weight * pose_mat.view<3, 3>();
    }
  }

  void accumulate_bbone(const bPoseChannel &pchan,
                        const int index,
                        const float3 &co,
                        const float weight,
                        const int bbone_index)
  {
    const Span<float4x4> pose_mats = Span<Mat4>(pchan.runtime.bbone_deform_mats,
                                                pchan.runtime.bbone_segments + 2)
                                         .cast<float4x4>();
    const float4x4 &pose_mat = pose_mats[bbone_index + 1];

    position_delta[index] += weight * (math::transform_point(pose_mat, co) - co);
    if constexpr (full_deform) {
      deform[index] += weight * pose_mat.view<3, 3>();
    }
  }

  void finalize(const int index, const float total, float3 &r_position) const
  {
    const float scale_factor = 1.0f / total;
    r_position = r_position + position_delta[index] * scale_factor;
  };

  void finalize(const int index, const float total, float3 &r_position, float3x3 &r_deform_mat)
  {
    BLI_assert(full_deform);
    const float scale_factor = 1.0f / total;
    r_position = r_position + position_delta[index] * scale_factor;
    r_deform_mat = (deform[index] * scale_factor) * r_deform_mat;
  };
};

/**
 * Utility class for accumulating dual quaternion bone deformation.
 * If full_deform is true the deformation matrix is also computed.
 */
template<bool full_deform> struct BoneDeformDualQuaternionMixer {
  Array<DualQuat> dq;

  BoneDeformDualQuaternionMixer(const int size) : dq(size, DualQuat{}) {}

  void accumulate(const bPoseChannel &pchan, const int index, const float3 &co, const float weight)
  {
    const DualQuat &deform_quat = pchan.runtime.deform_dual_quat;

    add_weighted_dq_dq_pivot(&dq[index], &deform_quat, co, weight, full_deform);
  }

  void accumulate_bbone(const bPoseChannel &pchan,
                        const int index,
                        const float3 &co,
                        const float weight,
                        const int bbone_index)
  {
    const Span<DualQuat> quats = {pchan.runtime.bbone_dual_quats,
                                  pchan.runtime.bbone_segments + 1};
    const DualQuat &deform_quat = quats[bbone_index];

    add_weighted_dq_dq_pivot(&dq[index], &deform_quat, co, weight, full_deform);
  }

  void finalize(const int index, const float total, float3 &r_position)
  {
    normalize_dq(&dq[index], total);
    float3x3 local_deform_mat;
    mul_v3m3_dq(r_position, local_deform_mat.ptr(), &dq[index]);
  }

  void finalize(const int index, const float total, float3 &r_position, float3x3 &r_deform_mat)
  {
    BLI_assert(full_deform);
    normalize_dq(&dq[index], total);
    float3x3 local_deform_mat;
    mul_v3m3_dq(r_position, local_deform_mat.ptr(), &dq[index]);
    r_deform_mat = local_deform_mat * r_deform_mat;
  }
};

using MixerVariant = std::variant<BoneDeformLinearMixer<true>,
                                  BoneDeformLinearMixer<false>,
                                  BoneDeformDualQuaternionMixer<true>,
                                  BoneDeformDualQuaternionMixer<false>>;

static void deform_single_group_with_mixer(const bPoseChannel &pchan,
                                           const ArmatureDeformGroup &deform_group,
                                           const Span<float3> positions,
                                           MixerVariant &mixer,
                                           MutableSpan<float> contrib)
{
  BLI_assert(pchan.bone != nullptr);
  const Bone &bone = *pchan.bone;
  BLI_assert(!(bone.flag & BONE_NO_DEFORM));

  constexpr GrainSize grain_size = GrainSize(1024);

  if (bone.segments > 1 && pchan.runtime.bbone_segments == bone.segments) {
    std::visit(
        [&](auto &mixer) {
          deform_group.mask.foreach_index(grain_size, [&](const int index, const int pos) {
            const float3 &position = positions[index];
            /* Weights are stored as compressed array. */
            const float weight = deform_group.weights[pos];

            /* Calculate the indices of the 2 affecting b_bone segments. */
            int bbone_index;
            float bbone_blend;
            BKE_pchan_bbone_deform_segment_index(&pchan, position, &bbone_index, &bbone_blend);
            mixer.accumulate_bbone(
                pchan, index, position, weight * (1.0f - bbone_blend), bbone_index);
            mixer.accumulate_bbone(pchan, index, position, weight * bbone_blend, bbone_index + 1);
            contrib[index] += weight;
          });
        },
        mixer);
  }
  else {
    std::visit(
        [&](auto &mixer) {
          deform_group.mask.foreach_index(grain_size, [&](const int index, const int pos) {
            const float3 &position = positions[index];
            /* Weights are stored as compressed array. */
            const float weight = deform_group.weights[pos];

            mixer.accumulate(pchan, index, position, weight);
            contrib[index] += weight;
          });
        },
        mixer);
  }
}

static void deform_groups_with_mixer(const Span<PoseChannelDeformGroup> pchan_groups,
                                     const Span<float3> positions,
                                     MixerVariant &mixer,
                                     MutableSpan<float> contrib,
                                     IndexMask *r_undeformed_mask,
                                     IndexMaskMemory &memory)
{
  for (const PoseChannelDeformGroup &group : pchan_groups) {
    if (!group.pose_channel) {
      continue;
    }

    deform_single_group_with_mixer(
        *group.pose_channel, group.deform_group, positions, mixer, contrib);

    if (r_undeformed_mask) {
      *r_undeformed_mask = IndexMask::from_difference(
          *r_undeformed_mask, group.deform_group.mask, memory);
    }
  }
}

static void deform_with_mixer(
    const ListBase &pose_channels,
    const blender::IndexMask &selection,
    const std::optional<Span<float>> point_weights,
    const Span<PoseChannelDeformGroup> custom_groups,
    const std::optional<ArmatureDeformVertexGroupParams> &vertex_group_params,
    const bool use_envelope,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    MixerVariant &mixer,
    const MutableSpan<float3> positions,
    const std::optional<MutableSpan<float3x3>> deform_mats)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  BLI_assert(positions.size() >= selection.min_array_size());

  IndexMaskMemory memory;

  Array<float> contrib(selection.min_array_size(), 0.0f);
  IndexMask undeformed_mask = selection;

  deform_groups_with_mixer(custom_groups, positions, mixer, contrib, &undeformed_mask, memory);

  if (vertex_group_params) {
    Vector<PoseChannelDeformGroup> pchan_groups = pose_channel_groups_from_vertex_groups(
        pose_channels,
        undeformed_mask,
        vertex_group_params,
        weight_threshold,
        use_envelope_multiply,
        positions,
        memory);
    deform_groups_with_mixer(pchan_groups, positions, mixer, contrib, &undeformed_mask, memory);
  }

  if (use_envelope) {
    Vector<PoseChannelDeformGroup> pchan_groups = pose_channel_groups_from_envelopes(
        pose_channels, undeformed_mask, weight_threshold, positions, memory);
    /* Note: undeformed_mask is not updated here, no further groups are added. */
    deform_groups_with_mixer(
        pchan_groups, positions, mixer, contrib, nullptr /*&undeformed_mask*/, memory);
  }

  if (deform_mats) {
    std::visit(
        [&](auto &mixer) {
          selection.foreach_index(grain_size, [&](const int index) {
            /* Overall influence. */
            const float mask_weight = point_weights ? (*point_weights)[index] : 1.0f;
            const float weight = contrib[index];
            const float total_weight = weight * mask_weight;
            if (weight_threshold && total_weight <= weight_threshold) {
              return;
            }

            float3 &position = positions[index];
            float3x3 &deform_mat = (*deform_mats)[index];
            mixer.finalize(index, weight, position, deform_mat);
          });
        },
        mixer);
  }
  else {
    std::visit(
        [&](auto &mixer) {
          // TODO can optimize based on full_deform case here.
          selection.foreach_index(grain_size, [&](const int index) {
            /* Overall influence. */
            const float mask_weight = point_weights ? (*point_weights)[index] : 1.0f;
            const float weight = contrib[index];
            const float total_weight = weight * mask_weight;
            if (weight_threshold && total_weight <= weight_threshold) {
              return;
            }

            float3 &position = positions[index];
            mixer.finalize(index, weight, position);
          });
        },
        mixer);
  }
}

/* Determine the correct mixer type. */
static MixerVariant get_deform_mixer(const ArmatureDeformSkinningMode skinning_mode,
                                     const bool has_deform_mats,
                                     const int size)
{
  if (has_deform_mats) {
    switch (skinning_mode) {
      case ArmatureDeformSkinningMode::Linear:
        return BoneDeformLinearMixer<true>(size);
      case ArmatureDeformSkinningMode::DualQuatenrion:
        return BoneDeformDualQuaternionMixer<true>(size);
    }
  }
  else {
    switch (skinning_mode) {
      case ArmatureDeformSkinningMode::Linear:
        return BoneDeformLinearMixer<false>(size);
      case ArmatureDeformSkinningMode::DualQuatenrion:
        return BoneDeformDualQuaternionMixer<false>(size);
    }
  }
  BLI_assert_unreachable();
  return BoneDeformLinearMixer<false>(0);
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

  MixerVariant mixer = get_deform_mixer(skinning_mode, false, selection.min_array_size());
  deform_with_mixer(ob_arm.pose->chanbase,
                    selection,
                    point_weights,
                    custom_groups,
                    vertex_group_params,
                    use_envelope,
                    0.0f,
                    true,
                    mixer,
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

  MixerVariant mixer = get_deform_mixer(skinning_mode, false, selection.min_array_size());
  deform_with_mixer(ob_arm.pose->chanbase,
                    selection,
                    point_weights,
                    custom_groups,
                    vertex_group_params,
                    use_envelope,
                    0.0f,
                    true,
                    mixer,
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
