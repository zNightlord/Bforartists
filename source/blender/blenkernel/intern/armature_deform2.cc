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

#include "CLG_log.h"

static CLG_LogRef LOG = {"geom.armature_deform"};

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
      if (!dw.def_nr != def_nr) {
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

  Vector<int> indices;
  indices.reserve(def_nr_count);
  Array<float> weights(def_nr_count);

  universe.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (!dw.def_nr == def_nr) {
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

template<typename MixerT>
static void deform_with_mixer(
    const ListBase &pose_channels,
    const blender::IndexMask &selection,
    const std::optional<ArmatureDeformVertexGroupParams> &vertex_group_params,
    const bool use_envelope,
    const std::optional<float> weight_threshold,
    const bool use_envelope_multiply,
    const Span<float3> positions)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  BLI_assert(positions.size() >= selection.min_array_size());

  IndexMaskMemory memory;

  Array<MixerT> mixers(selection.min_array_size(), {});
  Array<float> contrib(selection.min_array_size(), 0.0f);
  Array<bool> deformed(selection.min_array_size(), false);

  if (vertex_group_params) {
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
          selection,
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
    /* Apply envelope deformation on undeformed points. */
    const IndexMask undeformed = IndexMask::from_batch_predicate(
        selection,
        grain_size,
        memory,
        [&](const IndexMaskSegment universe_segment, IndexRangesBuilder<int16_t> &builder) {
          for (const int16_t local_index : universe_segment.base_span()) {
            if (!deformed[universe_segment.offset() + local_index]) {
              builder.add(local_index);
            }
            return universe_segment.offset();
          }
        });

    ArmatureDeformGroup deform_group = build_deform_group_for_envelope(
        undeformed, positions, *bone, weight_threshold, memory);
  }
}

void armature_deform_vectors(const Object &ob_arm,
                             const blender::float4x4 &target_to_world,
                             const blender::IndexMask &selection,
                             std::optional<ArmatureDeformVertexGroupParams> vertex_group_params,
                             bool use_envelope,
                             const ArmatureDeformSkinningMode skinning_mode,
                             blender::MutableSpan<blender::float3> vectors)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  const float4x4 target_to_armature = ob_arm.world_to_object() * target_to_world;
  const float4x4 armature_to_target = ob_arm.object_to_world() * math::invert(target_to_world);

  /* Input coordinates to start from. */
  Array<float3> positions = vectors.as_span();
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    float3 &position = positions[index];
    position = math::transform_point(target_to_armature, position);
  });

  switch (skinning_mode) {
    case ArmatureDeformSkinningMode::Linear:
      deform_with_mixer<BoneDeformLinearMixer<false>>(ob_arm.pose->chanbase,
                                                      selection,
                                                      vertex_group_params,
                                                      use_envelope,
                                                      0.0f,
                                                      true,
                                                      positions);
      break;
    case ArmatureDeformSkinningMode::DualQuatenrion:
      deform_with_mixer<BoneDeformDualQuaternionMixer<false>>(ob_arm.pose->chanbase,
                                                              selection,
                                                              vertex_group_params,
                                                              use_envelope,
                                                              0.0f,
                                                              true,
                                                              positions);
      break;
  }

  // Array<MixerT> mixers(positions.size(), MixerT{});
  // Array<float> contrib(positions.size(), 0.0f);
  // Array<bool> deformed(positions.size(), false);
  // /* Apply vertex group deformation if enabled. */
  // if (params.use_dverts && dverts) {
  //   /* Range of valid def_nr in MDeformWeight. */
  //   const IndexRange def_nr_range = params.pose_channel_by_vertex_group.index_range();
  //   IndexMaskMemory memory;
  //   /* Build index masks for all vertex groups.
  //    * TODO It may be worthwile to cache this - unless in the future the sparse attribute
  //    storage
  //    * will provide vertex weights in this order anyway, which would make index masks a lot
  //    * cheaper to create. */
  //   const Vector<VertexGroupMask> group_masks = vertex_groups_to_masks(
  //       *dverts, selection, def_nr_range, std::nullopt, true, memory);

  //   for (const VertexGroupMask &group_mask : group_masks) {
  //     BLI_assert(def_nr_range.contains(group_mask.def_nr));
  //     const bPoseChannel *pchan = params.pose_channel_by_vertex_group[group_mask.def_nr];
  //     if (pchan == nullptr) {
  //       continue;
  //     }
  //     /* Bone option to mix with envelope weight. */
  //     const Bone *bone = pchan->bone;
  //     const bool use_envelope_multiply = (bone && bone->flag & BONE_MULT_VG_ENV);

  //     group_mask.mask.foreach_index(grain_size, [&](const int index, const int pos) {
  //       const float3 &position = positions[index];
  //       /* Weights are stored as compressed array. */
  //       float weight = group_mask.weights[pos];
  //       /* Multiply with bone envelope falloff if enabled for this bone. */
  //       if (use_envelope_multiply) {
  //         weight *= distfactor_to_bone(position,
  //                                      float3(bone->arm_head),
  //                                      float3(bone->arm_tail),
  //                                      bone->rad_head,
  //                                      bone->rad_tail,
  //                                      bone->dist);
  //       }

  //       contrib[index] += pchan_bone_deform(*pchan, weight, position, mixers[index]);
  //       deformed[index] = true;
  //     });
  //   }
  // }
  // /* Use envelope if enabled and no bone deformed the vertex yet. */
  // if (params.use_envelope) {
  //   selection.foreach_index(grain_size, [&](const int index) {
  //     if (deformed[index]) {
  //       return;
  //     }

  //     const float3 &position = positions[index];
  //     for (const bPoseChannel *pchan : params.pose_channels) {
  //       if (!(pchan->bone->flag & BONE_NO_DEFORM)) {
  //         contrib[index] += dist_bone_deform(*pchan, position, mixers[index]);
  //       }
  //     }
  //   });
  // }

  // selection.foreach_index(grain_size, [&](const int index) {
  //   float3 &position = positions[index];

  //   /* Overall influence, can change by masking with a vertex group. */
  //   float armature_weight = 1.0f;
  //   float prevco_weight = 0.0f; /* weight for optional cached vertexcos */
  //   if (params.vert_influence) {
  //     const float mask_weight = (*params.vert_influence)[index];
  //     /* On multi-modifier the mask is used to blend with previous coordinates. */
  //     if (params.vert_coords_prev) {
  //       prevco_weight = params.invert_vgroup ? mask_weight : 1.0f - mask_weight;
  //       if (prevco_weight == 1.0f) {
  //         return;
  //       }
  //     }
  //     else {
  //       armature_weight = params.invert_vgroup ? 1.0f - mask_weight : mask_weight;
  //       if (armature_weight == 0.0f) {
  //         return;
  //       }
  //     }
  //   }

  //   /* TODO Actually should be EPSILON? Weight values and contrib can be like 10e-39 small. */
  //   constexpr float contrib_threshold = 0.0001f;
  //   if (contrib[index] > contrib_threshold) {
  //     float3 delta_co;
  //     float3x3 local_deform_mat;
  //     mixers[index].finalize(
  //         position, contrib[index], armature_weight, delta_co, local_deform_mat);

  //     position += delta_co;
  //     if (full_deform) {
  //       float3x3 &deform_mat = (*params.vert_deform_mats)[index];
  //       const float3x3 armature_to_target = params.armature_to_target.view<3, 3>();
  //       const float3x3 target_to_armature = params.target_to_armature.view<3, 3>();
  //       deform_mat = armature_to_target * local_deform_mat * target_to_armature * deform_mat;
  //     }
  //   }

  //   /* Transform back to target object space. */
  //   position = math::transform_point(params.armature_to_target, position);

  //   /* Multi-modifier: Interpolate with previous modifier position using the vertex group mask.
  //    */
  //   if (params.vert_coords_prev) {
  //     copy_v3_v3(params.vert_coords[index],
  //                math::interpolate(position, params.vert_coords[index], prevco_weight));
  //   }
  //   else {
  //     copy_v3_v3(params.vert_coords[index], position);
  //   }
  // });
}

}  // namespace blender::bke

/** \} */
