/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "BLI_index_mask.hh"
#include "BLI_listbase.h"
#include "BLI_listbase_wrapper.hh"
#include "BLI_math_matrix.h"
#include "BLI_math_quaternion.hh"
#include "BLI_math_rotation.h"
#include "BLI_math_vector.h"
#include "BLI_task.hh"

#include "DNA_armature_types.h"
#include "DNA_listBase.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BKE_armature.hh"
#include "BKE_armature_deform_fields.hh"
#include "BKE_curves.hh"
#include "BKE_grease_pencil.hh"
#include "BKE_mesh.hh"

namespace blender::bke::armature_deform {

/* Vertex group information:
 * dverts contain weight information and index the vertex_group list.
 * The vertex_groups list is matched against the pose channels of an armature by name, which maps
 * the deformation weight to a bone pose channel. */
struct VertexGroupData {
  ConstListBaseWrapper<bDeformGroup> vertex_groups;
  Span<MDeformVert> dverts;
};

/* Return vertex group information for supported geometry types (mesh, grease pencil).
 * This data is only available on original object geometry with user-defined vertex groups. */
static std::optional<VertexGroupData> vertex_groups_from_field_context(
    const bke::GeometryFieldContext &context)
{
  switch (context.type()) {
    case GeometryComponent::Type::Mesh:
      if (const Mesh *mesh = context.mesh()) {
        return VertexGroupData{{mesh->vertex_group_names}, mesh->deform_verts()};
      }
      break;
    case GeometryComponent::Type::GreasePencil: {
      const GreasePencil *grease_pencil = context.grease_pencil();
      const bke::greasepencil::Drawing *drawing = context.grease_pencil_layer_drawing();
      if (grease_pencil && drawing) {
        return VertexGroupData{{grease_pencil->vertex_group_names},
                               drawing->geometry.wrap().deform_verts()};
      }
      break;
    }
    default:
      break;
  }
  return std::nullopt;
}

static std::optional<int> find_vertex_group_def_nr(
    const ConstListBaseWrapper<bDeformGroup> &vertex_groups, StringRef name)
{
  int def_nr = 0;
  for (const bDeformGroup *def_grp : vertex_groups) {
    if (def_grp->name == name) {
      return def_nr;
    }
    ++def_nr;
  }
  return std::nullopt;
}

static void weights_from_deform_verts(const IndexMask &mask,
                                      const Span<MDeformVert> dverts,
                                      const int def_nr,
                                      MutableSpan<float> weights)
{
  constexpr GrainSize grain_size(1024);

  BLI_assert(dverts.size() >= mask.min_array_size());

  mask.foreach_index(grain_size, [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    weights[index] = 0.0f;
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr == def_nr) {
        weights[index] = dw.weight;
        break;
      }
    }
  });
}

static void selection_from_deform_verts(const IndexMask &mask,
                                        const Span<MDeformVert> dverts,
                                        const int def_nr,
                                        const float threshold,
                                        MutableSpan<bool> selection)
{
  constexpr GrainSize grain_size(1024);

  BLI_assert(dverts.size() >= mask.min_array_size());

  mask.foreach_index(grain_size, [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    selection[index] = false;
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr == def_nr && dw.weight > threshold) {
        selection[index] = true;
        break;
      }
    }
  });
}

VertexGroupWeightInput::VertexGroupWeightInput(const Bone &bone)
    : bke::GeometryFieldInput(CPPType::get<float>(), "Vertex Group Weight Input"), bone_(&bone)
{
}

GVArray VertexGroupWeightInput::get_varray_for_context(const bke::GeometryFieldContext &context,
                                                       const IndexMask &mask) const
{
  const auto default_weights = VArray<float>::from_single(0.0f, mask.min_array_size());
  std::optional<VertexGroupData> vgroup_data = vertex_groups_from_field_context(context);
  if (!vgroup_data) {
    return default_weights;
  }
  const std::optional<int> def_nr = find_vertex_group_def_nr(vgroup_data->vertex_groups,
                                                             bone_->name);
  if (!def_nr) {
    return default_weights;
  }

  Array<float> weights(mask.min_array_size());
  weights_from_deform_verts(mask, vgroup_data->dverts, *def_nr, weights);
  return VArray<float>::from_container(std::move(weights));
}

uint64_t VertexGroupWeightInput::hash() const
{
  return get_default_hash(bone_);
}

bool VertexGroupWeightInput::is_equal_to(const fn::FieldNode &other) const
{
  if (const auto *other_field = dynamic_cast<const VertexGroupWeightInput *>(&other)) {
    return other_field->bone_ == bone_;
  }
  return false;
}

std::optional<AttrDomain> VertexGroupWeightInput::preferred_domain(
    const GeometryComponent & /*component*/) const
{
  return AttrDomain::Point;
}

VertexGroupSelectionInput::VertexGroupSelectionInput(const Bone &bone, const float threshold)
    : bke::GeometryFieldInput(CPPType::get<bool>(), "Vertex Group Selection Input"),
      bone_(&bone),
      threshold_(threshold)
{
}

GVArray VertexGroupSelectionInput::get_varray_for_context(const bke::GeometryFieldContext &context,
                                                          const IndexMask &mask) const
{
  const auto default_selection = VArray<bool>::from_single(false, mask.min_array_size());
  std::optional<VertexGroupData> vgroup_data = vertex_groups_from_field_context(context);
  if (!vgroup_data) {
    return default_selection;
  }
  const std::optional<int> def_nr = find_vertex_group_def_nr(vgroup_data->vertex_groups,
                                                             bone_->name);
  if (!def_nr) {
    return default_selection;
  }

  Array<bool> selection(mask.min_array_size());
  selection_from_deform_verts(mask, vgroup_data->dverts, *def_nr, threshold_, selection);
  return VArray<bool>::from_container(std::move(selection));
}

uint64_t VertexGroupSelectionInput::hash() const
{
  return get_default_hash(bone_, threshold_);
}

bool VertexGroupSelectionInput::is_equal_to(const fn::FieldNode &other) const
{
  if (const auto *other_field = dynamic_cast<const VertexGroupSelectionInput *>(&other)) {
    return other_field->bone_ == bone_ && other_field->threshold_ == threshold_;
  }
  return false;
}

std::optional<AttrDomain> VertexGroupSelectionInput::preferred_domain(
    const GeometryComponent & /*component*/) const
{
  return AttrDomain::Point;
}

BoneEnvelopeMultiFunction::BoneEnvelopeMultiFunction(const float4x4 &target_to_armature,
                                                     const Bone &bone,
                                                     const float threshold)
    : target_to_armature_(target_to_armature), bone_(&bone), threshold_(threshold)
{
  static const mf::Signature signature = []() {
    mf::Signature signature;
    mf::SignatureBuilder builder{"Bone Envelope", signature};
    builder.single_input<float3>("Position");
    builder.single_output<float>("Weight");
    builder.single_output<bool>("Selection");
    return signature;
  }();
  this->set_signature(&signature);
}

void BoneEnvelopeMultiFunction::call(const IndexMask &mask,
                                     mf::Params params,
                                     mf::Context /*context*/) const
{
  const VArraySpan<float3> positions = params.readonly_single_input<float3>(0, "Position");
  MutableSpan<float> weights = params.uninitialized_single_output<float>(1, "Weight");
  MutableSpan<bool> selection = params.uninitialized_single_output<bool>(2, "Selection");
  mask.foreach_index([&](const int64_t i) {
    const float3 &pos = positions[i];
    const float3 arm_pos = math::transform_point(target_to_armature_, pos);
    weights[i] = distfactor_to_bone(
        arm_pos, bone_->head, bone_->tail, bone_->rad_head, bone_->rad_tail, bone_->dist);
    selection[i] = weights[i] > threshold_;
  });
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

static void deform_with_mixer(const IndexMask &selection,
                              const std::optional<Span<float>> point_weights,
                              const Span<PoseChannelDeformGroup> deform_groups,
                              const std::optional<float> weight_threshold,
                              MixerVariant &mixer,
                              const MutableSpan<float3> positions,
                              const std::optional<MutableSpan<float3x3>> deform_mats)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  BLI_assert(positions.size() >= selection.min_array_size());

  Array<float> contrib(selection.min_array_size(), 0.0f);

  for (const PoseChannelDeformGroup &group : deform_groups) {
    if (!group.pose_channel) {
      continue;
    }

    deform_single_group_with_mixer(
        *group.pose_channel, group.deform_group, positions, mixer, contrib);
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
static MixerVariant get_deform_mixer(const SkinningMode skinning_mode,
                                     const bool has_deform_mats,
                                     const int size)
{
  if (has_deform_mats) {
    switch (skinning_mode) {
      case SkinningMode::Linear:
        return BoneDeformLinearMixer<true>(size);
      case SkinningMode::DualQuaternion:
        return BoneDeformDualQuaternionMixer<true>(size);
    }
  }
  else {
    switch (skinning_mode) {
      case SkinningMode::Linear:
        return BoneDeformLinearMixer<false>(size);
      case SkinningMode::DualQuaternion:
        return BoneDeformDualQuaternionMixer<false>(size);
    }
  }
  BLI_assert_unreachable();
  return BoneDeformLinearMixer<false>(0);
}

void deform_positions(const float4x4 &target_to_armature,
                      const IndexMask &selection,
                      const std::optional<Span<float>> point_weights,
                      const Span<PoseChannelDeformGroup> deform_groups,
                      const SkinningMode skinning_mode,
                      MutableSpan<float3> positions)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  const float4x4 armature_to_target = math::invert(target_to_armature);

  /* Input coordinates to start from. */
  Array<float3> armature_positions(positions.size());
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    armature_positions[index] = math::transform_point(target_to_armature, positions[index]);
  });

  MixerVariant mixer = get_deform_mixer(skinning_mode, false, selection.min_array_size());
  deform_with_mixer(
      selection, point_weights, deform_groups, 0.0f, mixer, armature_positions, std::nullopt);

  /* Transform back to target object space. */
  selection.foreach_index(grain_size, [&](const int index) {
    positions[index] = math::transform_point(armature_to_target, armature_positions[index]);
  });
}

void deform_matrices(const float4x4 &target_to_armature,
                     const IndexMask &selection,
                     const std::optional<Span<float>> point_weights,
                     const Span<PoseChannelDeformGroup> deform_groups,
                     const SkinningMode skinning_mode,
                     MutableSpan<float4x4> matrices)
{
  constexpr GrainSize grain_size = GrainSize(1024);

  const float4x4 armature_to_target = math::invert(target_to_armature);

  /* Input coordinates to start from. */
  Array<float3> armature_positions(matrices.size());
  Array<float3x3> deform_mats(matrices.size());
  /* Transform to armature space. */
  selection.foreach_index(grain_size, [&](const int index) {
    const float4x4 armature_matrix = target_to_armature * matrices[index];
    armature_positions[index] = armature_matrix.location();
    deform_mats[index] = armature_matrix.view<3, 3>();
  });

  MixerVariant mixer = get_deform_mixer(skinning_mode, true, selection.min_array_size());
  deform_with_mixer(
      selection, point_weights, deform_groups, 0.0f, mixer, armature_positions, deform_mats);

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

}  // namespace blender::bke::armature_deform

/** \} */
