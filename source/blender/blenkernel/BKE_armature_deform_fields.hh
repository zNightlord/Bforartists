/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bke
 */

#include "BLI_math_matrix_types.hh"

#include "BKE_geometry_fields.hh"

struct bArmature;
struct Bone;
struct bPoseChannel;

namespace blender::bke::armature_deform {

/* XXX Getting vertex groups in field evaluation isn't possible because it relies on special vertex
 * group data arrays (MDeformVert) and the declared list of vertex group names (bDeformGroup) to
 * map to the pose channels of the armature.
 *
 * The only feasible way to get vertex group information is in fields, where the field context can
 * provide the necessary access to both the `dverts` array and the `vertex_group_names` list.
 *
 * Future implementation of sparse attributes should solve these issues. For now this should be
 * considered a proof-of-concept. */

/**
 * Outputs weights for a vertex group associated with the bone.
 * The vertex group data is based on the geometry context.
 */
class VertexGroupWeightInput final : public bke::GeometryFieldInput {
 private:
  const Bone *bone_;

 public:
  VertexGroupWeightInput(const Bone &bone);

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const;
  uint64_t hash() const;
  bool is_equal_to(const fn::FieldNode &other) const;
  std::optional<AttrDomain> preferred_domain(const GeometryComponent &component) const;
};

/**
 * Outputs a selection field for a vertex group associated with the bone.
 * The vertex group data is based on the geometry context.
 */
class VertexGroupSelectionInput final : public bke::GeometryFieldInput {
 private:
  const Bone *bone_;
  float threshold_;

 public:
  VertexGroupSelectionInput(const Bone &bone, const float threshold);

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const;
  uint64_t hash() const;
  bool is_equal_to(const fn::FieldNode &other) const;
  std::optional<AttrDomain> preferred_domain(const GeometryComponent &component) const;
};

/**
 * Computes deformation weight of the bone envelope falloff function based on position,
 * as well as a boolean value for selecting points within the falloff range.
 */
class BoneEnvelopeMultiFunction : public mf::MultiFunction {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  float threshold_;

 public:
  BoneEnvelopeMultiFunction(const float4x4 &target_to_armature,
                            const Bone &bone,
                            const float threshold);

  void call(const IndexMask &mask, mf::Params params, mf::Context context) const override;
};

/**
 * Group of deformation weights defined on a set of points.
 */
struct ArmatureDeformGroup {
  /**
   * Point indices included in this group.
   */
  IndexMask mask;
  /**
   * Compressed weights, matches positions in the index mask.
   */
  Array<float> weights;
};

/**
 * Armature deform group associated with a pose channel.
 */
struct PoseChannelDeformGroup {
  ArmatureDeformGroup deform_group;
  const bPoseChannel *pose_channel;
};

/**
 * The skinning mode used by armature deformation functions.
 */
enum class SkinningMode {
  /* Linear interpolation of the bone transformation. */
  Linear,
  /* Dual-quaternion interpolation, aka. "Preserve Volume" method. */
  DualQuaternion,
};

/**
 * Deform position vectors by blending one or more groups.
 *
 * \param target_to_armature Transforms points into the armature pose space.
 * \param selection Set of points to deform.
 * \param point_weights Optional global influence factor for each point.
 * \param deform_groups One or more groups whose influence is normalized for each point.
 * \param skinning_mode Skinning mode to use for the deformation.
 * \param positions Position data to deform.
 */
void deform_positions(const float4x4 &target_to_armature,
                      const IndexMask &selection,
                      std::optional<Span<float>> point_weights,
                      Span<PoseChannelDeformGroup> deform_groups,
                      const SkinningMode skinning_mode,
                      MutableSpan<float3> positions);

/**
 * Deform matrices by blending one or more groups.
 *
 * \param target_to_armature Transforms points into the armature pose space.
 * \param selection Set of points to deform.
 * \param point_weights Optional global influence factor for each point.
 * \param deform_groups One or more groups whose influence is normalized for each point.
 * \param skinning_mode Skinning mode to use for the deformation.
 * \param matrices Matrix data to deform.
 */
void deform_matrices(const float4x4 &target_to_armature,
                     const IndexMask &selection,
                     std::optional<Span<float>> point_weights,
                     Span<PoseChannelDeformGroup> deform_groups,
                     const SkinningMode skinning_mode,
                     MutableSpan<float4x4> matrices);

}  // namespace blender::bke::armature_deform
