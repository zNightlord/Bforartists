/* SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
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

namespace blender::bke {
//
// class VertexGroupWeightInput final : public bke::GeometryFieldInput {
// private:
//  float4x4 target_to_armature_;
//  int def_nr_;
//
// public:
//  VertexGroupWeightInput(const float4x4 &target_to_armature, const int def_nr);
//  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
//                                 const IndexMask &mask) const final;
//  uint64_t hash() const override;
//  bool is_equal_to(const fn::FieldNode &other) const override;
//  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
//  final;
//};
//
// class VertexGroupSelectionInput final : public bke::GeometryFieldInput {
// private:
//  float4x4 target_to_armature_;
//  int def_nr_;
//  float threshold_;
//
// public:
//  VertexGroupSelectionInput(const float4x4 &target_to_armature,
//                            const int def_nr,
//                            const float threshold);
//  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
//                                 const IndexMask &mask) const final;
//  uint64_t hash() const override;
//  bool is_equal_to(const fn::FieldNode &other) const override;
//  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
//  final;
//};
//
// class BoneEnvelopeWeightInput final : public bke::GeometryFieldInput {
// private:
//  float4x4 target_to_armature_;
//  const Bone *bone_;
//  fn::Field<float3> position_field_;
//
// public:
//  BoneEnvelopeWeightInput(const float4x4 &target_to_armature,
//                          const Bone &bone,
//                          fn::Field<float3> position_field);
//  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
//                                 const IndexMask &mask) const final;
//  uint64_t hash() const override;
//  bool is_equal_to(const fn::FieldNode &other) const override;
//  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
//  final;
//};
//
// class BoneEnvelopeSelectionInput final : public bke::GeometryFieldInput {
// private:
//  float4x4 target_to_armature_;
//  const Bone *bone_;
//  fn::Field<float3> position_field_;
//  float threshold_;
//
// public:
//  BoneEnvelopeSelectionInput(const float4x4 &target_to_armature,
//                             const Bone &bone,
//                             fn::Field<float3> position_field,
//                             const float threshold);
//  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
//                                 const IndexMask &mask) const final;
//  uint64_t hash() const override;
//  bool is_equal_to(const fn::FieldNode &other) const override;
//  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
//  final;
//};

class BoneEnvelopeMultiFunction : public mf::MultiFunction {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  float threshold_;

 public:
  BoneEnvelopeMultiFunction(const float4x4 &target_to_armature,
                            const Bone &bone,
                            const float threshold)
      : target_to_armature_(target_to_armature), bone_(&bone), threshold_(threshold)
  {
    static mf::Signature signature_;
    mf::SignatureBuilder builder{"Quaternion to Axis Angle", signature_};
    builder.single_input<math::Quaternion>("Quaternion");
    builder.single_output<float3>("Axis");
    builder.single_output<float>("Angle");
    this->set_signature(&signature_);
  }

  void call(const IndexMask &mask, mf::Params params, mf::Context /*context*/) const override
  {
    const VArraySpan<math::Quaternion> quaternions =
        params.readonly_single_input<math::Quaternion>(0, "Quaternion");
    MutableSpan<float3> axes = params.uninitialized_single_output<float3>(1, "Axis");
    MutableSpan<float> angles = params.uninitialized_single_output<float>(2, "Angle");
    mask.foreach_index([&](const int64_t i) {
      const math::AxisAngle axis_angle = math::to_axis_angle(quaternions[i]);
      axes[i] = axis_angle.axis();
      angles[i] = axis_angle.angle().radian();
    });
  }
};

}  // namespace blender::bke
