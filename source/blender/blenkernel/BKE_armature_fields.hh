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

class BoneEnvelopeWeightInput final : public bke::GeometryFieldInput {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  fn::Field<float3> position_field_;

 public:
  BoneEnvelopeWeightInput(const float4x4 &target_to_armature,
                          const Bone &bone,
                          fn::Field<float3> position_field);
  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const final;
  uint64_t hash() const override;
  bool is_equal_to(const fn::FieldNode &other) const override;
  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const final;
};

class BoneEnvelopeSelectionInput final : public bke::GeometryFieldInput {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  fn::Field<float3> position_field_;
  float threshold_;

 public:
  BoneEnvelopeSelectionInput(const float4x4 &target_to_armature,
                             const Bone &bone,
                             fn::Field<float3> position_field,
                             const float threshold);
  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const final;
  uint64_t hash() const override;
  bool is_equal_to(const fn::FieldNode &other) const override;
  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const final;
};

}  // namespace blender::bke
