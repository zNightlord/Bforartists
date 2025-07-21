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

/* XXX Getting vertex groups in field evaluation isn't possible because it relies on special vertex
 * group data arrays (MDeformVert) and the declared list of vertex group names (bDeformGroup) to
 * map to the pose channels of the armature.
 *
 * The only feasible way to get vertex group information is in fields, where the field context can
 * provide the necessary access to both the `dverts` array and the `vertex_group_names` list.
 *
 * Future implementation of sparse attributes should solve these issues. For now this should be
 * considered a proof-of-concept. */

class VertexGroupWeightInput final : public bke::GeometryFieldInput {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;

 public:
  VertexGroupWeightInput(const float4x4 &target_to_armature, const Bone *bone)
      : bke::GeometryFieldInput(CPPType::get<float>(), "Vertex Group Weight Input"),
        target_to_armature_(target_to_armature),
        bone_(bone)
  {
  }

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const
  {
    Array<float> weights(mask.min_array_size());
    devirtualize_varray(positions, [&](const auto positions) {
      mask.foreach_index(GrainSize(4096), [&](const int64_t index) {
        const float3 &pos = positions[index];
        const float3 arm_pos = math::transform_point(target_to_armature_, pos);
        weights[index] = distfactor_to_bone(
            arm_pos, bone_->head, bone_->tail, bone_->rad_head, bone_->rad_tail, bone_->dist);
      });
    });
    return VArray<float>::from_container(std::move(weights));
  }

  uint64_t hash() const
  {
    return get_default_hash(target_to_armature_, def_nr_);
  }

  bool is_equal_to(const fn::FieldNode &other) const
  {
    if (const auto *other_field = dynamic_cast<const VertexGroupWeightInput *>(&other)) {
      return other_field->target_to_armature_ == target_to_armature_ &&
             other_field->def_nr_ == def_nr_;
    }
    return false;
  }

  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
  {
    return AttrDomain::Point;
  }
};

class VertexGroupSelectionInput final : public bke::GeometryFieldInput {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  float threshold_;

 public:
  VertexGroupSelectionInput(const float4x4 &target_to_armature,
                            const Bone *bone,
                            const float threshold)
      : bke::GeometryFieldInput(CPPType::get<bool>(), "Vertex Group Selection Input"),
        target_to_armature_(target_to_armature),
        bone_(bone),
        threshold_(threshold)
  {
  }

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const
  {
    fn::FieldEvaluator evaluator{context, mask.min_array_size()};
    evaluator.add(position_field_);
    evaluator.evaluate();
    const VArray<float3> positions = evaluator.get_evaluated<float3>(0);

    Array<bool> selection(mask.min_array_size());
    devirtualize_varray(positions, [&](const auto positions) {
      mask.foreach_index(GrainSize(4096), [&](const int64_t index) {
        const float3 &pos = positions[index];
        const float3 arm_pos = math::transform_point(target_to_armature_, pos);
        /* TODO optimize with simplified "distfactor" function for threshold selection. */
        const float weight = distfactor_to_bone(
            arm_pos, bone_->head, bone_->tail, bone_->rad_head, bone_->rad_tail, bone_->dist);
        selection[index] = weight > threshold_;
      });
    });
    return VArray<bool>::from_container(std::move(selection));
  }

  uint64_t hash() const
  {
    return get_default_hash(target_to_armature_, def_nr_, threshold_);
  }

  bool is_equal_to(const fn::FieldNode &other) const
  {
    if (const auto *other_field = dynamic_cast<const VertexGroupSelectionInput *>(&other)) {
      return other_field->target_to_armature_ == target_to_armature_ &&
             other_field->def_nr_ == def_nr_ && other_field->threshold_ == threshold_;
    }
    return false;
  }

  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
  {
    return AttrDomain::Point;
  }
};

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
    mf::SignatureBuilder builder{"Bone Envelope", signature_};
    builder.single_input<float3>("Position");
    builder.single_output<float>("Weight");
    builder.single_output<bool>("Selection");
    this->set_signature(&signature_);
  }

  void call(const IndexMask &mask, mf::Params params, mf::Context /*context*/) const override
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
};

}  // namespace blender::bke
