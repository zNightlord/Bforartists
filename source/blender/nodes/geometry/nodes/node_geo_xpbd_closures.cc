/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_attribute.hh"
#include "BKE_geometry_set.hh"

#include "NOD_xpbd_constraints.hh"

#include "node_geometry_util.hh"

#include <fmt/format.h>

namespace blender::nodes::xpbd_constraints {

template<typename T>
static AttributeReader<T> lookup_or_warn(AttributeAccessor &attributes,
                                         const StringRef attribute_id,
                                         const AttrDomain domain,
                                         const T &default_value,
                                         ConstraintEvalParams::ErrorFn error_fn)
{
  if (!attributes.contains(attribute_id)) {
    error_fn(fmt::format("Missing \"{}\" attribute", attribute_id));
  }
  return attributes.lookup_or_default<T>(attribute_id, domain, default_value);
}

template<typename T>
static AttributeReader<T> lookup_or_warn(const GeoNodeExecParams &params,
                                         AttributeAccessor &attributes,
                                         const StringRef attribute_id,
                                         const AttrDomain domain,
                                         const T &default_value)
{
  return lookup_or_warn(params.error_message_add, attributes, attribute_id, domain, default_value);
}

/* Constraint attributes. */
constexpr StringRef ATTR_SOLVER_GROUP = "solver_group";
constexpr StringRef ATTR_ALPHA = "compliance";
constexpr StringRef ATTR_BETA = "damping";
constexpr StringRef ATTR_POINT1 = "point1";
constexpr StringRef ATTR_POINT2 = "point2";
constexpr StringRef ATTR_ACTIVE = "active";

ConstraintClosure::ConstraintClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
    : geometry_set(std::move(geometry_set))
{
  BLI_assert(this->geometry_set.has_component<PointCloudComponent>());
  const PointCloudComponent &component = *this->geometry_set.get_component<PointCloudComponent>();
  std::optional<bke::AttributeAccessor> attributes = component.attributes();
  this->solver_groups = *lookup_or_warn<int>(
      *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0, error_fn);
}

struct PositionGoalClosure : public ConstraintClosure {
  /* Include damping in the positional constraint update. */
  static constexpr bool use_damping = true;

  VArraySpan<int> points;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float3> goal_positions;

  bke::SpanAttributeWriter<float> position_lambda;

  PositionGoalClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
    this->points = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->goal_positions = *lookup_or_warn<float3>(
        *attributes, "goal_position", AttrDomain::Point, float3(0.0f), error_fn);

    this->position_lambda = attributes->lookup_or_add_for_write_span<float>(
        "position_lambda", AttrDomain::Point, bke::AttributeInitDefaultValue());
  }

  template<bool debug_check>
  void do_apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        params.positions.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      float &lambda = this->position_lambda.span[index];
      const float3 &goal = this->goal_positions[index];

      position_checker.claim_variable(point1);

      float3 &position1 = params.positions[point1];
      if constexpr (use_damping) {
        const float3 &old_position1 = params.old_positions[point1];
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        const float gamma = this->alpha[index] * this->beta[index] * params.inv_delta_time;
        xpbd_constraints::apply_position_goal(
            goal, alpha, gamma, old_position1, lambda, position1);
      }
      else {
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        xpbd_constraints::apply_position_goal(goal, alpha, lambda, position1);
      }
    });

    if (position_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, group_mask);
    }
  }

  void apply_to_velocities(ConstraintEvalParams & /*params*/,
                           const IndexMask & /*group_mask*/) override
  {
  }

  void reset_lambda() override
  {
    for (const int index : this->position_lambda.span.index_range()) {
      this->position_lambda.span[index] = 0.0f;
    }
  }

  void finish_attributes() override
  {
    this->position_lambda.finish();
  }
};

struct RotationGoalClosure : public ConstraintClosure {
  /* Include damping in the positional constraint update. */
  static constexpr bool use_damping = true;

  VArraySpan<int> points;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<math::Quaternion> goal_rotations;

  bke::SpanAttributeWriter<float> position_lambda_w;
  bke::SpanAttributeWriter<float3> position_lambda_xyz;

  RotationGoalClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
    this->points = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->goal_rotations = *lookup_or_warn<math::Quaternion>(
        *attributes, "goal_rotation", AttrDomain::Point, math::Quaternion::identity(), error_fn);

    this->position_lambda_w = attributes->lookup_or_add_for_write_span<float>(
        "position_lambda_w", AttrDomain::Point, bke::AttributeInitDefaultValue());
    this->position_lambda_xyz = attributes->lookup_or_add_for_write_span<float3>(
        "position_lambda_xyz", AttrDomain::Point, bke::AttributeInitDefaultValue());
  }

  template<bool debug_check>
  void do_apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      float &lambda_w = this->position_lambda_w.span[index];
      float3 &lambda_xyz = this->position_lambda_xyz.span[index];
      const math::Quaternion &goal = this->goal_rotations[index];

      rotation_checker.claim_variable(point1);

      math::Quaternion &rotation1 = params.rotations[point1];
      float4 lambda = float4(lambda_w, lambda_xyz);
      if constexpr (use_damping) {
        const math::Quaternion &old_rotation1 = params.old_rotations[point1];
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        const float gamma = this->alpha[index] * this->beta[index] * params.inv_delta_time;
        xpbd_constraints::apply_rotation_goal2<true>(
            goal, alpha, gamma, old_rotation1, lambda, rotation1);
      }
      else {
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        xpbd_constraints::apply_rotation_goal2<true>(goal, alpha, lambda, rotation1);
      }
      lambda_w = lambda.x;
      lambda_xyz = lambda.yzw();
    });

    if (rotation_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, group_mask);
    }
  }

  void apply_to_velocities(ConstraintEvalParams & /*params*/,
                           const IndexMask & /*group_mask*/) override
  {
  }

  void reset_lambda() override
  {
    for (const int index : this->position_lambda_w.span.index_range()) {
      this->position_lambda_w.span[index] = 0.0f;
      this->position_lambda_xyz.span[index] = float3(0.0f);
    }
  }

  void finish_attributes() override
  {
    this->position_lambda_w.finish();
    this->position_lambda_xyz.finish();
  }
};

struct StretchShearClosure : public ConstraintClosure {
  /* Include damping in the positional constraint update. */
  static constexpr bool use_damping = true;

  VArraySpan<int> points1;
  VArraySpan<int> points2;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float> edge_lengths;

  bke::SpanAttributeWriter<float3> position_lambda;

  StretchShearClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
    this->points1 = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->points2 = *lookup_or_warn<int>(*attributes, ATTR_POINT2, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->edge_lengths = *lookup_or_warn<float>(
        *attributes, "edge_length", AttrDomain::Point, 0.0f, error_fn);

    this->position_lambda = attributes->lookup_or_add_for_write_span<float3>(
        "position_lambda", AttrDomain::Point, bke::AttributeInitDefaultValue());
  }

  template<bool debug_check>
  void do_apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        params.positions.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      float3 &lambda = this->position_lambda.span[index];
      const float edge_length = this->edge_lengths[index];

      position_checker.claim_variable(point1);
      position_checker.claim_variable(point2);
      rotation_checker.claim_variable(point1);

      float3 &position1 = params.positions[point1];
      float3 &position2 = params.positions[point2];
      math::Quaternion &rotation = params.rotations[point1];
      const float weight_pos1 = params.position_weights[point1];
      const float weight_pos2 = params.position_weights[point2];
      const float weight_rot = params.rotation_weights[point1];
      if constexpr (use_damping) {
        const float3 old_position1 = params.old_positions[point1];
        const float3 old_position2 = params.old_positions[point2];
        const math::Quaternion &old_rotation = params.old_rotations[point1];
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        const float gamma = this->alpha[index] * this->beta[index] * params.inv_delta_time;
        xpbd_constraints::apply_position_stretch_shear<true>(weight_pos1,
                                                             weight_pos2,
                                                             weight_rot,
                                                             edge_length,
                                                             alpha,
                                                             gamma,
                                                             old_position1,
                                                             old_position2,
                                                             old_rotation,
                                                             lambda,
                                                             position1,
                                                             position2,
                                                             rotation);
      }
      else {
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        xpbd_constraints::apply_position_stretch_shear<true>(weight_pos1,
                                                             weight_pos2,
                                                             weight_rot,
                                                             edge_length,
                                                             alpha,
                                                             lambda,
                                                             position1,
                                                             position2,
                                                             rotation);
      }
    });

    if (position_checker.has_overlap() || rotation_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, group_mask);
    }
  }

  void apply_to_velocities(ConstraintEvalParams & /*params*/,
                           const IndexMask & /*group_mask*/) override
  {
  }

  void reset_lambda() override
  {
    for (const int index : this->position_lambda.span.index_range()) {
      this->position_lambda.span[index] = float3(0.0f);
    }
  }

  void finish_attributes() override
  {
    this->position_lambda.finish();
  }
};

struct BendTwistClosure : public ConstraintClosure {
  /* Include damping in the positional constraint update. */
  static constexpr bool use_damping = true;

  VArraySpan<int> points1;
  VArraySpan<int> points2;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float> edge_lengths;
  VArraySpan<float> darboux_w;
  VArraySpan<float3> darboux_xyz;

  bke::SpanAttributeWriter<float> position_lambda_w;
  bke::SpanAttributeWriter<float3> position_lambda_xyz;

  BendTwistClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
    this->points1 = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->points2 = *lookup_or_warn<int>(*attributes, ATTR_POINT2, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->edge_lengths = *lookup_or_warn<float>(
        *attributes, "edge_length", AttrDomain::Point, 0.0f, error_fn);
    this->darboux_w = *lookup_or_warn<float>(
        *attributes, "darboux_w", AttrDomain::Point, float(0.0f), error_fn);
    this->darboux_xyz = *lookup_or_warn<float3>(
        *attributes, "darboux_xyz", AttrDomain::Point, float3(0.0f), error_fn);

    /* XXX plain float4 attribute is not supported, have to store it as float + float3. */
    this->position_lambda_w = attributes->lookup_or_add_for_write_span<float>(
        "position_lambda_w", AttrDomain::Point, bke::AttributeInitDefaultValue());
    this->position_lambda_xyz = attributes->lookup_or_add_for_write_span<float3>(
        "position_lambda_xyz", AttrDomain::Point, bke::AttributeInitDefaultValue());
  }

  template<bool debug_check>
  void do_apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      float &lambda_w = this->position_lambda_w.span[index];
      float3 &lambda_xyz = this->position_lambda_xyz.span[index];
      const float edge_length = this->edge_lengths[index];
      const math::Quaternion darboux_vector = math::Quaternion(this->darboux_w[index],
                                                               this->darboux_xyz[index]);

      rotation_checker.claim_variable(point1);
      rotation_checker.claim_variable(point2);

      math::Quaternion &rotation1 = params.rotations[point1];
      math::Quaternion &rotation2 = params.rotations[point2];
      const float weight_rot1 = params.rotation_weights[point1];
      const float weight_rot2 = params.rotation_weights[point2];

      float4 lambda = float4(lambda_w, lambda_xyz);
      if constexpr (use_damping) {
        const math::Quaternion &old_rotation1 = params.old_rotations[point1];
        const math::Quaternion &old_rotation2 = params.old_rotations[point2];
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        const float gamma = this->alpha[index] * this->beta[index] * params.inv_delta_time;
        xpbd_constraints::apply_position_bend_twist<true>(weight_rot1,
                                                          weight_rot2,
                                                          edge_length,
                                                          darboux_vector,
                                                          alpha,
                                                          gamma,
                                                          old_rotation1,
                                                          old_rotation2,
                                                          lambda,
                                                          rotation1,
                                                          rotation2);
      }
      else {
        const float alpha = this->alpha[index] * params.inv_delta_time_squared;
        xpbd_constraints::apply_position_bend_twist<true>(weight_rot1,
                                                          weight_rot2,
                                                          edge_length,
                                                          darboux_vector,
                                                          alpha,
                                                          lambda,
                                                          rotation1,
                                                          rotation2);
      }
      lambda_w = lambda[0];
      lambda_xyz = float3(lambda[1], lambda[2], lambda[3]);
    });

    if (rotation_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, group_mask);
    }
  }

  void apply_to_velocities(ConstraintEvalParams & /*params*/,
                           const IndexMask & /*group_mask*/) override
  {
  }

  void reset_lambda() override
  {
    for (const int index : this->position_lambda_w.span.index_range()) {
      this->position_lambda_w.span[index] = 0.0f;
      this->position_lambda_xyz.span[index] = float3(0.0f);
    }
  }

  void finish_attributes() override
  {
    this->position_lambda_w.finish();
    this->position_lambda_xyz.finish();
  }
};

struct ContactClosure : public ConstraintClosure {
  VArraySpan<int> points;
  VArraySpan<int> collider_index;
  VArraySpan<float3> local_position1;
  VArraySpan<float3> local_position2;
  VArraySpan<float3> normal;

  VArraySpan<float> friction;
  VArraySpan<float> restitution;

  bke::SpanAttributeWriter<float> position_lambda;
  bke::SpanAttributeWriter<float> restitution_lambda;
  bke::SpanAttributeWriter<float> friction_lambda;

  /* Remember active contacts for later velocity update. */
  bke::SpanAttributeWriter<bool> active;

  ContactClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
    this->points = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->collider_index = *lookup_or_warn<int>(
        *attributes, "collider_index", AttrDomain::Point, 0, error_fn);
    this->local_position1 = *lookup_or_warn<float3>(
        *attributes, "local_position1", AttrDomain::Point, float3(0.0f), error_fn);
    this->local_position2 = *lookup_or_warn<float3>(
        *attributes, "local_position2", AttrDomain::Point, float3(0.0f), error_fn);
    this->normal = *lookup_or_warn<float3>(
        *attributes, "normal", AttrDomain::Point, float3(0.0f), error_fn);

    this->friction = *attributes->lookup_or_default<float>("friction", AttrDomain::Point, 0.0f);
    this->restitution = *attributes->lookup_or_default<float>(
        "restitution", AttrDomain::Point, 0.0f);

    const int num_constraints = attributes->domain_size(AttrDomain::Point);
    this->active = attributes->lookup_or_add_for_write_span<bool>(
        ATTR_ACTIVE,
        AttrDomain::Point,
        bke::AttributeInitVArray(VArray<bool>::ForSingle(false, num_constraints)));
    this->position_lambda = attributes->lookup_or_add_for_write_span<float>(
        "position_lambda", AttrDomain::Point, bke::AttributeInitDefaultValue());
    this->restitution_lambda = attributes->lookup_or_add_for_write_span<float>(
        "restitution_lambda", AttrDomain::Point, bke::AttributeInitDefaultValue());
    this->friction_lambda = attributes->lookup_or_add_for_write_span<float>(
        "friction_lambda", AttrDomain::Point, bke::AttributeInitDefaultValue());
  }

  template<bool debug_check>
  void do_apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        params.positions.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        params.rotations.index_range());
    [[maybe_unused]] std::atomic_bool error_unit_contact_normal = false;

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point = this->points[index];
      const int collider_index = this->collider_index[index];
      float &lambda = this->position_lambda.span[index];

      position_checker.claim_variable(point);
      rotation_checker.claim_variable(point);

      if (!params.collider_transforms.index_range().contains(collider_index)) {
        return;
      };
      const float4x4 collider_transform = params.collider_transforms[collider_index];
      float3 collider_position;
      math::Quaternion collider_rotation;
      float3 collider_scale;
      math::to_loc_rot_scale(
          collider_transform, collider_position, collider_rotation, collider_scale);

      float3 &position = params.positions[point];
      math::Quaternion &rotation = params.rotations[point];
      bool &active = this->active.span[index];

      const float3 &local_position1 = this->local_position1[index];
      const float3 &local_position2 = this->local_position2[index];
      const float3 &normal = this->normal[index];
      if constexpr (debug_check) {
        if (!math::is_unit(normal)) {
          error_unit_contact_normal.store(true, std::memory_order_relaxed);
        }
      }

      /* Contact constraints are stiff. */
      const float alpha = 0.0f * params.inv_delta_time_squared;

      /* Zero weights for the collider, only the point can move. */
      active = xpbd_constraints::apply_position_contact(1.0f,
                                                        0.0f,
                                                        1.0f,
                                                        0.0f,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        alpha,
                                                        lambda,
                                                        position,
                                                        collider_position,
                                                        rotation,
                                                        collider_rotation);
    });

    if (position_checker.has_overlap() || rotation_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
    if (error_unit_contact_normal) {
      params.error_message_add("Contact normal vector not normalized");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(ConstraintEvalParams &params, const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> velocity_checker(
        params.velocities.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> angular_velocity_checker(
        params.angular_velocities.index_range());
    [[maybe_unused]] std::atomic_bool error_unit_contact_normal = false;

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      /* Active status is determined by the position evaluation. */
      if (!this->active.span[index]) {
        return;
      }

      const int point = this->points[index];

      velocity_checker.claim_variable(point);
      angular_velocity_checker.claim_variable(point);

      const int collider_index = this->collider_index[index];
      if (!params.collider_transforms.index_range().contains(collider_index)) {
        return;
      };

      float &lambda_restitution = this->restitution_lambda.span[index];
      float &lambda_friction = this->friction_lambda.span[index];

      /* Local positions are relative to moving point and collider respectively. */
      const float3 &local_position1 = this->local_position1[index];
      const float3 &local_position2 = this->local_position2[index];
      /* Normal is a fixed shared direction for both participants. */
      const float3 &normal = this->normal[index];
      if constexpr (debug_check) {
        if (!math::is_unit(normal)) {
          error_unit_contact_normal.store(true, std::memory_order_relaxed);
        }
      }
      const float restitution = this->restitution[index];
      const float friction = this->friction[index];

      const float4x4 &collider_transform = params.collider_transforms[collider_index];
      const float4x4 &old_collider_transform = params.old_collider_transforms[collider_index];

      float3 &velocity = params.velocities[point];
      float3 &angular_velocity = params.angular_velocities[point];
      const float3 &orig_velocity = params.orig_velocities[point];
      const float3 &orig_angular_velocity = params.orig_angular_velocities[point];

      /* Compute velocity from old/new collider transforms. */
      float3 collider_loc, old_collider_loc;
      math::Quaternion collider_rot, old_collider_rot;
      float3 collider_scale, old_collider_scale;
      math::to_loc_rot_scale(collider_transform, collider_loc, collider_rot, collider_scale);
      math::to_loc_rot_scale(
          old_collider_transform, old_collider_loc, old_collider_rot, old_collider_scale);
      float3 collider_velocity = (collider_loc - old_collider_loc) * params.inv_delta_time;
      float3 collider_angular_velocity =
          2.0f * (math::invert_normalized(old_collider_rot) * collider_rot).imaginary_part() *
          params.inv_delta_time;
      /* No change in animated collider velocity. */
      const float3 orig_collider_velocity = collider_velocity;
      const float3 orig_collider_angular_velocity = collider_angular_velocity;

      xpbd_constraints::apply_velocity_contact(orig_velocity,
                                               orig_collider_velocity,
                                               orig_angular_velocity,
                                               orig_collider_angular_velocity,
                                               local_position1,
                                               local_position2,
                                               normal,
                                               restitution,
                                               friction,
                                               lambda_restitution,
                                               lambda_friction,
                                               velocity,
                                               collider_velocity,
                                               angular_velocity,
                                               collider_angular_velocity);
    });

    if (velocity_checker.has_overlap() || angular_velocity_checker.has_overlap()) {
      params.error_message_add("Overlapping constraint solver groups");
    }
    if (error_unit_contact_normal) {
      params.error_message_add("Contact normal vector not normalized");
    }
  }

  void apply_to_positions(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, group_mask);
    }
  }

  void apply_to_velocities(ConstraintEvalParams &params, const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, group_mask);
    }
  }

  void reset_lambda() override
  {
    for (const int index : this->position_lambda.span.index_range()) {
      this->position_lambda.span[index] = 0.0f;
      this->restitution_lambda.span[index] = 0.0f;
      this->friction_lambda.span[index] = 0.0f;
    }
  }
  void finish_attributes() override
  {
    this->position_lambda.finish();
    this->restitution_lambda.finish();
    this->friction_lambda.finish();
    this->active.finish();
  }
};

static Array<ConstraintTypeInfo> create_constraint_info()
{
  using ErrorFn = ConstraintClosure::ErrorFn;

  ConstraintTypeInfo position_goal_info = {
      "Position Goal Constraints",
      "Set position of a point to a target vector",
      [](GeometrySet &&geometry_set, ErrorFn error_fn) -> ConstraintClosure * {
        return new PositionGoalClosure(std::move(geometry_set), error_fn);
      }};
  ConstraintTypeInfo rotation_goal_info = {
      "Rotation Goal Constraints",
      "Set orientation of an edge to a target rotation",
      [](GeometrySet &&geometry_set, ErrorFn error_fn) -> ConstraintClosure * {
        return new RotationGoalClosure(std::move(geometry_set), error_fn);
      }};
  ConstraintTypeInfo stretch_shear_info = {
      "Stretch/Shear Constraints",
      "Enforces edge length and aligns forward direction with the edge vector",
      [](GeometrySet &&geometry_set, ErrorFn error_fn) -> ConstraintClosure * {
        return new StretchShearClosure(std::move(geometry_set), error_fn);
      }};
  ConstraintTypeInfo bend_twist_info = {
      "Bend/Twist Constraints",
      "Enforces angles between neighboring edges to their relative rest orientation",
      [](GeometrySet &&geometry_set, ErrorFn error_fn) -> ConstraintClosure * {
        return new BendTwistClosure(std::move(geometry_set), error_fn);
      }};
  ConstraintTypeInfo contact_info = {
      "Contact Constraints",
      "Keep contact points from penetrating",
      [](GeometrySet &&geometry_set, ErrorFn error_fn) -> ConstraintClosure * {
        return new ContactClosure(std::move(geometry_set), error_fn);
      }};

  /* Order of constraint passes is chosen by increasing "importance":
   * Later constraints have less residual error, and the last constraint type is solved exactly.
   */
  Array<ConstraintTypeInfo> constraint_info = {
      std::move(bend_twist_info),
      std::move(stretch_shear_info),
      std::move(rotation_goal_info),
      std::move(position_goal_info),
      std::move(contact_info),
  };

  return constraint_info;
}

Span<ConstraintTypeInfo> get_constraint_info()
{
  static const Array<ConstraintTypeInfo> constraint_info = create_constraint_info();
  return constraint_info;
}

Span<ConstraintTypeInfo> get_constraint_info_ordered()
{
  /* TODO currently relies on fixed order in get_constraint_info(),
   * could also re-order based on some priority value. */
  return get_constraint_info();
}

}  // namespace blender::nodes::xpbd_constraints
