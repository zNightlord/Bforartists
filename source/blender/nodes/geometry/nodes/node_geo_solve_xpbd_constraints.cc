/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_array_utils.hh"

#include "BKE_attribute.hh"
#include "BKE_geometry_set.hh"
#include "BKE_instances.hh"

#include "NOD_xpbd_constraints.hh"

#include "node_geometry_util.hh"

#include <fmt/format.h>

namespace blender::nodes::node_geo_solve_xpbd_constraints_cc {

template<typename T>
static AttributeReader<T> lookup_or_warn(
    AttributeAccessor &attributes,
    const StringRef attribute_id,
    const AttrDomain domain,
    const T &default_value,
    FunctionRef<void(const NodeWarningType type, const StringRef message)> error_fn)
{
  if (!attributes.contains(attribute_id)) {
    error_fn(geo_eval_log::NodeWarningType::Warning,
             fmt::format("Missing \"{}\" attribute", attribute_id));
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

constexpr float default_fps = 1.0f / 25.0f;

// enum class ConstraintGeometryType {
//   /* Set position of a point to a target vector. */
//   PositionGoal,
//   /* Set orientation of an edge to a target rotation. */
//   RotationGoal,
//   /* Enforces edge length and aligns forward direction with the edge vector. */
//   StretchShear,
//   /* Enforces angles between neighboring edges to their relative rest orientation. */
//   BendTwist,
//   /* Keep contact points from penetrating. */
//   Contact,
// };

// struct ConstraintEvalParams {
//   /* TODO split this by EvaluationTarget, only either (positions + rotations) or (velocities +
//    * angular_velocities) are used. */
//   Array<float3> positions;
//   Array<math::Quaternion> rotations;
//   VArraySpan<float3> old_positions;
//   VArraySpan<math::Quaternion> old_rotations;
//
//   Array<float3> velocities;
//   Array<float3> angular_velocities;
//   VArraySpan<float3> orig_velocities;
//   VArraySpan<float3> orig_angular_velocities;
//
//   VArraySpan<float> position_weights;
//   VArraySpan<float> rotation_weights;
//
//   Span<float4x4> collider_transforms;
//   /* Collider transforms of the previous frame for computing velocity constraints. */
//   Span<float4x4> old_collider_transforms;
//
//   struct {
//     VArray<int> solver_group;
//     VArraySpan<int> points;
//     VArraySpan<float> alpha;
//     VArraySpan<float> beta;
//     VArraySpan<float3> goal_positions;
//     VArraySpan<float3> goal_velocities;
//   } position_goal;
//   struct {
//     VArray<int> solver_group;
//     VArraySpan<int> points;
//     VArraySpan<float> alpha;
//     VArraySpan<float> beta;
//     VArraySpan<math::Quaternion> goal_rotations;
//     VArraySpan<float3> goal_angular_velocities;
//   } rotation_goal;
//   struct {
//     VArray<int> solver_group;
//     VArraySpan<int> points1;
//     VArraySpan<int> points2;
//     VArraySpan<float> alpha;
//     VArraySpan<float> beta;
//     VArraySpan<float> edge_lengths;
//   } stretch_shear;
//   struct {
//     VArray<int> solver_group;
//     VArraySpan<int> points1;
//     VArraySpan<int> points2;
//     VArraySpan<float> alpha;
//     VArraySpan<float> beta;
//     VArraySpan<float> edge_lengths;
//     VArraySpan<float> darboux_w;
//     VArraySpan<float3> darboux_xyz;
//   } bend_twist;
//   struct {
//     VArray<int> solver_group;
//     VArraySpan<int> points;
//     VArraySpan<int> collider_index;
//     VArraySpan<float3> local_position1;
//     VArraySpan<float3> local_position2;
//     VArraySpan<float3> normal;
//
//     VArraySpan<float> friction;
//     VArraySpan<float> restitution;
//
//     /* Remember active contacts for later velocity update. */
//     bke::SpanAttributeWriter<bool> active;
//   } contact;
// };

enum class EvaluationTarget {
  Positions,
  Velocities,
};

enum class SolverMethod {
  GaussSeidel,
  Jacobi,
};

struct SolverParams {
  float delta_time;
  float delta_time_squared;
  float inv_delta_time_squared;
  EvaluationTarget target;

  std::function<void(const NodeWarningType type, const StringRef message)> error_message_add;
  /* Perform debug checks on user inputs at runtime. This helps avoid common errors but has a
   * significant performance cost, so should be optional. */
  bool debug_check;
};

struct ConstraintEvalParams {
  Array<float3> positions;
  Array<math::Quaternion> rotations;
  VArraySpan<float3> old_positions;
  VArraySpan<math::Quaternion> old_rotations;

  /* Velocity inputs are not defined during the position update stage. */
  Array<float3> velocities;
  Array<float3> angular_velocities;
  VArraySpan<float3> orig_velocities;
  VArraySpan<float3> orig_angular_velocities;

  VArraySpan<float> position_weights;
  VArraySpan<float> rotation_weights;

  Span<float4x4> collider_transforms;
  /* Collider transforms of the previous frame for computing velocity constraints. */
  Span<float4x4> old_collider_transforms;
};

struct ConstraintClosure {
  using ErrorFn = FunctionRef<void(const NodeWarningType type, const StringRef message)>;

  GeometrySet geometry_set;
  VArray<int> solver_groups;

  ConstraintClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : geometry_set(std::move(geometry_set))
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::AttributeAccessor> attributes = component.attributes();
    this->solver_groups = *lookup_or_warn<int>(
        *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0, error_fn);
  }
  virtual ~ConstraintClosure() = default;

  virtual void apply_to_positions(const SolverParams &solver_params,
                                  ConstraintEvalParams &eval_params,
                                  const IndexMask &group_mask) = 0;
  virtual void apply_to_velocities(const SolverParams &solver_params,
                                   ConstraintEvalParams &eval_params,
                                   const IndexMask &group_mask) = 0;

  virtual void finish_attributes() = 0;
};

struct PositionGoalClosure : public ConstraintClosure {
  VArraySpan<int> points;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float3> goal_positions;
  VArraySpan<float3> goal_velocities;

  PositionGoalClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::AttributeAccessor> attributes = component.attributes();
    this->points = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->goal_positions = *lookup_or_warn<float3>(
        *attributes, "goal_position", AttrDomain::Point, float3(0.0f), error_fn);
    this->goal_velocities = *lookup_or_warn<float3>(
        *attributes, "goal_velocity", AttrDomain::Point, float3(0.0f), error_fn);
  }

  template<bool debug_check>
  void do_apply_to_positions(const SolverParams &solver_params,
                             ConstraintEvalParams &eval_params,
                             const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        eval_params.positions.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      const float alpha = this->alpha[index];
      const float3 &goal = this->goal_positions[index];

      position_checker.claim_variable(point1);

      float lambda = 0.0f;
      float3 &position1 = eval_params.positions[point1];
      xpbd_constraints::eval_position_goal(
          goal, alpha * solver_params.inv_delta_time_squared, lambda, position1);
    });

    if (position_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(const SolverParams &solver_params,
                              ConstraintEvalParams &eval_params,
                              const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> velocity_checker(
        eval_params.velocities.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      const float beta = this->beta[index];
      const float3 &goal = this->goal_velocities[index];

      velocity_checker.claim_variable(point1);

      float lambda = 0.0f;
      float3 &velocity1 = eval_params.velocities[point1];
      xpbd_constraints::eval_velocity_goal(
          goal, beta * solver_params.inv_delta_time_squared, lambda, velocity1);
    });

    if (velocity_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(const SolverParams &params,
                          ConstraintEvalParams &constraint_params,
                          const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, constraint_params, group_mask);
    }
  }

  void apply_to_velocities(const SolverParams &params,
                           ConstraintEvalParams &constraint_params,
                           const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, constraint_params, group_mask);
    }
  }

  void finish_attributes() override {}
};

struct RotationGoalClosure : public ConstraintClosure {
  VArraySpan<int> points;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<math::Quaternion> goal_rotations;
  VArraySpan<float3> goal_angular_velocities;

  RotationGoalClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::AttributeAccessor> attributes = component.attributes();
    this->points = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->goal_rotations = *lookup_or_warn<math::Quaternion>(
        *attributes, "goal_rotation", AttrDomain::Point, math::Quaternion::identity(), error_fn);
    this->goal_angular_velocities = *lookup_or_warn<float3>(
        *attributes, "goal_angular_velocity", AttrDomain::Point, float3(0.0f), error_fn);
  }

  template<bool debug_check>
  void do_apply_to_positions(const SolverParams &solver_params,
                             ConstraintEvalParams &eval_params,
                             const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        eval_params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      const float alpha = this->alpha[index];
      const math::Quaternion &goal = this->goal_rotations[index];

      rotation_checker.claim_variable(point1);

      float lambda = 0.0f;
      math::Quaternion &rotation1 = eval_params.rotations[point1];
      xpbd_constraints::eval_rotation_goal<true>(
          goal, alpha * solver_params.inv_delta_time_squared, lambda, rotation1);
    });

    if (rotation_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(const SolverParams &solver_params,
                              ConstraintEvalParams &eval_params,
                              const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> angular_velocity_checker(
        eval_params.angular_velocities.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points[index];
      const float beta = this->beta[index];
      const float3 &goal = this->goal_angular_velocities[index];

      angular_velocity_checker.claim_variable(point1);

      float lambda = 0.0f;
      float3 &angular_velocity1 = eval_params.angular_velocities[point1];
      xpbd_constraints::eval_angular_velocity_goal(
          goal, beta * solver_params.inv_delta_time_squared, lambda, angular_velocity1);
    });

    if (angular_velocity_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(const SolverParams &params,
                          ConstraintEvalParams &constraint_params,
                          const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, constraint_params, group_mask);
    }
  }

  void apply_to_velocities(const SolverParams &params,
                           ConstraintEvalParams &constraint_params,
                           const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, constraint_params, group_mask);
    }
  }

  void finish_attributes() override {}
};

struct StretchShearClosure : public ConstraintClosure {
  VArraySpan<int> points1;
  VArraySpan<int> points2;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float> edge_lengths;

  StretchShearClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::AttributeAccessor> attributes = component.attributes();
    this->points1 = *lookup_or_warn<int>(*attributes, ATTR_POINT1, AttrDomain::Point, 0, error_fn);
    this->points2 = *lookup_or_warn<int>(*attributes, ATTR_POINT2, AttrDomain::Point, 0, error_fn);
    this->alpha = *attributes->lookup_or_default<float>(ATTR_ALPHA, AttrDomain::Point, 0.0f);
    this->beta = *attributes->lookup_or_default<float>(ATTR_BETA, AttrDomain::Point, 0.0f);
    this->edge_lengths = *lookup_or_warn<float>(
        *attributes, "edge_length", AttrDomain::Point, 0.0f, error_fn);
  }

  template<bool debug_check>
  void do_apply_to_positions(const SolverParams &solver_params,
                             ConstraintEvalParams &eval_params,
                             const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        eval_params.positions.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        eval_params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      const float edge_length = this->edge_lengths[index];
      const float alpha = this->alpha[index];

      position_checker.claim_variable(point1);
      position_checker.claim_variable(point2);
      rotation_checker.claim_variable(point1);

      float3 lambda = float3(0.0f);
      float3 &position1 = eval_params.positions[point1];
      float3 &position2 = eval_params.positions[point2];
      math::Quaternion &rotation = eval_params.rotations[point1];
      const float weight_pos1 = eval_params.position_weights[point1];
      const float weight_pos2 = eval_params.position_weights[point2];
      const float weight_rot = eval_params.rotation_weights[point1];
      xpbd_constraints::eval_position_stretch_shear<true>(weight_pos1,
                                                          weight_pos2,
                                                          weight_rot,
                                                          edge_length,
                                                          alpha *
                                                              solver_params.inv_delta_time_squared,
                                                          lambda,
                                                          position1,
                                                          position2,
                                                          rotation);
    });

    if (position_checker.has_overlap() || rotation_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(const SolverParams &solver_params,
                              ConstraintEvalParams &eval_params,
                              const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> velocity_checker(
        eval_params.velocities.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> angular_velocity_checker(
        eval_params.angular_velocities.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      const float edge_length = this->edge_lengths[index];
      const float beta = this->beta[index];

      velocity_checker.claim_variable(point1);
      velocity_checker.claim_variable(point2);
      angular_velocity_checker.claim_variable(point1);

      float3 lambda = float3(0.0f);
      const math::Quaternion &rotation = eval_params.rotations[point1];
      float3 &velocity1 = eval_params.velocities[point1];
      float3 &velocity2 = eval_params.velocities[point2];
      float3 &angular_velocity = eval_params.angular_velocities[point1];
      const float weight_pos1 = eval_params.position_weights[point1];
      const float weight_pos2 = eval_params.position_weights[point2];
      const float weight_rot = eval_params.rotation_weights[point1];
      xpbd_constraints::eval_velocity_stretch_shear(rotation,
                                                    weight_pos1,
                                                    weight_pos2,
                                                    weight_rot,
                                                    edge_length,
                                                    beta * solver_params.inv_delta_time_squared,
                                                    lambda,
                                                    velocity1,
                                                    velocity2,
                                                    angular_velocity);
    });

    if (velocity_checker.has_overlap() || angular_velocity_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(const SolverParams &params,
                          ConstraintEvalParams &constraint_params,
                          const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, constraint_params, group_mask);
    }
  }

  void apply_to_velocities(const SolverParams &params,
                           ConstraintEvalParams &constraint_params,
                           const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, constraint_params, group_mask);
    }
  }

  void finish_attributes() override {}
};

struct BendTwistClosure : public ConstraintClosure {
  VArraySpan<int> points1;
  VArraySpan<int> points2;
  VArraySpan<float> alpha;
  VArraySpan<float> beta;
  VArraySpan<float> edge_lengths;
  VArraySpan<float> darboux_w;
  VArraySpan<float3> darboux_xyz;

  BendTwistClosure(GeometrySet &&geometry_set, ErrorFn error_fn)
      : ConstraintClosure(std::move(geometry_set), error_fn)
  {
    PointCloudComponent &component =
        this->geometry_set.get_component_for_write<PointCloudComponent>();
    std::optional<bke::AttributeAccessor> attributes = component.attributes();
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
  }

  template<bool debug_check>
  void do_apply_to_positions(const SolverParams &solver_params,
                             ConstraintEvalParams &eval_params,
                             const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        eval_params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      const float edge_length = this->edge_lengths[index];
      const math::Quaternion darboux_vector = math::Quaternion(this->darboux_w[index],
                                                               this->darboux_xyz[index]);
      const float alpha = this->alpha[index];

      rotation_checker.claim_variable(point1);
      rotation_checker.claim_variable(point2);

      float4 lambda = float4(0.0f);
      math::Quaternion &rotation1 = eval_params.rotations[point1];
      math::Quaternion &rotation2 = eval_params.rotations[point2];
      const float weight_rot1 = eval_params.rotation_weights[point1];
      const float weight_rot2 = eval_params.rotation_weights[point2];
      xpbd_constraints::eval_position_bend_twist<true>(weight_rot1,
                                                       weight_rot2,
                                                       edge_length,
                                                       darboux_vector,
                                                       alpha *
                                                           solver_params.inv_delta_time_squared,
                                                       lambda,
                                                       rotation1,
                                                       rotation2);
    });

    if (rotation_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(const SolverParams &solver_params,
                              ConstraintEvalParams &eval_params,
                              const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> angular_velocity_checker(
        eval_params.angular_velocities.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point1 = this->points1[index];
      const int point2 = this->points2[index];
      const float edge_length = this->edge_lengths[index];
      const float beta = this->beta[index];

      angular_velocity_checker.claim_variable(point1);
      angular_velocity_checker.claim_variable(point2);

      float3 lambda = float3(0.0f);
      float3 &angular_velocity1 = eval_params.angular_velocities[point1];
      float3 &angular_velocity2 = eval_params.angular_velocities[point2];
      const float weight_rot1 = eval_params.rotation_weights[point1];
      const float weight_rot2 = eval_params.rotation_weights[point2];
      xpbd_constraints::eval_velocity_bend_twist(weight_rot1,
                                                 weight_rot2,
                                                 edge_length,
                                                 beta * solver_params.inv_delta_time_squared,
                                                 lambda,
                                                 angular_velocity1,
                                                 angular_velocity2);
    });

    if (angular_velocity_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(const SolverParams &params,
                          ConstraintEvalParams &constraint_params,
                          const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, constraint_params, group_mask);
    }
  }

  void apply_to_velocities(const SolverParams &params,
                           ConstraintEvalParams &constraint_params,
                           const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, constraint_params, group_mask);
    }
  }

  void finish_attributes() override {}
};

struct ContactClosure : public ConstraintClosure {
  VArraySpan<int> points;
  VArraySpan<int> collider_index;
  VArraySpan<float3> local_position1;
  VArraySpan<float3> local_position2;
  VArraySpan<float3> normal;

  VArraySpan<float> friction;
  VArraySpan<float> restitution;

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
  }

  template<bool debug_check>
  void do_apply_to_positions(const SolverParams &solver_params,
                             ConstraintEvalParams &eval_params,
                             const IndexMask &group_mask)
  {
    xpbd_constraints::error_check::VariableChecker<debug_check> position_checker(
        eval_params.positions.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> rotation_checker(
        eval_params.rotations.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      const int point = this->points[index];
      const int collider_index = this->collider_index[index];

      position_checker.claim_variable(point);
      rotation_checker.claim_variable(point);

      if (!eval_params.collider_transforms.index_range().contains(collider_index)) {
        return;
      };
      const float4x4 collider_transform = eval_params.collider_transforms[collider_index];
      float3 collider_position;
      math::Quaternion collider_rotation;
      float3 collider_scale;
      math::to_loc_rot_scale(
          collider_transform, collider_position, collider_rotation, collider_scale);

      float3 &position = eval_params.positions[point];
      math::Quaternion &rotation = eval_params.rotations[point];
      bool &active = this->active.span[index];

      const float3 &local_position1 = this->local_position1[index];
      const float3 &local_position2 = this->local_position2[index];
      const float3 &normal = this->normal[index];
      if constexpr (debug_check) {
        if (!math::is_unit(normal)) {
          solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                          "Contact normal vector not normalized");
        }
      }

      /* Contact constraints are stiff. */
      const float alpha = 0.0f;

      /* TODO store this for warm-starting. */
      float lambda = 0.0f;
      /* Zero weights for the collider, only the point can move. */
      active = xpbd_constraints::eval_position_contact(1.0f,
                                                       0.0f,
                                                       1.0f,
                                                       0.0f,
                                                       local_position1,
                                                       local_position2,
                                                       normal,
                                                       alpha *
                                                           solver_params.inv_delta_time_squared,
                                                       lambda,
                                                       position,
                                                       collider_position,
                                                       rotation,
                                                       collider_rotation);
    });

    if (position_checker.has_overlap() || rotation_checker.has_overlap()) {
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  template<bool debug_check>
  void do_apply_to_velocities(const SolverParams &solver_params,
                              ConstraintEvalParams &eval_params,
                              const IndexMask &group_mask)
  {
    const float inv_delta_time = math::safe_rcp(solver_params.delta_time);

    xpbd_constraints::error_check::VariableChecker<debug_check> velocity_checker(
        eval_params.velocities.index_range());
    xpbd_constraints::error_check::VariableChecker<debug_check> angular_velocity_checker(
        eval_params.angular_velocities.index_range());

    group_mask.foreach_index(GrainSize(1024), [&](const int index) {
      /* Active status is determined by the position evaluation. */
      if (!this->active.span[index]) {
        return;
      }

      const int point = this->points[index];
      velocity_checker.claim_variable(point);
      angular_velocity_checker.claim_variable(point);
      const int collider_index = this->collider_index[index];
      if (!eval_params.collider_transforms.index_range().contains(collider_index)) {
        return;
      };

      /* Local positions are relative to moving point and collider respectively. */
      const float3 &local_position1 = this->local_position1[index];
      const float3 &local_position2 = this->local_position2[index];
      /* Normal is a fixed shared direction for both participants. */
      const float3 &normal = this->normal[index];
      if constexpr (debug_check) {
        if (!math::is_unit(normal)) {
          solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                          "Contact normal vector not normalized");
        }
      }
      const float restitution = this->restitution[index];
      const float friction = this->friction[index];

      const float4x4 &collider_transform = eval_params.collider_transforms[collider_index];
      const float4x4 &old_collider_transform = eval_params.old_collider_transforms[collider_index];

      float3 &velocity = eval_params.velocities[point];
      float3 &angular_velocity = eval_params.angular_velocities[point];
      const float3 &orig_velocity = eval_params.orig_velocities[point];
      const float3 &orig_angular_velocity = eval_params.orig_angular_velocities[point];

      /* Compute velocity from old/new collider transforms. */
      float3 collider_loc, old_collider_loc;
      math::Quaternion collider_rot, old_collider_rot;
      float3 collider_scale, old_collider_scale;
      math::to_loc_rot_scale(collider_transform, collider_loc, collider_rot, collider_scale);
      math::to_loc_rot_scale(
          old_collider_transform, old_collider_loc, old_collider_rot, old_collider_scale);
      float3 collider_velocity = (collider_loc - old_collider_loc) * inv_delta_time;
      float3 collider_angular_velocity =
          2.0f * (math::invert_normalized(old_collider_rot) * collider_rot).imaginary_part() *
          inv_delta_time;
      /* No change in animated collider velocity. */
      const float3 orig_collider_velocity = collider_velocity;
      const float3 orig_collider_angular_velocity = collider_angular_velocity;

      /* TODO store for warm-starting. */
      float lambda_restitution = 0.0f;
      float lambda_friction = 0.0f;
      xpbd_constraints::eval_velocity_contact(orig_velocity,
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
      solver_params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                      "Overlapping constraint solver groups");
    }
  }

  void apply_to_positions(const SolverParams &params,
                          ConstraintEvalParams &constraint_params,
                          const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_positions<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_positions<false>(params, constraint_params, group_mask);
    }
  }

  void apply_to_velocities(const SolverParams &params,
                           ConstraintEvalParams &constraint_params,
                           const IndexMask &group_mask) override
  {
    if (params.debug_check) {
      this->do_apply_to_velocities<true>(params, constraint_params, group_mask);
    }
    else {
      this->do_apply_to_velocities<false>(params, constraint_params, group_mask);
    }
  }

  void finish_attributes() override
  {
    this->active.finish();
  }
};

struct ConstraintTypeInfo {
  using ErrorFn = ConstraintClosure::ErrorFn;

  std::string ui_name;
  std::string ui_description;

  ConstraintClosure *(*get_closure)(GeometrySet &&geometry_set, ErrorFn error_fn);
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

static Span<ConstraintTypeInfo> get_constraint_info()
{
  static const Array<ConstraintTypeInfo> constraint_info = create_constraint_info();
  return constraint_info;
}

static Span<ConstraintTypeInfo> get_constraint_info_ordered()
{
  /* TODO currently relies on fixed order in get_constraint_info(),
   * could also re-order based on some priority value. */
  return get_constraint_info();
}

static void node_declare_positions(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();
  b.add_input<decl::Int>("Gauss-Seidel Steps").default_value(1).min(0);

  b.add_input<decl::Geometry>("Geometry");
  b.add_output<decl::Geometry>("Geometry").align_with_previous();

  b.add_input<decl::Vector>("Position").implicit_field_on(implicit_field_inputs::position, {2});
  b.add_output<decl::Vector>("Position").field_on({0}).align_with_previous();
  b.add_input<decl::Rotation>("Rotation").field_on({2}).hide_value();
  b.add_output<decl::Rotation>("Rotation").field_on({0}).align_with_previous();
  b.add_input<decl::Vector>("Old Position").field_on({2}).hide_value();
  b.add_input<decl::Rotation>("Old Rotation").field_on({2}).hide_value();

  b.add_input<decl::Float>("Position Weight")
      .default_value(1.0f)
      .field_on({2})
      .description("Influence weight of constraints for each point (inverse mass)");
  b.add_input<decl::Float>("Rotation Weight")
      .default_value(1.0f)
      .field_on({2})
      .description("Influence weight of constraints for each point (inverse moment of inertia)");

  for (const ConstraintTypeInfo &info : get_constraint_info()) {
    b.add_input<decl::Geometry>(info.ui_name)
        .supported_type(GeometryComponent::Type::PointCloud)
        .description(info.ui_description);
    b.add_output<decl::Geometry>(info.ui_name)
        .description(info.ui_description)
        .align_with_previous();
  }

  b.add_input<decl::Geometry>("Colliders")
      .only_instances()
      .description("Instances of colliders to evaluate contact transforms");

  b.add_input<decl::Bool>("Debug Checks")
      .default_value(false)
      .description("Perform checks on input data, which can impact performance");
}

static void node_declare_velocities(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();

  b.add_input<decl::Geometry>("Geometry");
  b.add_output<decl::Geometry>("Geometry").align_with_previous();

  b.add_input<decl::Vector>("Position").implicit_field_on(implicit_field_inputs::position, {2});
  b.add_input<decl::Rotation>("Rotation").field_on({2}).hide_value();
  b.add_input<decl::Vector>("Velocity").field_on({1}).hide_value();
  b.add_output<decl::Vector>("Velocity").field_on({0}).align_with_previous();
  b.add_input<decl::Vector>("Angular Velocity").field_on({1}).hide_value();
  b.add_output<decl::Vector>("Angular Velocity").field_on({0}).align_with_previous();
  b.add_input<decl::Vector>("Original Velocity").field_on({1}).hide_value();
  b.add_input<decl::Vector>("Original Angular Velocity").field_on({1}).hide_value();

  b.add_input<decl::Float>("Position Weight")
      .default_value(1.0f)
      .field_on({2})
      .description("Influence weight of constraints for each point (inverse mass)");
  b.add_input<decl::Float>("Rotation Weight")
      .default_value(1.0f)
      .field_on({2})
      .description("Influence weight of constraints for each point (inverse moment of inertia)");

  for (const ConstraintTypeInfo &info : get_constraint_info()) {
    b.add_input<decl::Geometry>(info.ui_name)
        .supported_type(GeometryComponent::Type::PointCloud)
        .description(info.ui_description);
    b.add_output<decl::Geometry>(info.ui_name)
        .description(info.ui_description)
        .align_with_previous();
  }

  b.add_input<decl::Geometry>("Colliders")
      .only_instances()
      .description("Instances of colliders to evaluate contact transforms");

  b.add_input<decl::Bool>("Debug Checks")
      .default_value(false)
      .description("Perform checks on input data, which can impact performance");
}

/* Evaluate a group of constraints in parallel.
 * It's important that none of the constraints write to the same variables. Solver groups must be
 * computed in such a way that each variable is only affected by one constraint in each group.
 */
static void evaluate_constraint_group(const SolverParams &solver_params,
                                      ConstraintEvalParams &eval_params,
                                      ConstraintClosure &closure,
                                      const IndexMask &group_mask)
{
  BLI_assert(!closure.solver_groups.is_empty());

  switch (solver_params.target) {
    case EvaluationTarget::Positions:
      closure.apply_to_positions(solver_params, eval_params, group_mask);
      break;
    case EvaluationTarget::Velocities:
      closure.apply_to_velocities(solver_params, eval_params, group_mask);
      break;
  }
}

static void do_single_constraint_passes(const SolverParams &solver_params,
                                        ConstraintEvalParams &eval_params,
                                        ConstraintClosure &closure)
{
  const VArray<int> &solver_groups = closure.solver_groups;
  if (solver_groups.is_empty()) {
    return;
  }

  const IndexRange constraints = solver_groups.index_range();
  IndexMaskMemory memory;
  VectorSet<int> unique_group_ids;
  Vector<IndexMask> group_index_masks = IndexMask::from_group_ids(
      constraints, solver_groups, memory, unique_group_ids);

  /* Sort group indices to ensure solver groups are executed in consistent order.
   * IndexMask::from_group_ids creates masks in the order they appear in the data:
   * whichever element comes first creates a mask and index, regardless of the actual group ID
   * values and their relative ordering.
   */
  Array<int> sorted_group_indices(unique_group_ids.size());
  array_utils::fill_index_range(sorted_group_indices.as_mutable_span());
  std::sort(sorted_group_indices.begin(),
            sorted_group_indices.end(),
            [&](const int index_a, const int index_b) {
              const int group_id_a = unique_group_ids[index_a];
              const int group_id_b = unique_group_ids[index_b];
              return group_id_a < group_id_b;
            });

  /* Solve in consistent order by using the sorted index set. */
  for (const int i_group : sorted_group_indices) {
    const IndexMask &group_mask = group_index_masks[i_group];
    /* Group ID is not really relevant at this point. */
    /* const int group_id = unique_group_ids[i_group]; */

    evaluate_constraint_group(solver_params, eval_params, closure, group_mask);
  }
}

static void do_gauss_seidel_step(const SolverParams &solver_params,
                                 ConstraintEvalParams &eval_params,
                                 const Span<ConstraintClosure *> closures)
{
  for (ConstraintClosure *closure : closures) {
    do_single_constraint_passes(solver_params, eval_params, *closure);
  }
}

static void do_jacobi_step(const SolverParams & /*solver_params*/,
                           ConstraintEvalParams & /*eval_params*/,
                           const Span<ConstraintClosure *> /*closures*/)
{
  /* TODO */
}

static void do_solver_steps(const SolverMethod method,
                            const int steps,
                            const SolverParams &solver_params,
                            ConstraintEvalParams &eval_params,
                            const Span<ConstraintClosure *> closures)
{
  if (steps < 1) {
    return;
  }

  for ([[maybe_unused]] const int step : IndexRange(steps)) {
    switch (method) {
      case SolverMethod::GaussSeidel:
        do_gauss_seidel_step(solver_params, eval_params, closures);
        break;
      case SolverMethod::Jacobi:
        do_jacobi_step(solver_params, eval_params, closures);
        break;
    }
  }
}

static void bind_constraint_closures(GeoNodeExecParams params,
                                     Vector<ConstraintClosure *> &closures)
{
  auto error_fn = [params](const NodeWarningType type, const StringRef message) {
    params.error_message_add(type, message);
  };

  const Span<ConstraintTypeInfo> constraint_infos = get_constraint_info_ordered();
  closures.reinitialize(constraint_infos.size());
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    GeometrySet geometry_set = params.extract_input<GeometrySet>(info.ui_name);
    if (geometry_set.has_component<PointCloudComponent>()) {
      // PointCloudComponent &constraint_component =
      //     geometry_set.get_component_for_write<PointCloudComponent>();
      closures[i] = info.get_closure(std::move(geometry_set), error_fn);
    }
  }

  // if (std::optional<MutableAttributeAccessor> attributes =
  //         constraint_attributes[int(ConstraintGeometryType::PositionGoal)])
  //{
  //   constraint_params.position_goal.solver_group = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
  //   constraint_params.position_goal.points = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
  //   constraint_params.position_goal.alpha = *attributes->lookup_or_default<float>(
  //       ATTR_ALPHA, AttrDomain::Point, 0.0f);
  //   constraint_params.position_goal.beta = *attributes->lookup_or_default<float>(
  //       ATTR_BETA, AttrDomain::Point, 0.0f);
  //   constraint_params.position_goal.goal_positions = *lookup_or_warn<float3>(
  //       params, *attributes, "goal_position", AttrDomain::Point, float3(0.0f));
  //   constraint_params.position_goal.goal_velocities = *lookup_or_warn<float3>(
  //       params, *attributes, "goal_velocity", AttrDomain::Point, float3(0.0f));
  // }
  // if (std::optional<MutableAttributeAccessor> attributes =
  //         constraint_attributes[int(ConstraintGeometryType::RotationGoal)])
  //{
  //   constraint_params.rotation_goal.solver_group = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
  //   constraint_params.rotation_goal.points = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
  //   constraint_params.rotation_goal.alpha = *attributes->lookup_or_default<float>(
  //       ATTR_ALPHA, AttrDomain::Point, 0.0f);
  //   constraint_params.rotation_goal.beta = *attributes->lookup_or_default<float>(
  //       ATTR_BETA, AttrDomain::Point, 0.0f);
  //   constraint_params.rotation_goal.goal_rotations = *lookup_or_warn<math::Quaternion>(
  //       params, *attributes, "goal_rotation", AttrDomain::Point, math::Quaternion::identity());
  //   constraint_params.rotation_goal.goal_angular_velocities = *lookup_or_warn<float3>(
  //       params, *attributes, "goal_angular_velocity", AttrDomain::Point, float3(0.0f));
  // }
  // if (std::optional<MutableAttributeAccessor> attributes =
  //         constraint_attributes[int(ConstraintGeometryType::StretchShear)])
  //{
  //   constraint_params.stretch_shear.solver_group = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
  //   constraint_params.stretch_shear.points1 = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
  //   constraint_params.stretch_shear.points2 = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT2, AttrDomain::Point, 0);
  //   constraint_params.stretch_shear.alpha = *attributes->lookup_or_default<float>(
  //       ATTR_ALPHA, AttrDomain::Point, 0.0f);
  //   constraint_params.stretch_shear.beta = *attributes->lookup_or_default<float>(
  //       ATTR_BETA, AttrDomain::Point, 0.0f);
  //   constraint_params.stretch_shear.edge_lengths = *lookup_or_warn<float>(
  //       params, *attributes, "edge_length", AttrDomain::Point, 0.0f);
  // }
  // if (std::optional<MutableAttributeAccessor> attributes =
  //         constraint_attributes[int(ConstraintGeometryType::BendTwist)])
  //{
  //   constraint_params.bend_twist.solver_group = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
  //   constraint_params.bend_twist.points1 = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
  //   constraint_params.bend_twist.points2 = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT2, AttrDomain::Point, 0);
  //   constraint_params.bend_twist.alpha = *attributes->lookup_or_default<float>(
  //       ATTR_ALPHA, AttrDomain::Point, 0.0f);
  //   constraint_params.bend_twist.beta = *attributes->lookup_or_default<float>(
  //       ATTR_BETA, AttrDomain::Point, 0.0f);
  //   constraint_params.bend_twist.edge_lengths = *lookup_or_warn<float>(
  //       params, *attributes, "edge_length", AttrDomain::Point, 0.0f);
  //   constraint_params.bend_twist.darboux_w = *lookup_or_warn<float>(
  //       params, *attributes, "darboux_w", AttrDomain::Point, float(0.0f));
  //   constraint_params.bend_twist.darboux_xyz = *lookup_or_warn<float3>(
  //       params, *attributes, "darboux_xyz", AttrDomain::Point, float3(0.0f));
  // }
  // if (std::optional<MutableAttributeAccessor> attributes =
  //         constraint_attributes[int(ConstraintGeometryType::Contact)])
  //{
  //   constraint_params.contact.solver_group = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
  //   constraint_params.contact.points = *lookup_or_warn<int>(
  //       params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
  //   constraint_params.contact.collider_index = *lookup_or_warn<int>(
  //       params, *attributes, "collider_index", AttrDomain::Point, 0);
  //   constraint_params.contact.local_position1 = *lookup_or_warn<float3>(
  //       params, *attributes, "local_position1", AttrDomain::Point, float3(0.0f));
  //   constraint_params.contact.local_position2 = *lookup_or_warn<float3>(
  //       params, *attributes, "local_position2", AttrDomain::Point, float3(0.0f));
  //   constraint_params.contact.normal = *lookup_or_warn<float3>(
  //       params, *attributes, "normal", AttrDomain::Point, float3(0.0f));

  //  constraint_params.contact.friction = *attributes->lookup_or_default<float>(
  //      "friction", AttrDomain::Point, 0.0f);
  //  constraint_params.contact.restitution = *attributes->lookup_or_default<float>(
  //      "restitution", AttrDomain::Point, 0.0f);

  //  colliders_geometry_set = params.extract_input<GeometrySet>("Colliders");
  //  constraint_params.collider_transforms =
  //      colliders_geometry_set.has_instances() ?
  //          colliders_geometry_set.get_instances()->transforms() :
  //          Span<float4x4>{};
  //  /* XXX Transforms of the previous frame are not currently available, these are always the
  //   * same as the current frame. Eventually this will allow transfer of velocity from animated
  //   * colliders. */
  //  constraint_params.old_collider_transforms = constraint_params.collider_transforms;

  //  const int num_constraints = attributes->domain_size(AttrDomain::Point);
  //  constraint_params.contact.active = attributes->lookup_or_add_for_write_span<bool>(
  //      ATTR_ACTIVE,
  //      AttrDomain::Point,
  //      bke::AttributeInitVArray(VArray<bool>::ForSingle(false, num_constraints)));
  //}
}

// static void write_constraint_attributes(ConstraintEvalParams &constraint_params)
//{
//   if (constraint_params.contact.active) {
//     constraint_params.contact.active.finish();
//   }
// }

static void node_geo_exec_positions(GeoNodeExecParams params)
{
  const int gauss_seidel_steps = params.extract_input<int>("Gauss-Seidel Steps");
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
  const float delta_time_squared = delta_time * delta_time;
  const float inv_delta_time_squared = math::safe_rcp(delta_time_squared);
  const bool debug_check = params.extract_input<bool>("Debug Checks");
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");
  Field<float3> position_field = params.extract_input<Field<float3>>("Position");
  Field<math::Quaternion> rotation_field = params.extract_input<Field<math::Quaternion>>(
      "Rotation");
  Field<float3> old_position_field = params.extract_input<Field<float3>>("Old Position");
  Field<math::Quaternion> old_rotation_field = params.extract_input<Field<math::Quaternion>>(
      "Old Rotation");
  std::optional<std::string> position_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Position");
  std::optional<std::string> rotation_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Rotation");
  Field<float> position_weight_field = params.extract_input<Field<float>>("Position Weight");
  Field<float> rotation_weight_field = params.extract_input<Field<float>>("Rotation Weight");

  GeometrySet colliders_geometry_set = params.extract_input<GeometrySet>("Colliders");
  Span<float4x4> collider_transforms = colliders_geometry_set.has_instances() ?
                                           colliders_geometry_set.get_instances()->transforms() :
                                           Span<float4x4>{};

  Vector<ConstraintClosure *> constraint_closures;
  bind_constraint_closures(params, constraint_closures);

  static const Array<GeometryComponent::Type> types = {bke::GeometryComponent::Type::Mesh,
                                                       bke::GeometryComponent::Type::PointCloud,
                                                       bke::GeometryComponent::Type::Curve,
                                                       bke::GeometryComponent::Type::GreasePencil};
  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    for (const bke::GeometryComponent::Type type : types) {
      if (geometry_set.has(type)) {
        bke::GeometryComponent &component = geometry_set.get_component_for_write(type);
        std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
        if (!attributes) {
          continue;
        }

        const int num_points = attributes->domain_size(AttrDomain::Point);
        ConstraintEvalParams eval_params;
        eval_params.positions.reinitialize(num_points);
        eval_params.rotations.reinitialize(num_points);

        const bke::GeometryFieldContext field_context{component, AttrDomain::Point};
        fn::FieldEvaluator evaluator{field_context, num_points};
        evaluator.add_with_destination(position_field, eval_params.positions.as_mutable_span());
        evaluator.add_with_destination(rotation_field, eval_params.rotations.as_mutable_span());
        evaluator.add(old_position_field);
        evaluator.add(old_rotation_field);
        evaluator.add(position_weight_field);
        evaluator.add(rotation_weight_field);
        evaluator.evaluate();
        eval_params.old_positions = evaluator.get_evaluated<float3>(2);
        eval_params.old_rotations = evaluator.get_evaluated<math::Quaternion>(3);
        eval_params.position_weights = evaluator.get_evaluated<float>(4);
        eval_params.rotation_weights = evaluator.get_evaluated<float>(5);

        eval_params.collider_transforms = collider_transforms;
        /* XXX Transforms of the previous frame are not currently available, these are always the
         * same as the current frame. Eventually this will allow transfer of velocity from animated
         * colliders. */
        eval_params.old_collider_transforms = eval_params.collider_transforms;

        SolverParams solver_params = {
            delta_time, delta_time_squared, inv_delta_time_squared, EvaluationTarget::Positions};
        solver_params.error_message_add = [params](const NodeWarningType type,
                                                   const StringRef message) {
          params.error_message_add(type, message);
        };
        solver_params.debug_check = debug_check;

        do_solver_steps(SolverMethod::GaussSeidel,
                        gauss_seidel_steps,
                        solver_params,
                        eval_params,
                        constraint_closures);

        if (position_output_id) {
          AttributeWriter<float3> positions_writer = attributes->lookup_or_add_for_write<float3>(
              *position_output_id, AttrDomain::Point);
          BLI_assert(eval_params.positions.size() == num_points);
          positions_writer.varray.set_all(eval_params.positions);
          positions_writer.finish();
        }
        if (rotation_output_id) {
          AttributeWriter<math::Quaternion> rotations_writer =
              attributes->lookup_or_add_for_write<math::Quaternion>(*rotation_output_id,
                                                                    AttrDomain::Point);
          BLI_assert(eval_params.rotations.size() == num_points);
          rotations_writer.varray.set_all(eval_params.rotations);
          rotations_writer.finish();
        }
      }
    }
  });

  for (ConstraintClosure *closure : constraint_closures) {
    closure->finish_attributes();
  }

  params.set_output("Geometry", geometry_set);

  const Span<ConstraintTypeInfo> constraint_infos = get_constraint_info();
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    params.set_output(info.ui_name, constraint_closures[i]->geometry_set);
  }
}

static void node_geo_exec_velocities(GeoNodeExecParams params)
{
  constexpr int gauss_seidel_steps = 1;
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
  const float delta_time_squared = delta_time * delta_time;
  const float inv_delta_time_squared = math::safe_rcp(delta_time_squared);
  const bool debug_check = params.extract_input<bool>("Debug Checks");
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");
  Field<float3> position_field = params.extract_input<Field<float3>>("Position");
  Field<math::Quaternion> rotation_field = params.extract_input<Field<math::Quaternion>>(
      "Rotation");
  Field<float3> velocity_field = params.extract_input<Field<float3>>("Velocity");
  Field<float3> angular_velocity_field = params.extract_input<Field<float3>>("Angular Velocity");
  std::optional<std::string> velocity_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Velocity");
  std::optional<std::string> angular_velocity_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Angular Velocity");
  Field<float3> orig_velocity_field = params.extract_input<Field<float3>>("Original Velocity");
  Field<float3> orig_angular_velocity_field = params.extract_input<Field<float3>>(
      "Original Angular Velocity");
  Field<float> position_weight_field = params.extract_input<Field<float>>("Position Weight");
  Field<float> rotation_weight_field = params.extract_input<Field<float>>("Rotation Weight");

  GeometrySet colliders_geometry_set = params.extract_input<GeometrySet>("Colliders");
  Span<float4x4> collider_transforms = colliders_geometry_set.has_instances() ?
                                           colliders_geometry_set.get_instances()->transforms() :
                                           Span<float4x4>{};

  Vector<ConstraintClosure *> constraint_closures;
  bind_constraint_closures(params, constraint_closures);

  static const Array<GeometryComponent::Type> types = {bke::GeometryComponent::Type::Mesh,
                                                       bke::GeometryComponent::Type::PointCloud,
                                                       bke::GeometryComponent::Type::Curve,
                                                       bke::GeometryComponent::Type::GreasePencil};
  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    for (const bke::GeometryComponent::Type type : types) {
      if (geometry_set.has(type)) {
        bke::GeometryComponent &component = geometry_set.get_component_for_write(type);
        std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
        if (!attributes) {
          continue;
        }

        const int num_points = attributes->domain_size(AttrDomain::Point);
        ConstraintEvalParams eval_params;
        eval_params.positions.reinitialize(num_points);
        eval_params.rotations.reinitialize(num_points);
        eval_params.velocities.reinitialize(num_points);
        eval_params.angular_velocities.reinitialize(num_points);

        const bke::GeometryFieldContext field_context{component, AttrDomain::Point};
        fn::FieldEvaluator evaluator{field_context, num_points};
        evaluator.add_with_destination(position_field, eval_params.positions.as_mutable_span());
        evaluator.add_with_destination(rotation_field, eval_params.rotations.as_mutable_span());
        evaluator.add_with_destination(velocity_field, eval_params.velocities.as_mutable_span());
        evaluator.add_with_destination(angular_velocity_field,
                                       eval_params.angular_velocities.as_mutable_span());
        evaluator.add(orig_velocity_field);
        evaluator.add(orig_angular_velocity_field);
        evaluator.add(position_weight_field);
        evaluator.add(rotation_weight_field);
        evaluator.evaluate();
        eval_params.orig_velocities = evaluator.get_evaluated<float3>(4);
        eval_params.orig_angular_velocities = evaluator.get_evaluated<float3>(5);
        eval_params.position_weights = evaluator.get_evaluated<float>(6);
        eval_params.rotation_weights = evaluator.get_evaluated<float>(7);

        eval_params.collider_transforms = collider_transforms;
        /* XXX Transforms of the previous frame are not currently available, these are always the
         * same as the current frame. Eventually this will allow transfer of velocity from animated
         * colliders. */
        eval_params.old_collider_transforms = eval_params.collider_transforms;

        SolverParams solver_params = {
            delta_time, delta_time_squared, inv_delta_time_squared, EvaluationTarget::Velocities};
        solver_params.error_message_add = [params](const NodeWarningType type,
                                                   const StringRef message) {
          params.error_message_add(type, message);
        };
        solver_params.debug_check = debug_check;

        do_solver_steps(SolverMethod::GaussSeidel,
                        gauss_seidel_steps,
                        solver_params,
                        eval_params,
                        constraint_closures);

        if (velocity_output_id) {
          AttributeWriter<float3> velocities_writer = attributes->lookup_or_add_for_write<float3>(
              *velocity_output_id, AttrDomain::Point);
          BLI_assert(eval_params.velocities.size() == num_points);
          velocities_writer.varray.set_all(eval_params.velocities);
          velocities_writer.finish();
        }
        if (angular_velocity_output_id) {
          AttributeWriter<float3> angular_velocities_writer =
              attributes->lookup_or_add_for_write<float3>(*angular_velocity_output_id,
                                                          AttrDomain::Point);
          BLI_assert(eval_params.angular_velocities.size() == num_points);
          angular_velocities_writer.varray.set_all(eval_params.angular_velocities);
          angular_velocities_writer.finish();
        }
      }
    }
  });

  for (ConstraintClosure *closure : constraint_closures) {
    closure->finish_attributes();
  }

  params.set_output("Geometry", geometry_set);

  const Span<ConstraintTypeInfo> constraint_infos = get_constraint_info();
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    params.set_output(info.ui_name, constraint_closures[i]->geometry_set);
  }
}

static void node_register_position_solve()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeSolvePositionConstraints");
  ntype.ui_name = "Solve XPBD Position Constraints";
  ntype.ui_description =
      "Solve position and rotation constraints on geometry using the XPBD framework";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  node_type_size(&ntype, 200, 120, 300);
  ntype.geometry_node_execute = node_geo_exec_positions;
  ntype.declare = node_declare_positions;
  blender::bke::node_register_type(&ntype);
}

static void node_register_velocity_solve()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeSolveVelocityConstraints");
  ntype.ui_name = "Solve XPBD Velocity Constraints";
  ntype.ui_description =
      "Solve linear and angular velocity constraints on geometry using the XPBD framework";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  node_type_size(&ntype, 200, 120, 300);
  ntype.geometry_node_execute = node_geo_exec_velocities;
  ntype.declare = node_declare_velocities;
  blender::bke::node_register_type(&ntype);
}

static void node_register()
{
  node_register_position_solve();
  node_register_velocity_solve();
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_solve_xpbd_constraints_cc
