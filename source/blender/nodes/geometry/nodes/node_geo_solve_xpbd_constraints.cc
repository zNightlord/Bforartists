/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_array_utils.hh"

#include "BKE_attribute.hh"
#include "BKE_geometry_set.hh"
#include "BKE_instances.hh"

#include "node_geometry_util.hh"

#include <fmt/format.h>

namespace blender::nodes::node_geo_solve_xpbd_constraints_cc {

enum class ConstraintType {
  /* Set position of a point to a target vector. */
  PositionGoal,
  /* Set orientation of an edge to a target rotation. */
  RotationGoal,
  /* Enforces edge length and aligns forward direction with the edge vector. */
  StretchShear,
  /* Enforces angles between neighboring edges to their relative rest orientation. */
  BendTwist,
  /* Keep contact points from penetrating. */
  Contact,
};
constexpr int NumConstraintTypes = 5;

enum class EvaluationTarget {
  Positions,
  Velocities,
};

enum class SolverMethod {
  GaussSeidel,
  Jacobi,
};

/* Constraint attributes. */
constexpr StringRef ATTR_SOLVER_GROUP = "solver_group";
constexpr StringRef ATTR_POINT1 = "point1";
constexpr StringRef ATTR_POINT2 = "point2";
constexpr StringRef ATTR_ACTIVE = "active";

static StringRef constraint_ui_geometry(const ConstraintType type)
{
  switch (type) {
    case ConstraintType::PositionGoal:
      return "Position Goal Constraints";
    case ConstraintType::RotationGoal:
      return "Rotation Goal Constraints";
    case ConstraintType::StretchShear:
      return "Stretch/Shear Constraints";
    case ConstraintType::BendTwist:
      return "Bend/Twist Constraints";
    case ConstraintType::Contact:
      return "Contact Constraints";
  }
  BLI_assert_unreachable();
  return "";
}

static StringRef constraint_ui_description(const ConstraintType type)
{
  switch (type) {
    case ConstraintType::PositionGoal:
      return "Set position of a point to a target vector";
    case ConstraintType::RotationGoal:
      return "Set orientation of an edge to a target rotation";
    case ConstraintType::StretchShear:
      return "Enforces edge length and aligns forward direction with the edge vector";
    case ConstraintType::BendTwist:
      return "Enforces angles between neighboring edges to their relative rest orientation";
    case ConstraintType::Contact:
      return "Keep contact points from penetrating";
  }
  BLI_assert_unreachable();
  return "";
}

constexpr float default_fps = 1.0f / 25.0f;

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

  auto add_constraint_input_output = [&](const ConstraintType type) {
    const StringRef name = constraint_ui_geometry(type);
    const StringRef description = constraint_ui_description(type);

    b.add_input<decl::Geometry>(name)
        .supported_type(GeometryComponent::Type::PointCloud)
        .description(description);
    b.add_output<decl::Geometry>(name).description(description).align_with_previous();
  };

  add_constraint_input_output(ConstraintType::PositionGoal);
  add_constraint_input_output(ConstraintType::RotationGoal);
  add_constraint_input_output(ConstraintType::StretchShear);
  add_constraint_input_output(ConstraintType::BendTwist);
  add_constraint_input_output(ConstraintType::Contact);
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

  b.add_input<decl::Vector>("Velocity").field_on({1}).hide_value();
  b.add_output<decl::Vector>("Velocity").field_on({0}).align_with_previous();
  b.add_input<decl::Vector>("Angular Velocity").field_on({1}).hide_value();
  b.add_output<decl::Vector>("Angular Velocity").field_on({0}).align_with_previous();
  b.add_input<decl::Vector>("Original Velocity").field_on({1}).hide_value();
  b.add_input<decl::Vector>("Original Angular Velocity").field_on({1}).hide_value();

  auto add_constraint_input_output = [&](const ConstraintType type) {
    const StringRef name = constraint_ui_geometry(type);
    const StringRef description = constraint_ui_description(type);

    b.add_input<decl::Geometry>(name)
        .supported_type(GeometryComponent::Type::PointCloud)
        .description(description);
    b.add_output<decl::Geometry>(name).description(description).align_with_previous();
  };

  add_constraint_input_output(ConstraintType::PositionGoal);
  add_constraint_input_output(ConstraintType::RotationGoal);
  add_constraint_input_output(ConstraintType::StretchShear);
  add_constraint_input_output(ConstraintType::BendTwist);
  add_constraint_input_output(ConstraintType::Contact);
  b.add_input<decl::Geometry>("Colliders")
      .only_instances()
      .description("Instances of colliders to evaluate contact transforms");

  b.add_input<decl::Bool>("Debug Checks")
      .default_value(false)
      .description("Perform checks on input data, which can impact performance");
}

template<typename T>
static AttributeReader<T> lookup_or_warn(const GeoNodeExecParams &params,
                                         AttributeAccessor &attributes,
                                         const StringRef attribute_id,
                                         const AttrDomain domain,
                                         const T &default_value)
{
  if (!attributes.contains(attribute_id)) {
    params.error_message_add(geo_eval_log::NodeWarningType::Warning,
                             fmt::format("Missing \"{}\" attribute", attribute_id));
  }
  return attributes.lookup_or_default<T>(attribute_id, domain, default_value);
}

struct ConstraintEvalInputs {
  /* TODO split this by EvaluationTarget, only either (positions + rotations) or (velocities +
   * angular_velocities) are used. */
  VArraySpan<float3> positions;
  VArraySpan<math::Quaternion> rotations;
  Array<float3> positions_buffer;
  Array<math::Quaternion> rotations_buffer;
  VArraySpan<float3> old_positions;
  VArraySpan<math::Quaternion> old_rotations;

  VArraySpan<float3> velocities;
  VArraySpan<float3> angular_velocities;
  Array<float3> velocities_buffer;
  Array<float3> angular_velocities_buffer;
  VArraySpan<float3> orig_velocities;
  VArraySpan<float3> orig_angular_velocities;

  Span<float4x4> collider_transforms;
  /* Collider transforms of the previous frame for computing velocity constraints. */
  Span<float4x4> old_collider_transforms;

  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<float3> goals;
  } position_goal;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<math::Quaternion> goals;
  } rotation_goal;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points1;
    VArraySpan<int> points2;
    VArraySpan<float> edge_length;
  } stretch_shear;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points1;
    VArraySpan<int> points2;
    VArraySpan<float3> darboux_vector;
  } bend_twist;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<int> collider_index;
    VArraySpan<float3> local_position1;
    VArraySpan<float3> local_position2;
    VArraySpan<float3> normal;

    VArraySpan<float> friction;
    VArraySpan<float> restitution;
  } contact;
};

struct ConstraintEvalOutputs {
  Array<float3> positions;
  Array<math::Quaternion> rotations;
  Array<float3> velocities;
  Array<float3> angular_velocities;

  struct {
    /* Remember active contacts for later velocity update. */
    SpanAttributeWriter<bool> active;
  } contact;
};

static const VArray<int> &constraints_solver_groups(const ConstraintEvalInputs &inputs,
                                                    const ConstraintType type)
{
  switch (type) {
    case ConstraintType::PositionGoal:
      return inputs.position_goal.solver_group;
    case ConstraintType::RotationGoal:
      return inputs.rotation_goal.solver_group;
    case ConstraintType::StretchShear:
      return inputs.stretch_shear.solver_group;
    case ConstraintType::BendTwist:
      return inputs.bend_twist.solver_group;
    case ConstraintType::Contact:
      return inputs.contact.solver_group;
  }
  BLI_assert_unreachable();
  static const VArray<int> dummy;
  return dummy;
}

struct SolverParams {
  float delta_time;
  EvaluationTarget target;
  ConstraintEvalInputs &inputs;
  ConstraintEvalOutputs &outputs;

  std::function<void(const NodeWarningType type, const StringRef message)> error_message_add;
  /* Perform debug checks on user inputs at runtime. This helps avoid common errors but has a
   * significant performance cost, so should be optional. */
  bool debug_check;
};

void evaluate_constraint_group_position_goal(const SolverParams &params,
                                             const IndexMask &group_mask)
{
  switch (params.target) {
    case EvaluationTarget::Positions:
      group_mask.foreach_index(GrainSize(1024), [&](const int index) {
        const int point1 = params.inputs.position_goal.points[index];
        const float3 position1 = params.inputs.positions[point1];
        const float3 goal = params.inputs.position_goal.goals[index];

        const float3 residual = position1 - goal;
        /* Gradient is identity transform. */
        const float3 delta_position = -residual;

        params.outputs.positions[point1] = position1 + delta_position;
      });
      break;
    case EvaluationTarget::Velocities:
      break;
  }
}

void evaluate_constraint_group_rotation_goal(const SolverParams &params,
                                             const IndexMask &group_mask)
{
  switch (params.target) {
    case EvaluationTarget::Positions:
      break;
    case EvaluationTarget::Velocities:
      break;
  }
}

void evaluate_constraint_group_stretch_shear(const SolverParams &params,
                                             const IndexMask &group_mask)
{
  switch (params.target) {
    case EvaluationTarget::Positions:
      break;
    case EvaluationTarget::Velocities:
      break;
  }
}

void evaluate_constraint_group_bend_twist(const SolverParams &params, const IndexMask &group_mask)
{
  switch (params.target) {
    case EvaluationTarget::Positions:
      break;
    case EvaluationTarget::Velocities:
      break;
  }
}

template<bool debug_check>
void evaluate_constraint_group_contact(const SolverParams &params, const IndexMask &group_mask)
{
  switch (params.target) {
    case EvaluationTarget::Positions: {
      group_mask.foreach_index(GrainSize(1024), [&](const int index) {
        const int point1 = params.inputs.contact.points[index];
        const int collider_index = params.inputs.contact.collider_index[index];
        if (!params.inputs.collider_transforms.index_range().contains(collider_index)) {
          return;
        };

        /* Local positions are relative to moving point and collider respectively. */
        const float3 &local_position1 = params.inputs.contact.local_position1[index];
        const float3 &local_position2 = params.inputs.contact.local_position2[index];
        /* Normal is a fixed shared direction for both participants. */
        const float3 &normal = params.inputs.contact.normal[index];
        if constexpr (debug_check) {
          if (!math::is_unit(normal)) {
            params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                     "Contact normal vector not normalized");
          }
        }

        /* Construct transforms for both the point (translation only for now) as well as the
         * collider. */
        const float3 position1 = params.inputs.positions[point1];
        const float4x4 point_transform = math::from_location<float4x4>(position1);
        const float4x4 &collider_transform = params.inputs.collider_transforms[collider_index];
        /* Contact points are computed by applying the transforms to relative local positions. */
        const float3 contact_point1 = math::transform_point(point_transform, local_position1);
        const float3 contact_point2 = math::transform_point(collider_transform, local_position2);

        /* Positional constraint for penetration depth along the normal. */
        const float residual_depth = math::dot(contact_point1 - contact_point2, normal);
        /* Only act on contact. */
        const bool active = residual_depth < 0.0f;
        params.outputs.contact.active.span[index] = active;
        if (!active) {
          return;
        }

        /* Gradient is normal, length is 1, no need to compute gradient norm. */
        const float3 delta_position = -residual_depth * normal;
        params.outputs.positions[point1] += delta_position;
      });
      break;
    }
    case EvaluationTarget::Velocities: {
      group_mask.foreach_index(GrainSize(1024), [&](const int index) {
        /* Active status is determined by the position evaluation. */
        if (!params.outputs.contact.active.span[index]) {
          return;
        }

        const int point1 = params.inputs.contact.points[index];
        const int collider_index = params.inputs.contact.collider_index[index];
        if (!params.inputs.collider_transforms.index_range().contains(collider_index)) {
          return;
        };

        /* Local positions are relative to moving point and collider respectively. */
        const float3 &local_position2 = params.inputs.contact.local_position2[index];
        /* Normal is a fixed shared direction for both participants. */
        const float3 &normal = params.inputs.contact.normal[index];
        if constexpr (debug_check) {
          if (!math::is_unit(normal)) {
            params.error_message_add(geo_eval_log::NodeWarningType::Error,
                                     "Contact normal vector not normalized");
          }
        }
        const float restitution = params.inputs.contact.restitution[index];
        const float friction = params.inputs.contact.friction[index];

        const float4x4 &collider_transform = params.inputs.collider_transforms[collider_index];
        const float4x4 &old_collider_transform =
            params.inputs.old_collider_transforms[collider_index];
        /* Contact points are computed by applying the transforms to relative local positions. */
        const float3 contact_point2 = math::transform_point(collider_transform, local_position2);
        const float3 old_contact_point2 = math::transform_point(old_collider_transform,
                                                                local_position2);

        const float3 velocity1 = params.inputs.orig_velocities[point1];
        /* Compute velocity of the collider contact point. */
        const float3 velocity2 = contact_point2 - old_contact_point2;
        const float3 relative_velocity = velocity1 - velocity2;
        const float normal_relative_velocity = math::dot(relative_velocity, normal);
        const float3 surface_relative_velocity = relative_velocity -
                                                 normal * normal_relative_velocity;

        const float residual_restitution = normal_relative_velocity;
        /* Folded into delta. */
        /* const float residual_friction = math::length(surface_relative_velocity); */
        /* Gradients are normalized, no need to compute gradient norm. */
        const float3 delta_velocity_restitution = (residual_restitution < 0.0f ?
                                                       -(1.0f + restitution) *
                                                           residual_restitution * normal :
                                                       float3(0.0f));
        const float3 delta_velocity_friction = -friction * surface_relative_velocity;

        params.outputs.velocities[point1] += delta_velocity_restitution + delta_velocity_friction;
      });
      break;
    }
  }
}

/* Evaluate a group of constraints in parallel.
 * It's important that none of the constraints write to the same variables. Solver groups must be
 * computed in such a way that each variable is only affected by one constraint in each group.
 */
template<bool debug_check>
static void evaluate_constraint_group(const SolverParams &params,
                                      const ConstraintType type,
                                      const IndexMask &group_mask)
{
  BLI_assert(!constraints_solver_groups(params.inputs, type).is_empty());

  switch (type) {
    case ConstraintType::PositionGoal:
      evaluate_constraint_group_position_goal(params, group_mask);
      break;
    case ConstraintType::RotationGoal:
      evaluate_constraint_group_rotation_goal(params, group_mask);
      break;
    case ConstraintType::StretchShear:
      evaluate_constraint_group_stretch_shear(params, group_mask);
      break;
    case ConstraintType::BendTwist:
      evaluate_constraint_group_bend_twist(params, group_mask);
      break;
    case ConstraintType::Contact:
      evaluate_constraint_group_contact<debug_check>(params, group_mask);
      break;
  }

  switch (params.target) {
    case EvaluationTarget::Positions:
      /* Move output array to input buffer for next update. */
      params.inputs.positions_buffer = params.outputs.positions;
      params.inputs.rotations_buffer = params.outputs.rotations;
      params.inputs.positions = VArray<float3>::ForContainer(params.inputs.positions_buffer);
      params.inputs.rotations = VArray<math::Quaternion>::ForContainer(
          params.inputs.rotations_buffer);
      break;
    case EvaluationTarget::Velocities:
      /* Move output array to input buffer for next update. */
      params.inputs.velocities_buffer = params.outputs.velocities;
      params.inputs.angular_velocities_buffer = params.outputs.angular_velocities;
      params.inputs.velocities = VArray<float3>::ForContainer(params.inputs.velocities_buffer);
      params.inputs.angular_velocities = VArray<float3>::ForContainer(
          params.inputs.angular_velocities_buffer);
      break;
  }
}

static void do_single_constraint_passes(const SolverParams &params, const ConstraintType type)
{
  const VArray<int> solver_groups = constraints_solver_groups(params.inputs, type);
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

    if (params.debug_check) {
      evaluate_constraint_group<true>(params, type, group_mask);
    }
    else {
      evaluate_constraint_group<false>(params, type, group_mask);
    }
  }
}

static void do_gauss_seidel_step(const SolverParams &params)
{
  /* Order of constraint passes is chosen by increasing "importance":
   * Later constraints have less residual error, and the last constraint type is solved exactly.
   */
  do_single_constraint_passes(params, ConstraintType::BendTwist);
  do_single_constraint_passes(params, ConstraintType::StretchShear);
  do_single_constraint_passes(params, ConstraintType::RotationGoal);
  do_single_constraint_passes(params, ConstraintType::PositionGoal);
  do_single_constraint_passes(params, ConstraintType::Contact);
}

static void do_jacobi_step(const SolverParams & /*params*/)
{
  /* TODO */
}

static void do_solver_steps(const SolverMethod method, const int steps, const SolverParams &params)
{
  if (steps < 1) {
    return;
  }

  for ([[maybe_unused]] const int step : IndexRange(steps)) {
    switch (method) {
      case SolverMethod::GaussSeidel:
        do_gauss_seidel_step(params);
        break;
      case SolverMethod::Jacobi:
        do_jacobi_step(params);
        break;
    }
  }
}

static void prepare_constraint_data(GeoNodeExecParams params,
                                    Vector<GeometrySet> &constraint_geometry_sets,
                                    ConstraintEvalInputs &constraint_inputs,
                                    ConstraintEvalOutputs &constraint_outputs)
{
  constraint_geometry_sets.resize(NumConstraintTypes);
  Array<std::optional<MutableAttributeAccessor>> constraint_attributes(NumConstraintTypes,
                                                                       std::nullopt);
  auto extract_constraint_attributes = [&](const ConstraintType type) {
    BLI_assert(constraint_geometry_sets.index_range().contains(int(type)));
    BLI_assert(constraint_attributes.index_range().contains(int(type)));
    constraint_geometry_sets[int(type)] = params.extract_input<GeometrySet>(
        constraint_ui_geometry(type));
    GeometrySet &geometry_set = constraint_geometry_sets[int(type)];
    if (geometry_set.has_component<PointCloudComponent>()) {
      PointCloudComponent &constraint_component =
          geometry_set.get_component_for_write<PointCloudComponent>();
      constraint_attributes[int(type)] = constraint_component.attributes_for_write();
    }
  };
  extract_constraint_attributes(ConstraintType::PositionGoal);
  extract_constraint_attributes(ConstraintType::RotationGoal);
  extract_constraint_attributes(ConstraintType::StretchShear);
  extract_constraint_attributes(ConstraintType::BendTwist);
  extract_constraint_attributes(ConstraintType::Contact);

  if (std::optional<MutableAttributeAccessor> attributes =
          constraint_attributes[int(ConstraintType::PositionGoal)])
  {
    constraint_inputs.position_goal.solver_group = *lookup_or_warn<int>(
        params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
    constraint_inputs.position_goal.points = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
    constraint_inputs.position_goal.goals = *lookup_or_warn<float3>(
        params, *attributes, "goal", AttrDomain::Point, float3(0.0f));
  }
  if (std::optional<MutableAttributeAccessor> attributes =
          constraint_attributes[int(ConstraintType::RotationGoal)])
  {
    constraint_inputs.rotation_goal.solver_group = *lookup_or_warn<int>(
        params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
    constraint_inputs.rotation_goal.points = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
    constraint_inputs.rotation_goal.goals = *lookup_or_warn<math::Quaternion>(
        params, *attributes, "goal", AttrDomain::Point, math::Quaternion::identity());
  }
  if (std::optional<MutableAttributeAccessor> attributes =
          constraint_attributes[int(ConstraintType::StretchShear)])
  {
    constraint_inputs.stretch_shear.solver_group = *lookup_or_warn<int>(
        params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
    constraint_inputs.stretch_shear.points1 = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
    constraint_inputs.stretch_shear.points2 = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT2, AttrDomain::Point, 0);
    constraint_inputs.stretch_shear.edge_length = *lookup_or_warn<float>(
        params, *attributes, "edge_length", AttrDomain::Point, 0.0f);
  }
  if (std::optional<MutableAttributeAccessor> attributes =
          constraint_attributes[int(ConstraintType::BendTwist)])
  {
    constraint_inputs.bend_twist.solver_group = *lookup_or_warn<int>(
        params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
    constraint_inputs.bend_twist.points1 = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
    constraint_inputs.bend_twist.points2 = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT2, AttrDomain::Point, 0);
    constraint_inputs.bend_twist.darboux_vector = *lookup_or_warn<float3>(
        params, *attributes, "darboux_vector", AttrDomain::Point, float3(0.0f));
  }
  if (std::optional<MutableAttributeAccessor> attributes =
          constraint_attributes[int(ConstraintType::Contact)])
  {
    constraint_inputs.contact.solver_group = *lookup_or_warn<int>(
        params, *attributes, ATTR_SOLVER_GROUP, AttrDomain::Point, 0);
    constraint_inputs.contact.points = *lookup_or_warn<int>(
        params, *attributes, ATTR_POINT1, AttrDomain::Point, 0);
    constraint_inputs.contact.collider_index = *lookup_or_warn<int>(
        params, *attributes, "collider_index", AttrDomain::Point, 0);
    constraint_inputs.contact.local_position1 = *lookup_or_warn<float3>(
        params, *attributes, "local_position1", AttrDomain::Point, float3(0.0f));
    constraint_inputs.contact.local_position2 = *lookup_or_warn<float3>(
        params, *attributes, "local_position2", AttrDomain::Point, float3(0.0f));
    constraint_inputs.contact.normal = *lookup_or_warn<float3>(
        params, *attributes, "normal", AttrDomain::Point, float3(0.0f));

    constraint_inputs.contact.friction = *attributes->lookup_or_default<float>(
        "friction", AttrDomain::Point, 0.0f);
    constraint_inputs.contact.restitution = *attributes->lookup_or_default<float>(
        "restitution", AttrDomain::Point, 0.0f);

    const GeometrySet collider_geometry = params.extract_input<GeometrySet>("Colliders");
    constraint_inputs.collider_transforms = collider_geometry.has_instances() ?
                                                collider_geometry.get_instances()->transforms() :
                                                Span<float4x4>{};
    /* XXX Transforms of the previous frame are not currently available, these are always the
     * same as the current frame. Eventually this will allow transfer of velocity from animated
     * colliders. */
    constraint_inputs.old_collider_transforms = constraint_inputs.collider_transforms;

    const int num_constraints = attributes->domain_size(AttrDomain::Point);
    constraint_outputs.contact.active = attributes->lookup_or_add_for_write_span<bool>(
        ATTR_ACTIVE,
        AttrDomain::Point,
        bke::AttributeInitVArray(VArray<bool>::ForSingle(false, num_constraints)));
  }
}

static void write_constraint_attributes(ConstraintEvalOutputs &constraint_outputs)
{
  if (constraint_outputs.contact.active) {
    constraint_outputs.contact.active.finish();
  }
}

static void node_geo_exec_positions(GeoNodeExecParams params)
{
  const int gauss_seidel_steps = params.extract_input<int>("Gauss-Seidel Steps");
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
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

  Vector<GeometrySet> constraint_geometry_sets;
  ConstraintEvalInputs constraint_inputs;
  ConstraintEvalOutputs constraint_outputs;
  prepare_constraint_data(params, constraint_geometry_sets, constraint_inputs, constraint_outputs);

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

        const bke::GeometryFieldContext field_context{component, AttrDomain::Point};
        fn::FieldEvaluator evaluator{field_context, num_points};
        evaluator.add(position_field);
        evaluator.add(rotation_field);
        evaluator.add(old_position_field);
        evaluator.add(old_rotation_field);
        evaluator.evaluate();
        constraint_inputs.positions = evaluator.get_evaluated<float3>(0);
        constraint_inputs.rotations = evaluator.get_evaluated<math::Quaternion>(1);
        constraint_inputs.old_positions = evaluator.get_evaluated<float3>(2);
        constraint_inputs.old_rotations = evaluator.get_evaluated<math::Quaternion>(3);

        /* Full input copy to initialize outputs, since not all variables may be changed. */
        constraint_outputs.positions.reinitialize(num_points);
        constraint_outputs.rotations.reinitialize(num_points);
        array_utils::copy(constraint_inputs.positions,
                          constraint_outputs.positions.as_mutable_span());
        array_utils::copy(constraint_inputs.rotations,
                          constraint_outputs.rotations.as_mutable_span());

        SolverParams solver_params = {
            delta_time, EvaluationTarget::Positions, constraint_inputs, constraint_outputs};
        solver_params.error_message_add = [params](const NodeWarningType type,
                                                   const StringRef message) {
          params.error_message_add(type, message);
        };
        solver_params.debug_check = debug_check;

        do_solver_steps(SolverMethod::GaussSeidel, gauss_seidel_steps, solver_params);

        if (position_output_id) {
          AttributeWriter<float3> positions_writer = attributes->lookup_or_add_for_write<float3>(
              *position_output_id, AttrDomain::Point);
          BLI_assert(constraint_inputs.positions.size() == num_points);
          positions_writer.varray.set_all(constraint_inputs.positions);
          positions_writer.finish();
        }
        if (rotation_output_id) {
          AttributeWriter<math::Quaternion> rotations_writer =
              attributes->lookup_or_add_for_write<math::Quaternion>(*rotation_output_id,
                                                                    AttrDomain::Point);
          BLI_assert(constraint_inputs.rotations.size() == num_points);
          rotations_writer.varray.set_all(constraint_inputs.rotations);
          rotations_writer.finish();
        }
      }
    }
  });

  write_constraint_attributes(constraint_outputs);

  params.set_output("Geometry", geometry_set);

  auto set_constraint_outputs = [&](const ConstraintType type) {
    params.set_output(constraint_ui_geometry(type), constraint_geometry_sets[int(type)]);
  };
  set_constraint_outputs(ConstraintType::PositionGoal);
  set_constraint_outputs(ConstraintType::RotationGoal);
  set_constraint_outputs(ConstraintType::StretchShear);
  set_constraint_outputs(ConstraintType::BendTwist);
  set_constraint_outputs(ConstraintType::Contact);
}

static void node_geo_exec_velocities(GeoNodeExecParams params)
{
  constexpr int gauss_seidel_steps = 1;
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
  const bool debug_check = params.extract_input<bool>("Debug Checks");
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");
  Field<float3> velocity_field = params.extract_input<Field<float3>>("Velocity");
  Field<float3> angular_velocity_field = params.extract_input<Field<float3>>("Angular Velocity");
  std::optional<std::string> velocity_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Velocity");
  std::optional<std::string> angular_velocity_output_id =
      params.get_output_anonymous_attribute_id_if_needed("Angular Velocity");
  Field<float3> orig_velocity_field = params.extract_input<Field<float3>>("Original Velocity");
  Field<float3> orig_angular_velocity_field = params.extract_input<Field<float3>>(
      "Original Angular Velocity");

  Vector<GeometrySet> constraint_geometry_sets;
  ConstraintEvalInputs constraint_inputs;
  ConstraintEvalOutputs constraint_outputs;
  prepare_constraint_data(params, constraint_geometry_sets, constraint_inputs, constraint_outputs);

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

        const bke::GeometryFieldContext field_context{component, AttrDomain::Point};
        fn::FieldEvaluator evaluator{field_context, num_points};
        evaluator.add(velocity_field);
        evaluator.add(angular_velocity_field);
        evaluator.add(orig_velocity_field);
        evaluator.add(orig_angular_velocity_field);
        evaluator.evaluate();
        constraint_inputs.velocities = evaluator.get_evaluated<float3>(0);
        constraint_inputs.angular_velocities = evaluator.get_evaluated<float3>(1);
        constraint_inputs.orig_velocities = evaluator.get_evaluated<float3>(2);
        constraint_inputs.orig_angular_velocities = evaluator.get_evaluated<float3>(3);

        /* Full input copy to initialize outputs, since not all variables may be changed. */
        constraint_outputs.velocities.reinitialize(num_points);
        constraint_outputs.angular_velocities.reinitialize(num_points);
        array_utils::copy(constraint_inputs.velocities,
                          constraint_outputs.velocities.as_mutable_span());
        array_utils::copy(constraint_inputs.angular_velocities,
                          constraint_outputs.angular_velocities.as_mutable_span());

        SolverParams solver_params = {
            delta_time, EvaluationTarget::Velocities, constraint_inputs, constraint_outputs};
        solver_params.error_message_add = [params](const NodeWarningType type,
                                                   const StringRef message) {
          params.error_message_add(type, message);
        };
        solver_params.debug_check = debug_check;

        do_solver_steps(SolverMethod::GaussSeidel, gauss_seidel_steps, solver_params);

        if (velocity_output_id) {
          AttributeWriter<float3> velocities_writer = attributes->lookup_or_add_for_write<float3>(
              *velocity_output_id, AttrDomain::Point);
          BLI_assert(constraint_inputs.velocities.size() == num_points);
          velocities_writer.varray.set_all(constraint_inputs.velocities);
          velocities_writer.finish();
        }
        if (angular_velocity_output_id && !constraint_inputs.angular_velocities_buffer.is_empty())
        {
          AttributeWriter<float3> angular_velocities_writer =
              attributes->lookup_or_add_for_write<float3>(*angular_velocity_output_id,
                                                          AttrDomain::Point);
          BLI_assert(constraint_inputs.angular_velocities.size() == num_points);
          angular_velocities_writer.varray.set_all(constraint_inputs.angular_velocities);
          angular_velocities_writer.finish();
        }
      }
    }
  });

  params.set_output("Geometry", geometry_set);

  auto set_constraint_outputs = [&](const ConstraintType type) {
    params.set_output(constraint_ui_geometry(type), constraint_geometry_sets[int(type)]);
  };
  set_constraint_outputs(ConstraintType::PositionGoal);
  set_constraint_outputs(ConstraintType::RotationGoal);
  set_constraint_outputs(ConstraintType::StretchShear);
  set_constraint_outputs(ConstraintType::BendTwist);
  set_constraint_outputs(ConstraintType::Contact);
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
