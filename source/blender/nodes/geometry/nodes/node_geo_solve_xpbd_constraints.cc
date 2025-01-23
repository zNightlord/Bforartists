/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_array_utils.hh"

#include "BKE_attribute.hh"
#include "BKE_geometry_set.hh"

#include "node_geometry_util.hh"

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

enum class SolverMethod {
  GaussSeidel,
  Jacobi,
};

const char *ATTR_SOLVER_GROUP = "solver_group";

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

static void node_declare(NodeDeclarationBuilder &b)
{
  constexpr float default_fps = 1.0f / 25.0f;

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();
  b.add_input<decl::Int>("Gauss-Seidel Steps").default_value(1).min(0);

  b.add_input<decl::Geometry>("Geometry");
  b.add_output<decl::Geometry>("Geometry").align_with_previous();

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
}

struct SolverParams {
  float delta_time;
  bke::MutableAttributeAccessor geometry_attributes;
  Vector<std::optional<bke::MutableAttributeAccessor>> constraint_attributes;

  std::function<void(const NodeWarningType type, const StringRef message)> error_message_add;
};

/* Evaluate a group of constraints in parallel.
 * It's important that none of the constraints write to the same variables. Solver groups must be
 * computed in such a way that each variable is only affected by one constraint in each group.
 */
static void evaluate_constraint_group(const IndexMask &group_mask)
{
  group_mask.foreach_index(GrainSize(1024), [&](const int index) {

  });
}

static void do_single_constraint_passes(const SolverParams &params, const ConstraintType type)
{
  if (!params.constraint_attributes[int(type)]) {
    return;
  }

  MutableAttributeAccessor constraint_attributes = *params.constraint_attributes[int(type)];
  const IndexRange constraints = IndexRange(constraint_attributes.domain_size(AttrDomain::Point));
  if (!constraint_attributes.contains(ATTR_SOLVER_GROUP)) {
    params.error_message_add(geo_eval_log::NodeWarningType::Warning,
                             "Missing 'solver_group' attribute");
  }
  VArray<int> solver_groups = *constraint_attributes.lookup_or_default<int>(
      ATTR_SOLVER_GROUP, AttrDomain::Point, 0);

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

    evaluate_constraint_group(group_mask);
  }
}

static void do_gauss_seidel_step(const SolverParams &params)
{
  /* Order of constraint passes is chosen by increasing "importance":
   * Later constraints have less residual error, and the last constraint type is solved exactly. */
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

static void node_geo_exec(GeoNodeExecParams params)
{
  const int gauss_seidel_steps = params.extract_input<int>("Gauss-Seidel Steps");
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");

  Array<PointCloudComponent *> constraint_components(5, nullptr);
  auto extract_constraint_component = [&](const ConstraintType type) {
    BLI_assert(constraint_components.index_range().contains(int(type)));
    GeometrySet constraint_geometry_set = params.extract_input<GeometrySet>(
        constraint_ui_geometry(type));
    if (constraint_geometry_set.is_empty()) {
      return;
    }
    constraint_components[int(type)] =
        &constraint_geometry_set.get_component_for_write<PointCloudComponent>();
  };
  extract_constraint_component(ConstraintType::PositionGoal);
  extract_constraint_component(ConstraintType::RotationGoal);
  extract_constraint_component(ConstraintType::StretchShear);
  extract_constraint_component(ConstraintType::BendTwist);
  extract_constraint_component(ConstraintType::Contact);

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

        SolverParams solver_params = {delta_time, *attributes};
        solver_params.constraint_attributes.reinitialize(constraint_components.size());
        for (const int i : constraint_components.index_range()) {
          solver_params.constraint_attributes[i] =
              constraint_components[i] ? constraint_components[i]->attributes_for_write() :
                                         std::nullopt;
        }
        solver_params.error_message_add = [params](const NodeWarningType type,
                                                   const StringRef message) {
          params.error_message_add(type, message);
        };

        do_solver_steps(SolverMethod::GaussSeidel, gauss_seidel_steps, solver_params);
      }
    }
  });

  params.set_default_remaining_outputs();
}

static void node_register()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeSolveXPBDConstraints");
  ntype.ui_name = "Solve XPBD Constraints";
  ntype.ui_description =
      "Solve position and rotation constraints on geometry using the XPBD framework";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.declare = node_declare;
  blender::bke::node_register_type(&ntype);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_solve_xpbd_constraints_cc
