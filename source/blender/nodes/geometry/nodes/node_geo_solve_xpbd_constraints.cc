/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_attribute.hh"
#include "BKE_geometry_set.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_solve_xpbd_constraints_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  constexpr float default_fps = 1.0f / 25.0f;

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();
  b.add_input<decl::Int>("Gauss-Seidel Steps").default_value(1).min(0);

  b.add_input<decl::Geometry>("Geometry");
  b.add_input<decl::Geometry>("Stretch/Shear Constraints");
  b.add_input<decl::Geometry>("Bend/Twist Constraints");
  b.add_input<decl::Geometry>("Goal Constraints");
  b.add_input<decl::Geometry>("Contact Constraints");
  b.add_output<decl::Geometry>("Geometry").propagate_all();
  b.add_output<decl::Geometry>("Bend/Twist Constraints").propagate_all();
  b.add_output<decl::Geometry>("Goal Constraints").propagate_all();
  b.add_input<decl::Geometry>("Contact Constraints");
}

enum class SolverMethod {
  GaussSeidel,
  Jacobi,
};

struct Params {
  float delta_time;
  bke::MutableAttributeAccessor geometry_attributes;
  std::optional<bke::MutableAttributeAccessor> stretch_shear_constraint_attributes;
  std::optional<bke::MutableAttributeAccessor> bend_twist_constraint_attributes;
  std::optional<bke::MutableAttributeAccessor> goal_constraint_attributes;
  std::optional<bke::MutableAttributeAccessor> contact_constraint_attributes;
};

static void do_gauss_seidel_step(const Params &params) {}

static void do_jacobi_step(const Params &params)
{
  /* TODO */
}

static void do_solver_steps(const SolverMethod method, const int steps, const Params &params)
{
  if (steps < 1) {
    return;
  }

  for (const int step : IndexRange(steps)) {
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
  GeometrySet stretch_shear_geometry = params.extract_input<GeometrySet>(
      "Stretch/Shear Constraints");
  GeometrySet bend_twist_geometry = params.extract_input<GeometrySet>("Bend/Twist Constraints");
  GeometrySet goal_geometry = params.extract_input<GeometrySet>("Goal Constraints");
  GeometrySet contact_geometry = params.extract_input<GeometrySet>("Contact Constraints");

  auto get_constraint_attributes =
      [&](GeometrySet &constraint_geometry) -> std::optional<bke::MutableAttributeAccessor> {
    if (constraint_geometry.has(bke::GeometryComponent::Type::PointCloud)) {
      return constraint_geometry.get_component_for_write(bke::GeometryComponent::Type::PointCloud)
          .attributes_for_write();
    }
    return std::nullopt;
  };

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

        const Params params = {
            delta_time,
            *attributes,
            get_constraint_attributes(stretch_shear_geometry),
            get_constraint_attributes(bend_twist_geometry),
            get_constraint_attributes(goal_geometry),
            get_constraint_attributes(contact_geometry),
        };
        do_solver_steps(SolverMethod::GaussSeidel, gauss_seidel_steps, params);
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
