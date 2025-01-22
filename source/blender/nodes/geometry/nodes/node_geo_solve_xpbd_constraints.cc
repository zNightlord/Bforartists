/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_solve_xpbd_constraints_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  constexpr float default_fps = 1.0f / 25.0f;

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).hide_value();
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

static void node_geo_exec(GeoNodeExecParams params)
{
  // const float delta_time =

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
