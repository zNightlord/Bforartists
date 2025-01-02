/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_mesh.hh"

#include "UI_interface.hh"
#include "UI_resources.hh"

#include "NOD_rna_define.hh"

#include "RNA_enum_types.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_set_mesh_normal_cc {

enum class Mode {
  Free = 0,
  SmoothFanSpace = 1,
};

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Geometry>("Mesh").supported_type({GeometryComponent::Type::Mesh});
  b.add_input<decl::Bool>("Selection").default_value(true).hide_value().field_on_all();
  b.add_input<decl::Vector>("Normal")
      .subtype(PROP_XYZ)
      .implicit_field(nodes::implicit_field_inputs::normal)
      .hide_value();
  b.add_output<decl::Geometry>("Mesh").propagate_all();
}

static void node_layout(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  const bNode &node = *static_cast<const bNode *>(ptr->data);
  uiItemR(layout, ptr, "mode", UI_ITEM_NONE, "", ICON_NONE);
}

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->custom1 = int16_t(Mode::Free);
}

static void set_mesh_normal(Mesh &mesh, const Mode mode, const Field<float3> &custom_normal) {}

static void node_geo_exec(GeoNodeExecParams params)
{
  const bNode &node = params.node();
  const Mode mode = static_cast<Mode>(node.custom1);
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Mesh");
  fn::Field<float3> custom_normal = params.extract_input<fn::Field<float3>>("Normal");

  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    if (Mesh *mesh = geometry_set.get_mesh_for_write()) {
      switch (mode) {
        case Mode::Free: {
          const bke::AttrDomain domain = bke::AttrDomain(node.custom2);
          bke::try_capture_field_on_geometry(mesh->attributes_for_write(),
                                             bke::MeshFieldContext(*mesh, domain),
                                             "custom_normal",
                                             domain,
                                             fn::make_constant_field(true),
                                             custom_normal);
          break;
        }
        case Mode::SmoothFanSpace: {
          BKE_mesh_set_custom_normals_normalized(mesh, ) break;
        }
      }
    }
  });

  params.set_output("Curve", std::move(geometry_set));
}

static void node_rna(StructRNA *srna)
{
  RNA_def_node_enum(srna,
                    "mode",
                    "Mode",
                    "Mode for curve normal evaluation",
                    rna_enum_curve_normal_mode_items,
                    NOD_inline_enum_accessors(custom1));
}

static void node_register()
{
  static blender::bke::bNodeType ntype;
  geo_node_type_base(&ntype, GEO_NODE_SET_MESH_NORMAL, "Set Mesh Normal", NODE_CLASS_GEOMETRY);
  ntype.declare = node_declare;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.initfunc = node_init;
  ntype.draw_buttons = node_layout;

  blender::bke::node_register_type(&ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_set_mesh_normal_cc
