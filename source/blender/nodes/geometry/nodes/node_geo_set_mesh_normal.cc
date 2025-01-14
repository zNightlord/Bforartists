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
  CornerFanSpace = 1,
};

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Geometry>("Mesh").supported_type(GeometryComponent::Type::Mesh);
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
  if (Mode(node.custom1) == Mode::Free) {
    uiItemR(layout, ptr, "domain", UI_ITEM_NONE, "", ICON_NONE);
  }
}

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->custom1 = int16_t(Mode::Free);
  node->custom2 = int16_t(bke::AttrDomain::Point);
}

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
        case Mode::CornerFanSpace: {
          const bke::MeshFieldContext context(*mesh, bke::AttrDomain::Corner);
          fn::FieldEvaluator evaluator(context, mesh->corners_num);
          Array<float3> corner_normals(mesh->corners_num);
          evaluator.add_with_destination<float3>(custom_normal, corner_normals);
          evaluator.evaluate();
          bke::mesh_set_custom_normals(*mesh, corner_normals);
          break;
        }
      }
    }
  });

  params.set_output("Mesh", std::move(geometry_set));
}

static void node_rna(StructRNA *srna)
{
  static const EnumPropertyItem mode_items[] = {
      {int(Mode::Free),
       "Free",
       0,
       "Free",
       "Store custom normals as simple vectors in the local space of the mesh. Values are not "
       "necessarily updated automatically later on as the mesh is deformed."},
      {int(Mode::CornerFanSpace),
       "CORNER_FAN_SPACE",
       0,
       "Corner Fan Space",
       "Store normals in a deformation-independent custom transformation space. This method is "
       "slower, but can be better when subsequent operations change the mesh without handling "
       "normals specifically."},
      {0, nullptr, 0, nullptr, nullptr},
  };

  RNA_def_node_enum(srna,
                    "mode",
                    "Mode",
                    "Storage mode for custom normal data",
                    mode_items,
                    NOD_inline_enum_accessors(custom1));
  RNA_def_node_enum(srna,
                    "domain",
                    "Domain",
                    "Attribute domain to store free custom normals",
                    rna_enum_attribute_domain_only_mesh_no_edge_items,
                    NOD_inline_enum_accessors(custom2));
}

static void node_register()
{
  static blender::bke::bNodeType ntype;
  geo_node_type_base(&ntype, "GeometryNodeSetMeshNormal");
  ntype.ui_name = "Set Mesh Normal";
  ntype.ui_description = "Store a normal vector for each mesh element";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  ntype.declare = node_declare;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.initfunc = node_init;
  ntype.draw_buttons = node_layout;

  blender::bke::node_register_type(&ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_set_mesh_normal_cc
