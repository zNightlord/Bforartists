/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "DNA_mesh_types.h"

#include "NOD_rna_define.hh"

#include "GEO_mesh_bevel.hh"

#include "UI_interface.hh"
#include "UI_resources.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_bevel_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Geometry>("Mesh").supported_type(GeometryComponent::Type::Mesh);
  b.add_input<decl::Bool>("Selection").default_value(true).field_on_all().hide_value();
  b.add_input<decl::Float>("Offset").default_value(0.0f).min(0.0f);
  b.add_input<decl::Int>("Segments").default_value(1);
  b.add_output<decl::Geometry>("Mesh").propagate_all();
}

static void node_layout(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiItemR(layout, ptr, "offset", UI_ITEM_NONE, "", ICON_NONE);
  uiItemR(layout, ptr, "affect", UI_ITEM_NONE, "", ICON_NONE);
}

static void geo_bevel_init(bNodeTree * /*tree*/, bNode *node)
{
  /* Use node->custom1 for affect. */
  node->custom1 = int(geometry::BevelAffect::Edges);
}

static void node_geo_exec(GeoNodeExecParams params)
{
  GeometrySet geometry_set = params.extract_input<GeometrySet>("Mesh");
  Field<bool> selection_field = params.extract_input<Field<bool>>("Selection");
  const AttributeFilter &attribute_filter = params.get_attribute_filter("Mesh");
  float offset = params.extract_input<float>("Offset");
  int segments = params.extract_input<int>("Segments");

  geometry::BevelAffect affect = geometry::BevelAffect(params.node().custom1);

  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    const Mesh *src_mesh = geometry_set.get_mesh();
    if (!src_mesh) {
      return;
    }

    const bke::MeshFieldContext context(*src_mesh, AttrDomain::Edge);
    FieldEvaluator evaluator{context, src_mesh->edges_num};
    evaluator.add(selection_field);
    evaluator.evaluate();
    const IndexMask selection = evaluator.get_evaluated_as_mask(0);
    if (selection.is_empty()) {
      return;
    }

    geometry::BevelParameters bevel_params;
    bevel_params.offset = offset;
    bevel_params.affect_type = affect;
    bevel_params.segments = segments;
    std::optional<Mesh *> mesh = geometry::mesh_bevel(
        *src_mesh, selection, bevel_params, attribute_filter);
    if (!mesh) {
      return;
    }

    geometry_set.replace_mesh(*mesh);
  });

  params.set_output("Mesh", std::move(geometry_set));
}

static void node_rna(StructRNA *srna)
{
  static const EnumPropertyItem rna_node_geometry_bevel_affect_items[] = {
      {int(geometry::BevelAffect::Vertices), "VERTICES", 0, "Vertices", "Bevel selected vertices"},
      {int(geometry::BevelAffect::Edges), "EDGES", 0, "Edges", "Bevel selected edges"},
      {0, nullptr, 0, nullptr, nullptr},
  };

  RNA_def_node_enum(srna,
                    "affect",
                    "Affect Kind",
                    "Which mesh element type to bevel",
                    rna_node_geometry_bevel_affect_items,
                    NOD_inline_enum_accessors(custom3),
                    std::optional<int>(int(geometry::BevelAffect::Edges)),
                    nullptr,
                    true);

  RNA_def_float(srna, "offset", 0.0f, 0.0f, 1e10f, "Offset", "How much to bevel", 0.0f, 10.0f);
}

static void node_register()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeMeshBevel", std::nullopt);
  ntype.ui_name = "Mesh Bevel";
  ntype.ui_description = "Bevel selected edges or vertices";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  ntype.draw_buttons = node_layout;
  ntype.declare = node_declare;
  ntype.initfunc = geo_bevel_init;
  ntype.geometry_node_execute = node_geo_exec;
  blender::bke::node_register_type(ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_bevel_cc
