/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_armature.hh"
#include "BKE_curves.hh"
#include "BKE_deform.hh"
#include "BKE_grease_pencil.hh"
#include "BKE_mesh.hh"

#include "NOD_geometry_nodes_closure.hh"
#include "NOD_geometry_nodes_closure_eval.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_armature_deform_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_input<decl::Geometry>("Geometry");
  b.add_output<decl::Geometry>("Geometry").align_with_previous();

  b.add_input<decl::Object>("Armature Object").hide_label();

  b.add_input<decl::Bool>("Preserve Volume")
      .default_value(false)
      .description(
          "Use dual quaternion deformation to better preserve the initial volume of the geometry");
  b.add_input<decl::Bool>("Use Envelope")
      .default_value(false)
      .description("Use bone envelope distance to compute deformation weights");
  b.add_input<decl::Bool>("Use Vertex Groups")
      .default_value(true)
      .description("Use vertex groups as deformation weights");
  b.add_input<decl::Bool>("Invert Vertex Groups")
      .default_value(false)
      .description("Invert vertex group weights");
  b.add_input<decl::Closure>("Custom Weights");

  b.add_input<decl::Float>("Mask").default_value(1.0f).hide_value().field_on_all().description(
      "Influence of the deformation for each point");

  b.add_output<decl::Matrix>("Deform Matrix")
      .field_on_all()
      .description("Local deformation gradient for each point");
}

static void node_geo_exec(GeoNodeExecParams params)
{
  const Object *armature_object = params.extract_input<Object *>("Armature Object");
  if (!armature_object || armature_object->type != OB_ARMATURE) {
    params.set_output("Geometry", params.extract_input<GeometrySet>("Geometry"));
    params.set_default_remaining_outputs();
    return;
  }
  const Object *self_object = params.self_object();
  BLI_assert(self_object != nullptr);
  const float4x4 &target_to_world = self_object->object_to_world();
  const ListBase *vertex_groups = self_object->data ?
                                      BKE_id_defgroup_list_get(
                                          static_cast<const ID *>(self_object->data)) :
                                      nullptr;

  const bool use_quaternion = params.extract_input<bool>("Preserve Volume");
  const bool use_envelope = params.extract_input<bool>("Use Envelope");
  const bool use_vertex_groups = params.extract_input<bool>("Use Vertex Groups");
  const bool invert_vertex_groups = params.extract_input<bool>("Invert Vertex Groups");
  const ClosurePtr custom_weights = params.extract_input<ClosurePtr>("Custom Weights");
  const Field<float> mask_field = params.extract_input<Field<float>>("Mask");

  GeometrySet geometry_set = params.extract_input<GeometrySet>("Geometry");
  /* Deform matrix is only computed if necessary. */
  const std::optional<std::string> deform_matrix_id =
      params.get_output_anonymous_attribute_id_if_needed("Deform Matrix");

  static const Array<GeometryComponent::Type> types = {GeometryComponent::Type::Mesh,
                                                       GeometryComponent::Type::PointCloud,
                                                       GeometryComponent::Type::Curve,
                                                       GeometryComponent::Type::GreasePencil};
  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    if (geometry_set.has_mesh()) {
      Mesh &mesh = *geometry_set.get_mesh_for_write();
      std::optional<ArmatureDeformVertexGroupParams> vgroup_params;
      if (use_vertex_groups && vertex_groups) {
        vgroup_params = {*vertex_groups, mesh.deform_verts(), invert_vertex_groups};
      }

      std::optional<Span<float>> vert_influence;
      Array<float> vert_influence_data;
      if (mask_field) {
        bke::MeshFieldContext field_context(mesh, AttrDomain::Point);
        FieldEvaluator evaluator(field_context, mesh.verts_num);
        vert_influence_data.reinitialize(mesh.verts_num);
        evaluator.add_with_destination(mask_field, vert_influence_data.as_mutable_span());
        evaluator.evaluate();
        vert_influence = vert_influence_data;
      }

      SpanAttributeWriter<float4x4> deform_matrix_writer;
      std::optional<Array<float3x3>> deform_matrix_buffer;
      if (deform_matrix_id) {
        deform_matrix_writer = mesh.attributes_for_write().lookup_or_add_for_write_span<float4x4>(
            *deform_matrix_id, AttrDomain::Point);
        // TODO armature deform should support outputting float4x4 directly instead of only
        // float3x3. Then this extra buffer and conversion wouldn't be necessary.
        deform_matrix_buffer = Array<float3x3>(mesh.verts_num, float3x3::identity());
      }

      BKE_armature_deform_coords(*armature_object,
                                 target_to_world,
                                 use_envelope,
                                 use_quaternion,
                                 vgroup_params,
                                 vert_influence,
                                 mesh.vert_positions_for_write(),
                                 deform_matrix_writer ? deform_matrix_buffer : std::nullopt);

      mesh.tag_positions_changed();

      if (deform_matrix_writer) {
        const Span<float3x3> deform_matrix3 = *deform_matrix_buffer;
        threading::parallel_for(IndexRange(mesh.verts_num), 1024, [&](const IndexRange range) {
          for (const int i : range) {
            deform_matrix_writer.span[i] = float4x4(deform_matrix3[i]);
          }
        });
        deform_matrix_writer.finish();
      }
    }
  });

  params.set_output("Geometry", std::move(geometry_set));
}

static void node_register()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeArmatureDeform");
  ntype.ui_name = "Armature Deform";
  ntype.ui_description = "Deform geometry using armature bones";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.declare = node_declare;
  blender::bke::node_register_type(ntype);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_armature_deform_cc
