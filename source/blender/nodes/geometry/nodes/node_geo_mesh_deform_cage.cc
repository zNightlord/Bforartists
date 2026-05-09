/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup geo_nodes
 *
 * GEO_NODE_MESH_CAGE_DEFORM
 *
 * Deforms the input geometry by tracking how a "cage" mesh has moved away from
 * its rest pose and propagating those displacements to nearby points via
 * inverse-distance-weighted (IDW) interpolation.
 *
 * Registration
 * ------------
 * 1.  Add  GEO_NODE_MESH_CAGE_DEFORM  to  DNA_node_types.h  (NodeGeometryType enum).
 * 2.  Add  NODE_GEO_MESH_CAGE_DEFORM  entry in  NOD_static_types.h  DefNode table.
 * 3.  Include this file in  CMakeLists.txt  under the geometry nodes sources.
 * 4.  Add the node to the geometry node add menu / search list.
 *
 * API notes (Bforartists-specific)
 * ---------------------------------
 * - Socket names : add_input / add_output / extract_input / set_output  ->  "name"_ustr
 * - .description()                                                       ->  plain "text" (std::string)
 * - geo_node_type_base : (bNodeType*, UString name, optional<int16_t> enum_id)
 * - node_register_type : takes reference, not pointer
 */

#include "BLI_array.hh"
#include "BLI_index_range.hh"
#include "BLI_math_vector.hh"
#include "BLI_task.hh"
#include "BLI_ustring.hh"

#include "DNA_mesh_types.h"

#include "BKE_attribute.hh"
#include "BKE_curves.hh"
#include "BKE_geometry_set.hh"
#include "BKE_mesh.hh"
#include "BKE_pointcloud.hh"

#include "node_geometry_util.hh"

/* ------------------------------------------------------------------ */
/* Tuning constants                                                     */
/* ------------------------------------------------------------------ */

#define MESHDEFORM_GEO_MIN_DIST 1e-6f
#define MESHDEFORM_GEO_IDW_POWER 2.0f
#define MESHDEFORM_GEO_MAX_INFLUENCE 64

namespace blender::nodes::node_geo_mesh_deform_cage_cc {

/* ------------------------------------------------------------------ */
/* Node declaration                                                     */
/* ------------------------------------------------------------------ */

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Geometry>("Geometry"_ustr)
      .supported_type({GeometryComponent::Type::Mesh,
                       GeometryComponent::Type::PointCloud,
                       GeometryComponent::Type::Curve});

  b.add_input<decl::Geometry>("Cage"_ustr)
      .supported_type(GeometryComponent::Type::Mesh)
      .description("Deformed cage mesh for the current frame");

  b.add_input<decl::Geometry>("Rest Cage"_ustr)
      .supported_type(GeometryComponent::Type::Mesh)
      .description("Original rest-pose cage mesh that matches the bind frame");

  /* Selection: bool field, masks which points are deformed. */
  b.add_input<decl::Bool>("Selection"_ustr)
      .default_value(true)
      .field_on_all()
      .description("Which points are affected. Unselected points keep their original positions");

  b.add_input<decl::Int>("Influence Count"_ustr)
      .default_value(8)
      .min(1)
      .max(MESHDEFORM_GEO_MAX_INFLUENCE)
      .description("Number of nearest cage vertices used per point. Higher values are smoother but slower");

  b.add_input<decl::Float>("Strength"_ustr)
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .description("Global blend factor between original and deformed positions");

  b.add_input<decl::Float>("Factor"_ustr)
      .default_value(1.0f)
      .min(0.0f)
      .max(1.0f)
      .subtype(PROP_FACTOR)
      .field_on_all()
      .description("Per-point weight. Connect a Named Attribute or vertex group here");

  b.add_output<decl::Geometry>("Geometry"_ustr).propagate_all();
}

/* ------------------------------------------------------------------ */
/* K-nearest via brute-force partial sort                               */
/* ------------------------------------------------------------------ */

static int find_k_nearest(const float3 &co,
                          const Span<float3> rest_positions,
                          const int k,
                          MutableSpan<int> r_indices,
                          MutableSpan<float> r_dist2)
{
  const int cage_verts_num = int(rest_positions.size());
  const int actual_k = min_ii(k, cage_verts_num);

  Array<std::pair<float, int>> buf(cage_verts_num);
  for (int i = 0; i < cage_verts_num; i++) {
    buf[i] = {math::distance_squared(co, rest_positions[i]), i};
  }
  std::partial_sort(buf.begin(), buf.begin() + actual_k, buf.end());

  for (int i = 0; i < actual_k; i++) {
    r_dist2[i] = buf[i].first;
    r_indices[i] = buf[i].second;
  }
  return actual_k;
}

/* ------------------------------------------------------------------ */
/* Per-point IDW displacement                                            */
/* ------------------------------------------------------------------ */

static float3 compute_idw_displacement(const float3 &co,
                                       const Span<float3> rest_positions,
                                       const Span<float3> cage_positions,
                                       const int k)
{
  const int actual_k = min_ii(k, int(rest_positions.size()));
  if (actual_k == 0) {
    return float3(0.0f);
  }

  Array<int> indices(actual_k);
  Array<float> dist2(actual_k);
  const int found = find_k_nearest(co, rest_positions, actual_k, indices, dist2);

  float total_weight = 0.0f;
  float3 delta(0.0f);

  for (int qi = 0; qi < found; qi++) {
    const int vi = indices[qi];
    const float3 disp = cage_positions[vi] - rest_positions[vi];

    if (dist2[qi] < (MESHDEFORM_GEO_MIN_DIST * MESHDEFORM_GEO_MIN_DIST)) {
      return disp; /* Exactly on cage vertex — full weight. */
    }

    const float w = 1.0f / powf(dist2[qi], MESHDEFORM_GEO_IDW_POWER * 0.5f);
    delta += disp * w;
    total_weight += w;
  }

  if (total_weight > 0.0f) {
    delta /= total_weight;
  }
  return delta;
}

/* ------------------------------------------------------------------ */
/* Deform a flat position span                                          */
/* ------------------------------------------------------------------ */

static void deform_positions(MutableSpan<float3> positions,
                             const VArray<bool> &selection,
                             const VArray<float> &factors,
                             const Span<float3> rest_positions,
                             const Span<float3> cage_positions,
                             const int influence_k,
                             const float strength)
{
  threading::parallel_for(positions.index_range(), 512, [&](const IndexRange range) {
    for (const int i : range) {
      if (!selection[i]) {
        continue; /* Masked out — keep original position. */
      }
      const float per_point_fac = factors[i] * strength;
      if (per_point_fac <= 0.0f) {
        continue;
      }
      const float3 disp = compute_idw_displacement(
          positions[i], rest_positions, cage_positions, influence_k);
      positions[i] += disp * per_point_fac;
    }
  });
}

/* ------------------------------------------------------------------ */
/* Per-component dispatch                                               */
/* ------------------------------------------------------------------ */

static void deform_component(GeometryComponent &component,
                             const Field<bool> &selection_field,
                             const Field<float> &factor_field,
                             const Span<float3> rest_cage_positions,
                             const Span<float3> cage_positions,
                             const int influence_k,
                             const float strength)
{
  const int domain_size = component.attribute_domain_size(AttrDomain::Point);
  if (domain_size == 0) {
    return;
  }

  const bke::GeometryFieldContext field_ctx(component, AttrDomain::Point);
  fn::FieldEvaluator evaluator(field_ctx, domain_size);
  evaluator.add(selection_field);
  evaluator.add(factor_field);
  evaluator.evaluate();
  const VArray<bool> selection = evaluator.get_evaluated<bool>(0);
  const VArray<float> factors = evaluator.get_evaluated<float>(1);

  switch (component.type()) {
    case GeometryComponent::Type::Mesh: {
      MeshComponent &mesh_comp = static_cast<MeshComponent &>(component);
      Mesh *mesh = mesh_comp.get_for_write();
      if (!mesh) {
        return;
      }
      MutableSpan<float3> positions = mesh->vert_positions_for_write();
      deform_positions(
          positions, selection, factors, rest_cage_positions, cage_positions, influence_k, strength);
      mesh->tag_positions_changed();
      break;
    }
    case GeometryComponent::Type::PointCloud: {
      PointCloudComponent &pc_comp = static_cast<PointCloudComponent &>(component);
      PointCloud *pointcloud = pc_comp.get_for_write();
      if (!pointcloud) {
        return;
      }
      bke::MutableAttributeAccessor attrs = pointcloud->attributes_for_write();
      bke::SpanAttributeWriter<float3> positions =
          attrs.lookup_or_add_for_write_span<float3>("position", AttrDomain::Point);
      deform_positions(positions.span,
                       selection,
                       factors,
                       rest_cage_positions,
                       cage_positions,
                       influence_k,
                       strength);
      positions.finish();
      break;
    }
    case GeometryComponent::Type::Curve: {
      CurveComponent &curve_comp = static_cast<CurveComponent &>(component);
      Curves *curves_id = curve_comp.get_for_write();
      if (!curves_id) {
        return;
      }
      bke::CurvesGeometry &curves = curves_id->geometry.wrap();
      MutableSpan<float3> positions = curves.positions_for_write();
      deform_positions(positions,
                       selection,
                       factors,
                       rest_cage_positions,
                       cage_positions,
                       influence_k,
                       strength);
      curves.tag_positions_changed();
      break;
    }
    default:
      break;
  }
}

/* ------------------------------------------------------------------ */
/* Node execution                                                       */
/* ------------------------------------------------------------------ */

static void node_geo_exec(GeoNodeExecParams params)
{
  GeometrySet geometry = params.extract_input<GeometrySet>("Geometry"_ustr);
  const GeometrySet cage_set = params.extract_input<GeometrySet>("Cage"_ustr);
  const GeometrySet rest_cage_set = params.extract_input<GeometrySet>("Rest Cage"_ustr);

  Field<bool> selection_field = params.extract_input<Field<bool>>("Selection"_ustr);
  const int influence_k = clamp_i(
      params.extract_input<int>("Influence Count"_ustr), 1, MESHDEFORM_GEO_MAX_INFLUENCE);
  const float strength = params.extract_input<float>("Strength"_ustr);
  Field<float> factor_field = params.extract_input<Field<float>>("Factor"_ustr);

  const Mesh *cage_mesh = cage_set.get_mesh();
  const Mesh *rest_cage_mesh = rest_cage_set.get_mesh();

  if (cage_mesh == nullptr) {
    params.error_message_add(NodeWarningType::Error, TIP_("Cage input is not a mesh"));
    params.set_output("Geometry"_ustr, std::move(geometry));
    return;
  }
  if (rest_cage_mesh == nullptr) {
    params.error_message_add(NodeWarningType::Error, TIP_("Rest Cage input is not a mesh"));
    params.set_output("Geometry"_ustr, std::move(geometry));
    return;
  }

  const Span<float3> cage_positions = cage_mesh->vert_positions();
  const Span<float3> rest_positions = rest_cage_mesh->vert_positions();

  if (cage_positions.size() != rest_positions.size()) {
    params.error_message_add(
        NodeWarningType::Error,
        TIP_("Cage and Rest Cage must have the same number of vertices"));
    params.set_output("Geometry"_ustr, std::move(geometry));
    return;
  }

  if (cage_positions.is_empty()) {
    params.set_output("Geometry"_ustr, std::move(geometry));
    return;
  }

  for (const GeometryComponent::Type type : {GeometryComponent::Type::Mesh,
                                             GeometryComponent::Type::PointCloud,
                                             GeometryComponent::Type::Curve})
  {
    if (!geometry.has(type)) {
      continue;
    }
    GeometryComponent &component = geometry.get_component_for_write(type);
    deform_component(component,
                     selection_field,
                     factor_field,
                     rest_positions,
                     cage_positions,
                     influence_k,
                     strength);
  }

  params.set_output("Geometry"_ustr, std::move(geometry));
}

/* ------------------------------------------------------------------ */
/* Node type registration                                               */
/* ------------------------------------------------------------------ */

static void node_register()
{
  static blender::bke::bNodeType ntype;
  /* BFA: geo_node_type_base(bNodeType*, UString name, optional<int16_t> enum_id)
   * Name is arg 2 as UString; enum is optional arg 3.               */
  geo_node_type_base(&ntype, "GeometryNodeMeshDeformCage"_ustr, GEO_NODE_MESH_DEFORM_CAGE);
  ntype.declare = node_declare;
  ntype.geometry_node_execute = node_geo_exec;
  blender::bke::node_register_type(ntype); /* reference, not pointer */
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_mesh_deform_cage_cc
