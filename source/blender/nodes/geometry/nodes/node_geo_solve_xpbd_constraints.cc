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

namespace blender::nodes::xpbd_constraints {

static void append_instance_item(GeometrySet &container,
                                 GeometrySet item,
                                 const StringRef name,
                                 const float4x4 &transform = float4x4::identity())
{
  if (!container.has_instances()) {
    container.replace_instances(new bke::Instances);
  }
  bke::Instances &instances =
      *container.get_component_for_write<InstancesComponent>().get_for_write();

  item.name = name;
  const int handle = instances.add_new_reference(std::move(item));
  instances.add_instance(handle, transform);
}

DebugRecorder::DebugRecorder(const bke::GeometrySet &debug_steps) : debug_steps_(debug_steps) {}

void DebugRecorder::set_geometry(const GeometrySet &geometry_set,
                                 GeometryComponent::Type component_type)
{
  geometry_set_ = geometry_set;
  component_type_ = component_type;
}

void DebugRecorder::record_step(const StringRef label,
                                ConstraintEvalParams &eval_params,
                                const ConstraintClosure *closure,
                                const IndexMask &group_mask)
{
  GeometrySet step_geometry;

  {
    GeometrySet updated_geometry = geometry_set_;
    GeometryComponent &component = updated_geometry.get_component_for_write(component_type_);
    MutableAttributeAccessor attributes = *component.attributes_for_write();
    if (!eval_params.positions.is_empty()) {
      AttributeWriter<float3> positions_writer = attributes.lookup_or_add_for_write<float3>(
          "position", AttrDomain::Point);
      positions_writer.varray.set_all(eval_params.positions);
      positions_writer.finish();
    }
    if (!eval_params.rotations.is_empty()) {
      AttributeWriter<math::Quaternion> rotations_writer =
          attributes.lookup_or_add_for_write<math::Quaternion>("rotation", AttrDomain::Point);
      rotations_writer.varray.set_all(eval_params.rotations);
      rotations_writer.finish();
    }
    if (!eval_params.velocities.is_empty()) {
      AttributeWriter<float3> velocities_writer = attributes.lookup_or_add_for_write<float3>(
          "velocity", AttrDomain::Point);
      velocities_writer.varray.set_all(eval_params.velocities);
      velocities_writer.finish();
    }
    if (!eval_params.angular_velocities.is_empty()) {
      AttributeWriter<float3> angular_velocities_writer =
          attributes.lookup_or_add_for_write<float3>("angular_velocity", AttrDomain::Point);
      angular_velocities_writer.varray.set_all(eval_params.angular_velocities);
      angular_velocities_writer.finish();
    }

    append_instance_item(step_geometry, updated_geometry, "Geometry");
  }

  if (closure) {
    GeometrySet constraint_geometry = closure->geometry_set;
    PointCloudComponent &constraint_component =
        constraint_geometry.get_component_for_write<PointCloudComponent>();
    MutableAttributeAccessor attributes = *constraint_component.attributes_for_write();
    attributes.remove("active");
    SpanAttributeWriter<bool> group_active_writer =
        attributes.lookup_or_add_for_write_only_span<bool>("group_active", AttrDomain::Point);
    group_mask.foreach_index(GrainSize(4096),
                             [&](const int index) { group_active_writer.span[index] = true; });
    group_active_writer.finish();

    append_instance_item(step_geometry, constraint_geometry, "Constraints");
  }

  append_instance_item(debug_steps_, step_geometry, label);
}

const bke::GeometrySet &DebugRecorder::debug_steps() const
{
  return debug_steps_;
}

}  // namespace blender::nodes::xpbd_constraints

namespace blender::nodes::node_geo_solve_xpbd_constraints_cc {

using xpbd_constraints::ConstraintClosure;
using xpbd_constraints::ConstraintEvalParams;
using xpbd_constraints::ConstraintTypeInfo;

constexpr float default_fps = 1.0f / 25.0f;

enum class EvaluationTarget {
  Positions,
  Velocities,
};

enum class SolverMethod {
  GaussSeidel,
  Jacobi,
};

enum class ConstraintInit {
  /* Zero initialize lambda and compute change from current positions. */
  ZeroInit,
  /* Warm start: Use previous lambda as initial step. */
  WarmStart,
};

static void node_declare_positions(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();
  b.add_input<decl::Int>("Gauss-Seidel Steps").default_value(1).min(0);
  b.add_input<decl::Bool>("Initialize")
      .default_value(true)
      .description("Initialize lambda values and positions from zero or warm start");
  b.add_input<decl::Bool>("Warm Start")
      .default_value(false)
      .description("Use previous lambda value when initializing instead of starting from zero");

  const int geometry_in = b.add_input<decl::Geometry>("Geometry").index();
  const int geometry_out = b.add_output<decl::Geometry>("Geometry").align_with_previous().index();

  b.add_input<decl::Vector>("Position")
      .implicit_field_on(implicit_field_inputs::position, {geometry_in});
  b.add_output<decl::Vector>("Position").field_on({geometry_out}).align_with_previous();
  b.add_input<decl::Rotation>("Rotation").field_on({geometry_in}).hide_value();
  b.add_output<decl::Rotation>("Rotation").field_on({geometry_out}).align_with_previous();
  b.add_input<decl::Vector>("Old Position").field_on({geometry_in}).hide_value();
  b.add_input<decl::Rotation>("Old Rotation").field_on({geometry_in}).hide_value();

  b.add_input<decl::Float>("Position Weight")
      .default_value(1.0f)
      .field_on({geometry_in})
      .description("Influence weight of constraints for each point (inverse mass)");
  b.add_input<decl::Float>("Rotation Weight")
      .default_value(1.0f)
      .field_on({geometry_in})
      .description("Influence weight of constraints for each point (inverse moment of inertia)");

  for (const ConstraintTypeInfo &info : xpbd_constraints::get_constraint_info()) {
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
  b.add_input<decl::Geometry>("Debug Steps")
      .description("Complete constraint and geometry information for each solver iteration");
  b.add_output<decl::Geometry>("Debug Steps")
      .description("Complete constraint and geometry information for each solver iteration")
      .align_with_previous();
}

static void node_declare_velocities(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_input<decl::Float>("Delta Time").default_value(default_fps).min(0.0f).hide_value();
  b.add_input<decl::Int>("Gauss-Seidel Steps").default_value(1).min(0);
  b.add_input<decl::Bool>("Initialize")
      .default_value(false)
      .description("Initialize lambda values and positions from zero or warm start");
  b.add_input<decl::Bool>("Warm Start")
      .default_value(true)
      .description("Use previous lambda value when initializing instead of starting from zero");

  const int geometry_in = b.add_input<decl::Geometry>("Geometry").index();
  const int geometry_out = b.add_output<decl::Geometry>("Geometry").align_with_previous().index();

  b.add_input<decl::Vector>("Position")
      .implicit_field_on(implicit_field_inputs::position, {geometry_in});
  b.add_input<decl::Rotation>("Rotation").field_on({geometry_in}).hide_value();
  b.add_input<decl::Vector>("Velocity").field_on({geometry_in}).hide_value();
  b.add_output<decl::Vector>("Velocity").field_on({geometry_out}).align_with_previous();
  b.add_input<decl::Vector>("Angular Velocity").field_on({geometry_in}).hide_value();
  b.add_output<decl::Vector>("Angular Velocity").field_on({geometry_out}).align_with_previous();
  b.add_input<decl::Vector>("Original Velocity").field_on({geometry_in}).hide_value();
  b.add_input<decl::Vector>("Original Angular Velocity").field_on({geometry_in}).hide_value();

  b.add_input<decl::Float>("Position Weight")
      .default_value(1.0f)
      .field_on({geometry_in})
      .description("Influence weight of constraints for each point (inverse mass)");
  b.add_input<decl::Float>("Rotation Weight")
      .default_value(1.0f)
      .field_on({geometry_in})
      .description("Influence weight of constraints for each point (inverse moment of inertia)");

  for (const ConstraintTypeInfo &info : xpbd_constraints::get_constraint_info()) {
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
  b.add_input<decl::Geometry>("Debug Steps")
      .description("Complete constraint and geometry information for each solver iteration");
  b.add_output<decl::Geometry>("Debug Steps")
      .description("Complete constraint and geometry information for each solver iteration")
      .align_with_previous();
}

/* Evaluate a group of constraints in parallel.
 * It's important that none of the constraints write to the same variables. Solver groups must be
 * computed in such a way that each variable is only affected by one constraint in each group.
 */
static void evaluate_constraint_group(EvaluationTarget target,
                                      ConstraintEvalParams &eval_params,
                                      const ConstraintTypeInfo &constraint_info,
                                      ConstraintClosure &closure,
                                      const IndexMask &group_mask)
{
  switch (target) {
    case EvaluationTarget::Positions:
      closure.apply_to_positions(eval_params, group_mask);
      break;
    case EvaluationTarget::Velocities:
      closure.apply_to_velocities(eval_params, group_mask);
      break;
  }

  if (eval_params.debug_recorder) {
    const std::string label = fmt::format("Evaluate: {}", constraint_info.ui_name);
    eval_params.debug_recorder->record_step(label, eval_params, &closure, group_mask);
  }
}

static void do_single_constraint_passes(EvaluationTarget target,
                                        ConstraintEvalParams &eval_params,
                                        const ConstraintTypeInfo &constraint_info,
                                        ConstraintClosure &closure,
                                        const Span<IndexMask> group_masks)
{
  /* Solve in consistent order by using the sorted index set. */
  for (const IndexMask &group_mask : group_masks) {
    evaluate_constraint_group(target, eval_params, constraint_info, closure, group_mask);
  }
}

static Vector<IndexMask> build_group_masks(ConstraintClosure &closure, IndexMaskMemory &memory)
{
  const VArray<int> &solver_groups = closure.solver_groups;
  if (solver_groups.is_empty()) {
    return {};
  }

  const IndexRange constraints = solver_groups.index_range();
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

  Vector<IndexMask> group_masks;
  group_masks.resize(sorted_group_indices.size());
  for (const int i : sorted_group_indices.index_range()) {
    /* Group ID is not really relevant at this point. */
    /* const int group_id = unique_group_ids[i_group]; */

    group_masks[i] = std::move(group_index_masks[sorted_group_indices[i]]);
  }
  return group_masks;
}

/* A closure and associated solver group masks. */
struct ClosureEvalInfo {
  ConstraintClosure *closure;
  Vector<IndexMask> group_masks;
};

static void do_gauss_seidel_step(EvaluationTarget target,
                                 ConstraintEvalParams &eval_params,
                                 const Span<ClosureEvalInfo> closures)
{

  const Span<ConstraintTypeInfo> constraint_infos =
      xpbd_constraints::get_constraint_info_ordered();
  for (const int i : closures.index_range()) {
    const ClosureEvalInfo &info = closures[i];
    if (info.closure) {
      do_single_constraint_passes(
          target, eval_params, constraint_infos[i], *info.closure, info.group_masks);
    }
  }
}

static void do_jacobi_step(EvaluationTarget /*target*/,
                           ConstraintEvalParams & /*eval_params*/,
                           const Span<ClosureEvalInfo> /*closures*/)
{
  /* TODO */
}

static void do_solver_steps(const SolverMethod method,
                            const int steps,
                            EvaluationTarget target,
                            std::optional<ConstraintInit> init_mode,
                            ConstraintEvalParams &eval_params,
                            const Span<ClosureEvalInfo> closures)
{
  if (init_mode) {
    for (const ClosureEvalInfo &info : closures) {
      if (!info.closure) {
        continue;
      }
      switch (*init_mode) {
        case ConstraintInit::ZeroInit:
          info.closure->reset_lambda();
          break;
        case ConstraintInit::WarmStart:
          /* TODO */
          eval_params.error_message_add("Warm starting not yet implemented");
          info.closure->reset_lambda();
          // for (const IndexMask &group_mask : info.group_masks) {
          //   switch (target) {
          //     case EvaluationTarget::Positions:
          //       info.closure->init_positions(eval_params, group_mask);
          //       break;
          //     case EvaluationTarget::Velocities:
          //       info.closure->init_velocities(eval_params, group_mask);
          //       break;
          //   }
          // }
          break;
      }
    }
  }

  std::string label;
  switch (target) {
    case EvaluationTarget::Positions:
      label = fmt::format("Init position target, ");
      break;
    case EvaluationTarget::Velocities:
      label = fmt::format("Init velocity target, ");
      break;
  }
  switch (method) {
    case SolverMethod::GaussSeidel:
      label = fmt::format("{}, Gauss-Seidel steps", label);
      break;
    case SolverMethod::Jacobi:
      label = fmt::format("{}, Jacobi steps", label);
      break;
  }
  if (eval_params.debug_recorder) {
    eval_params.debug_recorder->record_step(label, eval_params, nullptr, {});
  }

  for ([[maybe_unused]] const int step : IndexRange(steps)) {
    switch (method) {
      case SolverMethod::GaussSeidel:
        do_gauss_seidel_step(target, eval_params, closures);
        break;
      case SolverMethod::Jacobi:
        do_jacobi_step(target, eval_params, closures);
        break;
    }
  }
}

static const std::optional<ConstraintInit> extract_init_mode(GeoNodeExecParams params)
{
  return params.extract_input<bool>("Initialize") ?
             std::make_optional(params.extract_input<bool>("Warm Start") ?
                                    ConstraintInit::WarmStart :
                                    ConstraintInit::ZeroInit) :
             std::nullopt;
}

static ConstraintEvalParams extract_eval_params(GeoNodeExecParams params)
{
  const float delta_time = std::max(params.extract_input<float>("Delta Time"), 0.0f);
  const float delta_time_squared = delta_time * delta_time;
  const float inv_delta_time = math::safe_rcp(delta_time);
  const float inv_delta_time_squared = math::safe_rcp(delta_time_squared);
  const bool debug_check = params.extract_input<bool>("Debug Checks");
  const bool use_debug_steps = params.output_is_required("Debug Steps");

  ConstraintEvalParams eval_params;
  eval_params.delta_time = delta_time;
  eval_params.delta_time_squared = delta_time_squared;
  eval_params.inv_delta_time = inv_delta_time;
  eval_params.inv_delta_time_squared = inv_delta_time_squared;
  eval_params.error_message_add = [params](const StringRef message) {
    params.error_message_add(geo_eval_log::NodeWarningType::Warning, message);
  };
  eval_params.debug_check = debug_check;
  if (use_debug_steps) {
    eval_params.debug_recorder = xpbd_constraints::DebugRecorder(
        params.extract_input<GeometrySet>("Debug Steps"));
  }

  return eval_params;
}

static void bind_constraint_closures(GeoNodeExecParams params,
                                     Vector<ClosureEvalInfo> &closures,
                                     IndexMaskMemory &memory)
{
  auto error_fn = [params](const StringRef message) {
    params.error_message_add(geo_eval_log::NodeWarningType::Warning, message);
  };

  const Span<ConstraintTypeInfo> constraint_infos =
      xpbd_constraints::get_constraint_info_ordered();
  closures.reinitialize(constraint_infos.size());
  closures.fill({});
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    GeometrySet geometry_set = params.extract_input<GeometrySet>(info.ui_name);
    if (geometry_set.has_component<PointCloudComponent>()) {
      ConstraintClosure *closure = info.get_closure(std::move(geometry_set), error_fn);
      Vector<IndexMask> group_masks = build_group_masks(*closure, memory);
      closures[i] = {closure, std::move(group_masks)};
    }
  }
}

static void node_geo_exec_positions(GeoNodeExecParams params)
{
  std::optional<ConstraintInit> init_mode = extract_init_mode(params);
  const int gauss_seidel_steps = std::max(params.extract_input<int>("Gauss-Seidel Steps"), 0);
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

  Vector<ClosureEvalInfo> constraint_closures;
  IndexMaskMemory memory;
  bind_constraint_closures(params, constraint_closures, memory);

  ConstraintEvalParams eval_params = extract_eval_params(params);

  static const Array<GeometryComponent::Type> types = {bke::GeometryComponent::Type::Mesh,
                                                       bke::GeometryComponent::Type::PointCloud,
                                                       bke::GeometryComponent::Type::Curve,
                                                       bke::GeometryComponent::Type::GreasePencil};
  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    for (const bke::GeometryComponent::Type component_type : types) {
      if (geometry_set.has(component_type)) {
        bke::GeometryComponent &component = geometry_set.get_component_for_write(component_type);
        std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
        if (!attributes) {
          continue;
        }

        if (eval_params.debug_recorder) {
          eval_params.debug_recorder->set_geometry(geometry_set, component_type);
        }

        const int num_points = attributes->domain_size(AttrDomain::Point);
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

        do_solver_steps(SolverMethod::GaussSeidel,
                        gauss_seidel_steps,
                        EvaluationTarget::Positions,
                        init_mode,
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

  for (const ClosureEvalInfo &info : constraint_closures) {
    if (info.closure) {
      info.closure->finish_attributes();
    }
  }

  params.set_output("Geometry", geometry_set);

  const Span<ConstraintTypeInfo> constraint_infos = xpbd_constraints::get_constraint_info();
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    if (constraint_closures[i].closure) {
      params.set_output(info.ui_name, constraint_closures[i].closure->geometry_set);

      delete constraint_closures[i].closure;
    }
    else {
      params.set_output(info.ui_name, GeometrySet{});
    }
  }

  if (eval_params.debug_recorder) {
    params.set_output("Debug Steps", eval_params.debug_recorder->debug_steps());
  }
}

static void node_geo_exec_velocities(GeoNodeExecParams params)
{
  std::optional<ConstraintInit> init_mode = extract_init_mode(params);
  const int gauss_seidel_steps = std::max(params.extract_input<int>("Gauss-Seidel Steps"), 0);
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

  Vector<ClosureEvalInfo> constraint_closures;
  IndexMaskMemory memory;
  bind_constraint_closures(params, constraint_closures, memory);

  ConstraintEvalParams eval_params = extract_eval_params(params);

  static const Array<GeometryComponent::Type> types = {bke::GeometryComponent::Type::Mesh,
                                                       bke::GeometryComponent::Type::PointCloud,
                                                       bke::GeometryComponent::Type::Curve,
                                                       bke::GeometryComponent::Type::GreasePencil};
  geometry_set.modify_geometry_sets([&](GeometrySet &geometry_set) {
    for (const bke::GeometryComponent::Type component_type : types) {
      if (geometry_set.has(component_type)) {
        bke::GeometryComponent &component = geometry_set.get_component_for_write(component_type);
        std::optional<bke::MutableAttributeAccessor> attributes = component.attributes_for_write();
        if (!attributes) {
          continue;
        }

        if (eval_params.debug_recorder) {
          eval_params.debug_recorder->set_geometry(geometry_set, component_type);
        }

        const int num_points = attributes->domain_size(AttrDomain::Point);
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

        do_solver_steps(SolverMethod::GaussSeidel,
                        gauss_seidel_steps,
                        EvaluationTarget::Velocities,
                        init_mode,
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

  for (const ClosureEvalInfo &info : constraint_closures) {
    if (info.closure) {
      info.closure->finish_attributes();
    }
  }

  params.set_output("Geometry", geometry_set);

  const Span<ConstraintTypeInfo> constraint_infos = xpbd_constraints::get_constraint_info();
  for (const int i : constraint_infos.index_range()) {
    const ConstraintTypeInfo &info = constraint_infos[i];
    if (constraint_closures[i].closure) {
      params.set_output(info.ui_name, constraint_closures[i].closure->geometry_set);

      delete constraint_closures[i].closure;
    }
    else {
      params.set_output(info.ui_name, GeometrySet{});
    }
  }

  if (eval_params.debug_recorder) {
    params.set_output("Debug Steps", eval_params.debug_recorder->debug_steps());
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
  node_type_size(ntype, 200, 120, 300);
  ntype.geometry_node_execute = node_geo_exec_positions;
  ntype.declare = node_declare_positions;
  blender::bke::node_register_type(ntype);
}

static void node_register_velocity_solve()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeSolveVelocityConstraints");
  ntype.ui_name = "Solve XPBD Velocity Constraints";
  ntype.ui_description =
      "Solve linear and angular velocity constraints on geometry using the XPBD framework";
  ntype.nclass = NODE_CLASS_GEOMETRY;
  node_type_size(ntype, 200, 120, 300);
  ntype.geometry_node_execute = node_geo_exec_velocities;
  ntype.declare = node_declare_velocities;
  blender::bke::node_register_type(ntype);
}

static void node_register()
{
  node_register_position_solve();
  node_register_velocity_solve();
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_solve_xpbd_constraints_cc
