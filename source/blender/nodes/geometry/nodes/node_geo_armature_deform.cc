/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_listbase_wrapper.hh"

#include "BKE_armature.hh"
#include "BKE_curves.hh"
#include "BKE_deform.hh"
#include "BKE_grease_pencil.hh"
#include "BKE_mesh.hh"

#include "NOD_geometry_nodes_closure.hh"
#include "NOD_geometry_nodes_closure_eval.hh"
#include "NOD_rna_define.hh"

#include "UI_interface_layout.hh"
#include "UI_resources.hh"

#include "RNA_access.hh"
#include "RNA_enum_types.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_armature_deform_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.use_custom_socket_order();
  b.allow_any_socket_order();

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

  const bNode *node = b.node_or_null();
  if (node == nullptr) {
    return;
  }
  const eNodeSocketDatatype data_type = eNodeSocketDatatype(node->custom1);

  b.add_default_layout();

  switch (data_type) {
    case SOCK_VECTOR:
      b.add_input<decl::Vector>("Position", "Value")
          .implicit_field_on_all(NodeDefaultInputType::NODE_DEFAULT_INPUT_POSITION_FIELD);
      b.add_output<decl::Vector>("Position", "Value").field_on_all().align_with_previous();
      break;
    case SOCK_MATRIX:
      b.add_input<decl::Matrix>("Transform", "Value")
          .field_on_all()
          .description("Local deformation gradient for each point");
      b.add_output<decl::Matrix>("Transform", "Value")
          .field_on_all()
          .description("Local deformation gradient for each point")
          .align_with_previous();
      break;
    default:
      BLI_assert_unreachable();
      break;
  }
}

static void node_layout(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  layout->use_property_split_set(true);
  layout->use_property_decorate_set(false);
  layout->prop(ptr, "data_type", UI_ITEM_NONE, "", ICON_NONE);
}

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->custom1 = SOCK_VECTOR;
}

struct CustomDeformGroupFields {
  const bPoseChannel *pose_channel;
  Field<float> weights_field;
  Field<bool> selection_field;
};

class ArmatureDeformField final : public bke::GeometryFieldInput {
 private:
  const Object &armature_object_;
  float4x4 target_to_world_;
  const GField value_field_;
  const Field<float> mask_field_;
  bool use_envelope_;
  bool use_vertex_groups_;
  bool invert_vertex_groups_;
  Array<CustomDeformGroupFields> custom_group_fields_;
  bke::ArmatureDeformSkinningMode skinning_mode_;

 public:
  ArmatureDeformField(const Object &armature_object,
                      const float4x4 &target_to_world,
                      GField value_field,
                      Field<float> mask_field,
                      const bool use_envelope,
                      const bool use_vertex_groups,
                      const bool invert_vertex_groups,
                      const Span<CustomDeformGroupFields> custom_group_fields,
                      const bke::ArmatureDeformSkinningMode skinning_mode)
      : bke::GeometryFieldInput(value_field.cpp_type(), "Armature Deform"),
        armature_object_(armature_object),
        target_to_world_(target_to_world),
        value_field_(std::move(value_field)),
        mask_field_(std::move(mask_field)),
        use_envelope_(use_envelope),
        use_vertex_groups_(use_vertex_groups),
        invert_vertex_groups_(invert_vertex_groups),
        custom_group_fields_(custom_group_fields),
        skinning_mode_(skinning_mode)
  {
  }

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const final
  {
    const int64_t domain_size = context.attributes()->domain_size(context.domain());
    if (context.domain() != AttrDomain::Point) {
      FieldEvaluator evaluator(context, domain_size);
      evaluator.add(value_field_);
      evaluator.evaluate();
      return evaluator.get_evaluated(0);
    }

    GArray<> value_buffer(*type_, domain_size);
    Array<float> mask_buffer(domain_size);
    Array<bke::PoseChannelDeformGroup> custom_groups(custom_group_fields_.size());

    FieldEvaluator evaluator(context, domain_size);
    evaluator.add_with_destination(mask_field_, mask_buffer.as_mutable_span());
    evaluator.add_with_destination(value_field_, value_buffer.as_mutable_span());
    for (const int i : custom_group_fields_.index_range()) {
      const CustomDeformGroupFields &group_fields = custom_group_fields_[i];
      custom_groups[i].pose_channel = group_fields.pose_channel;
      custom_groups[i].deform_group.weights.reinitialize(domain_size);
      evaluator.add(group_fields.selection_field);
      MutableSpan<float> weights = custom_groups[i].deform_group.weights;
      evaluator.add_with_destination(group_fields.weights_field, weights);
    }
    /* Output indices of custom group selection masks. */
    evaluator.evaluate();
    const std::optional<Span<float>> vert_influence = mask_buffer;
    for (const int i : custom_groups.index_range()) {
      const int mask_field_output = 2 + i * 2;
      custom_groups[i].deform_group.mask = evaluator.get_evaluated_as_mask(mask_field_output);
    }

    std::optional<bke::ArmatureDeformVertexGroupParams> vgroup_params;
    if (use_vertex_groups_) {
      switch (context.type()) {
        case GeometryComponent::Type::Mesh:
          if (const Mesh *mesh = context.mesh()) {
            vgroup_params = {
                mesh->vertex_group_names, mesh->deform_verts(), invert_vertex_groups_};
          }
          break;
        case GeometryComponent::Type::Curve:
          if (const bke::CurvesGeometry *curves = context.curves()) {
            vgroup_params = {
                curves->vertex_group_names, curves->deform_verts(), invert_vertex_groups_};
          }
          break;
        case GeometryComponent::Type::GreasePencil: {
          const GreasePencil *grease_pencil = context.grease_pencil();
          const bke::greasepencil::Drawing *drawing = context.grease_pencil_layer_drawing();
          if (grease_pencil && drawing) {
            vgroup_params = {grease_pencil->vertex_group_names,
                             drawing->geometry.wrap().deform_verts(),
                             invert_vertex_groups_};
          }
          break;
        }
        default:
          break;
      }
    }

    if (type_ == &CPPType::get<float3>()) {
      bke::armature_deform_positions(armature_object_,
                                     target_to_world_,
                                     mask,
                                     vert_influence,
                                     custom_groups,
                                     vgroup_params,
                                     use_envelope_,
                                     skinning_mode_,
                                     value_buffer.as_mutable_span().typed<float3>());
    }
    else if (type_ == &CPPType::get<float4x4>()) {
      bke::armature_deform_matrices(armature_object_,
                                    target_to_world_,
                                    mask,
                                    vert_influence,
                                    custom_groups,
                                    vgroup_params,
                                    use_envelope_,
                                    skinning_mode_,
                                    value_buffer.as_mutable_span().typed<float4x4>());
    }
    else {
      /* Unsupported field type for armature deformation. */
      BLI_assert_unreachable();
    }

    return GVArray::from_garray(std::move(value_buffer));
  }

  void for_each_field_input_recursive(FunctionRef<void(const FieldInput &)> fn) const override
  {
    value_field_.node().for_each_field_input_recursive(fn);
    mask_field_.node().for_each_field_input_recursive(fn);
  }

  uint64_t hash() const override
  {
    uint64_t custom_group_fields_hash = 0;
    for (const CustomDeformGroupFields &group_fields : custom_group_fields_) {
      custom_group_fields_hash = get_default_hash(custom_group_fields_hash,
                                                  get_default_hash(group_fields.pose_channel,
                                                                   group_fields.weights_field,
                                                                   group_fields.selection_field));
    }
    return get_default_hash(
        get_default_hash(&armature_object_, target_to_world_, value_field_),
        get_default_hash(mask_field_, use_envelope_, use_vertex_groups_),
        get_default_hash(invert_vertex_groups_, skinning_mode_, custom_group_fields_hash));
  }

  bool is_equal_to(const fn::FieldNode &other) const override
  {
    if (const ArmatureDeformField *other_deform = dynamic_cast<const ArmatureDeformField *>(
            &other))
    {
      if (custom_group_fields_.size() != other_deform->custom_group_fields_.size()) {
        return false;
      }
      for (const int i : custom_group_fields_.index_range()) {
        if (custom_group_fields_[i].pose_channel !=
                other_deform->custom_group_fields_[i].pose_channel ||
            custom_group_fields_[i].weights_field !=
                other_deform->custom_group_fields_[i].weights_field ||
            custom_group_fields_[i].selection_field !=
                other_deform->custom_group_fields_[i].selection_field)
        {
          return false;
        }
      }
      return &armature_object_ == &other_deform->armature_object_ &&
             target_to_world_ == other_deform->target_to_world_ &&
             value_field_ == other_deform->value_field_ &&
             mask_field_ == other_deform->mask_field_ &&
             use_envelope_ == other_deform->use_envelope_ &&
             use_vertex_groups_ == other_deform->use_vertex_groups_ &&
             invert_vertex_groups_ == other_deform->invert_vertex_groups_ &&
             skinning_mode_ == other_deform->skinning_mode_;
    }
    return false;
  }

  std::optional<AttrDomain> preferred_domain(const GeometryComponent &component) const override
  {
    const std::optional<AttrDomain> domain = bke::try_detect_field_domain(component, value_field_);
    if (domain.has_value() && *domain == AttrDomain::Corner) {
      return AttrDomain::Point;
    }
    return domain;
  }
};

static Vector<CustomDeformGroupFields> evaluate_custom_deform_groups(
    const Object &armature_object,
    GeoNodesUserData *user_data,
    const Closure &custom_groups_closure)
{
  const ConstListBaseWrapper<bPoseChannel> pose_channels(armature_object.pose->chanbase);
  const int pose_channels_num = BLI_listbase_count(&armature_object.pose->chanbase);
  const bke::bNodeSocketType *stype_string = bke::node_socket_type_find_static(SOCK_STRING);
  const bke::bNodeSocketType *stype_float = bke::node_socket_type_find_static(SOCK_FLOAT);
  const bke::bNodeSocketType *stype_bool = bke::node_socket_type_find_static(SOCK_BOOLEAN);

  Vector<CustomDeformGroupFields> custom_group_fields;
  custom_group_fields.reserve(pose_channels_num);

  for (const bPoseChannel *pose_channel : pose_channels) {
    if (!pose_channel->bone || bool(pose_channel->bone->flag & BONE_NO_DEFORM)) {
      continue;
    }

    SocketValueVariant bone_variant(std::string(pose_channel->name));
    SocketValueVariant weights_variant;
    SocketValueVariant selection_variant;
    ClosureEagerEvalParams eval_params = {
        {{SocketInterfaceKey("Bone"), stype_string, &bone_variant}},
        {{SocketInterfaceKey("Weight"), stype_float, &weights_variant},
         {SocketInterfaceKey("Selection"), stype_bool, &selection_variant}},
        user_data};
    evaluate_closure_eagerly(custom_groups_closure, eval_params);

    Field<float> weights_field = weights_variant.extract<Field<float>>();
    Field<bool> selection_field = selection_variant.extract<Field<bool>>();
    custom_group_fields.append_unchecked(
        {pose_channel, std::move(weights_field), std::move(selection_field)});
  }
  return custom_group_fields;
}

static void node_geo_exec(GeoNodeExecParams params)
{
  GField value_field = params.extract_input<GField>("Value");

  const Object *armature_object = params.extract_input<Object *>("Armature Object");
  if (!armature_object || armature_object->type != OB_ARMATURE) {
    params.set_output("Value", value_field);
    return;
  }
  const Object *self_object = params.self_object();
  BLI_assert(self_object != nullptr);
  const float4x4 &target_to_world = self_object->object_to_world();

  const bke::ArmatureDeformSkinningMode skinning_mode =
      params.extract_input<bool>("Preserve Volume") ?
          bke::ArmatureDeformSkinningMode::DualQuatenrion :
          bke::ArmatureDeformSkinningMode::Linear;
  const bool use_envelope = params.extract_input<bool>("Use Envelope");
  const bool use_vertex_groups = params.extract_input<bool>("Use Vertex Groups");
  const bool invert_vertex_groups = params.extract_input<bool>("Invert Vertex Groups");
  const ClosurePtr custom_weights = params.extract_input<ClosurePtr>("Custom Weights");
  Field<float> mask_field = params.extract_input<Field<float>>("Mask");

  Vector<CustomDeformGroupFields> custom_group_fields;
  if (custom_weights) {
    custom_group_fields = evaluate_custom_deform_groups(
        *armature_object, params.user_data(), *custom_weights);
  }

  GField output_field{std::make_shared<ArmatureDeformField>(*armature_object,
                                                            target_to_world,
                                                            std::move(value_field),
                                                            std::move(mask_field),
                                                            use_envelope,
                                                            use_vertex_groups,
                                                            invert_vertex_groups,
                                                            custom_group_fields,
                                                            skinning_mode)};
  params.set_output<GField>("Value", std::move(output_field));
}

static void node_rna(StructRNA *srna)
{
  RNA_def_node_enum(
      srna,
      "data_type",
      "Data Type",
      "Type of grid data",
      rna_enum_node_socket_data_type_items,
      NOD_inline_enum_accessors(custom1),
      SOCK_FLOAT,
      [](bContext * /*C*/, PointerRNA * /*ptr*/, PropertyRNA * /*prop*/, bool *r_free) {
        *r_free = true;
        return enum_items_filter(rna_enum_node_socket_data_type_items,
                                 [](const EnumPropertyItem &item) -> bool {
                                   return ELEM(item.value, SOCK_VECTOR, SOCK_MATRIX);
                                 });
      });
}

static void node_register()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeArmatureDeform");
  ntype.ui_name = "Armature Deform";
  ntype.ui_description = "Deformation using armature bones";
  ntype.nclass = NODE_CLASS_CONVERTER;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.declare = node_declare;
  ntype.draw_buttons = node_layout;
  ntype.initfunc = node_init;
  blender::bke::node_register_type(ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_armature_deform_cc
