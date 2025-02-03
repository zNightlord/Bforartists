/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_cpp_types.hh"
#include "BLI_parameter_pack_utils.hh"

#include "FN_multi_function_builder.hh"

#include "NOD_rna_define.hh"
#include "NOD_socket.hh"
#include "NOD_socket_search_link.hh"
#include "NOD_xpbd_constraints.hh"

#include "RNA_define.hh"
#include "RNA_enum_types.hh"

#include "node_function_util.hh"

#include "UI_interface.hh"
#include "UI_resources.hh"

namespace blender::fn::multi_function::build {

/** Build multi-function with 13 single-input and 5 single output parameters. */
template<typename In1,
         typename In2,
         typename In3,
         typename In4,
         typename In5,
         typename In6,
         typename In7,
         typename In8,
         typename In9,
         typename In10,
         typename In11,
         typename In12,
         typename In13,
         typename Out1,
         typename Out2,
         typename Out3,
         typename Out4,
         typename Out5,
         typename ElementFn,
         typename ExecPreset = exec_presets::Materialized>
inline auto SI13_SO5(const char *name,
                     const ElementFn element_fn,
                     const ExecPreset exec_preset = exec_presets::Materialized())
{
  constexpr auto param_tags = TypeSequence<ParamTag<ParamCategory::SingleInput, In1>,
                                           ParamTag<ParamCategory::SingleInput, In2>,
                                           ParamTag<ParamCategory::SingleInput, In3>,
                                           ParamTag<ParamCategory::SingleInput, In4>,
                                           ParamTag<ParamCategory::SingleInput, In5>,
                                           ParamTag<ParamCategory::SingleInput, In6>,
                                           ParamTag<ParamCategory::SingleInput, In7>,
                                           ParamTag<ParamCategory::SingleInput, In8>,
                                           ParamTag<ParamCategory::SingleInput, In9>,
                                           ParamTag<ParamCategory::SingleInput, In10>,
                                           ParamTag<ParamCategory::SingleInput, In11>,
                                           ParamTag<ParamCategory::SingleInput, In12>,
                                           ParamTag<ParamCategory::SingleInput, In13>,
                                           ParamTag<ParamCategory::SingleOutput, Out1>,
                                           ParamTag<ParamCategory::SingleOutput, Out2>,
                                           ParamTag<ParamCategory::SingleOutput, Out3>,
                                           ParamTag<ParamCategory::SingleOutput, Out4>,
                                           ParamTag<ParamCategory::SingleOutput, Out5>>();
  auto call_fn = detail::build_multi_function_call_from_element_fn(
      element_fn, exec_preset, param_tags);
  return detail::CustomMF(name, call_fn, param_tags);
}

}  // namespace blender::fn::multi_function::build

namespace blender::nodes::node_fn_constraint_cc {

using xpbd_constraints::ConstraintType;

static const EnumPropertyItem rna_enum_constraint_type_items[] = {
    {int(xpbd_constraints::ConstraintType::PositionGoal), "POSITION_GOAL", 0, "Position Goal", ""},
    {int(xpbd_constraints::ConstraintType::RotationGoal), "ROTATION_GOAL", 0, "Rotation Goal", ""},
    {int(xpbd_constraints::ConstraintType::StretchShear), "STRETCH_SHEAR", 0, "Stretch/Shear", ""},
    {int(xpbd_constraints::ConstraintType::BendTwist), "BEND_TWIST", 0, "Bend/Twist", ""},
    {int(xpbd_constraints::ConstraintType::Contact), "CONTACT", 0, "Contact", ""},
    {0, nullptr, 0, nullptr, nullptr},
};

static void node_declare(NodeDeclarationBuilder &b)
{
  b.is_function_node();
  b.use_custom_socket_order();
  b.allow_any_socket_order();

  b.add_default_layout();

  const bNode *node = b.node_or_null();
  if (!node) {
    return;
  }

  const ConstraintType constraint_type = ConstraintType(node->custom1);
  switch (constraint_type) {
    case ConstraintType::PositionGoal:
      break;
    case ConstraintType::RotationGoal:
      break;
    case ConstraintType::StretchShear:
      break;
    case ConstraintType::BendTwist:
      break;
    case ConstraintType::Contact:
      b.add_input<decl::Float>("Lambda");
      b.add_output<decl::Float>("Lambda").align_with_previous();
      b.add_input<decl::Vector>("Position 1").hide_value();
      b.add_output<decl::Vector>("Position 1").align_with_previous();
      b.add_input<decl::Vector>("Position 2").hide_value();
      b.add_output<decl::Vector>("Position 2").align_with_previous();
      b.add_input<decl::Rotation>("Rotation 1").hide_value();
      b.add_output<decl::Rotation>("Rotation 1").align_with_previous();
      b.add_input<decl::Rotation>("Rotation 2").hide_value();
      b.add_output<decl::Rotation>("Rotation 2").align_with_previous();
      b.add_separator();
      b.add_input<decl::Float>("Position Weight 1").default_value(1.0f);
      b.add_input<decl::Float>("Position Weight 2").default_value(1.0f);
      b.add_input<decl::Float>("Rotation Weight 1").default_value(1.0f);
      b.add_input<decl::Float>("Rotation Weight 2").default_value(1.0f);
      b.add_separator();
      b.add_input<decl::Vector>("Local Position 1");
      b.add_input<decl::Vector>("Local Position 2");
      b.add_input<decl::Vector>("Normal");
      b.add_input<decl::Float>("Alpha");
      break;
  }
}

static void node_layout(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiItemR(layout, ptr, "constraint_type", UI_ITEM_NONE, "", ICON_NONE);
}

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->custom1 = int(ConstraintType::Contact);
}

static const mf::MultiFunction *get_multi_function(const bNode &bnode)
{
  namespace mf = fn::multi_function;

  const ConstraintType constraint_type = ConstraintType(bnode.custom1);
  static auto exec_preset = mf::build::exec_presets::AllSpanOrSingle();

  static auto fn_contact = mf::build::SI13_SO5<float,
                                               float3,
                                               float3,
                                               math::Quaternion,
                                               math::Quaternion,
                                               float,
                                               float,
                                               float,
                                               float,
                                               float3,
                                               float3,
                                               float3,
                                               float,
                                               float,
                                               float3,
                                               float3,
                                               math::Quaternion,
                                               math::Quaternion>(
      "XPBD Contact",
      [](float lambda,
         float3 position1,
         float3 position2,
         math::Quaternion rotation1,
         math::Quaternion rotation2,
         const float weight_pos1,
         const float weight_pos2,
         const float weight_rot1,
         const float weight_rot2,
         const float3 &local_position1,
         const float3 &local_position2,
         const float3 &normal,
         const float alpha,
         float &lambda_out,
         float3 &position1_out,
         float3 &position2_out,
         math::Quaternion &rotation1_out,
         math::Quaternion &rotation2_out) -> void {
        xpbd_constraints::eval_position_contact(weight_pos1,
                                                weight_pos2,
                                                weight_rot1,
                                                weight_rot2,
                                                local_position1,
                                                local_position2,
                                                normal,
                                                alpha,
                                                lambda,
                                                position1,
                                                position2,
                                                rotation1,
                                                rotation2);
        lambda = lambda_out;
        position1_out = position1;
        position2_out = position2;
        rotation1_out = rotation1;
        rotation2_out = rotation2;
      },
      exec_preset);

  switch (constraint_type) {
    case ConstraintType::PositionGoal:
    case ConstraintType::RotationGoal:
    case ConstraintType::StretchShear:
    case ConstraintType::BendTwist:
      BLI_assert_unreachable();
      return nullptr;

    case ConstraintType::Contact:
      return &fn_contact;
    default:
      BLI_assert_unreachable();
      return nullptr;
  }
}

static void node_build_multi_function(NodeMultiFunctionBuilder &builder)
{
  const mf::MultiFunction *fn = get_multi_function(builder.node());
  builder.set_matching_fn(fn);
}

static void node_rna(StructRNA *srna)
{
  RNA_def_node_enum(srna,
                    "constraint_type",
                    "Constraint Type",
                    "",
                    rna_enum_constraint_type_items,
                    NOD_inline_enum_accessors(custom1),
                    int(ConstraintType::PositionGoal));
}

static void node_register()
{
  static blender::bke::bNodeType ntype;
  fn_node_type_base(&ntype, "FunctionNodeConstraint", FN_NODE_CONSTRAINT);
  ntype.ui_name = "Constraint";
  ntype.enum_name_legacy = "CONSTRAINT";
  ntype.nclass = NODE_CLASS_CONVERTER;
  ntype.declare = node_declare;
  ntype.initfunc = node_init;
  ntype.build_multi_function = node_build_multi_function;
  ntype.draw_buttons = node_layout;
  bke::node_type_size(&ntype, 200, 100, 300);
  blender::bke::node_register_type(&ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_fn_constraint_cc
