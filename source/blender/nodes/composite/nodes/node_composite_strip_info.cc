/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup cmpnodes
 */

#include "BLI_math_vector_types.hh"

#include "DNA_node_types.h"

#include "COM_node_operation.hh"

#include "node_composite_util.hh"

namespace blender::nodes::node_composite_strip_info_cc {

static void cmp_node_strip_info_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Int>("Start Frame").default_value(0);
  b.add_output<decl::Int>("End Frame").default_value(0);
  b.add_output<decl::Vector>("Position").default_value({0, 0});
  b.add_output<decl::Float>("Rotation").default_value(0.0);
  b.add_output<decl::Vector>("Scale").default_value({0, 0});
}

using namespace blender::compositor;

class StripInfoOperation : public NodeOperation {
 public:
  using NodeOperation::NodeOperation;

  void execute() override
  {
    Result &result_start = this->get_result("Start Frame");
    Result &result_end = this->get_result("End Frame");
    Result &result_position = this->get_result("Position");
    Result &result_rotation = this->get_result("Rotation");
    Result &result_scale = this->get_result("Scale");

    result_start.allocate_single_value();
    result_end.allocate_single_value();
    result_position.allocate_single_value();
    result_rotation.allocate_single_value();
    result_scale.allocate_single_value();

    Strip strip = context().get_strip();

    result_start.set_single_value(static_cast<int>(strip.start + strip.startofs));
    result_end.set_single_value(strip.enddisp);

    result_position.set_single_value(
        float2(strip.data->transform->xofs, strip.data->transform->yofs));
    result_rotation.set_single_value(strip.data->transform->rotation);
    result_scale.set_single_value(
        float2(strip.data->transform->scale_x, strip.data->transform->scale_y));
  }
};

static NodeOperation *get_compositor_operation(Context &context, DNode node)
{
  return new StripInfoOperation(context, node);
}

}  // namespace blender::nodes::node_composite_strip_info_cc

static void register_node_type_cmp_strip_info()
{
  namespace file_ns = blender::nodes::node_composite_strip_info_cc;

  static blender::bke::bNodeType ntype;

  cmp_node_type_base(&ntype, "CompositorNodeStripInfo");
  ntype.ui_name = "Strip Info";
  ntype.ui_description = "Retrieve information about the strip";
  ntype.enum_name_legacy = "Strip Info";
  ntype.nclass = NODE_CLASS_INPUT;
  ntype.declare = file_ns::cmp_node_strip_info_declare;
  blender::bke::node_type_size_preset(ntype, blender::bke::eNodeSizePreset::Default);
  ntype.get_compositor_operation = file_ns::get_compositor_operation;

  blender::bke::node_register_type(ntype);
}
NOD_REGISTER_NODE(register_node_type_cmp_strip_info)
