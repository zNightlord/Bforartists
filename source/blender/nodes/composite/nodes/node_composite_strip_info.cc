/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup cmpnodes
 */

#include "BLI_math_vector_types.hh"

#include "DNA_node_types.h"

#include "COM_node_operation.hh"

#include "SEQ_time.hh"

#include "node_composite_util.hh"

namespace blender::nodes::node_composite_strip_info_cc {

static void cmp_node_strip_info_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Int>("Start Frame").default_value(0);
  b.add_output<decl::Int>("End Frame").default_value(0);
  b.add_output<decl::Vector>("Location").default_value({0, 0});
  b.add_output<decl::Float>("Rotation").default_value(0.0);
  b.add_output<decl::Vector>("Scale").default_value({0, 0});
}

using namespace blender::compositor;

class StripInfoOperation : public NodeOperation {
 public:
  using NodeOperation::NodeOperation;

  void execute() override
  {
    Strip strip = context().get_strip();

    Result &start_frame_result = this->get_result("Start Frame");
    if (start_frame_result.should_compute()) {
      start_frame_result.allocate_single_value();
      start_frame_result.set_single_value(
          seq::time_left_handle_frame_get(&context().get_scene(), &strip));
    }

    Result &end_frame_result = this->get_result("End Frame");
    if (end_frame_result.should_compute()) {
      end_frame_result.allocate_single_value();
      end_frame_result.set_single_value(
          seq::time_right_handle_frame_get(&context().get_scene(), &strip));
    }

    Result &location_result = this->get_result("Location");
    if (location_result.should_compute()) {
      location_result.allocate_single_value();
      location_result.set_single_value(
          float2(strip.data->transform->xofs, strip.data->transform->yofs));
    }

    Result &rotation_result = this->get_result("Rotation");
    if (rotation_result.should_compute()) {
      rotation_result.allocate_single_value();
      rotation_result.set_single_value(strip.data->transform->rotation);
    }

    Result &scale_result = this->get_result("Scale");
    if (scale_result.should_compute()) {
      scale_result.allocate_single_value();
      scale_result.set_single_value(
          float2(strip.data->transform->scale_x, strip.data->transform->scale_y));
    }
  };
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
  ntype.nclass = NODE_CLASS_INPUT;
  ntype.declare = file_ns::cmp_node_strip_info_declare;
  ntype.get_compositor_operation = file_ns::get_compositor_operation;

  blender::bke::node_register_type(ntype);
}
NOD_REGISTER_NODE(register_node_type_cmp_strip_info)
