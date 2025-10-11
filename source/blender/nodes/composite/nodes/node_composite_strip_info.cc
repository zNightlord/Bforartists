/* SPDX-FileCopyrightText: 2006 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup cmpnodes
 */

#include "BLI_math_vector_types.hh"

#include "DNA_node_types.h"

#include "COM_node_operation.hh"

#include "node_composite_util.hh"

/* **************** Strip Info ******************** */

namespace blender::nodes::node_composite_strip_info_cc {

static void cmp_node_strip_info_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Int>("Start").default_value(0);
  b.add_output<decl::Int>("End").default_value(0);
}

using namespace blender::compositor;

class StripInfoOperation : public NodeOperation {
 public:
  using NodeOperation::NodeOperation;

  void execute() override
  {
    Result &result_start = get_result("Start");
    Result &result_end = get_result("End");

    result_start.allocate_single_value();
    result_end.allocate_single_value();

    Strip strip = context().get_strip();
    result_start.set_single_value(static_cast<int>(strip.start + strip.startofs));
    result_end.set_single_value(strip.enddisp);
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

  cmp_node_type_base(&ntype, "CompositorNodeStripInfo", CMP_NODE_STRIP_INFO);
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
