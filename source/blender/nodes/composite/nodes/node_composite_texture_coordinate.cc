/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup cmpnodes
 */

#include "COM_node_operation.hh"
#include "COM_utilities.hh"

#include "node_composite_util.hh"

/* **************** TEXTURE ******************** */

namespace blender::nodes::node_composite_texture_coordinate_cc {

static void cmp_node_texture_coordinate_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Vector>("Generated");
}

using namespace blender::realtime_compositor;

class TextureCoordinateOperation : public NodeOperation {
 public:
  using NodeOperation::NodeOperation;

  void execute() override {}

  Domain compute_domain() override
  {
    return Domain(context().get_compositing_region_size());
  }
};

static NodeOperation *get_compositor_operation(Context &context, DNode node)
{
  return new TextureCoordinateOperation(context, node);
}

}  // namespace blender::nodes::node_composite_texture_coordinate_cc

void register_node_type_cmp_texture_coordinate()
{
  namespace file_ns = blender::nodes::node_composite_texture_coordinate_cc;

  static bNodeType ntype;

  cmp_node_type_base(&ntype, CMP_NODE_TEXTURE_COORDINATE, "Texture Coordinate", NODE_CLASS_INPUT);
  ntype.declare = file_ns::cmp_node_texture_coordinate_declare;
  ntype.get_compositor_operation = file_ns::get_compositor_operation;

  nodeRegisterType(&ntype);
}
