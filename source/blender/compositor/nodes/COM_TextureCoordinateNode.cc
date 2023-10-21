/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "COM_TextureCoordinateNode.h"

#include "DNA_node_types.h"
#include "DNA_scene_types.h"

#include "COM_TextureCoordinateOperation.h"

namespace blender::compositor {

void TextureCoordinateNode::convert_to_operations(NodeConverter &converter,
                                                  const CompositorContext &context) const
{
  TextureCoordinateOperation *operation = new TextureCoordinateOperation();
  NodeOutput *output = this->get_output_socket(0);

  operation->set_render_data(context.get_render_data());
  converter.add_operation(operation);

  converter.map_output_socket(output, operation->get_output_socket());
}

}  // namespace blender::composito
