/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "COM_TextureCoordinateOperation.h"

#include "BLI_math_base.hh"

#include "BKE_scene.h"

namespace blender::compositor {

TextureCoordinateOperation::TextureCoordinateOperation()
{
  this->add_output_socket(DataType::Vector);
}

void TextureCoordinateOperation::execute_pixel_sampled(float output[4],
                                                       const float x,
                                                       const float y,
                                                       const PixelSampler /*sampler*/)
{
  const float width = this->get_width();
  const float height = this->get_height();
  const float side = math::max(width, height);

  output[0] = x / side;
  output[1] = y / side;
  output[2] = 0;
}

void TextureCoordinateOperation::update_memory_buffer_partial(MemoryBuffer *output,
                                                              const rcti &area,
                                                              Span<MemoryBuffer *> inputs)
{
  // XXX
}

void TextureCoordinateOperation::determine_canvas(const rcti &preferred_area, rcti &r_area)
{
  r_area = preferred_area;
  if (BLI_rcti_is_empty(&preferred_area)) {
    int width, height;
    BKE_render_resolution(rd_, false, &width, &height);
    r_area.xmax = preferred_area.xmin + width;
    r_area.ymax = preferred_area.ymin + height;
  }

  if (execution_model_ == eExecutionModel::FullFrame) {
    /* Determine inputs. */
    rcti temp = COM_AREA_NONE;
    NodeOperation::determine_canvas(r_area, temp);
  }
}

}  // namespace blender::compositor
