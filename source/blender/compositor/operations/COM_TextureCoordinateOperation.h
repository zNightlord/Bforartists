/* SPDX-FileCopyrightText: 2032 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "COM_ConstantOperation.h"

namespace blender::compositor {

class TextureCoordinateOperation : public MultiThreadedOperation {
  const RenderData *rd_;

 public:
  TextureCoordinateOperation();

  void set_render_data(const RenderData *rd)
  {
    rd_ = rd;
  }

  void execute_pixel_sampled(float output[4], float x, float y, PixelSampler sampler) override;
  void update_memory_buffer_partial(MemoryBuffer *output,
                                    const rcti &area,
                                    Span<MemoryBuffer *> inputs) override;

  void determine_canvas(const rcti &preferred_area, rcti &r_area) override;
};

}  // namespace blender::compositor
