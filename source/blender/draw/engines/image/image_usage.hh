/* SPDX-FileCopyrightText: 2022 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup draw_engine
 */

#pragma once

#include <cstring>

#include "DNA_color_types.h"
#include "DNA_image_types.h"

namespace blender::image_engine {

/**
 * ImageUsage contains data of the image and image user to identify changes that require a rebuild
 * the texture slots.
 */
struct ImageUsage {
  /**
   * Layer/pass selection by name (the authoritative selection), and the view
   * selection inputs (selected view, stereo eye, stereo-display flag). Storing
   * the selection inputs rather than the resolved indices means any change that
   * would change the displayed buffer is detected without depending on resolved
   * runtime state.
   */
  char layer_name[/*MAX_NAME*/ 64] = "";
  char pass_name[/*MAX_NAME*/ 64] = "";
  short view = 0;
  /** Render slot of the image that is displayed. */
  short render_slot = 0;
  char multiview_eye = 0;
  bool show_stereo = false;

  /** Colorspace name of the image (#ColorManagedColorspaceSettings::name). */
  char colorspace_name[/*MAX_COLORSPACE_NAME*/ 64] = "";
  /** IMA_ALPHA_* */
  char alpha_mode = 0;
  bool last_tile_drawing = false;

  const void *last_image = nullptr;
  const void *last_scene = nullptr;

  ImageUsage() = default;
  ImageUsage(const blender::Image *image,
             const blender::ImageUser *image_user,
             bool do_tile_drawing)
  {
    if (image_user) {
      memcpy(layer_name, image_user->layer_name, sizeof(layer_name));
      memcpy(pass_name, image_user->pass_name, sizeof(pass_name));
      view = image_user->view;
      render_slot = image->render_slot;
      multiview_eye = image_user->multiview_eye;
      show_stereo = (image_user->flag & IMA_SHOW_STEREO) != 0;
    }
    memcpy(colorspace_name, image->colorspace_settings.name, sizeof(colorspace_name));
    alpha_mode = image->alpha_mode;
    last_image = static_cast<const void *>(image);
    last_scene = image_user ? static_cast<const void *>(image_user->scene) : nullptr;
    last_tile_drawing = do_tile_drawing;
  }

  bool operator==(const ImageUsage &other) const = default;
};

}  // namespace blender::image_engine
