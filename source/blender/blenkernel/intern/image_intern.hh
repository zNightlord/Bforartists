/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 *
 * Internal interface between image code in blenkernel.
 */

#pragma once

namespace blender {

struct Image;
struct ImBuf;
struct ImagePass;

enum class ImageUDIMTexture {
  /** Not a UDIM aggregate. */
  None = 0,
  /** All tiles packed into a single 2D array texture. */
  Atlas,
  /** 1D-array lookup texture mapping tile number to atlas layer + UV offset/scale. */
  TileMapping,
};

#define IMA_NO_VIEW 0x7FEFEFEF

struct ImageCacheKey {
  /** UDIM buffer type. */
  ImageUDIMTexture udim_type = ImageUDIMTexture::None;
  /** Image layer/pass. */
  const ImagePass *pass = nullptr;
  /* UDIM tile number. */
  int tile_number = 0;
  /** Image sequence frame. */
  int frame = 0;
  /** Multiview. */
  int view = IMA_NO_VIEW;

  friend bool operator==(const ImageCacheKey &a, const ImageCacheKey &b)
  {
    return a.udim_type == b.udim_type && a.pass == b.pass && a.frame == b.frame &&
           a.tile_number == b.tile_number && a.view == b.view;
  }
};

/** Insert image buffer into image cache. */
void imagecache_put(Image *image, ImageCacheKey key, ImBuf *ibuf);

}  // namespace blender
