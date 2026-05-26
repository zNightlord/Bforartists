/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "scene/image_loader.h"

#include "util/cache_limiter.h"
#include "util/image.h"
#include "util/progress.h"
#include "util/string.h"

CCL_NAMESPACE_BEGIN

class OIIOImageLoader : public ImageLoader {
 public:
  /* The optional subimage_name selects one OpenImageIO subimage of the source,
   * used to read a single layer/pass of a multi-layer EXR. Resolved to an
   * integer index in load_metadata and cached on the ImageMetaData. */
  OIIOImageLoader(const string &filepath, const string &subimage_name = "");
  ~OIIOImageLoader() override;

  bool load_metadata(ImageMetaData &metadata,
                     const ImageLoaderParams &params,
                     Progress &progress) override;

  bool load_pixels(const ImageMetaData &metadata, void *pixels) override;

  bool load_pixels_tile(const ImageMetaData &metadata,
                        int miplevel,
                        int64_t x,
                        int64_t y,
                        int64_t w,
                        int64_t h,
                        int64_t x_stride,
                        int64_t y_stride,
                        int64_t padding,
                        ExtensionType extension,
                        uint8_t *pixels) override;

  string name() const override;

  bool equals(const ImageLoader &other) const override;

 protected:
  const string &get_filepath() const;

  string original_filepath_;
  /* Name of the OIIO subimage to read. Empty for plain (single-image) files;
   * non-empty selects a layer/pass of a multi-layer EXR. The name applies only
   * to the original file; a .tx cache is always single-subimage so the
   * resolved index is 0 once #texture_cache_filepath_ is in use. */
  string subimage_name_;
  string texture_cache_filepath_;
  CacheHandle<ImageInput> texture_cache_file_handle;
  bool texture_cache_file_handle_failed = false;
};

CCL_NAMESPACE_END
