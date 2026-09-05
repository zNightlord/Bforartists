/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include <string>

namespace blender {

struct Image;

/* Create python module _cycles used by addon. */
void *CCL_python_module_init();

void CCL_log_init();
void CCL_implicit_sharing_init();

/* Texture cache generation. */

/* The optional subimage_name selects a layer/pass of a multi-layer EXR; the
 * matching OpenImageIO subimage is then folded into the .tx cache hash so each
 * pass gets its own cache file. An empty string (default) preserves the legacy
 * single-buffer behavior, with a bit-identical hash for non-multilayer images. */
bool CCL_resolve_texture_cache(const Image *image,
                               const char *filepath,
                               const char *texture_cache_directory,
                               const char *subimage_name,
                               std::string &r_tx_filepath);

bool CCL_generate_texture_cache(const Image *image,
                                const char *filepath,
                                const char *texture_cache_directory = "",
                                const char *subimage_name = "");

}  // namespace blender
