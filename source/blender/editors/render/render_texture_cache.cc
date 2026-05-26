/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edrend
 */

#ifdef WITH_CYCLES

#  include "render_intern.hh" /* own include */

#  include <fmt/format.h>
#  include <mutex>
#  include <string>
#  include <utility>

#  include "MEM_guardedalloc.h"

#  include "DNA_image_types.h"
#  include "DNA_node_types.h"
#  include "DNA_userdef_types.h"
#  include "DNA_windowmanager_enums.h"

#  include "BLI_fileops.hh"
#  include "BLI_hash.hh"
#  include "BLI_listbase.hh"
#  include "BLI_path_utils.hh"
#  include "BLI_set.hh"
#  include "BLI_string.hh"
#  include "BLI_task.hh"
#  include "BLI_utildefines.hh"
#  include "BLI_vector.hh"

#  include "BKE_bpath.hh"
#  include "BKE_context.hh"
#  include "BKE_global.hh"
#  include "BKE_image.hh"
#  include "BKE_lib_query.hh"
#  include "BKE_main.hh"
#  include "BKE_node.hh"
#  include "BKE_node_runtime.hh"
#  include "BKE_report.hh"

#  include "BLT_translation.hh"

#  include "RNA_access.hh"
#  include "RNA_define.hh"

#  include "UI_interface_icons.hh"
#  include "UI_interface_layout.hh"
#  include "UI_resources.hh"

#  include "WM_api.hh"

#  include "CCL_api.h"

namespace blender {

/* -------------------------------------------------------------------- */
/** \name Texture Cache Image Gathering
 * \{ */

/* An image plus the multi-layer EXR subimage (layer/pass) selected by a shader
 * node. For non-multilayer images the subimage name stays empty and the legacy
 * single-buffer .tx is used. For multilayer EXRs we only cache the (layer, pass)
 * combinations that some shader node selected, not the cartesian product of the
 * whole catalog. The OpenImageIO subimage name is composed as "<layer>.<pass>"
 * to match Blender's multi-part EXR writer. */
struct ImageCacheRef {
  const Image *image;
  std::string subimage_name;

  bool operator==(const ImageCacheRef &other) const
  {
    return image == other.image && subimage_name == other.subimage_name;
  }
  uint64_t hash() const
  {
    return get_default_hash(image, subimage_name);
  }
};

/* A source image file to cache, tagged with its selected subimage. */
struct ImageCacheFile {
  const Image *image;
  std::string filepath;
  std::string subimage_name;
};

/* Compose the OpenImageIO subimage name "<layer>.<pass>" from a resolved
 * ImageUser, or just one half when the other is empty, or empty for a
 * non-multilayer image. */
static std::string compose_subimage_name(const ImageUser &iuser)
{
  const bool has_layer = iuser.layer_name[0] != '\0';
  const bool has_pass = iuser.pass_name[0] != '\0';
  if (!has_layer && !has_pass) {
    return std::string();
  }
  if (has_layer && has_pass) {
    return std::string(iuser.layer_name) + "." + iuser.pass_name;
  }
  return std::string(has_layer ? iuser.layer_name : iuser.pass_name);
}

/* Gather (image, subimage) pairs actually referenced by shader node trees. */
static Set<ImageCacheRef> gather_cache_images(Main *bmain)
{
  Set<ImageCacheRef> image_refs;

  FOREACH_NODETREE_BEGIN (bmain, ntree, id) {
    if (ntree->type != NTREE_SHADER) {
      continue;
    }
    for (bNode *node : ntree->all_nodes()) {
      Image *image = nullptr;
      ImageUser *iuser = nullptr;
      if (node->is_type("ShaderNodeTexImage"_ustr)) {
        image = blender::id_cast<Image *>(node->id);
        iuser = &static_cast<NodeTexImage *>(node->storage)->iuser;
      }
      else if (node->is_type("ShaderNodeTexEnvironment"_ustr)) {
        image = blender::id_cast<Image *>(node->id);
        iuser = &static_cast<NodeTexEnvironment *>(node->storage)->iuser;
      }
      if (!image || !iuser) {
        continue;
      }

      ImageUser local_iuser = *iuser;
      BKE_image_user_resolve_from_names(image, &local_iuser);
      image_refs.add({image, compose_subimage_name(local_iuser)});
    }
  }
  FOREACH_NODETREE_END;

  return image_refs;
}

/* Gather all source image files to the texture cache, expanding UDIM tiles and
 * optionally image sequences. */
static Vector<ImageCacheFile> gather_cache_filepaths(Main *bmain, const bool include_sequences)
{
  Vector<ImageCacheFile> filepaths;

  for (const ImageCacheRef &ref : gather_cache_images(bmain)) {
    const Image *image = ref.image;
    const std::string &subimage_name = ref.subimage_name;
    if (ELEM(image->source, IMA_SRC_MOVIE, IMA_SRC_GENERATED, IMA_SRC_VIEWER)) {
      continue;
    }
    if (BKE_image_has_packedfile(image)) {
      continue;
    }
    if (image->filepath[0] == '\0') {
      continue;
    }

    /* Get regular absolute path. */
    char filepath[FILE_MAX];
    STRNCPY(filepath, image->filepath);
    BLI_path_abs(filepath, ID_BLEND_PATH_FROM_GLOBAL(&image->id));
    BLI_path_normalize(filepath);

    /* Handle each UDIM tile. */
    if (image->source == IMA_SRC_TILED) {
      eUDIM_TILE_FORMAT tile_format = UDIM_TILE_FORMAT_NONE;
      char *udim_pattern = BKE_image_get_tile_strformat(filepath, &tile_format);

      if (tile_format != UDIM_TILE_FORMAT_NONE) {
        for (const ImageTile &tile : image->tiles) {
          char tile_filepath[FILE_MAX];
          BKE_image_set_filepath_from_tile_number(
              tile_filepath, udim_pattern, tile_format, tile.tile_number);
          if (BLI_is_file(tile_filepath)) {
            filepaths.append({image, tile_filepath, subimage_name});
          }
        }
        MEM_delete(udim_pattern);
        continue;
      }
    }

    /* Handle each image sequence frame. */
    if (image->source == IMA_SRC_SEQUENCE) {
      if (include_sequences) {
        BKE_bpath_sequence_filepaths_foreach(filepath, [&](StringRef frame_filepath) {
          filepaths.append({image, std::string(frame_filepath), subimage_name});
        });
      }
      continue;
    }

    /* Handle regular image. */
    if (BLI_is_file(filepath)) {
      filepaths.append({image, filepath, subimage_name});
    }
  }

  return filepaths;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Generate Texture Cache Operator
 * \{ */

struct GenerateTextureCacheJob {
  Main *bmain;
  bool include_sequences;
};

struct GenerateTextureCacheUI {
  int sequence_generate_num = 0;
};

static void generate_texture_cache(Main *bmain,
                                   ReportList *reports,
                                   const bool include_sequences,
                                   wmJobWorkerStatus *worker_status = nullptr)
{
  const Vector<ImageCacheFile> filepaths = gather_cache_filepaths(bmain, include_sequences);

  /* Determine which texture cache files need to be generated. */
  Set<std::string> found_tx;
  int uptodate_num = 0;
  Vector<ImageCacheFile> outdated;
  for (const ImageCacheFile &item : filepaths) {
    std::string tx_filepath;
    const bool up_to_date = CCL_resolve_texture_cache(item.image,
                                                      item.filepath.c_str(),
                                                      U.texture_cachedir,
                                                      item.subimage_name.c_str(),
                                                      tx_filepath);
    if (tx_filepath.empty() || !found_tx.add(tx_filepath)) {
      continue;
    }
    if (up_to_date) {
      uptodate_num++;
    }
    else {
      outdated.append(item);
    }
  }

  /* Generate texture cache. */
  std::atomic<int> completed_num = 0;
  std::atomic<int> failed_num = 0;
  std::mutex reports_mutex;

  blender::threading::parallel_for_each(outdated, [&](const ImageCacheFile &item) {
    if (worker_status) {
      if (G.is_break || worker_status->stop) {
        return;
      }
      worker_status->progress = (completed_num + failed_num) / float(outdated.size());
      worker_status->do_update = true;
    }

    if (CCL_generate_texture_cache(
            item.image, item.filepath.c_str(), U.texture_cachedir, item.subimage_name.c_str()))
    {
      completed_num++;
    }
    else {
      std::scoped_lock lock(reports_mutex);
      BKE_reportf(
          reports, RPT_ERROR, "Failed to generate texture cache for: %s", item.filepath.c_str());
      failed_num++;
    }
  });

  /* Report stats. */
  if (found_tx.is_empty()) {
    BKE_report(reports, RPT_INFO, "No image files found to generate tx files");
  }
  else {
    BKE_reportf(reports,
                (failed_num) ? RPT_WARNING : RPT_INFO,
                "Generated %d tx files, %d failed, %d up to date",
                completed_num.load(),
                failed_num.load(),
                uptodate_num);
  }
}

static void generate_texture_cache_startjob(void *customdata, wmJobWorkerStatus *worker_status)
{
  GenerateTextureCacheJob *job = static_cast<GenerateTextureCacheJob *>(customdata);
  generate_texture_cache(
      job->bmain, worker_status->reports, job->include_sequences, worker_status);
}

static wmOperatorStatus generate_texture_cache_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);
  const bool include_sequences = RNA_boolean_get(op->ptr, "generate_sequences");

  /* No dialog, execute immediately without job. */
  if (op->customdata == nullptr) {
    generate_texture_cache(bmain, op->reports, include_sequences);
    return OPERATOR_FINISHED;
  }

  MEM_delete(static_cast<GenerateTextureCacheUI *>(op->customdata));
  op->customdata = nullptr;

  wmWindowManager *wm = CTX_wm_manager(C);
  wmJob *wm_job = WM_jobs_get(wm,
                              CTX_wm_window(C),
                              bmain,
                              RPT_("Generating texture cache..."),
                              WM_JOB_PRIORITY | WM_JOB_PROGRESS,
                              WM_JOB_TYPE_GENERATE_TEXTURE_CACHE);

  GenerateTextureCacheJob *job = MEM_new<GenerateTextureCacheJob>(__func__);
  job->bmain = bmain;
  job->include_sequences = include_sequences;
  WM_jobs_customdata_set(wm_job, job, [](void *customdata) {
    MEM_delete(static_cast<GenerateTextureCacheJob *>(customdata));
  });

  WM_jobs_timer(wm_job, 0.2, NC_WM | ND_JOB, 0);
  WM_jobs_callbacks(wm_job, generate_texture_cache_startjob, nullptr, nullptr, nullptr);

  G.is_break = false;
  WM_jobs_start(wm, wm_job);

  return OPERATOR_FINISHED;
}

static wmOperatorStatus generate_texture_cache_invoke(bContext *C,
                                                      wmOperator *op,
                                                      const wmEvent * /*event*/)
{
  Main *bmain = CTX_data_main(C);

  /* Count files. */
  Set<std::string> found_tx;
  int regular_generate_num = 0;
  int sequence_generate_num = 0;
  int uptodate_num = 0;

  const bool include_sequences = true;
  const Vector<ImageCacheFile> filepaths = gather_cache_filepaths(bmain, include_sequences);
  for (const ImageCacheFile &item : filepaths) {
    const bool is_sequence = item.image->source == IMA_SRC_SEQUENCE;
    std::string tx_filepath;
    const bool up_to_date = CCL_resolve_texture_cache(item.image,
                                                      item.filepath.c_str(),
                                                      U.texture_cachedir,
                                                      item.subimage_name.c_str(),
                                                      tx_filepath);
    if (tx_filepath.empty() || !found_tx.add(tx_filepath)) {
      continue;
    }
    if (up_to_date) {
      uptodate_num++;
    }
    else if (is_sequence) {
      sequence_generate_num++;
    }
    else {
      regular_generate_num++;
    }
  }

  /* Nothing to generate, just execute to show the message. */
  if (regular_generate_num == 0 && sequence_generate_num == 0) {
    return generate_texture_cache_exec(C, op);
  }

  /* Show confirmation dialog. */
  std::string message;
  if (uptodate_num > 0) {
    message += fmt::format(fmt::runtime(IFACE_("{} tx files up to date")), uptodate_num);
  }
  if (regular_generate_num > 0) {
    if (!message.empty()) {
      message += '\n';
    }
    message += fmt::format(fmt::runtime(IFACE_("{} tx files to be generated")),
                           regular_generate_num);
  }

  GenerateTextureCacheUI *data = MEM_new<GenerateTextureCacheUI>(__func__);
  data->sequence_generate_num = sequence_generate_num;
  op->customdata = data;

  return WM_operator_props_dialog_popup(C,
                                        op,
                                        300,
                                        IFACE_("Generate Texture Cache?"),
                                        IFACE_("Generate"),
                                        false,
                                        std::move(message),
                                        true);
}

static void generate_texture_cache_ui(bContext * /*C*/, wmOperator *op)
{
  ui::Layout &layout = *op->layout;
  const GenerateTextureCacheUI *data = static_cast<GenerateTextureCacheUI *>(op->customdata);

  if (data) {
    if (data->sequence_generate_num > 0) {
      const std::string label = fmt::format(
          fmt::runtime(IFACE_("Include {} image sequence tx files")), data->sequence_generate_num);
      layout.prop(op->ptr, "generate_sequences", UI_ITEM_NONE, label, ICON_NONE);
    }
  }
  else {
    /* Show unconditionally outside the confirmation dialog, like a redo panel. */
    layout.prop(op->ptr, "generate_sequences", UI_ITEM_NONE, nullptr, ICON_NONE);
  }
}

static void generate_texture_cache_cancel(bContext * /*C*/, wmOperator *op)
{
  if (op->customdata) {
    MEM_delete(static_cast<GenerateTextureCacheUI *>(op->customdata));
    op->customdata = nullptr;
  }
}

void RENDER_OT_generate_texture_cache(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Generate Texture Cache";
  ot->idname = "RENDER_OT_generate_texture_cache";
  ot->description = "Generate Cycles texture cache files for all images used in shader nodes";

  /* API callbacks. */
  ot->exec = generate_texture_cache_exec;
  ot->invoke = generate_texture_cache_invoke;
  ot->ui = generate_texture_cache_ui;
  ot->cancel = generate_texture_cache_cancel;

  /* flags */
  ot->flag = OPTYPE_REGISTER;

  /* properties */
  RNA_def_boolean(ot->srna,
                  "generate_sequences",
                  true,
                  "Image Sequences",
                  "Generate texture cache files for all frames of image sequences");
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Clear Texture Cache Operator
 * \{ */

static Set<std::string> gather_clear_tx_files(Main *bmain)
{
  Set<std::string> tx_files;

  const bool include_sequences = true;
  const Vector<ImageCacheFile> filepaths = gather_cache_filepaths(bmain, include_sequences);

  Set<std::string> source_filepaths;
  for (const ImageCacheFile &item : filepaths) {
    source_filepaths.add(item.filepath);
  }

  for (const ImageCacheFile &item : filepaths) {
    BKE_image_texture_cache_filepaths_foreach(
        item.filepath.c_str(), U.texture_cachedir, [&](StringRef cache_filepath) {
          /* Defensive check just in case a tx file was referenced directly even
           * though this should not be done. */
          if (source_filepaths.contains(cache_filepath)) {
            return;
          }
          tx_files.add(cache_filepath);
        });
  }

  return tx_files;
}

static wmOperatorStatus clear_texture_cache_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);

  const Set<std::string> tx_files = gather_clear_tx_files(bmain);

  int deleted_num = 0;
  for (const std::string &tx_filepath : tx_files) {
    if (BLI_delete(tx_filepath.c_str(), false, false) == 0) {
      deleted_num++;
    }
    else {
      BKE_reportf(op->reports, RPT_ERROR, "Failed to delete tx file: %s", tx_filepath.c_str());
    }
  }

  BKE_reportf(op->reports, RPT_INFO, "Deleted %d tx files", deleted_num);
  return OPERATOR_FINISHED;
}

static wmOperatorStatus clear_texture_cache_invoke(bContext *C,
                                                   wmOperator *op,
                                                   const wmEvent * /*event*/)
{
  Main *bmain = CTX_data_main(C);
  const int used_num = gather_clear_tx_files(bmain).size();

  if (used_num == 0) {
    BKE_report(op->reports, RPT_INFO, "No tx files to delete");
    return OPERATOR_CANCELLED;
  }

  const std::string message = fmt::format(fmt::runtime(IFACE_("{} tx files to be deleted")),
                                          used_num);

  return WM_operator_confirm_ex(C,
                                op,
                                IFACE_("Clear Texture Cache?"),
                                message.c_str(),
                                IFACE_("Clear Cache"),
                                ui::AlertIcon::Info,
                                false);
}

void RENDER_OT_clear_texture_cache(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Clear Texture Cache";
  ot->idname = "RENDER_OT_clear_texture_cache";
  ot->description = "Delete Cycles texture cache files from disk";

  /* API callbacks. */
  ot->exec = clear_texture_cache_exec;
  ot->invoke = clear_texture_cache_invoke;

  /* flags */
  ot->flag = OPTYPE_REGISTER;
}

/** \} */

}  // namespace blender

#endif
