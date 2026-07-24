/* SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup DNA
 */

#pragma once

#include "DNA_ID.h"
#include "DNA_color_types.h" /* for color management */
#include "DNA_defs.h"
#include "DNA_scene_enums.h"

#include "BLI_enum_flags.hh"

namespace blender {

namespace bke {
struct ImageRuntime;
}  // namespace bke

struct ImBufCache;
struct MovieReader;
struct PackedFile;
struct RenderResult;
struct Scene;

/** #ImageUser::flag */
enum eImageUser_Flag : short {
  IMA_ANIM_ALWAYS = 1 << 0,
  IMA_SHOW_SEQUENCER_SCENE = 1 << 1,
  // IMA_UNUSED_2 = 1 << 2,
  IMA_NEED_FRAME_RECALC = 1 << 3,
  IMA_SHOW_STEREO = 1 << 4,
  // IMA_UNUSED_5 = 1 << 5,
  IMA_USER_FRAME_IN_RANGE = (1 << 10),
};
ENUM_OPERATORS(eImageUser_Flag)

/** #Image.flag */
enum eImage_Flag : int {
  IMA_HIGH_BITDEPTH = (1 << 0),
  IMA_FLAG_UNUSED_1 = (1 << 1), /* cleared */
#ifdef DNA_DEPRECATED_ALLOW
  IMA_DO_PREMUL = (1 << 2),
#endif
  IMA_FLAG_UNUSED_4 = (1 << 4), /* cleared */
  IMA_NOCOLLECT = (1 << 5),
  IMA_FLAG_UNUSED_6 = (1 << 6), /* cleared */
  IMA_OLD_PREMUL = (1 << 7),
  IMA_FLAG_UNUSED_8 = (1 << 8), /* cleared */
  IMA_USED_FOR_RENDER = (1 << 9),
  /** For image user, but these flags are mixed. */
  // IMA_USER_FRAME_IN_RANGE = (1 << 10),
  IMA_VIEW_AS_RENDER = (1 << 11),
  IMA_FLAG_UNUSED_12 = (1 << 12), /* cleared */
  IMA_DEINTERLACE = (1 << 13),
  IMA_USE_VIEWS = (1 << 14),
  IMA_FLAG_UNUSED_15 = (1 << 15), /* cleared */
  IMA_FLAG_UNUSED_16 = (1 << 16), /* cleared */
  /** Indicates that the image has autosave information */
  IMA_AUTOSAVE_TEMPPACK = (1 << 17),
};
ENUM_OPERATORS(eImage_Flag)

/* Image.source, where the image comes from */
enum eImageSource : short {
  /* IMA_SRC_CHECK = 0, */ /* UNUSED */
  IMA_SRC_FILE = 1,
  IMA_SRC_SEQUENCE = 2,
  IMA_SRC_MOVIE = 3,
  IMA_SRC_GENERATED = 4,
  IMA_SRC_VIEWER = 5,
  IMA_SRC_TILED = 6,
};

/* Image.type, how to handle or generate the image */
enum eImageType : short {
  IMA_TYPE_IMAGE = 0,
  IMA_TYPE_MULTILAYER = 1,
  /* generated */
  IMA_TYPE_UV_TEST = 2,
  /* viewers */
  IMA_TYPE_R_RESULT = 4,
  IMA_TYPE_COMPOSITE = 5,
};

/** #Image.gen_type */
enum eImageGenType : char {
  IMA_GENTYPE_BLANK = 0,
  IMA_GENTYPE_GRID = 1,
  IMA_GENTYPE_GRID_COLOR = 2,
};

/** Size of allocated string #RenderResult::text. */
#define IMA_MAX_RENDER_TEXT_SIZE 512

/** #Image.gen_flag */
enum eImage_GenFlag : char {
  IMA_GEN_FLOAT = (1 << 0),
  IMA_GEN_TILE = (1 << 1),
};
ENUM_OPERATORS(eImage_GenFlag)

/** #Image.alpha_mode */
enum eImageAlphaMode : char {
  IMA_ALPHA_STRAIGHT = 0,
  IMA_ALPHA_PREMUL = 1,
  IMA_ALPHA_CHANNEL_PACKED = 2,
  IMA_ALPHA_IGNORE = 3,
};

/**
 * ImageUser determines which specific image buffer to use for an Image,
 * selecting a specific layer, pass, view, frame and tile.
 *
 * It is saved persistently in DNA along side Image ID pointers, for the
 * image editor, image nodes, background image, etc.
 *
 * It is also used at runtime for acquiring a more specific image, for example
 * a UDIM tile or animation frame based on context.
 */
struct ImageUser {
  /********************************* Saved Data ******************************/

  /** Layer and pass selection.
   *
   * The names have priority if they are non-empty. The indices exist for backwards
   * compatibility in blend files and APIs, and temporary usage in UI menus. */
  char layer_name[/*MAX_NAME*/ 64] = "";
  char pass_name[/*MAX_NAME*/ 64] = "";
  short layer = 0;
  short pass = 0;

  /** Multi-view selection. */
  short view = 0;
  short _pad0 = 0;

  /** UDIM tile selection. */
  int tile = 0;

  /* Animation settings. */

  /** Total amount of frames of the sequence/movie to use. */
  int frames = 0;
  /** Frame offset, and start frame in global (scene) time. */
  int offset = 0, sfra = 0;
  /** Cyclic: loop the frame range. */
  char cycl = 0;
  char _pad1 = 0;

  /* General Settings. */
  eImageUser_Flag flag = {};

  /****************************** Runtime Data ********************************/

  /** Resolved frame number, computed from animation settings and scene current frame. */
  int framenr = 0;

  /** Current scene for acquiring a render result or viewer. */
  struct Scene *scene = nullptr;

  /** Multi-view eye to draw, for stereo drawing. */
  char multiview_eye = 0;
  char _pad2[7] = {};
};

struct ImageAnim {
  struct ImageAnim *next = nullptr, *prev = nullptr;
  struct MovieReader *anim = nullptr;
};

struct ImageView {
  struct ImageView *next = nullptr, *prev = nullptr;
  char name[/*MAX_NAME*/ 64] = "";
  char filepath[/*FILE_MAX*/ 1024] = "";
};

struct ImagePackedFile {
  struct ImagePackedFile *next = nullptr, *prev = nullptr;
  struct PackedFile *packedfile = nullptr;

  /* Which view and tile this ImagePackedFile represents. Normal images will use 0 and 1001
   * respectively when creating their ImagePackedFile. Must be provided for each packed image. */
  int view = 0;
  int tile_number = 0;
  char filepath[/*FILE_MAX*/ 1024] = "";
};

struct RenderSlot {
  struct RenderSlot *next = nullptr, *prev = nullptr;
  char name[/*MAX_NAME*/ 64] = "";
  struct RenderResult *render = nullptr;
};

struct ImageTile_Runtime {
  int tilearray_layer = 0;
  int _pad = {};
  int tilearray_offset[2] = {};
  int tilearray_size[2] = {};
};

struct ImageTile {
  struct ImageTile *next = nullptr, *prev = nullptr;

  struct ImageTile_Runtime runtime;

  int tile_number = 0;

  /* for generated images */
  int gen_x = 0, gen_y = 0;
  eImageGenType gen_type = IMA_GENTYPE_BLANK;
  eImage_GenFlag gen_flag = {};
  short gen_depth = 0;
  float gen_color[4] = {};

  char label[64] = "";
};

struct ImagePass {
  struct ImagePass *next = nullptr, *prev = nullptr;

  /** Pass name, e.g. "Combined", "Depth", "AO", "Roughness". */
  char name[/*MAX_NAME*/ 64] = "";
  /** Channels IDs (like RGBA or XYZ). */
  char chan_id[24] = "";
  /** Number of channels. */
  int channels_num = 0;
  char _pad[4] = {};
  /** Solid fill of a generated pass buffer (authored catalogs only), stored in
   * the image's colorspace (no sRGB conversion on fill). */
  float gen_color[4] = {};
};

struct ImageLayer {
  struct ImageLayer *next = nullptr, *prev = nullptr;

  /** Layer name. */
  char name[/*MAX_NAME*/ 64] = "";

  ListBaseT<ImagePass> passes = {nullptr, nullptr};
};

struct Image {
#ifdef __cplusplus
  /** See #ID_Type comment for why this is here. */
  static constexpr ID_Type id_type = ID_IM;
#endif

  ID id;
  struct AnimData *adt = nullptr;

  /** File path. */
  char filepath[/*FILE_MAX*/ 1024] = "";

  /* sources from: */
  ListBaseT<ImageAnim> anims = {nullptr, nullptr};

  ListBaseT<RenderSlot> renderslots = {nullptr, nullptr};
  short render_slot = 0, last_render_slot = 0;

  eImage_Flag flag = {};
  eImageSource source = {};
  eImageType type = IMA_TYPE_IMAGE;
  int lastframe = 0;

  /* Number of iterations to perform when extracting mask for uv seam fixing. */
  short seam_margin = 8;

  char _pad2[6] = {};

  /** Deprecated. */
  DNA_DEPRECATED struct PackedFile *packedfile = nullptr;
  ListBaseT<ImagePackedFile> packedfiles = {nullptr, nullptr};
  struct PreviewImage *preview = nullptr;

  ListBaseT<ImagePackedFile> autosave_packedfiles = {nullptr, nullptr};

  char _pad3[4] = {};

  /* for generated images */
  DNA_DEPRECATED int gen_x = 1024;
  DNA_DEPRECATED int gen_y = 1024;
  DNA_DEPRECATED eImageGenType gen_type = IMA_GENTYPE_GRID;
  DNA_DEPRECATED eImage_GenFlag gen_flag = {};
  DNA_DEPRECATED short gen_depth = 0;
  DNA_DEPRECATED float gen_color[4] = {};

  /* display aspect - for UV editing images resized for faster openGL display */
  float aspx = 1.0, aspy = 1.0;

  /* color management */
  ColorManagedColorspaceSettings colorspace_settings;
  eImageAlphaMode alpha_mode = IMA_ALPHA_STRAIGHT;

  char _pad = {};

  /* Multiview */
  /** For viewer node stereoscopy. */
  char eye = 0;
  eImageFormat_ViewsFormat views_format = {};

  /* ImageTile list for UDIMs. */
  int active_tile_index = 0;
  ListBaseT<ImageTile> tiles = {nullptr, nullptr};

  ListBaseT<ImageLayer> layers = {nullptr, nullptr};

  ListBaseT<ImageView> views = {nullptr, nullptr};
  struct Stereo3dFormat *stereo3d_format = nullptr;

  bke::ImageRuntime *runtime = nullptr;
};

}  // namespace blender
