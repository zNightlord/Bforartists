/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 * \brief Low-level operations for the 'Sequence' data-block.
 */

#include "BKE_idtype.hh"
#include "BKE_sequence.hh"

#include "BLT_translation.hh"

IDTypeInfo IDType_ID_SEQ = {
    /*id_code*/ ID_SEQ,
    /*id_filter*/ FILTER_ID_SEQ,
    /*dependencies_id_types*/ 0,
    /*main_listbase_index*/ INDEX_ID_SEQ,
    /*struct_size*/ sizeof(Sequence),
    /*name*/ "Sequence",
    /*name_plural*/ N_("sequences"),
    /*translation_context*/ "",
    /*flags*/ IDTYPE_FLAGS_NEVER_UNUSED,
    /*asset_type_info*/ nullptr,

    /*init_data*/ nullptr,
    /*copy_data*/ nullptr,
    /*free_data*/ nullptr,
    /*make_local*/ nullptr,
    /*foreach_id*/ nullptr,
    /*foreach_cache*/ nullptr,
    /*foreach_path*/ nullptr,
    /*owner_pointer_get*/ nullptr,

    /*blend_write*/ nullptr,
    /*blend_read_data*/ nullptr,
    /*blend_read_after_liblink*/ nullptr,

    /*blend_read_undo_preserve*/ nullptr,

    /*lib_override_apply_post*/ nullptr,
};
