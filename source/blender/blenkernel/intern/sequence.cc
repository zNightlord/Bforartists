/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 * \brief Low-level operations for the 'Sequence' data-block.
 */

#include "BLI_listbase.h"

#include "BKE_collection.hh"
#include "BKE_idtype.hh"
#include "BKE_layer.hh"
#include "BKE_lib_id.hh"
#include "BKE_main.hh"
#include "BKE_sequence.hh"

#include "BLT_translation.hh"

#include "SEQ_sequencer.hh"

static void sequence_data_init(ID *id)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  /* FIXME: Only here for compatibility with some scene APIs to avoid crashes. */
  Scene *new_scene = static_cast<Scene *>(BKE_id_new_nomain(ID_SCE, nullptr));
  sequence->legacy_scene_data = *new_scene;
  MEM_freeN(new_scene);

  SEQ_editing_ensure(&sequence->legacy_scene_data);
}

static void sequence_data_free(ID *id)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  BKE_libblock_free_datablock(&sequence->legacy_scene_data.id, 0);
}

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

    /*init_data*/ sequence_data_init,
    /*copy_data*/ nullptr,
    /*free_data*/ sequence_data_free,
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

Sequence *BKE_sequence_add(Main &bmain, const char *name)
{
  Sequence *sequence = static_cast<Sequence *>(BKE_id_new(&bmain, ID_SEQ, name));
  id_us_min(&sequence->id);
  id_us_ensure_real(&sequence->id);

  return sequence;
}

bool BKE_sequence_can_be_removed(const Main &bmain, const Sequence &sequence)
{
  /* Linked sequences can always be removed. */
  if (ID_IS_LINKED(&sequence)) {
    return true;
  }
  /* Local sequences can only be removed, when there is at least one local sequences left. */
  LISTBASE_FOREACH (Sequence *, other_sequence, &bmain.sequences) {
    if (other_sequence != &sequence && !ID_IS_LINKED(other_sequence)) {
      return true;
    }
  }
  return false;
}
