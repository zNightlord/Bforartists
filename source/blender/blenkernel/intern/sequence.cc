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
#include "BKE_scene.hh"
#include "BKE_sequence.hh"

#include "BLT_translation.hh"

#include "BLO_read_write.hh"
#include "BLO_readfile.hh"

#include "SEQ_sequencer.hh"

static void sequence_data_init(ID *id)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  /* FIXME: Only here for compatibility with some scene APIs to avoid crashes. */
  Scene *new_scene = static_cast<Scene *>(BKE_id_new_nomain(ID_SCE, nullptr));
  sequence->legacy_scene_data = *new_scene;
  MEM_freeN(new_scene);

  sequence->legacy_scene_data.id.flag |= ID_FLAG_EMBEDDED_DATA;

  SEQ_editing_ensure(&sequence->legacy_scene_data);
}

static void sequence_data_free(ID *id)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  BKE_libblock_free_datablock(&sequence->legacy_scene_data.id, 0);
}

static void sequence_foreach_id(ID *id, LibraryForeachIDData *data)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);
  ID *scene_id = &sequence->legacy_scene_data.id;
  // BKE_library_foreach_ID_embedded(data, &scene_id);
}

static void sequence_blend_write(BlendWriter *writer, ID *id, const void *id_address)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  BLO_write_id_struct(writer, Sequence, id_address, &sequence->id);
  BKE_id_blend_write(writer, &sequence->id);

  Scene *scene = &sequence->legacy_scene_data;
  BLO_write_struct(writer, Scene, scene);

  Editing *ed = scene->ed;
  BLI_assert(ed);
  BLO_write_struct(writer, Editing, ed);

  SEQ_blend_write(writer, &ed->seqbase);
  LISTBASE_FOREACH (SeqTimelineChannel *, channel, &ed->channels) {
    BLO_write_struct(writer, SeqTimelineChannel, channel);
  }
  /* new; meta stack too, even when its nasty restore code */
  LISTBASE_FOREACH (MetaStack *, ms, &ed->metastack) {
    BLO_write_struct(writer, MetaStack, ms);
  }

  // BKE_image_format_blend_write(writer, &scene->r.im_format);
}

static void link_recurs_seq(BlendDataReader *reader, ListBase *lb)
{
  BLO_read_struct_list(reader, Strip, lb);

  LISTBASE_FOREACH_MUTABLE (Strip *, seq, lb) {
    /* Sanity check. */
    if (!SEQ_is_valid_strip_channel(seq)) {
      BLI_freelinkN(lb, seq);
      BLO_read_data_reports(reader)->count.sequence_strips_skipped++;
    }
    else if (seq->seqbase.first) {
      link_recurs_seq(reader, &seq->seqbase);
    }
  }
}

static void sequence_read_data(BlendDataReader *reader, Scene *sce)
{
  Editing *ed = sce->ed;

  ListBase *old_seqbasep = &sce->ed->seqbase;
  ListBase *old_displayed_channels = &sce->ed->channels;

  ed->act_seq = static_cast<Strip *>(
      BLO_read_get_new_data_address_no_us(reader, ed->act_seq, sizeof(Strip)));
  ed->cache = nullptr;
  ed->prefetch_job = nullptr;
  ed->runtime.strip_lookup = nullptr;
  ed->runtime.media_presence = nullptr;
  ed->runtime.thumbnail_cache = nullptr;

  /* recursive link sequences, lb will be correctly initialized */
  link_recurs_seq(reader, &ed->seqbase);

  /* Read in sequence member data. */
  SEQ_blend_read(reader, &ed->seqbase);
  BLO_read_struct_list(reader, SeqTimelineChannel, &ed->channels);

  /* link metastack, slight abuse of structs here,
   * have to restore pointer to internal part in struct */
  {
    void *seqbase_poin;
    void *channels_poin;
    /* This whole thing with seqbasep offsets is really not good
     * and prevents changes to the Sequence struct. A more correct approach
     * would be to calculate offset using sDNA from the file (NOT from the
     * current Blender). Even better would be having some sort of dedicated
     * map of seqbase pointers to avoid this offset magic. */
    constexpr intptr_t seqbase_offset = offsetof(Strip, seqbase);
    constexpr intptr_t channels_offset = offsetof(Strip, channels);
#if ARCH_CPU_64_BITS
    static_assert(seqbase_offset == 264, "Sequence seqbase member offset cannot be changed");
    static_assert(channels_offset == 280, "Sequence channels member offset cannot be changed");
#else
    static_assert(seqbase_offset == 204, "Sequence seqbase member offset cannot be changed");
    static_assert(channels_offset == 212, "Sequence channels member offset cannot be changed");
#endif

    /* seqbase root pointer */
    if (ed->seqbasep == old_seqbasep) {
      ed->seqbasep = &ed->seqbase;
    }
    else {
      seqbase_poin = POINTER_OFFSET(ed->seqbasep, -seqbase_offset);

      seqbase_poin = BLO_read_get_new_data_address_no_us(reader, seqbase_poin, sizeof(Strip));

      if (seqbase_poin) {
        ed->seqbasep = (ListBase *)POINTER_OFFSET(seqbase_poin, seqbase_offset);
      }
      else {
        ed->seqbasep = &ed->seqbase;
      }
    }

    /* Active channels root pointer. */
    if (ELEM(ed->displayed_channels, old_displayed_channels, nullptr)) {
      ed->displayed_channels = &ed->channels;
    }
    else {
      channels_poin = POINTER_OFFSET(ed->displayed_channels, -channels_offset);
      channels_poin = BLO_read_get_new_data_address_no_us(
          reader, channels_poin, sizeof(SeqTimelineChannel));

      if (channels_poin) {
        ed->displayed_channels = (ListBase *)POINTER_OFFSET(channels_poin, channels_offset);
      }
      else {
        ed->displayed_channels = &ed->channels;
      }
    }

    /* stack */
    BLO_read_struct_list(reader, MetaStack, &(ed->metastack));

    LISTBASE_FOREACH (MetaStack *, ms, &ed->metastack) {
      BLO_read_struct(reader, Strip, &ms->parseq);

      if (ms->oldbasep == old_seqbasep) {
        ms->oldbasep = &ed->seqbase;
      }
      else {
        seqbase_poin = POINTER_OFFSET(ms->oldbasep, -seqbase_offset);
        seqbase_poin = BLO_read_get_new_data_address_no_us(reader, seqbase_poin, sizeof(Strip));
        if (seqbase_poin) {
          ms->oldbasep = (ListBase *)POINTER_OFFSET(seqbase_poin, seqbase_offset);
        }
        else {
          ms->oldbasep = &ed->seqbase;
        }
      }

      if (ELEM(ms->old_channels, old_displayed_channels, nullptr)) {
        ms->old_channels = &ed->channels;
      }
      else {
        channels_poin = POINTER_OFFSET(ms->old_channels, -channels_offset);
        channels_poin = BLO_read_get_new_data_address_no_us(
            reader, channels_poin, sizeof(SeqTimelineChannel));

        if (channels_poin) {
          ms->old_channels = (ListBase *)POINTER_OFFSET(channels_poin, channels_offset);
        }
        else {
          ms->old_channels = &ed->channels;
        }
      }
    }
  }
}

static void sequence_blend_read_data(BlendDataReader *reader, ID *id)
{
  Sequence *sequence = reinterpret_cast<Sequence *>(id);

  BLO_read_struct(reader, Scene, &sequence->legacy_scene_data);
  Scene *scene = &sequence->legacy_scene_data;

  BLO_read_struct(reader, Editing, scene->ed);
  Editing *ed = scene->ed;
  BLI_assert(ed);
  sequence_read_data(reader, scene);

  // BKE_image_format_blend_read_data(reader, &scene->r.im_format);
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
    /*foreach_id*/ sequence_foreach_id,
    /*foreach_cache*/ nullptr,
    /*foreach_path*/ nullptr,
    /*owner_pointer_get*/ nullptr,

    /*blend_write*/ sequence_blend_write,
    /*blend_read_data*/ sequence_blend_read_data,
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
