/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edsequence
 */

#include "BLI_listbase.h"

#include "BKE_context.hh"
#include "BKE_lib_id.hh"
#include "BKE_main.hh"

#include "ED_sequence.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"

#include "BLT_translation.hh"

#include "WM_api.hh"
#include "WM_types.hh"

/* -------------------------------------------------------------------- */
/** \name Sequence Utilities
 * \{ */

static Sequence *sequence_add(Main &bmain, bContext &C, wmWindow &win, const AddSequenceMode mode)
{
  // Sequence &sequence_old = *WM_window_get_active_sequence(&win);
  Sequence *sequence_new = [&]() -> Sequence * {
    switch (mode) {
      case AddSequenceMode::Blank:
        return BKE_sequence_add(bmain, DATA_("Sequence"));
      case AddSequenceMode::CopySettings:
        /* FIXME: Copy Settings. */
        return BKE_sequence_add(bmain, DATA_("Sequence"));
      case AddSequenceMode::FullCopy:
        /* FIXME: Full Copy. */
        return BKE_sequence_add(bmain, DATA_("Sequence"));
      default:
        BLI_assert_unreachable();
    }
    return nullptr;
  }();

  WM_window_set_active_sequence(&bmain, &C, &win, sequence_new);

  /* TODO: Correct notifier! */
  WM_event_add_notifier(&C, NC_WINDOW, nullptr);

  return sequence_new;
}

static bool sequence_delete(bContext &C, Main &bmain, Sequence &sequence)
{
  Sequence *sequence_new = [&]() -> Sequence * {
    if (sequence.id.prev) {
      return static_cast<Sequence *>(sequence.id.prev);
    }
    if (sequence.id.next) {
      return static_cast<Sequence *>(sequence.id.next);
    }
    return nullptr;
  }();

  if (!sequence_new) {
    return false;
  }

  wmWindowManager *wm = static_cast<wmWindowManager *>(bmain.wm.first);

  LISTBASE_FOREACH (wmWindow *, win, &wm->windows) {
    if (win->parent != nullptr) { /* We only care about main windows here... */
      continue;
    }
    if (win->sequence == &sequence) {
      WM_window_set_active_sequence(&bmain, &C, win, sequence_new);
    }
  }

  BKE_id_delete(&bmain, &sequence);

  return true;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Sequence New Operator
 * \{ */

static int sequence_new_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);
  wmWindow *win = CTX_wm_window(C);
  const AddSequenceMode mode = AddSequenceMode(RNA_enum_get(op->ptr, "mode"));

  sequence_add(*bmain, *C, *win, mode);

  return OPERATOR_FINISHED;
}

static EnumPropertyItem sequence_new_items[] = {
    {int(AddSequenceMode::Blank),
     "BLANK",
     0,
     "Blank",
     "Add a new, empty sequence with default settings"},
    {int(AddSequenceMode::CopySettings),
     "EMPTY",
     0,
     "Copy Settings",
     "Add a new, empty sequence, and copy settings from the current sequence"},
    {int(AddSequenceMode::FullCopy),
     "FULL_COPY",
     0,
     "Full Copy",
     "Make a full copy of the current sequence"},
    {0, nullptr, 0, nullptr, nullptr},
};

static void SEQUENCE_OT_new(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "New Sequence";
  ot->description = "Add new sequence";
  ot->idname = "SEQUENCE_OT_new";

  /* api callbacks */
  ot->exec = sequence_new_exec;
  ot->invoke = WM_menu_invoke;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  /* properties */
  ot->prop = RNA_def_enum(
      ot->srna, "mode", sequence_new_items, int(AddSequenceMode::CopySettings), "Mode", "");
  RNA_def_property_translation_context(ot->prop, BLT_I18NCONTEXT_ID_SEQUENCE);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Scene Delete Operator
 * \{ */

static bool sequence_delete_poll(bContext *C)
{
  Main *bmain = CTX_data_main(C);
  Sequence *sequence = CTX_data_sequence(C);
  if (bmain && sequence) {
    return BKE_sequence_can_be_removed(*bmain, *sequence);
  }
  return false;
}

static int sequence_delete_exec(bContext *C, wmOperator * /*op*/)
{
  Main *bmain = CTX_data_main(C);
  Sequence *sequence = CTX_data_sequence(C);

  if (!sequence_delete(*C, *bmain, *sequence)) {
    return OPERATOR_CANCELLED;
  }

  /* TODO: Correct notifiers. */
  WM_event_add_notifier(C, NC_WINDOW, nullptr);

  return OPERATOR_FINISHED;
}

static void SEQUENCE_OT_delete(wmOperatorType *ot)
{
  /* identifiers */
  ot->name = "Delete Sequence";
  ot->description = "Delete active sequence";
  ot->idname = "SEQUENCE_OT_delete";

  /* api callbacks */
  ot->exec = sequence_delete_exec;
  ot->poll = sequence_delete_poll;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Registration
 * \{ */

void ED_operatortypes_sequence()
{
  WM_operatortype_append(SEQUENCE_OT_new);
  WM_operatortype_append(SEQUENCE_OT_delete);
}

/** \} */
