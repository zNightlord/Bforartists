/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edsequence
 */

#include "BLI_listbase.h"
#include "BLI_vector_set.hh"

#include "BKE_context.hh"
#include "BKE_lib_id.hh"
#include "BKE_main.hh"

#include "DNA_space_types.h"

#include "DEG_depsgraph.hh"

#include "ED_sequence.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"
#include "RNA_prototypes.hh"

#include "SEQ_channels.hh"
#include "SEQ_iterator.hh"
#include "SEQ_sequencer.hh"
#include "SEQ_time.hh"

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

static void sequence_delete(bContext &C, Main &bmain, Sequence &sequence)
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
}

namespace blender::ed::sequence {

void sync_scene_strip(bContext *C, Scene *sequence_scene)
{
  using namespace blender;
  wmWindow *win = CTX_wm_window(C);
  Scene *active_scene = WM_window_get_active_scene(win);

  Editing *ed = seq::editing_get(sequence_scene);
  ListBase *seqbase = seq::active_seqbase_get(ed);
  ListBase *channels = seq::channels_displayed_get(ed);
  VectorSet<Strip *> render_strips = seq::query_rendered_strips(
      sequence_scene, channels, seqbase, sequence_scene->r.cfra, 0);
  Vector<Strip *> strips = render_strips.extract_vector();
  /* Sort strips by channel. */
  std::sort(strips.begin(), strips.end(), [](const Strip *a, const Strip *b) {
    return a->machine > b->machine;
  });
  const Strip *scene_strip = [&]() -> const Strip * {
    for (const Strip *strip : strips) {
      if (strip->type == STRIP_TYPE_SCENE) {
        return strip;
      }
    }
    return nullptr;
  }();
  if (scene_strip && scene_strip->scene) {
    if (active_scene != scene_strip->scene) {
      /* Change active scene in window. */
      Main *bmain = CTX_data_main(C);
      WM_window_set_active_scene(bmain, C, win, scene_strip->scene);
      active_scene = scene_strip->scene;
    }
    if (scene_strip->scene_camera) {
      /* Update the camera in any 3D view that uses camera view. */
      PointerRNA camera_ptr = RNA_id_pointer_create(&scene_strip->scene_camera->id);
      bScreen *screen = WM_window_get_active_screen(win);
      LISTBASE_FOREACH (ScrArea *, area, &screen->areabase) {
        LISTBASE_FOREACH (SpaceLink *, sl, &area->spacedata) {
          if (sl->spacetype != SPACE_VIEW3D) {
            continue;
          }
          View3D *view3d = reinterpret_cast<View3D *>(sl);
          if (view3d->camera == scene_strip->scene_camera) {
            continue;
          }
          PointerRNA view3d_ptr = RNA_pointer_create_discrete(
              &screen->id, &RNA_SpaceView3D, view3d);
          RNA_pointer_set(&view3d_ptr, "camera", camera_ptr);
        }
      }
    }

    float frame_index = seq::give_frame_index(sequence_scene, scene_strip, sequence_scene->r.cfra);
    if (active_scene->r.flag & SCER_SHOW_SUBFRAME) {
      active_scene->r.cfra = int(frame_index);
      active_scene->r.subframe = frame_index - int(frame_index);
    }
    else {
      active_scene->r.cfra = round_fl_to_int(frame_index);
      active_scene->r.subframe = 0.0f;
    }
    FRAMENUMBER_MIN_CLAMP(active_scene->r.cfra);

    DEG_id_tag_update(&active_scene->id, ID_RECALC_FRAME_CHANGE);
    WM_event_add_notifier(C, NC_SCENE | ND_FRAME, active_scene);
  }
}

}  // namespace blender::ed::sequence

/** \} */

/* -------------------------------------------------------------------- */
/** \name Sequence New Operator
 * \{ */

static wmOperatorStatus sequence_new_exec(bContext *C, wmOperator *op)
{
  Main *bmain = CTX_data_main(C);
  wmWindow *win = CTX_wm_window(C);
  const AddSequenceMode mode = AddSequenceMode(RNA_enum_get(op->ptr, "mode"));

  sequence_add(*bmain, *C, *win, mode);

  return OPERATOR_FINISHED;
}

static wmOperatorStatus sequence_new_invoke(bContext *C, wmOperator *op, const wmEvent *event)
{
  if (!WM_window_get_active_sequence(CTX_wm_window(C))) {
    /* Always create a blank sequence if there is none in the window. */
    RNA_enum_set(op->ptr, "mode", int(AddSequenceMode::Blank));
    return sequence_new_exec(C, op);
  }
  return WM_menu_invoke(C, op, event);
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
  ot->invoke = sequence_new_invoke;

  /* flags */
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  /* properties */
  ot->prop = RNA_def_enum(
      ot->srna, "mode", sequence_new_items, int(AddSequenceMode::CopySettings), "Mode", "");
  RNA_def_property_translation_context(ot->prop, BLT_I18NCONTEXT_ID_SEQUENCE);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Sequence Delete Operator
 * \{ */

static bool sequence_delete_poll(bContext *C)
{
  Main *bmain = CTX_data_main(C);
  if (bmain) {
    return true;
  }
  return false;
}

static wmOperatorStatus sequence_delete_exec(bContext *C, wmOperator * /*op*/)
{
  Main *bmain = CTX_data_main(C);
  Sequence *sequence = CTX_data_sequence(C);

  sequence_delete(*C, *bmain, *sequence);

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
