# SPDX-License-Identifier: GPL-3.0-or-later
# Thanks to Znight and Spa Studios for the work of making this real

from typing import Optional

import bpy

from ..preferences import get_addon_prefs
from ..scene.core import (
    adjust_shot_duration,
    delete_scene,
    duplicate_scene,
    get_valid_shot_scenes,
    rename_scene,
    slip_shot_content,
)
from ..scene.naming import shot_naming, ShotNamingProperty
from ..sync.core import (
    get_sync_master_strip,
    get_sync_settings,
    remap_frame_value,
)
from ..utils import register_classes, unregister_classes


def get_last_sequence(
    strips: list[bpy.types.Strip],
) -> Optional[bpy.types.Strip]:
    """Get the last sequence, i.e. the one with the greatest final frame number."""
    return max(strips, key=lambda x: x.frame_final_end) if strips else None


def get_last_used_frame(
    strips: list[bpy.types.Strip], scene: bpy.types.Scene
) -> int:
    """
    Get the last used internal frame of `scene` from the given list of `strips`.
    """
    scene_strips = [
        s
        for s in strips
        if isinstance(s, bpy.types.SceneStrip) and s.scene == scene
    ]

    if not scene_strips:
        return scene.frame_start - 1

    return max(remap_frame_value(s.frame_final_end - 1, s) for s in scene_strips)


def get_selected_scene_sequences(
    strips: list[bpy.types.Strip],
) -> list[bpy.types.SceneStrip]:
    """
    :param strips: The strips to consider.
    :return: The list of selected scene sequence strips.
    """
    return [s for s in strips if isinstance(s, bpy.types.SceneStrip) and s.select]


def ensure_sequencer_frame_visible(context: bpy.types.Context, frame: int):
    """Ensure `frame` is visible within context's sequencer area.

    :param context: The context - context.area has to be a SequenceEditor.
    :param frame: The frame to consider.
    """
    if not context.area.type == "SEQUENCE_EDITOR":
        return
    frame_coord = context.region.view2d.view_to_region(frame, 0, clip=False)[0]
    if frame_coord < 0 or frame_coord > context.region.width:
        # Temp override current frame value and move view to frame.
        frame_old = context.scene.frame_current
        context.scene.frame_current = frame
        bpy.ops.sequencer.view_frame()
        context.scene.frame_current = frame_old


class SEQUENCER_OT_shot_timing_adjust(bpy.types.Operator):
    bl_idname = "sequencer.shot_timing_adjust"
    bl_label = "Adjust Timing"
    bl_description = "Adjust the timing of the active scene interactively"
    bl_options = {"GRAB_CURSOR_X", "BLOCKING", "UNDO"}

    offset: bpy.props.IntProperty(
        name="Offset",
    )

    mode: bpy.props.EnumProperty(
        name="Mode",
        items=(
            ("DURATION", "Duration", "Adjust strip duration"),
            ("SLIP", "Slip", "Slip strip content"),
        ),
        default="DURATION",
        options={"SKIP_SAVE"},
    )

    strip_handle: bpy.props.EnumProperty(
        name="Strip Adjustment Handle",
        items=(
            ("LEFT", "Left", "Left"),
            ("RIGHT", "Right", "Right"),
        ),
        default="RIGHT",
        options={"SKIP_SAVE"},
    )

    @classmethod
    def poll(cls, context: bpy.types.Context):
        return cls.get_active_strip(context) is not None

    @staticmethod
    def get_active_strip(
        context: bpy.types.Context,
    ) -> Optional[bpy.types.SceneStrip]:
        if context.area.type == "DOPESHEET_EDITOR":
            strip = get_sync_master_strip(use_cache=True)[0]
            return strip if strip and strip.scene == context.window.scene else None
        elif context.scene.sequence_editor and isinstance(
            context.scene.sequence_editor.active_strip, bpy.types.SceneStrip
        ):
            return context.scene.sequence_editor.active_strip

        return None

    def setup(self, context: bpy.types.Context):
        self.strip = self.get_active_strip(context)
        if not self.strip:
            self.report({"ERROR"}, "No current Scene Strip")
            return False
        return True

    def invoke(self, context: bpy.types.Context, event: bpy.types.Event):
        if not self.setup(context):
            return {"CANCELLED"}

        context.window.cursor_modal_set("SCROLL_X")

        self.start_mouse_coords = context.region.view2d.region_to_view(
            x=event.mouse_region_x, y=event.mouse_region_y
        )

        if context.area.type == "SEQUENCE_EDITOR":
            self.strip_handle = "LEFT" if self.strip.select_left_handle else "RIGHT"

        self.original_strip_duration = self.strip.frame_final_duration
        self.original_strip_scene_end = self.strip.scene.frame_end
        self.original_strip_offset_start = self.strip.frame_offset_start
        self.original_edit_frame_end = get_sync_settings().master_scene.frame_end

        context.window_manager.modal_handler_add(self)
        return {"RUNNING_MODAL"}

    def update_header_text(self, context, event):
        text = (
            f"Offset: {self.offset}"
            f" | New Scene Duration: {self.strip.frame_final_duration}"
        )
        context.area.header_text_set(text)

    def modal(self, context: bpy.types.Context, event: bpy.types.Event):
        self.update_header_text(context, event)
        # Cancel
        if event.type in {"RIGHTMOUSE", "ESC"}:
            self.cancel(context)
            return {"CANCELLED"}
        # Validate
        elif (event.type in {"LEFTMOUSE"} and event.value in {"PRESS", "RELEASE"}) or (
            event.type in {"RET", "NUMPAD_ENTER"} and event.value in {"PRESS"}
        ):
            if self.offset == 0:
                self.cancel(context)
                return {"CANCELLED"}
            self.restore_ui(context)
            return {"FINISHED"}
        # Update
        elif event.type in {"MOUSEMOVE"}:
            mouse_coords = context.region.view2d.region_to_view(
                x=event.mouse_region_x, y=event.mouse_region_y
            )
            offset = int(mouse_coords[0] - self.start_mouse_coords[0])
            if offset != self.offset:
                self.offset = offset
                self.execute(context)

        return {"RUNNING_MODAL"}

    def execute(self, context: bpy.types.Context):
        if not self.options.is_invoke:
            if not self.setup(context):
                return {"CANCELLED"}

        from_frame_start = self.strip_handle == "LEFT"
        select = self.mode == "DURATION"
        self.strip.select_left_handle = select and from_frame_start
        self.strip.select_right_handle = select and not from_frame_start
        # Adjust offset sign to match direction.
        # For instance, a positive offset (going to the right in modal):
        #  - SHRINKS the strip if using left handle (from frame start)
        #  - EXTENDS the strip otherwise (from frame end)
        offset = -self.offset if from_frame_start else self.offset
        # Compute current absolute offset from original duration
        if self.mode == "SLIP":
            delta = self.strip.frame_offset_start - self.original_strip_offset_start
            slip_shot_content(self.strip, offset - delta, clamp_start=True)
        else:
            delta = self.strip.frame_final_duration - self.original_strip_duration
            adjust_shot_duration(self.strip, offset - delta, from_frame_start)

        edit_scene = get_sync_settings().master_scene
        if from_frame_start or self.mode == "SLIP":
            # NOTE: When adjusting from frame start, the current frame does not change.
            #       Set time to a non meaningful value before re-setting the correct frame
            #       value, to trigger time-dependent updates (e.g: synchronization).
            edit_scene.frame_set(-1)
            update_frame = self.strip.frame_final_start
        else:
            update_frame = self.strip.frame_final_end - 1

        # Set sequencer's frame to strip's new end frame
        edit_scene.frame_set(update_frame)

        # Update both edit and internal scene's end frame if going past original ones
        frame_end = self.strip.frame_final_end - 1
        edit_scene.frame_end = max(frame_end, self.original_edit_frame_end)
        self.strip.scene.frame_end = max(
            remap_frame_value(frame_end, self.strip),
            self.original_strip_scene_end,
        )

        return {"FINISHED"}

    def restore_ui(self, context: bpy.types.Context):
        context.area.header_text_set(None)
        context.window.cursor_modal_restore()

    def cancel(self, context: bpy.types.Context):
        if self.offset:
            self.offset = 0
            self.execute(context)
        # Restore scenes' original end frames
        context.scene.frame_end = self.original_edit_frame_end
        self.strip.scene.frame_end = self.original_strip_scene_end
        self.restore_ui(context)


classes = (
    SEQUENCER_OT_shot_timing_adjust,
)


def register():
    register_classes(classes)


def unregister():
    unregister_classes(classes)
