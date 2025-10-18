# SPDX-License-Identifier: GPL-3.0-or-later
# Thanks to Znight and Spa Studios for the work of making this real

import bpy

from ..utils import register_classes, unregister_classes


class SEQUENCER_MT_shot(bpy.types.Menu):
    """Scene operators menu"""

    bl_idname = "SEQUENCER_MT_shot"
    bl_label = "Scene"

    def draw(self, context):
        layout = self.layout

        # Operator to slip and adjust time of active scene strip
        layout.operator("sequencer.shot_timing_adjust", icon="TIME")

        # Operator to update to active scene strip
        layout.operator('sequencer.change_3d_view_scene', text='Toggle Active Scene Strip', icon="FILE_REFRESH")

        # Operator to playback the master scene
        layout.operator('wm.timeline_sync_play_master', icon="PLAY")


def draw_MT_shot(self, context):
    layout = self.layout
    layout.menu(SEQUENCER_MT_shot.bl_idname)


classes = (
    SEQUENCER_MT_shot,
)


def register():
    register_classes(classes)

    # Add SEQUENCER_MT_shot to the sequencer editor menus
    bpy.types.SEQUENCER_MT_editor_menus.append(draw_MT_shot)


def unregister():
    unregister_classes(classes)

    # Remove SEQUENCER_MT_shot from the sequencer editor menus
    bpy.types.SEQUENCER_MT_editor_menus.remove(draw_MT_shot)
