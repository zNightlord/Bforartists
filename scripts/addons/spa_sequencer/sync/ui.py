# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2023, The SPA Studios. All rights reserved.

import bpy

from spa_sequencer.sync.core import get_sync_settings
from spa_sequencer.utils import register_classes, unregister_classes


class SEQUENCER_PT_SyncPanel(bpy.types.Panel):
    """Timeline Synchronization Panel."""

    bl_label = "Timeline Sync"
    bl_space_type = "SEQUENCE_EDITOR"
    bl_region_type = "UI"
    bl_region_type = 'HEADER'

    def draw(self, context):
        layout = self.layout
        col = layout.column()
        col.use_property_split = True
        col.use_property_decorate = False
        settings = get_sync_settings()
        col.operator(
            "wm.timeline_sync_toggle",
            text="Synchronize",
            icon="UV_SYNC_SELECT",
            depress=settings.enabled,
        )
        col.prop(settings, "master_scene")
        
        col = layout.column()
        box = col.box()
        col = box.column()
        col.enabled = bool(settings.master_scene)
        self.draw_settings(col)

    def draw_settings(self, layout):
        settings = get_sync_settings()
        layout.prop(settings, "keep_gpencil_tool_settings")
        layout.prop(settings, "bidirectional")
        layout.prop(settings, "use_preview_range")
        layout.prop(settings, "sync_all_windows")
        layout.prop(settings, "active_follows_playhead")

classes = (
    SEQUENCER_PT_SyncPanel,
)


def register():
    register_classes(classes)


def unregister():
    unregister_classes(classes)
