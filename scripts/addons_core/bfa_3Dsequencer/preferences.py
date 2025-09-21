# SPDX-License-Identifier: GPL-3.0-or-later


"""
Addon preferences management.
"""

import bpy

from ..utils import register_classes, unregister_classes


class SPASequencerAddonPreferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    shot_template_prefix: bpy.props.StringProperty(
        name="Shot Template Prefix",
        description="Scene name prefix that identifies Scene Templates",
        default="TEMPLATE_SHOT",
    )

    debug_mode: bpy.props.BoolProperty(
        name="Debug mode",
        description="Use for debugging, for other 3D Sequencer feature require restart",
        default=False,
    )

    def draw(self, context):
        self.layout.prop(self, "shot_template_prefix")


def get_addon_prefs() -> SPASequencerAddonPreferences:
    """Get the Addon Preferences instance."""
    return bpy.context.preferences.addons[__package__].preferences


classes = (SPASequencerAddonPreferences,)


def register():
    register_classes(classes)


def unregister():
    unregister_classes(classes)
