# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2023, The SPA Studios. All rights reserved.

"""
Addon preferences management.
"""

import bpy

from spa_sequencer.utils import register_classes, unregister_classes


class SPASequencerAddonPreferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    shot_template_prefix: bpy.props.StringProperty(
        name="Shot Template Prefix",
        description="Scene name prefix that identifies Shot Templates",
        default="TEMPLATE_SHOT",
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "shot_template_prefix")
        box = layout.box()
        col = box.column()
        has_otio = False
        if not has_otio:
            try:
                import opentimelineio
                has_otio = True
            except:
                pass
            col.operator("sequencer.pip_install", text="Install OpenTimelineIO").package = "opentimelineio"
        else:
            col.label(text="OpenTimelineIO Installed!")
            col.label(text="You might need to restart Blender.")

def python_exec():
    import os
    import bpy
    try:
        # 2.92 and older
        path = bpy.app.binary_path_python
    except AttributeError:
        # 2.93 and later
        import sys
        path = sys.executable
    return os.path.abspath(path)

class SEQ_OT_PipInstall(bpy.types.Operator):
    """Install the package by calling pip install"""
    bl_idname = 'sequencer.pip_install'
    bl_label = "Install the package"
    bl_options = {'REGISTER', 'INTERNAL'}

    package : bpy.props.StringProperty(name = "Package names")

    def execute(self, context):
        # https://github.com/robertguetzkow/blender-python-examples/tree/master/add_ons/install_dependencies
        environ_copy = dict(os.environ)
        environ_copy["PYTHONNOUSERSITE"] = "1"  # is set to disallow pip from checking the user site-packages
        cmd = [python_exec(), '-m', 'pip', 'install', '--upgrade'] + self.package.split(" ")
        ok = subprocess.call(cmd, env=environ_copy) == 0
        if ok:
            first_install = self.package in sv_dependencies and sv_dependencies[self.package] is None
            if first_install:
                self.report({'INFO'}, f"{self.package} installed successfully. Please restart Blender to see effect.")
            else:
                self.report({'INFO'}, f"{self.package} upgraded successfully. Please restart Blender to see effect.")
            return {'FINISHED'}
        else:
            self.report({'ERROR'}, "Cannot install %s, see console output for details" % self.package)
            return {'CANCELLED'}



def get_addon_prefs() -> SPASequencerAddonPreferences:
    """Get the Addon Preferences instance."""
    return bpy.context.preferences.addons[__package__].preferences


classes = (SPASequencerAddonPreferences, SEQ_OT_PipInstall)


def register():
    register_classes(classes)


def unregister():
    unregister_classes(classes)
