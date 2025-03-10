# SPDX-FileCopyrightText: 2016-2022 Blender Foundation
#
# SPDX-License-Identifier: GPL-2.0-or-later

# ----------------------------------------------------------
# Author: Antonio Vazquez (antonioya)
# ----------------------------------------------------------

# ----------------------------------------------
# Import modules
# ----------------------------------------------
if "bpy" in locals():
    import importlib

    importlib.reload(measureit_main)
    # print("measureit: Reloaded multifiles")
else:
    from . import measureit_main
    # print("measureit: Imported multifiles")

# noinspection PyUnresolvedReferences
import bpy
from bpy.types import (
        AddonPreferences,
        Scene,
        WindowManager,
        )
from bpy.props import (
        FloatVectorProperty,
        IntProperty,
        BoolProperty,
        StringProperty,
        FloatProperty,
        EnumProperty,
        )

from .measureit_sum_edge_lengths import *

# --------------------------------------------------------------
# Register all operators and panels
# --------------------------------------------------------------

# Add-ons Preferences Update Panel

# Define Panel classes for updating
panels = (
        measureit_main.MEASUREIT_PT_Edit,
        measureit_main.MEASUREIT_PT_Main,
        measureit_main.MEASUREIT_PT_Conf,
        measureit_main.MEASUREIT_PT_Render,
        )


def update_panel(self, context):
    message = "MeasureIt: Updating Panel locations has failed"
    try:
        for panel in panels:
            if "bl_rna" in panel.__dict__:
                bpy.utils.unregister_class(panel)

        for panel in panels:
            panel.bl_category = context.preferences.addons[__name__].preferences.category
            bpy.utils.register_class(panel)

    except Exception as e:
        print("\n[{}]\n{}\n\nError:\n{}".format(__name__, message, e))
        pass


class Measure_Pref(AddonPreferences):
    # this must match the addon name, use '__package__'
    # when defining this in a submodule of a python package.
    bl_idname = __name__

    category: StringProperty(
            name="Tab Category",
            description="Choose a name for the category of the panel",
            default="Display",
            update=update_panel
            )

    def draw(self, context):
        layout = self.layout

        row = layout.row()
        col = row.column()
        col.label(text="Tab Category:")
        col.prop(self, "category", text="")


# Define menu
# noinspection PyUnusedLocal
classes = (
    measureit_main.MEASUREIT_OT_RunHintDisplay,
    measureit_main.MEASUREIT_OT_AddSegment,
    measureit_main.MEASUREIT_OT_AddArea,
    measureit_main.MEASUREIT_OT_AddSegmentOrto,
    measureit_main.MEASUREIT_OT_AddAngle,
    measureit_main.MEASUREIT_OT_AddArc,
    measureit_main.MEASUREIT_OT_AddLabel,
    measureit_main.MEASUREIT_OT_AddNote,
    measureit_main.MEASUREIT_OT_AddLink,
    measureit_main.MEASUREIT_OT_AddOrigin,
    measureit_main.MEASUREIT_OT_DeleteSegment,
    measureit_main.MEASUREIT_OT_DeleteAllSegment,
    measureit_main.MEASUREIT_OT_DeleteAllSum,
    measureit_main.MEASUREIT_OT_RenderSegment,
    measureit_main.MEASUREIT_OT_ExpandAllSegment,
    measureit_main.MEASUREIT_OT_CollapseAllSegment,
    measureit_main.MEASUREIT_PT_Main,
    measureit_main.MEASUREIT_PT_Edit,
    measureit_main.MEASUREIT_PT_Conf,
    measureit_main.MEASUREIT_PT_Render,
#    Measure_Pref,
)

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)

    # Define properties
    Scene.measureit_default_color = FloatVectorProperty(
        name="Default color",
        description="Default Color",
        default=(0.173, 0.545, 1.0, 1.0),
        min=0.1,
        max=1,
        subtype='COLOR',
        size=4)
    Scene.measureit_font_size = IntProperty(name="Text Size",
                                            description="Default text size",
                                            default=14, min=10, max=150)
    Scene.measureit_hint_space = FloatProperty(name='Separation', min=0, max=100, default=0.1,
                                               precision=3,
                                               description="Default distance to display measure")
    Scene.measureit_gl_ghost = BoolProperty(name="All",
                                            description="Display measures for all objects,"
                                                        " not only selected",
                                            default=True)
    Scene.measureit_gl_txt = StringProperty(name="Text", maxlen=256,
                                            description="Short description (use | for line break)")

    Scene.measureit_gl_precision = IntProperty(name='Precision', min=0, max=5, default=2,
                                               description="Number of decimal precision")
    Scene.measureit_gl_show_d = BoolProperty(name="ShowDist",
                                             description="Display distances",
                                             default=True)
    Scene.measureit_gl_show_n = BoolProperty(name="ShowName",
                                             description="Display texts",
                                             default=False)
    Scene.measureit_scale = BoolProperty(name="Scale",
                                         description="Use scale factor",
                                         default=False)
    Scene.measureit_scale_factor = FloatProperty(name='Factor', min=0.001, max=9999999,
                                                 default=1.0,
                                                 precision=3,
                                                 description="Scale factor 1:x")
    Scene.measureit_scale_color = FloatVectorProperty(name="Scale color",
                                                      description="Scale Color",
                                                      default=(1, 1, 0, 1.0),
                                                      min=0.1,
                                                      max=1,
                                                      subtype='COLOR',
                                                      size=4)
    Scene.measureit_scale_font = IntProperty(name="Font",
                                             description="Text size",
                                             default=14, min=10, max=150)
    Scene.measureit_scale_pos_x = IntProperty(name="X Position",
                                              description="Margin on the X axis",
                                              default=5,
                                              min=0,
                                              max=100)
    Scene.measureit_scale_pos_y = IntProperty(name="Y Position",
                                              description="Margin on the Y axis",
                                              default=5,
                                              min=0,
                                              max=100)
    Scene.measureit_gl_scaletxt = StringProperty(name="ScaleText", maxlen=48,
                                                 description="Scale title",
                                                 default="Scale:")
    Scene.measureit_scale_precision = IntProperty(name='Precision', min=0, max=5, default=0,
                                                  description="Number of decimal precision")
    Scene.measureit_ovr = BoolProperty(name="Override",
                                       description="Override colors and fonts",
                                       default=False)
    Scene.measureit_ovr_font = IntProperty(name="Font",
                                           description="Override text size",
                                           default=14, min=10, max=150)
    Scene.measureit_ovr_color = FloatVectorProperty(name="Override color",
                                                    description="Override Color",
                                                    default=(1, 0, 0, 1.0),
                                                    min=0.1,
                                                    max=1,
                                                    subtype='COLOR',
                                                    size=4)
    Scene.measureit_ovr_width = IntProperty(name='Override width', min=1, max=10, default=1,
                                            description='override line width')

    Scene.measureit_ovr_font_rotation = IntProperty(name='Rotate', min=0, max=360, default=0,
                                                    description="Text rotation in degrees")
    Scene.measureit_ovr_font_align = EnumProperty(items=(('L', "Left Align", "Use current render"),
                                                         ('C', "Center Align", ""),
                                                         ('R', "Right Align", "")),
                                                  name="Align Font",
                                                  description="Set Font Alignment")
    Scene.measureit_units = EnumProperty(items=(('1', "Automatic", "Use scene units"),
                                                ('2', "Meters", ""),
                                                ('3', "Centimeters", ""),
                                                ('4', "Millimeters", ""),
                                                ('5', "Feet", ""),
                                                ('6', "Inches", "")),
                                         name="Units",
                                         default="2",
                                         description="Units")
    Scene.measureit_hide_units = BoolProperty(name="hide_units",
                                              description="Do not display unit of measurement on viewport",
                                              default=False)
    Scene.measureit_render = BoolProperty(name="Render",
                                          description="Save an image with measures over"
                                                      " render image",
                                          default=False)
    Scene.measureit_render_type = EnumProperty(items=(('1', "Frame", "Render current frame"),
                                                      ('2', "Animation", "")),
                                               name="Render type",
                                               description="Type of render image")
    Scene.measureit_sum = EnumProperty(items=(('99', "-", "Select a group for sum"),
                                              ('0', "A", ""),
                                              ('1', "B", ""),
                                              ('2', "C", ""),
                                              ('3', "D", ""),
                                              ('4', "E", ""),
                                              ('5', "F", ""),
                                              ('6', "G", ""),
                                              ('7', "H", ""),
                                              ('8', "I", ""),
                                              ('9', "J", ""),
                                              ('10', "K", ""),
                                              ('11', "L", ""),
                                              ('12', "M", ""),
                                              ('13', "N", ""),
                                              ('14', "O", ""),
                                              ('15', "P", ""),
                                              ('16', "Q", ""),
                                              ('17', "R", ""),
                                              ('18', "S", ""),
                                              ('19', "T", ""),
                                              ('20', "U", ""),
                                              ('21', "V", ""),
                                              ('22', "W", ""),
                                              ('23', "X", ""),
                                              ('24', "Y", ""),
                                              ('25', "Z", "")),
                                       name="Sum in Group",
                                       description="Add segment length in selected group")

    Scene.measureit_rf = BoolProperty(name="render_frame",
                                      description="Add a frame in render output",
                                      default=False)
    Scene.measureit_rf_color = FloatVectorProperty(name="Fcolor",
                                                   description="Frame Color",
                                                   default=(0.9, 0.9, 0.9, 1.0),
                                                   min=0.1,
                                                   max=1,
                                                   subtype='COLOR',
                                                   size=4)
    Scene.measureit_rf_border = IntProperty(name='fborder ', min=1, max=1000, default=10,
                                            description='Frame space from border')
    Scene.measureit_rf_line = IntProperty(name='fline', min=1, max=10, default=1,
                                          description='Line width for border')

    Scene.measureit_glarrow_a = EnumProperty(items=(('99', "--", "No arrow"),
                                                    ('1', "Line",
                                                     "The point of the arrow are lines"),
                                                    ('2', "Triangle",
                                                     "The point of the arrow is triangle"),
                                                    ('3', "TShape",
                                                     "The point of the arrow is a T")),
                                             name="A end",
                                             description="Add arrows to point A")
    Scene.measureit_glarrow_b = EnumProperty(items=(('99', "--", "No arrow"),
                                                    ('1', "Line",
                                                     "The point of the arrow are lines"),
                                                    ('2', "Triangle",
                                                     "The point of the arrow is triangle"),
                                                    ('3', "TShape",
                                                     "The point of the arrow is a T")),
                                             name="B end",
                                             description="Add arrows to point B")
    Scene.measureit_glarrow_s = IntProperty(name="Size",
                                            description="Arrow size",
                                            default=15, min=6, max=500)

    Scene.measureit_debug = BoolProperty(name="Debug",
                                         description="Display information for debugging"
                                                     " (expand/collapse for enabling or disabling)"
                                                     " this information is only rendered for "
                                                     "selected objects",
                                         default=False)
    Scene.measureit_debug_select = BoolProperty(name="Selected",
                                                description="Display information "
                                                            "only for selected items",
                                                default=False)
    Scene.measureit_debug_vertices = BoolProperty(name="Vertices",
                                                  description="Display vertex index number",
                                                  default=True)
    Scene.measureit_debug_objects = BoolProperty(name="Objects",
                                                 description="Display object scene index number",
                                                 default=False)
    Scene.measureit_debug_vert_loc = BoolProperty(name="Location",
                                                  description="Display vertex location",
                                                  default=False)
    Scene.measureit_debug_object_loc = BoolProperty(name="Location",
                                                    description="Display object location",
                                                    default=False)
    Scene.measureit_debug_edges = BoolProperty(name="Edges",
                                               description="Display edge index number",
                                               default=False)
    Scene.measureit_debug_faces = BoolProperty(name="Faces",
                                               description="Display face index number",
                                               default=False)
    Scene.measureit_debug_normals = BoolProperty(name="Normals",
                                                 description="Display face normal "
                                                             "vector and creation order",
                                                 default=False)
    Scene.measureit_debug_normal_details = BoolProperty(name="Details",
                                                        description="Display face normal details",
                                                        default=True)
    Scene.measureit_debug_font = IntProperty(name="Font",
                                             description="Debug text size",
                                             default=14, min=10, max=150)
    Scene.measureit_debug_vert_color = FloatVectorProperty(name="Debug color",
                                                           description="Debug Color",
                                                           default=(1, 0, 0, 1.0),
                                                           min=0.1,
                                                           max=1,
                                                           subtype='COLOR',
                                                           size=4)
    Scene.measureit_debug_face_color = FloatVectorProperty(name="Debug face color",
                                                           description="Debug face Color",
                                                           default=(0, 1, 0, 1.0),
                                                           min=0.1,
                                                           max=1,
                                                           subtype='COLOR',
                                                           size=4)
    Scene.measureit_debug_norm_color = FloatVectorProperty(name="Debug vector color",
                                                           description="Debug vector Color",
                                                           default=(1.0, 1.0, 0.1, 1.0),
                                                           min=0.1,
                                                           max=1,
                                                           subtype='COLOR',
                                                           size=4)
    Scene.measureit_debug_edge_color = FloatVectorProperty(name="Debug vector color",
                                                           description="Debug vector Color",
                                                           default=(0.1, 1.0, 1.0, 1.0),
                                                           min=0.1,
                                                           max=1,
                                                           subtype='COLOR',
                                                           size=4)
    Scene.measureit_debug_obj_color = FloatVectorProperty(name="Debug vector color",
                                                          description="Debug vector Color",
                                                          default=(1.0, 1.0, 1.0, 1.0),
                                                          min=0.1,
                                                          max=1,
                                                          subtype='COLOR',
                                                          size=4)
    Scene.measureit_debug_normal_size = FloatProperty(name='Len', min=0.001, max=9,
                                                      default=0.5,
                                                      precision=2,
                                                      description="Normal arrow size")
    Scene.measureit_debug_width = IntProperty(name='Debug width', min=1, max=10, default=2,
                                              description='Vector line thickness')
    Scene.measureit_debug_precision = IntProperty(name='Precision', min=0, max=5, default=1,
                                                  description="Number of decimal precision")
    Scene.measureit_debug_vert_loc_toggle = EnumProperty(items=(('1', "Local",
                                                                 "Uses local coordinates"),
                                                                ('2', "Global",
                                                                 "Uses global coordinates")),
                                                         name="Coordinates",
                                                         description="Choose coordinate system")
    Scene.measureit_font_rotation = IntProperty(name='Rotate', min=0, max=360, default=0,
                                                description="Default text rotation in degrees")
    Scene.measureit_font_align = EnumProperty(items=(('L', "Left Align", "Use current render"),
                                                     ('C', "Center Align", ""),
                                                     ('R', "Right Align", "")),
                                              name="Align Font",
                                              description="Set Font Alignment")

    # OpenGL flag
    wm = WindowManager
    # register internal property
    wm.measureit_run_opengl = BoolProperty(default=False)

    #Register Sum Edge Lengths
    bpy.utils.register_class(Sum_Edge_Lengths)
    bpy.types.VIEW3D_MT_edit_mesh_edges.append(measureit_sum_edge_lengths.menu_func)

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)

    # Remove properties
    del Scene.measureit_default_color
    del Scene.measureit_font_size
    del Scene.measureit_hint_space
    del Scene.measureit_gl_ghost
    del Scene.measureit_gl_txt
    del Scene.measureit_gl_precision
    del Scene.measureit_gl_show_d
    del Scene.measureit_gl_show_n
    del Scene.measureit_scale
    del Scene.measureit_scale_factor
    del Scene.measureit_scale_color
    del Scene.measureit_scale_font
    del Scene.measureit_scale_pos_x
    del Scene.measureit_scale_pos_y
    del Scene.measureit_gl_scaletxt
    del Scene.measureit_scale_precision
    del Scene.measureit_ovr
    del Scene.measureit_ovr_font
    del Scene.measureit_ovr_color
    del Scene.measureit_ovr_width
    del Scene.measureit_ovr_font_rotation
    del Scene.measureit_ovr_font_align
    del Scene.measureit_units
    del Scene.measureit_hide_units
    del Scene.measureit_render
    del Scene.measureit_render_type
    del Scene.measureit_sum
    del Scene.measureit_rf
    del Scene.measureit_rf_color
    del Scene.measureit_rf_border
    del Scene.measureit_rf_line
    del Scene.measureit_glarrow_a
    del Scene.measureit_glarrow_b
    del Scene.measureit_glarrow_s
    del Scene.measureit_debug
    del Scene.measureit_debug_select
    del Scene.measureit_debug_vertices
    del Scene.measureit_debug_objects
    del Scene.measureit_debug_edges
    del Scene.measureit_debug_faces
    del Scene.measureit_debug_normals
    del Scene.measureit_debug_normal_details
    del Scene.measureit_debug_font
    del Scene.measureit_debug_vert_color
    del Scene.measureit_debug_face_color
    del Scene.measureit_debug_norm_color
    del Scene.measureit_debug_edge_color
    del Scene.measureit_debug_obj_color
    del Scene.measureit_debug_normal_size
    del Scene.measureit_debug_width
    del Scene.measureit_debug_precision
    del Scene.measureit_debug_vert_loc
    del Scene.measureit_debug_object_loc
    del Scene.measureit_debug_vert_loc_toggle
    del Scene.measureit_font_rotation
    del Scene.measureit_font_align

    # remove OpenGL data
    measureit_main.MEASUREIT_OT_RunHintDisplay.handle_remove(measureit_main.MEASUREIT_OT_RunHintDisplay, bpy.context)
    wm = bpy.context.window_manager
    p = 'measureit_run_opengl'
    if p in wm:
        del wm[p]

    #Remove Sum Edge Lengths
    bpy.utils.unregister_class(Sum_Edge_Lengths)
    bpy.types.VIEW3D_MT_edit_mesh_edges.remove(measureit_sum_edge_lengths.menu_func)

if __name__ == '__main__':
    register()
