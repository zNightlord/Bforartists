# SPDX-FileCopyrightText: 2009-2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

import bpy
from bpy.types import Panel
from rna_prop_ui import PropertyPanel
from bl_ui.space_properties import PropertiesAnimationMixin

from bpy.types import Curve, SurfaceCurve, TextCurve


class CurveButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        return (context.curve is not None)


class CurveButtonsPanelCurve(CurveButtonsPanel):
    @classmethod
    def poll(cls, context):
        return (type(context.curve) is Curve)


class CurveButtonsPanelSurface(CurveButtonsPanel):
    @classmethod
    def poll(cls, context):
        return (type(context.curve) is SurfaceCurve)


class CurveButtonsPanelText(CurveButtonsPanel):
    @classmethod
    def poll(cls, context):
        return (type(context.curve) is TextCurve)


class CurveButtonsPanelActive(CurveButtonsPanel):
    """Same as above but for curves only"""

    @classmethod
    def poll(cls, context):
        curve = context.curve
        return (curve and type(curve) is not TextCurve and curve.splines.active)


class DATA_PT_context_curve(CurveButtonsPanel, Panel):
    bl_label = ""
    bl_options = {'HIDE_HEADER'}

    def draw(self, context):
        layout = self.layout

        obj = context.object
        curve = context.curve
        space = context.space_data

        if obj:
            layout.template_ID(obj, "data")
        elif curve:
            layout.template_ID(space, "pin_id")


class DATA_PT_shape_curve(CurveButtonsPanel, Panel):
    bl_label = "Shape"

    def draw(self, context):
        layout = self.layout

        curve = context.curve
        is_surf = type(curve) is SurfaceCurve
        is_curve = type(curve) is Curve
        is_text = type(curve) is TextCurve

        if is_curve:
            row = layout.row()
            row.prop(curve, "dimensions", expand=True)

        layout.use_property_split = True

        col = layout.column()
        sub = col.column(align=True)
        sub.prop(curve, "resolution_u", text="Resolution Preview U")
        if is_surf:
            sub.prop(curve, "resolution_v", text="V")

        sub = col.column(align=True)
        sub.prop(curve, "render_resolution_u", text="Render U")
        if is_surf:
            sub.prop(curve, "render_resolution_v", text="V")
        col.separator()

        if is_curve:
            col.prop(curve, "twist_mode")
            col.prop(curve, "twist_smooth", text="Smooth")
        elif is_text:
            row = layout.row()
            row.use_property_split = False
            row.prop(curve, "use_fast_edit", text="Fast Editing")
            row.prop_decorator(curve, "use_fast_edit")

        if is_curve or is_text:
            col = layout.column()
            col.separator()

            sub = col.column()
            sub.active = (curve.dimensions == '2D' or (curve.bevel_mode != 'OBJECT' and curve.dimensions == '3D'))
            sub.prop(curve, "fill_mode")

        if is_curve:
            col.label(text = "Curve Deform") # BFA

            col.use_property_split = False
            row = col.row()
            row.separator()
            row.prop(curve, "use_radius")
            row.prop_decorator(curve, "use_radius")
            row = col.row()
            row.separator()
            row.prop(curve, "use_stretch")
            row.prop_decorator(curve, "use_stretch")
            row = col.row()
            row.separator()
            row.prop(curve, "use_deform_bounds")
            row.prop_decorator(curve, "use_deform_bounds")


class DATA_PT_curve_texture_space(CurveButtonsPanel, Panel):
    bl_label = "Texture Space"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        curve = context.curve

        row = layout.row()
        row.use_property_split = False
        row.prop(curve, "use_auto_texspace")
        row.prop_decorator(curve, "use_auto_texspace")

        layout.use_property_split = True

        col = layout.column()
        col.prop(curve, "texspace_location")
        col.prop(curve, "texspace_size")

        layout.operator("curve.match_texture_space")


class DATA_PT_geometry_curve(CurveButtonsPanelCurve, Panel):
    bl_label = "Geometry"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        return (type(context.curve) in {Curve, TextCurve})

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        curve = context.curve

        col = layout.column()
        col.prop(curve, "offset")

        sub = col.column()
        sub.active = (curve.bevel_mode != 'OBJECT')
        sub.prop(curve, "extrude")

        col.prop(curve, "taper_object")
        if curve.taper_object is not None: # BFA
            row = col.row()
            row.separator()
            row.prop(curve, "taper_radius_mode")

        if type(curve) is not TextCurve:
            # This setting makes no sense for texts, since we have no control over start/end of the bevel object curve.
            row = col.row() # BFA
            if curve.taper_object is not None:
                row.use_property_split = False
                row.separator()
                row.prop(curve, "use_map_taper")
                row.prop_decorator(curve, "use_map_taper")


class DATA_PT_geometry_curve_bevel(CurveButtonsPanelCurve, Panel):
    bl_label = "Bevel"
    bl_parent_id = "DATA_PT_geometry_curve"

    @classmethod
    def poll(cls, context):
        return (type(context.curve) in {Curve, TextCurve})

    def draw(self, context):
        layout = self.layout

        curve = context.curve
        layout.prop(curve, "bevel_mode", expand=True)

        layout.use_property_split = True

        col = layout.column()
        if curve.bevel_mode == 'OBJECT':
            col.prop(curve, "bevel_object", text="Object")
        else:
            col.prop(curve, "bevel_depth", text="Depth")
            col.prop(curve, "bevel_resolution", text="Resolution")

        row = col.row() # BFA
        row.use_property_split = False
        row.prop(curve, "use_fill_caps")
        row.prop_decorator(curve, "use_fill_caps")

        if curve.bevel_mode == 'PROFILE':
            col.template_curveprofile(curve, "bevel_profile")


class DATA_PT_curve_animation(CurveButtonsPanel, PropertiesAnimationMixin, PropertyPanel, Panel):
    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        # MeshButtonsPanel.poll ensures this is not None.
        curve = context.curve

        col = layout.column(align=True)
        col.label(text=curve.bl_rna.name)  # "Surface Curve" or "Curve".
        self.draw_action_and_slot_selector(context, col, curve)

        if shape_keys := curve.shape_keys:
            col = layout.column(align=True)
            col.label(text="Shape Keys")
            self.draw_action_and_slot_selector(context, col, shape_keys)


class DATA_PT_geometry_curve_start_end(CurveButtonsPanelCurve, Panel):
    bl_label = "Start & End Mapping"
    bl_parent_id = "DATA_PT_geometry_curve"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        # Text objects don't support these properties
        return (type(context.curve) is Curve)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        curve = context.curve

        col = layout.column()

        col.active = (
            ((curve.bevel_depth > 0.0) or (curve.extrude > 0.0)) or
            ((curve.bevel_mode == 'OBJECT') and curve.bevel_object is not None)
        )
        sub = col.column(align=True)
        sub.prop(curve, "bevel_factor_start", text="Factor Start")
        sub.prop(curve, "bevel_factor_end", text="End")

        sub = col.column(align=True)
        sub.prop(curve, "bevel_factor_mapping_start", text="Mapping Start")
        sub.prop(curve, "bevel_factor_mapping_end", text="End")


class DATA_PT_pathanim(CurveButtonsPanelCurve, Panel):
    bl_label = "Path Animation"
    bl_options = {'DEFAULT_CLOSED'}

    def draw_header(self, context):
        curve = context.curve

        self.layout.prop(curve, "use_path", text="")

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        curve = context.curve

        layout.active = curve.use_path

        col = layout.column()
        col.prop(curve, "path_duration", text="Frames")
        col.prop(curve, "eval_time")

        # these are for paths only
        row = layout.row() # BFA
        row.use_property_split = False
        row.prop(curve, "use_path_clamp")
        row.prop_decorator(curve, "use_path_clamp")

        row = layout.row() # BFA
        row.use_property_split = False
        row.prop(curve, "use_path_follow")
        row.prop_decorator(curve, "use_path_follow")


class DATA_PT_font(CurveButtonsPanelText, Panel):
    bl_label = "Font"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout

        text = context.curve
        char = text.edit_format
        mode = context.mode

        row = layout.split(factor=0.25)
        row.label(text="Regular")
        row.template_ID(text, "font", open="font.open", unlink="font.unlink")
        row = layout.split(factor=0.25)
        row.label(text="Bold")
        row.template_ID(text, "font_bold", open="font.open", unlink="font.unlink")
        row = layout.split(factor=0.25)
        row.label(text="Italic")
        row.template_ID(text, "font_italic", open="font.open", unlink="font.unlink")
        row = layout.split(factor=0.25)
        row.label(text="Bold & Italic")
        row.template_ID(text, "font_bold_italic", open="font.open", unlink="font.unlink")

        if mode == 'EDIT_TEXT':
            layout.separator()

            if not text.has_selection:
                row = layout.row(align=True)
                row.prop(char, "use_bold", toggle=True)
                row.prop(char, "use_italic", toggle=True)
                row.prop(char, "use_underline", toggle=True)
                row.prop(char, "use_small_caps", toggle=True)
            else:
                row = layout.row(align=True)
                row.operator(
                    "font.style_toggle", text="Bold", icon='BOLD', depress=text.is_select_bold,
                ).style = 'BOLD'
                row.operator(
                    "font.style_toggle", text="Italic", icon='ITALIC', depress=text.is_select_italic,
                ).style = 'ITALIC'
                row.operator(
                    "font.style_toggle", text="Underline", icon='UNDERLINE', depress=text.is_select_underline,
                ).style = 'UNDERLINE'
                row.operator(
                    "font.style_toggle", text="Small Caps", icon='SMALL_CAPS', depress=text.is_select_smallcaps,
                ).style = 'SMALL_CAPS'


class DATA_PT_font_transform(CurveButtonsPanelText, Panel):
    bl_label = "Transform"
    bl_parent_id = "DATA_PT_font"

    def draw(self, context):
        layout = self.layout

        text = context.curve

        layout.use_property_split = True

        col = layout.column()

        col.separator()

        col.prop(text, "size", text="Size")
        col.prop(text, "shear")

        col.separator()

        col.prop(text, "family")
        col.prop(text, "follow_curve")

        col.separator()

        sub = col.column(align=True)
        sub.prop(text, "underline_position", text="Underline Position")
        sub.prop(text, "underline_height", text="Underline Thickness")

        col.prop(text, "small_caps_scale", text="Small Caps Scale")


class DATA_PT_paragraph(CurveButtonsPanelText, Panel):
    bl_label = "Paragraph"

    def draw(self, context):
        # Parent panel
        pass


class DATA_PT_paragraph_alignment(CurveButtonsPanelText, Panel):
    bl_parent_id = "DATA_PT_paragraph"
    bl_label = "Alignment"

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        text = context.curve

        col = layout.column()
        col.prop(text, "align_x", text="Horizontal")
        col.prop(text, "align_y", text="Vertical")


class DATA_PT_paragraph_spacing(CurveButtonsPanelText, Panel):
    bl_parent_id = "DATA_PT_paragraph"
    bl_label = "Spacing"

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True

        text = context.curve

        col = layout.column(align=True)
        col.prop(text, "space_character", text="Character Spacing")
        col.prop(text, "space_word", text="Word Spacing")
        col.prop(text, "space_line", text="Line Spacing")

        layout.separator()

        col = layout.column(align=True)
        col.prop(text, "offset_x", text="Offset X")
        col.prop(text, "offset_y", text="Y")


class DATA_PT_text_boxes(CurveButtonsPanelText, Panel):
    bl_label = "Text Boxes"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        layout = self.layout

        text = context.curve

        layout.operator("font.textbox_add", icon='ADD')
        layout.prop(text, "overflow", text="Overflow")

        for i, box in enumerate(text.text_boxes):

            boxy = layout.box()

            row = boxy.row()

            col = row.column()
            col.use_property_split = True

            sub = col.column(align=True)
            sub.prop(box, "width", text="Size X")
            sub.prop(box, "height", text="Y")

            sub = col.column(align=True)
            sub.prop(box, "x", text="Offset X")
            sub.prop(box, "y", text="Y")

            row.operator("font.textbox_remove", text="", icon='X', emboss=False).index = i


class DATA_PT_custom_props_curve(CurveButtonsPanel, PropertyPanel, Panel):
    _context_path = "object.data"
    _property_type = bpy.types.Curve


classes = (
    DATA_PT_context_curve,
    DATA_PT_shape_curve,
    DATA_PT_curve_texture_space,
    DATA_PT_geometry_curve,
    DATA_PT_geometry_curve_bevel,
    DATA_PT_geometry_curve_start_end,
    DATA_PT_pathanim,
    DATA_PT_font,
    DATA_PT_font_transform,
    DATA_PT_paragraph,
    DATA_PT_paragraph_alignment,
    DATA_PT_paragraph_spacing,
    DATA_PT_text_boxes,
    DATA_PT_curve_animation,
    DATA_PT_custom_props_curve,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
