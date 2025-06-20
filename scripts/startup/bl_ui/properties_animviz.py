# SPDX-FileCopyrightText: 2010-2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

# Generic Panels (Independent of DataType)

# NOTE:
# The specialized panel types are derived in their respective UI modules
# don't register these classes since they are only helpers.


class MotionPathButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = "Motion Paths"
    # bl_options = {'DEFAULT_CLOSED'}

    def draw_settings(self, _context, avs, mpath, bones=False):
        layout = self.layout

        mps = avs.motion_path

        layout.use_property_split = True
        layout.use_property_decorate = False

        col = layout.column(align=True)
        col.prop(mps, "type")
        range_group = col.column(align=True)
        # BFA - float left
        if mps.type == 'RANGE':
            row = range_group.row()
            row.separator()
            row.prop(mps, "range", text="Calculation Range")

        if mps.type == 'CURRENT_FRAME':
            col = layout.column(align=True)
            row = col.row()
            row.separator()
            row.prop(mps, "frame_before", text="Frame Range Before")
            row = col.row()
            row.separator()
            row.prop(mps, "frame_after", text="After")
            row = col.row()
            row.separator()
            row.prop(mps, "frame_step", text="Step")
        elif mps.type == 'RANGE':
            col = layout.column(align=True)
            start_end_group = col.column(align=True)
            if mps.range == 'MANUAL':
                row = start_end_group.row()
                row.separator()
                row.prop(mps, "frame_start", text="Frame Range Start")
                row = start_end_group.row()
                row.separator()
                row.prop(mps, "frame_end", text="End")
            row = col.row()
            row.separator()
            row.prop(mps, "frame_step", text="Step")

        row = col.row()
        row.use_property_split = False
        row.prop(mps, "use_camera_space_bake", text="Bake to Active Camera")

        if bones:
            op_category = "pose"
            icon = 'BONE_DATA'
        else:
            op_category = "object"
            icon = 'OBJECT_DATA'

        if mpath:
            col = layout.column(align=True)
            row = col.row(align=True)
            row.enabled = False
            row.prop(mpath, "frame_start", text="Cached Range (Info)") # BFA
            row.prop(mpath, "frame_end", text="")

            # Update Selected.
            col = layout.column(align=True)
            row = col.row(align=True)
            row.operator(op_category + ".paths_update", text="Update Path", icon=icon)
            row.operator(op_category + ".paths_clear", text="", icon='X').only_selected = True
        else:
            # Calculate.
            col = layout.column(align=True)
            col.label(text="Nothing to show yet...", icon='ERROR')
            col.operator(op_category + ".paths_calculate", text="Calculate...", icon=icon)

        # Update All & Clear All.
        # Note that `col` is from inside the preceding `if` or `else` block.
        row = col.row(align=True)
        row.operator("object.paths_update_visible", text="Update All Paths", icon='WORLD')
        row.operator(op_category + ".paths_clear", text="", icon='X').only_selected = False


class MotionPathButtonsPanel_display:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = "Display"

    def draw_settings(self, _context, avs, mpath, bones=False):
        layout = self.layout

        mps = avs.motion_path

        layout.use_property_split = False # BFA
        layout.use_property_decorate = False

        flow = layout.grid_flow(row_major=False, columns=0, even_columns=False, even_rows=False, align=True)

        flow.prop(mps, "show_frame_numbers", text="Frame Numbers")
        flow.prop(mps, "show_keyframe_highlight", text="Keyframes")
        sub = flow.column()
        sub.enabled = mps.show_keyframe_highlight
        if bones:
            sub.prop(mps, "show_keyframe_action_all", text="+ Non-Grouped Keyframes")
        sub.prop(mps, "show_keyframe_numbers", text="Keyframe Numbers")

        # Customize path
        if mpath is not None:
            flow.prop(mpath, "lines", text="Lines")

            col = layout.column()
            col.prop(mpath, "line_thickness", text="Thickness")

            # BFA - Float Left
            row = layout.row()
            row.use_property_split = False
            split = row.split(factor = 0.5)
            row = split.row()
            row.prop(mpath, "use_custom_color", text="Custom Color")
            row = split.row()
            if mpath.use_custom_color:
                row.label(icon='DISCLOSURE_TRI_DOWN')
            else:
                row.label(icon='DISCLOSURE_TRI_RIGHT')

            if mpath.use_custom_color:
                row = layout.row()
                row.separator()
                row.prop(mpath, "color", text="Before")
                row = layout.row()
                row.separator()
                row.prop(mpath, "color_post", text="After")


classes = (
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
