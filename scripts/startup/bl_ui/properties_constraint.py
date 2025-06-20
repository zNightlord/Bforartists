# SPDX-FileCopyrightText: 2009-2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

from bpy.types import Panel, Menu, Operator # BFA
from bpy.app.translations import contexts as i18n_contexts
from bl_ui.generic_column_menu import GenericColumnMenu, fetch_op_data, InvokeMenuOperator # BFA

class ObjectConstraintPanel:
    bl_context = "constraint"

    @classmethod
    def poll(cls, context):
        return (context.object)


class BoneConstraintPanel:
    bl_context = "bone_constraint"

    @classmethod
    def poll(cls, context):
        return (context.pose_bone)


class OBJECT_PT_constraints(ObjectConstraintPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = "Object Constraints"
    bl_options = {'HIDE_HEADER'}

    def draw(self, _context):
        layout = self.layout
        layout.operator("object.add_constraints_menu", icon='ADD') 

        layout.template_constraints(use_bone_constraints=False)

# BFA menu
class OBJECT_MT_constraint_add(GenericColumnMenu, Menu):
    bl_description = "Add a constraint to the active object"

    op_id = "object.constraint_add"
    OPERATOR_DATA, TRANSLATION_CONTEXT = fetch_op_data(class_name="Constraint")
    search_header = "Object Constraint"

    def draw(self, _context):
        layout = self.layout.row()

        self.draw_operator_column(layout, header="Motion Tracking",
            types=('CAMERA_SOLVER', 'FOLLOW_TRACK', 'OBJECT_SOLVER'))
        self.draw_operator_column(layout, header="Transform",
            types=('COPY_LOCATION', 'COPY_ROTATION', 'COPY_SCALE', 'COPY_TRANSFORMS', 'LIMIT_DISTANCE', 'LIMIT_LOCATION', 'LIMIT_ROTATION', 'LIMIT_SCALE', 'MAINTAIN_VOLUME', 'TRANSFORM', 'TRANSFORM_CACHE'))
        self.draw_operator_column(layout, header="Tracking",
            types=('CLAMP_TO', 'DAMPED_TRACK', 'LOCKED_TRACK', 'STRETCH_TO', 'TRACK_TO'))
        self.draw_operator_column(layout, header="Relationship",
            types=('ACTION', 'ARMATURE', 'CHILD_OF', 'FLOOR', 'FOLLOW_PATH', 'PIVOT', 'SHRINKWRAP'))

# BFA menu
class OBJECT_OT_add_constraints_menu(InvokeMenuOperator, Operator):
    bl_idname = "object.add_constraints_menu"
    bl_label = "Add Object Constraint"
    bl_description = "Add a constraint to the active object"

    menu_id = "OBJECT_MT_constraint_add"
    space_type = 'PROPERTIES'
    space_context = 'CONSTRAINT'


class BONE_PT_constraints(BoneConstraintPanel, Panel):
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = "Bone Constraints"
    bl_options = {'HIDE_HEADER'}

    def draw(self, _context):
        layout = self.layout
        layout.operator("bone.add_constraints_menu", icon='ADD')

        layout.template_constraints(use_bone_constraints=True)

# BFA - menu
class BONE_MT_constraint_add(GenericColumnMenu, Menu):
    bl_description = "Add a constraint to the active bone"

    op_id = "pose.constraint_add"
    OPERATOR_DATA, TRANSLATION_CONTEXT = fetch_op_data(class_name="Constraint")
    search_header = "Bone Constraint"

    def draw(self, _context):
        layout = self.layout.row()

        self.draw_operator_column(layout, header="Motion Tracking",
            types=('CAMERA_SOLVER', 'FOLLOW_TRACK', 'OBJECT_SOLVER'))
        self.draw_operator_column(layout, header="Transform",
            types=('COPY_LOCATION', 'COPY_ROTATION', 'COPY_SCALE', 'COPY_TRANSFORMS', 'LIMIT_DISTANCE', 'LIMIT_LOCATION', 'LIMIT_ROTATION', 'LIMIT_SCALE', 'MAINTAIN_VOLUME', 'TRANSFORM', 'TRANSFORM_CACHE'))
        self.draw_operator_column(layout, header="Tracking",
            types=('CLAMP_TO', 'DAMPED_TRACK', 'IK', 'LOCKED_TRACK', 'SPLINE_IK', 'STRETCH_TO', 'TRACK_TO'))
        self.draw_operator_column(layout, header="Relationship",
            types=('ACTION', 'ARMATURE', 'CHILD_OF', 'FLOOR', 'FOLLOW_PATH', 'PIVOT', 'SHRINKWRAP'))

# BFA - menu
class BONE_OT_add_constraints_menu(InvokeMenuOperator, Operator):
    bl_idname = "bone.add_constraints_menu"
    bl_label = "Add Bone Constraint"
    bl_description = "Add a constraint to the active bone"

    menu_id = "BONE_MT_constraint_add"
    space_type = 'PROPERTIES'
    space_context = 'BONE_CONSTRAINT'


# Parent class for constraint panels, with templates and drawing methods
# shared between the bone and object constraint panels
class ConstraintButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = ""
    bl_options = {'INSTANCED', 'HEADER_LAYOUT_EXPAND'}

    @staticmethod
    def draw_influence(layout, con):
        layout.separator()
        if con.type in {'IK', 'SPLINE_IK'}:
            # constraint.disable_keep_transform doesn't work well
            # for these constraints.
            layout.prop(con, "influence")
        else:
            row = layout.row(align=True)
            row.prop(con, "influence")
            row.operator("constraint.disable_keep_transform", text="", icon='CANCEL')

    @staticmethod
    def space_template(layout, con, target=True, owner=True, separator=True):
        if target or owner:
            if separator:
                layout.separator()
            if target:
                layout.prop(con, "target_space", text="Target")
            if owner:
                layout.prop(con, "owner_space", text="Owner")

            if con.target_space == 'CUSTOM' or con.owner_space == 'CUSTOM':
                col = layout.column()
                col.prop(con, "space_object")
                if space_object := con.space_object:
                    match space_object.type:
                        case 'ARMATURE':
                            col.prop_search(con, "space_subtarget", con.space_object.data, "bones", text="Bone")
                        case 'MESH', 'LATTICE':
                            col.prop_search(
                                con, "space_subtarget", con.space_object,
                                "vertex_groups", text="Vertex Group",
                            )

    @staticmethod
    def target_template(layout, con, subtargets=True):
        col = layout.column()
        col.prop(con, "target")  # XXX: limiting settings for only `curves` or some type of object.

        if con.target and subtargets:
            if con.target.type == 'ARMATURE':
                col.prop_search(con, "subtarget", con.target.data, "bones", text="Bone")

                if con.subtarget and hasattr(con, "head_tail"):
                    row = col.row(align=True)
                    row.use_property_decorate = False
                    sub = row.row(align=True)
                    sub.prop(con, "head_tail")
                    # XXX icon, and only when bone has segments?
                    sub.prop(con, "use_bbone_shape", text="", icon='IPO_BEZIER')
                    row.prop_decorator(con, "head_tail")
            elif con.target.type in {'MESH', 'LATTICE'}:
                col.prop_search(con, "subtarget", con.target, "vertex_groups", text="Vertex Group")

    def get_constraint(self, _context):
        con = self.custom_data
        self.layout.context_pointer_set("constraint", con)
        return con

    def draw_header(self, context):
        layout = self.layout
        con = self.get_constraint(context)

        layout.template_constraint_header(con)

    # Drawing methods for specific constraints. (Shared by object and bone constraint panels)

    def draw_childof(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row(heading="Location")
        row.use_property_decorate = False
        row.prop(con, "use_location_x", text="X", toggle=True)
        row.prop(con, "use_location_y", text="Y", toggle=True)
        row.prop(con, "use_location_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row(heading="Rotation")
        row.use_property_decorate = False
        row.prop(con, "use_rotation_x", text="X", toggle=True)
        row.prop(con, "use_rotation_y", text="Y", toggle=True)
        row.prop(con, "use_rotation_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row(heading="Scale")
        row.use_property_decorate = False
        row.prop(con, "use_scale_x", text="X", toggle=True)
        row.prop(con, "use_scale_y", text="Y", toggle=True)
        row.prop(con, "use_scale_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row()
        row.operator("constraint.childof_set_inverse")
        row.operator("constraint.childof_clear_inverse")

        self.draw_influence(layout, con)

    def draw_trackto(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "track_axis", expand=True)
        layout.prop(con, "up_axis", text="Up", expand=True)

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_target_z")
        row.prop_decorator(con, "use_target_z")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_follow_path(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        if con.use_fixed_location:
            layout.prop(con, "offset_factor", text="Offset Factor")
        else:
            layout.prop(con, "offset")

        layout.prop(con, "forward_axis", expand=True)
        layout.prop(con, "up_axis", expand=True)

        col = layout.column(align = True)
        col.use_property_split = False
        row = col.row()
        row.prop(con, "use_fixed_location")
        row.prop_decorator(con, "use_fixed_location")
        row = col.row()
        row.prop(con, "use_curve_radius")
        row.prop_decorator(con, "use_curve_radius")
        row = col.row()
        row.prop(con, "use_curve_follow")
        row.prop_decorator(con, "use_curve_follow")

        layout.operator("constraint.followpath_path_animate", text="Animate Path", icon='ANIM_DATA')

        self.draw_influence(layout, con)

    def draw_rot_limit(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        #########################################

        col = layout.column()
        split = col.split(factor = 0.38)
        split.use_property_split = False

        col = split.column()
        col.prop(con, "use_limit_x", text = "Limit X")

        col = split.column(align = True)
        if con.use_limit_x:
            col.use_property_decorate = False
            row = col.row(align = True)
            sub = row.column(align=True)
            sub.prop(con, "min_x", text="Min")
            sub.prop(con, "max_x", text="Max")
            row.label(icon='BLANK1')
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        #########################################

        col = layout.column()
        split = col.split(factor = 0.38)
        split.use_property_split = False

        col = split.column()
        col.prop(con, "use_limit_y", text = "Y")

        col = split.column(align = True)
        if con.use_limit_y:
            col.use_property_decorate = False
            row = col.row(align = True)
            sub = row.column(align=True)
            sub.prop(con, "min_y", text="Min")
            sub.prop(con, "max_y", text="Max")
            row.label(icon='BLANK1')
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        #########################################

        col = layout.column()
        split = col.split(factor = 0.38)
        split.use_property_split = False

        col = split.column(align = True)
        col.prop(con, "use_limit_z", text = "Z")

        col = split.column()
        if con.use_limit_z:
            col.use_property_decorate = False
            row = col.row(align = True)
            sub = row.column(align=True)
            sub.prop(con, "min_z", text="Min")
            sub.prop(con, "max_z", text="Max")
            row.label(icon='BLANK1')
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        ###########################################

        layout.prop(con, "euler_order", text="Order")

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_transform_limit")
        row.prop_decorator(con, "use_transform_limit")
        row.prop_decorator(con, "use_legacy_behavior")

        self.space_template(layout, con, target=False, owner=True)

        self.draw_influence(layout, con)

    def draw_loc_limit(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        #########################################

        col = layout.column(align = True)
        split = col.split(factor = 0.38)

        col = split.column(align = True)
        col.use_property_split = False
        col.prop(con, "use_min_x", text = "Minimum X")
        col.prop(con, "use_min_y", text = "Y")
        col.prop(con, "use_min_z", text = "Z")

        col = split.column(align = True)
        if con.use_min_x:
            row = col.row(align = True)
            row.prop(con, "min_x", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')
        if con.use_min_y:
            row = col.row(align = True)
            row.prop(con, "min_y", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')
        if con.use_min_z:
            row = col.row(align = True)
            row.prop(con, "min_z", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        ###########################################

        col = layout.column(align = True)
        split = col.split(factor = 0.38)

        col = split.column(align = True)
        col.use_property_split = False
        col.prop(con, "use_max_x", text = "Maximum X")
        col.prop(con, "use_max_y", text = "Y")
        col.prop(con, "use_max_z", text = "Z")

        col = split.column(align = True)
        if con.use_max_x:
            row = col.row(align = True)
            row.prop(con, "max_x", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')
        if con.use_max_y:
            row = col.row(align = True)
            row.prop(con, "max_y", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')
        if con.use_max_z:
            row = col.row(align = True)
            row.prop(con, "max_z", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        ###########################################

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_transform_limit")
        row.prop_decorator(con, "use_transform_limit")

        self.space_template(layout, con, target=False, owner=True)

        self.draw_influence(layout, con)

    def draw_size_limit(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        col = layout.column()

        row = col.row(heading="Minimum X", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_min_x", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_min_x
        subsub.prop(con, "min_x", text="")
        row.prop_decorator(con, "min_x")

        row = col.row(heading="Y", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_min_y", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_min_y
        subsub.prop(con, "min_y", text="")
        row.prop_decorator(con, "min_y")

        row = col.row(heading="Z", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_min_z", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_min_z
        subsub.prop(con, "min_z", text="")
        row.prop_decorator(con, "min_z")

        col.separator()

        row = col.row(heading="Maximum X", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_max_x", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_max_x
        subsub.prop(con, "max_x", text="")
        row.prop_decorator(con, "max_x")

        row = col.row(heading="Y", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_max_y", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_max_y
        subsub.prop(con, "max_y", text="")
        row.prop_decorator(con, "max_y")

        row = col.row(heading="Z", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_max_z", text="")
        subsub = sub.row(align=True)
        subsub.active = con.use_max_z
        subsub.prop(con, "max_z", text="")
        row.prop_decorator(con, "max_z")

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_transform_limit")
        row.prop_decorator(con, "use_transform_limit")

        self.space_template(layout, con, target=False, owner=True)

        self.draw_influence(layout, con)

    def draw_rotate_like(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "euler_order", text="Order")

        row = layout.row(heading="Axis", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_x", text="X", toggle=True)
        sub.prop(con, "use_y", text="Y", toggle=True)
        sub.prop(con, "use_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row(heading="Invert", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "invert_x", text="X", toggle=True)
        sub.prop(con, "invert_y", text="Y", toggle=True)
        sub.prop(con, "invert_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        # bfa - goo engine - copy rotation constraint invert rotation patch
        row = layout.row()
        row.use_property_split = False
        row.prop(con, "invert_all", text="Invert Rotation")
        row.prop_decorator(con, "invert_all")
        # end bfa - goo engine

        layout.prop(con, "mix_mode", text="Mix", text_ctxt=i18n_contexts.constraint)

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_locate_like(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row(heading="Axis", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_x", text="X", toggle=True)
        sub.prop(con, "use_y", text="Y", toggle=True)
        sub.prop(con, "use_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row(heading="Invert", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "invert_x", text="X", toggle=True)
        sub.prop(con, "invert_y", text="Y", toggle=True)
        sub.prop(con, "invert_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_offset")
        row.prop_decorator(con, "use_offset")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_size_like(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row(heading="Axis", align=True)
        row.use_property_decorate = False
        sub = row.row(align=True)
        sub.prop(con, "use_x", text="X", toggle=True)
        sub.prop(con, "use_y", text="Y", toggle=True)
        sub.prop(con, "use_z", text="Z", toggle=True)
        row.label(icon='BLANK1')

        col = layout.column()
        col.prop(con, "power")

        col = layout.column(align = True)
        col.use_property_split = False
        row = col.row()
        row.prop(con, "use_make_uniform")
        row.prop_decorator(con, "use_make_uniform")

        row = col.row()
        row.prop(con, "use_offset")
        row.prop_decorator(con, "use_offset")

        row = col.row()
        row.active = con.use_offset
        row.prop(con, "use_add")
        row.prop_decorator(con, "use_add")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_same_volume(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        layout.prop(con, "mode")

        row = layout.row(heading="Free Axis")
        row.prop(con, "free_axis", expand=True)

        layout.prop(con, "volume")

        self.space_template(layout, con, target=False, owner=True)

        self.draw_influence(layout, con)

    def draw_trans_like(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        # BFA        
        row = layout.row()
        row.use_property_split = False
        row.prop(con, "remove_target_shear")
        row.prop_decorator(con, "remove_target_shear")

        layout.prop(con, "mix_mode", text="Mix", text_ctxt=i18n_contexts.constraint)

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_action(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        target_row = layout.row(align=True)
        target_row.active = not con.use_eval_time
        self.target_template(target_row, con)

        # BFA ###########################################

        split = layout.split(factor = 0.38)
        col = split.column(align = True)
        col.use_property_split = False
        col.prop(con, "use_eval_time", text = "Evaluation Time")
        col = split.column()
        if con.use_eval_time:
            row = col.row()
            row.prop(con, "eval_time", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        # BFA ##########################################

        layout.prop(con, "mix_mode", text="Mix", text_ctxt=i18n_contexts.constraint)

        self.draw_influence(layout, con)

    def draw_lock_track(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "track_axis", expand=True)
        layout.prop(con, "lock_axis", expand=True)

        self.draw_influence(layout, con)

    def draw_dist_limit(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row()
        row.prop(con, "distance")
        row.operator("constraint.limitdistance_reset", text="", icon='X')

        layout.prop(con, "limit_mode", text="Clamp Region")

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_transform_limit")
        row.prop_decorator(con, "use_transform_limit")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_stretch_to(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row()
        row.prop(con, "rest_length")
        row.operator("constraint.stretchto_reset", text="", icon='X')

        layout.separator()

        col = layout.column()
        col.prop(con, "bulge", text="Volume Variation")

        ##########################################

        split = layout.split(factor = 0.38)
        col = split.column(align = True)
        col.use_property_split = False
        col.prop(con, "use_bulge_min", text = "Volume Min")
        col = split.column()
        if con.use_bulge_min:
            row = col.row()
            row.prop(con, "bulge_min", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        split = layout.split(factor = 0.38)
        col = split.column()
        col.use_property_split = False
        col.prop(con, "use_bulge_max", text = "Volume Max")
        col = split.column()
        if con.use_bulge_max:
            row = col.row()
            row.prop(con, "bulge_max", text="")
        else:
            col.label(icon='DISCLOSURE_TRI_RIGHT')

        if con.use_bulge_min or con.use_bulge_max:
            row = layout.row()
            row.separator()
            row.prop(con, "bulge_smooth", text="Smooth")

        ##########################################

        layout.separator()

        layout.prop(con, "volume", expand=True)
        layout.prop(con, "keep_axis", text="Rotation", expand=True)

        self.draw_influence(layout, con)

    def draw_min_max(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "offset")
        layout.prop(con, "floor_location", expand=True, text="Min/Max")

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_rotation")
        row.prop_decorator(con, "use_rotation")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_clamp_to(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "main_axis", expand=True)

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_cyclic")
        row.prop_decorator(con, "use_cyclic")

        self.draw_influence(layout, con)

    def draw_transform(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_motion_extrapolate", text="Extrapolate")
        row.prop_decorator(con, "use_motion_extrapolate")

        self.space_template(layout, con)

        self.draw_influence(layout, con)

    def draw_shrinkwrap(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con, False)

        layout.prop(con, "distance")
        layout.prop(con, "shrinkwrap_type", text="Mode")

        layout.separator()

        if con.shrinkwrap_type == 'PROJECT':
            layout.prop(con, "project_axis", expand=True, text="Project Axis")
            layout.prop(con, "project_axis_space", text="Space")

            if con.project_axis_space == 'CUSTOM':
                col = layout.column()
                col.prop(con, "space_object")
                if space_object := con.space_object:
                    match space_object.type:
                        case 'ARMATURE':
                            col.prop_search(
                                con, "space_subtarget",
                                con.space_object.data, "bones",
                                text="Bone",
                            )
                        case 'MESH', 'LATTICE':
                            col.prop_search(
                                con, "space_subtarget",
                                con.space_object, "vertex_groups",
                                text="Vertex Group",
                            )

            layout.prop(con, "project_limit", text="Distance")

            row = layout.row()
            row.use_property_split = False
            row.prop(con, "use_project_opposite")
            row.prop_decorator(con, "use_project_opposite")

            layout.separator()

            col = layout.column()
            row = col.row()
            row.prop(con, "cull_face", expand=True)
            row = col.row()
            row.active = con.use_project_opposite and con.cull_face != 'OFF'


            row = col.row()
            row.use_property_split = False
            row.prop(con, "use_invert_cull")
            row.prop_decorator(con, "use_invert_cull")

            layout.separator()

        if con.shrinkwrap_type in {'PROJECT', 'NEAREST_SURFACE', 'TARGET_PROJECT'}:
            layout.prop(con, "wrap_mode", text="Snap Mode")

            ###########################################

            split = layout.split(factor = 0.38)
            col = split.column(align = True)
            col.use_property_split = False
            col.prop(con, "use_track_normal", text = "Align to Normal")
            col = split.column()
            if con.use_track_normal:
                row = col.row()
                row.prop(con, "track_axis", text="")
            else:
                col.label(icon='DISCLOSURE_TRI_RIGHT')

            ##########################################

        self.draw_influence(layout, con)

    def draw_damp_track(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        layout.prop(con, "track_axis", expand=True)

        self.draw_influence(layout, con)

    def draw_spline_ik(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        self.draw_influence(layout, con)

    def draw_pivot(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        if con.target:
            layout.prop(con, "offset", text="Pivot Offset")
        else:
            row = layout.row()
            row.use_property_split = False
            row.prop(con, "use_relative_location")
            row.prop_decorator(con, "use_relative_location")
            if con.use_relative_location:
                layout.prop(con, "offset", text="Pivot Point")
            else:
                layout.prop(con, "offset", text="Pivot Point")

        col = layout.column()
        col.prop(con, "rotation_range", text="Rotation Range")

        self.draw_influence(layout, con)

    def draw_follow_track(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        clip = None
        if con.use_active_clip:
            clip = context.scene.active_clip
        else:
            clip = con.clip

        col = layout.column(align = True)
        col.use_property_split = False
        row = col.row()
        row.prop(con, "use_active_clip")
        row.prop_decorator(con, "use_active_clip")

        row = col.row()
        row.prop(con, "use_3d_position")
        row.prop_decorator(con, "use_3d_position")

        row = col.row()
        row.active = not con.use_3d_position
        row.prop(con, "use_undistorted_position")
        row.prop_decorator(con, "use_undistorted_position")

        if not con.use_active_clip:
            layout.prop(con, "clip")

        layout.prop(con, "frame_method")

        if clip:
            tracking = clip.tracking

            layout.prop_search(con, "object", tracking, "objects", icon='OBJECT_DATA')

            tracking_object = tracking.objects.get(con.object, tracking.objects[0])

            layout.prop_search(con, "track", tracking_object, "tracks", icon='ANIM_DATA')

        layout.prop(con, "camera")

        row = layout.row()
        row.active = not con.use_3d_position
        row.prop(con, "depth_object")

        layout.operator("clip.constraint_to_fcurve")

        self.draw_influence(layout, con)

    def draw_camera_solver(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = False
        layout.use_property_decorate = True

        row = layout.row()
        row.prop(con, "use_active_clip")
        row.prop_decorator(con, "use_active_clip")

        if not con.use_active_clip:
            layout.prop(con, "clip")

        layout.use_property_split = True
        layout.operator("clip.constraint_to_fcurve")

        self.draw_influence(layout, con)

    def draw_object_solver(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        clip = None
        if con.use_active_clip:
            clip = context.scene.active_clip
        else:
            clip = con.clip

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_active_clip")
        row.prop_decorator(con, "use_active_clip")

        if not con.use_active_clip:
            layout.prop(con, "clip")

        if clip:
            layout.prop_search(con, "object", clip.tracking, "objects", icon='OBJECT_DATA')

        layout.prop(con, "camera")

        row = layout.row()
        row.operator("constraint.objectsolver_set_inverse")
        row.operator("constraint.objectsolver_clear_inverse")

        layout.operator("clip.constraint_to_fcurve")

        self.draw_influence(layout, con)

    def draw_transform_cache(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        layout.template_cache_file(con, "cache_file")

        cache_file = con.cache_file

        if cache_file is not None:
            layout.prop_search(con, "object_path", cache_file, "object_paths")

        self.draw_influence(layout, con)

    def draw_armature(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        col = layout.column()

        row = col.row()
        row.use_property_split = False
        row.prop(con, "use_deform_preserve_volume")
        row.prop_decorator(con, "use_deform_preserve_volume")

        row = col.row()
        row.use_property_split = False
        row.prop(con, "use_bone_envelopes")
        row.prop_decorator(con, "use_bone_envelopes")

        if context.pose_bone:
            col.prop(con, "use_current_location")

        layout.operator("constraint.add_target", text="Add Target Bone")

        layout.operator("constraint.normalize_target_weights")

        self.draw_influence(layout, con)

        if not con.targets:
            layout.label(text="No target bones added", icon='ERROR')

    def draw_kinematic(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        self.target_template(layout, con)

        if context.object.pose.ik_solver == 'ITASC':
            layout.prop(con, "ik_type")

            # This button gives itself too much padding, so put it in a column with the subtarget
            col = layout.column()
            col.prop(con, "pole_target")

            if con.pole_target and con.pole_target.type == 'ARMATURE':
                col.prop_search(con, "pole_subtarget", con.pole_target.data, "bones", text="Bone")

            col = layout.column()
            if con.pole_target:
                col.prop(con, "pole_angle")
            col.prop(con, "use_tail")
            col.prop(con, "use_stretch")
            col.prop(con, "chain_count")

            if con.ik_type == 'COPY_POSE':
                layout.prop(con, "reference_axis", expand=True)

                # Use separate rows and columns here to avoid an alignment issue with the lock buttons
                loc_col = layout.column()
                loc_col.prop(con, "use_location")

                row = loc_col.row()
                row.active = con.use_location
                row.prop(con, "weight", text="Weight", slider=True)

                row = loc_col.row(heading="Lock", align=True)
                row.use_property_decorate = False
                row.active = con.use_location
                sub = row.row(align=True)
                sub.prop(con, "lock_location_x", text="X", toggle=True)
                sub.prop(con, "lock_location_y", text="Y", toggle=True)
                sub.prop(con, "lock_location_z", text="Z", toggle=True)
                row.label(icon='BLANK1')

                rot_col = layout.column()
                rot_col.prop(con, "use_rotation")

                row = rot_col.row()
                row.active = con.use_rotation
                row.prop(con, "orient_weight", text="Weight", slider=True)

                row = rot_col.row(heading="Lock", align=True)
                row.use_property_decorate = False
                row.active = con.use_rotation
                sub = row.row(align=True)
                sub.prop(con, "lock_rotation_x", text="X", toggle=True)
                sub.prop(con, "lock_rotation_y", text="Y", toggle=True)
                sub.prop(con, "lock_rotation_z", text="Z", toggle=True)
                row.label(icon='BLANK1')

            elif con.ik_type == 'DISTANCE':
                layout.prop(con, "limit_mode")

                col = layout.column()
                col.prop(con, "weight", text="Weight", slider=True)
                col.prop(con, "distance", text="Distance", slider=True)
        else:
            # Standard IK constraint
            col = layout.column()
            col.prop(con, "pole_target")

            if con.pole_target and con.pole_target.type == 'ARMATURE':
                col.prop_search(con, "pole_subtarget", con.pole_target.data, "bones", text="Bone")

            col = layout.column()
            if con.pole_target:
                col.prop(con, "pole_angle")
            col.prop(con, "iterations")
            col.prop(con, "chain_count")

            row = col.row()
            row.use_property_split = False
            row.prop(con, "use_tail")
            row.prop_decorator(con, "use_tail")

            row = col.row()
            row.use_property_split = False
            row.prop(con, "use_stretch")
            row.prop_decorator(con, "use_stretch")

            split = layout.split(factor = 0.38)
            col = split.column()
            col.use_property_split = False
            row = col.row()
            row.prop(con, "use_location", text = "Weight Position")
            row.prop_decorator(con, "use_location")
            col = split.column()
            if con.use_location:
                row = col.row()
                row.prop(con, "weight", text="")
            else:
                col.label(icon='DISCLOSURE_TRI_RIGHT')

            split = layout.split(factor = 0.38)
            col = split.column()
            col.use_property_split = False
            row = col.row()
            row.prop(con, "use_rotation", text = "Rotation")
            row.prop_decorator(con, "use_rotation")
            col = split.column()
            if con.use_rotation:
                row = col.row()
                row.prop(con, "orient_weight", text="")
            else:
                col.label(icon='DISCLOSURE_TRI_RIGHT')

        self.draw_influence(layout, con)


# Parent class for constraint sub-panels.
class ConstraintButtonsSubPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_label = ""

    def get_constraint(self, _context):
        con = self.custom_data
        self.layout.context_pointer_set("constraint", con)
        return con

    def draw_transform_from(self, context):
        layout = self.layout
        con = self.get_constraint(context)

        layout.prop(con, "map_from", expand=True)

        layout.use_property_split = True
        layout.use_property_decorate = True

        from_axes = [con.map_to_x_from, con.map_to_y_from, con.map_to_z_from]

        if con.map_from == 'ROTATION':
            layout.prop(con, "from_rotation_mode", text="Mode")

        ext = "" if con.map_from == 'LOCATION' else "_rot" if con.map_from == 'ROTATION' else "_scale"

        col = layout.column(align=True)
        col.active = "X" in from_axes
        col.prop(con, "from_min_x" + ext, text="X Min")
        col.prop(con, "from_max_x" + ext, text="Max")

        col = layout.column(align=True)
        col.active = "Y" in from_axes
        col.prop(con, "from_min_y" + ext, text="Y Min")
        col.prop(con, "from_max_y" + ext, text="Max")

        col = layout.column(align=True)
        col.active = "Z" in from_axes
        col.prop(con, "from_min_z" + ext, text="Z Min")
        col.prop(con, "from_max_z" + ext, text="Max")

    def draw_transform_to(self, context):
        layout = self.layout
        con = self.get_constraint(context)

        layout.prop(con, "map_to", expand=True)

        layout.use_property_split = True
        layout.use_property_decorate = True

        if con.map_to == 'ROTATION':
            layout.prop(con, "to_euler_order", text="Order")

        ext = "" if con.map_to == 'LOCATION' else "_rot" if con.map_to == 'ROTATION' else "_scale"

        col = layout.column(align=True)
        col.prop(con, "map_to_x_from", expand=False, text="X Source Axis")
        col.prop(con, "to_min_x" + ext, text="Min")
        col.prop(con, "to_max_x" + ext, text="Max")

        col = layout.column(align=True)
        col.prop(con, "map_to_y_from", expand=False, text="Y Source Axis")
        col.prop(con, "to_min_y" + ext, text="Min")
        col.prop(con, "to_max_y" + ext, text="Max")

        col = layout.column(align=True)
        col.prop(con, "map_to_z_from", expand=False, text="Z Source Axis")
        col.prop(con, "to_min_z" + ext, text="Min")
        col.prop(con, "to_max_z" + ext, text="Max")

        layout.prop(con, "mix_mode" + ext, text="Mix", text_ctxt=i18n_contexts.constraint)

    def draw_armature_bones(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        for i, tgt in enumerate(con.targets):
            has_target = tgt.target is not None

            box = layout.box()
            header = box.row()
            header.use_property_split = False

            split = header.split(factor=0.45, align=True)
            split.prop(tgt, "target", text="")

            row = split.row(align=True)
            row.active = has_target
            if has_target:
                row.prop_search(tgt, "subtarget", tgt.target.data, "bones", text="")
            else:
                row.prop(tgt, "subtarget", text="", icon='BONE_DATA')

            header.operator("constraint.remove_target", text="", icon='X').index = i

            row = box.row()
            row.active = has_target and tgt.subtarget != ""
            row.prop(tgt, "weight", slider=True, text="Weight")

    def draw_spline_ik_fitting(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        col = layout.column()
        col.prop(con, "chain_count")

        row = col.row()
        row.use_property_split = False
        row.prop(con, "use_even_divisions")
        row.prop_decorator(con, "use_even_divisions")

        row = col.row()
        row.use_property_split = False
        row.prop(con, "use_chain_offset")
        row.prop_decorator(con, "use_chain_offset")

    def draw_spline_ik_chain_scaling(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_curve_radius")
        row.prop_decorator(con, "use_curve_radius")

        layout.prop(con, "y_scale_mode")
        layout.prop(con, "xz_scale_mode")

        if con.xz_scale_mode in {'INVERSE_PRESERVE', 'VOLUME_PRESERVE'}:

            row = layout.row()
            row.use_property_split = False
            row.prop(con, "use_original_scale")
            row.prop_decorator(con, "use_original_scale")

        if con.xz_scale_mode == 'VOLUME_PRESERVE':
            col = layout.column()
            col.prop(con, "bulge", text="Volume Variation")

            split = layout.split(factor = 0.38)
            col = split.column()
            col.use_property_split = False
            row = col.row()
            row.prop(con, "use_bulge_min", text = "Volume Min")
            row.prop_decorator(con, "use_bulge_min")
            col = split.column()
            if con.use_bulge_min:
                row = col.row()
                row.prop(con, "bulge_min", text="")
            else:
                col.label(icon='DISCLOSURE_TRI_RIGHT')

            split = layout.split(factor = 0.38)
            col = split.column()
            col.use_property_split = False
            row = col.row()
            row.prop(con, "use_bulge_max", text = "Volume Max")
            row.prop_decorator(con, "use_bulge_max")
            col = split.column()
            if con.use_bulge_max:
                row = col.row()
                row.prop(con, "bulge_max", text="")
            else:
                col.label(icon='DISCLOSURE_TRI_RIGHT')

            if con.use_bulge_min or con.use_bulge_max:
                row = layout.row()
                row.prop(con, "bulge_smooth", text="Smooth")

    def draw_action_target(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        col = layout.column()
        col.active = not con.use_eval_time
        col.prop(con, "transform_channel", text="Channel")
        ConstraintButtonsPanel.space_template(col, con, target=True, owner=False, separator=False)

        sub = col.column(align=True)
        sub.prop(con, "min", text="Range Min")
        sub.prop(con, "max", text="Max")

    def draw_action_action(self, context):
        layout = self.layout
        con = self.get_constraint(context)
        layout.use_property_split = True
        layout.use_property_decorate = True

        col = layout.column(align=True)
        col.prop(con, "action")
        if con.action and con.action.is_action_layered:
            col.context_pointer_set("animated_id", con.id_data)
            col.template_search(
                con, "action_slot",
                con, "action_suitable_slots",
                new="",  # No use in making a new slot here.
                unlink="anim.slot_unassign_from_constraint",
                text="Slot",
            )


        row = layout.row()
        row.use_property_split = False
        row.prop(con, "use_bone_object_action")
        row.prop_decorator(con, "use_bone_object_action")

        col = layout.column(align=True)
        col.prop(con, "frame_start", text="Frame Start")
        col.prop(con, "frame_end", text="End")

    def draw_transform_cache_velocity(self, context):
        self.draw_transform_cache_subpanel(
            context, self.layout.template_cache_file_velocity
        )

    def draw_transform_cache_procedural(self, context):
        self.draw_transform_cache_subpanel(
            context, self.layout.template_cache_file_procedural
        )

    def draw_transform_cache_time(self, context):
        self.draw_transform_cache_subpanel(
            context, self.layout.template_cache_file_time_settings
        )

    def draw_transform_cache_layers(self, context):
        self.draw_transform_cache_subpanel(
            context, self.layout.template_cache_file_layers
        )

    def draw_transform_cache_subpanel(self, context, template_func):
        con = self.get_constraint(context)
        if con.cache_file is None:
            return

        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = True
        template_func(con, "cache_file")

# Child Of Constraint


class OBJECT_PT_bChildOfConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_childof(context)


class BONE_PT_bChildOfConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_childof(context)

# Track To Constraint


class OBJECT_PT_bTrackToConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_trackto(context)


class BONE_PT_bTrackToConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_trackto(context)

# Follow Path Constraint


class OBJECT_PT_bFollowPathConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_follow_path(context)


class BONE_PT_bFollowPathConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_follow_path(context)


# Rotation Limit Constraint

class OBJECT_PT_bRotLimitConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_rot_limit(context)


class BONE_PT_bRotLimitConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_rot_limit(context)


# Location Limit Constraint

class OBJECT_PT_bLocLimitConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_loc_limit(context)


class BONE_PT_bLocLimitConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_loc_limit(context)


# Size Limit Constraint

class OBJECT_PT_bSizeLimitConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_size_limit(context)


class BONE_PT_bSizeLimitConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_size_limit(context)


# Rotate Like Constraint

class OBJECT_PT_bRotateLikeConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_rotate_like(context)


class BONE_PT_bRotateLikeConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_rotate_like(context)


# Locate Like Constraint

class OBJECT_PT_bLocateLikeConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_locate_like(context)


class BONE_PT_bLocateLikeConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_locate_like(context)


# Size Like Constraint

class OBJECT_PT_bSizeLikeConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_size_like(context)


class BONE_PT_bSizeLikeConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_size_like(context)


# Same Volume Constraint

class OBJECT_PT_bSameVolumeConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_same_volume(context)


class BONE_PT_bSameVolumeConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_same_volume(context)


# Trans Like Constraint

class OBJECT_PT_bTransLikeConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_trans_like(context)


class BONE_PT_bTransLikeConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_trans_like(context)


# Action Constraint

class OBJECT_PT_bActionConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_action(context)


class BONE_PT_bActionConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_action(context)


class OBJECT_PT_bActionConstraint_target(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bActionConstraint"
    bl_label = "Target"

    def draw(self, context):
        self.draw_action_target(context)


class BONE_PT_bActionConstraint_target(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bActionConstraint"
    bl_label = "Target"

    def draw(self, context):
        self.draw_action_target(context)


class OBJECT_PT_bActionConstraint_action(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bActionConstraint"
    bl_label = "Action"

    def draw(self, context):
        self.draw_action_action(context)


class BONE_PT_bActionConstraint_action(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bActionConstraint"
    bl_label = "Action"

    def draw(self, context):
        self.draw_action_action(context)


# Lock Track Constraint

class OBJECT_PT_bLockTrackConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_lock_track(context)


class BONE_PT_bLockTrackConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_lock_track(context)


# Distance Limit Constraint

class OBJECT_PT_bDistLimitConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_dist_limit(context)


class BONE_PT_bDistLimitConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_dist_limit(context)


# Stretch To Constraint

class OBJECT_PT_bStretchToConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_stretch_to(context)


class BONE_PT_bStretchToConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_stretch_to(context)


# Min Max Constraint

class OBJECT_PT_bMinMaxConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_min_max(context)


class BONE_PT_bMinMaxConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_min_max(context)


# Clamp To Constraint

class OBJECT_PT_bClampToConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_clamp_to(context)


class BONE_PT_bClampToConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_clamp_to(context)


# Transform Constraint

class OBJECT_PT_bTransformConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_transform(context)


class BONE_PT_bTransformConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_transform(context)


class OBJECT_PT_bTransformConstraint_source(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformConstraint"
    bl_label = "Map From"

    def draw(self, context):
        self.draw_transform_from(context)


class BONE_PT_bTransformConstraint_from(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformConstraint"
    bl_label = "Map From"

    def draw(self, context):
        self.draw_transform_from(context)


class OBJECT_PT_bTransformConstraint_destination(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformConstraint"
    bl_label = "Map To"

    def draw(self, context):
        self.draw_transform_to(context)


class BONE_PT_bTransformConstraint_to(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformConstraint"
    bl_label = "Map To"

    def draw(self, context):
        self.draw_transform_to(context)


# Shrink-wrap Constraint.

class OBJECT_PT_bShrinkwrapConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_shrinkwrap(context)


class BONE_PT_bShrinkwrapConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_shrinkwrap(context)


# Damp Track Constraint

class OBJECT_PT_bDampTrackConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_damp_track(context)


class BONE_PT_bDampTrackConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_damp_track(context)


# Spline IK Constraint

class BONE_PT_bSplineIKConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_spline_ik(context)


class BONE_PT_bSplineIKConstraint_fitting(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bSplineIKConstraint"
    bl_label = "Fitting"

    def draw(self, context):
        self.draw_spline_ik_fitting(context)


class BONE_PT_bSplineIKConstraint_chain_scaling(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bSplineIKConstraint"
    bl_label = "Chain Scaling"

    def draw(self, context):
        self.draw_spline_ik_chain_scaling(context)


# Pivot Constraint

class OBJECT_PT_bPivotConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_pivot(context)


class BONE_PT_bPivotConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_pivot(context)


# Follow Track Constraint

class OBJECT_PT_bFollowTrackConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_follow_track(context)


class BONE_PT_bFollowTrackConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_follow_track(context)


# Camera Solver Constraint

class OBJECT_PT_bCameraSolverConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_camera_solver(context)


class BONE_PT_bCameraSolverConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_camera_solver(context)


# Object Solver Constraint

class OBJECT_PT_bObjectSolverConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_object_solver(context)


class BONE_PT_bObjectSolverConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_object_solver(context)


# Transform Cache Constraint

class OBJECT_PT_bTransformCacheConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_transform_cache(context)


class BONE_PT_bTransformCacheConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_transform_cache(context)


class OBJECT_PT_bTransformCacheConstraint_velocity(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformCacheConstraint"
    bl_label = "Velocity"

    def draw(self, context):
        self.draw_transform_cache_velocity(context)


class BONE_PT_bTransformCacheConstraint_velocity(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformCacheConstraint"
    bl_label = "Velocity"

    def draw(self, context):
        self.draw_transform_cache_velocity(context)


class OBJECT_PT_bTransformCacheConstraint_layers(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformCacheConstraint"
    bl_label = "Override Layers"

    def draw(self, context):
        self.draw_transform_cache_layers(context)


class BONE_PT_bTransformCacheConstraint_layers(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformCacheConstraint"
    bl_label = "Override Layers"

    def draw(self, context):
        self.draw_transform_cache_layers(context)


class OBJECT_PT_bTransformCacheConstraint_procedural(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformCacheConstraint"
    bl_label = "Render Procedural"

    def draw(self, context):
        self.draw_transform_cache_procedural(context)


class BONE_PT_bTransformCacheConstraint_procedural(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformCacheConstraint"
    bl_label = "Render Procedural"

    def draw(self, context):
        self.draw_transform_cache_procedural(context)


class OBJECT_PT_bTransformCacheConstraint_time(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bTransformCacheConstraint"
    bl_label = "Time"

    def draw(self, context):
        self.draw_transform_cache_time(context)


class BONE_PT_bTransformCacheConstraint_time(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bTransformCacheConstraint"
    bl_label = "Time"

    def draw(self, context):
        self.draw_transform_cache_time(context)


# Armature Constraint

class OBJECT_PT_bArmatureConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_armature(context)


class BONE_PT_bArmatureConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_armature(context)


class OBJECT_PT_bArmatureConstraint_bones(ObjectConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "OBJECT_PT_bArmatureConstraint"
    bl_label = "Bones"

    def draw(self, context):
        self.draw_armature_bones(context)


class BONE_PT_bArmatureConstraint_bones(BoneConstraintPanel, ConstraintButtonsSubPanel, Panel):
    bl_parent_id = "BONE_PT_bArmatureConstraint"
    bl_label = "Bones"

    def draw(self, context):
        self.draw_armature_bones(context)


# Inverse Kinematic Constraint

class OBJECT_PT_bKinematicConstraint(ObjectConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_kinematic(context)


class BONE_PT_bKinematicConstraint(BoneConstraintPanel, ConstraintButtonsPanel, Panel):
    def draw(self, context):
        self.draw_kinematic(context)


classes = (
    # Object Panels
    OBJECT_PT_constraints,
    OBJECT_MT_constraint_add, # BFA menu
    OBJECT_OT_add_constraints_menu, # BFA menu
    BONE_PT_constraints,
    BONE_MT_constraint_add, # BFA menu
    BONE_OT_add_constraints_menu, # BFA menu
    OBJECT_PT_bChildOfConstraint,
    OBJECT_PT_bTrackToConstraint,
    OBJECT_PT_bKinematicConstraint,
    OBJECT_PT_bFollowPathConstraint,
    OBJECT_PT_bRotLimitConstraint,
    OBJECT_PT_bLocLimitConstraint,
    OBJECT_PT_bSizeLimitConstraint,
    OBJECT_PT_bRotateLikeConstraint,
    OBJECT_PT_bLocateLikeConstraint,
    OBJECT_PT_bSizeLikeConstraint,
    OBJECT_PT_bSameVolumeConstraint,
    OBJECT_PT_bTransLikeConstraint,
    OBJECT_PT_bActionConstraint,
    OBJECT_PT_bActionConstraint_target,
    OBJECT_PT_bActionConstraint_action,
    OBJECT_PT_bLockTrackConstraint,
    OBJECT_PT_bDistLimitConstraint,
    OBJECT_PT_bStretchToConstraint,
    OBJECT_PT_bMinMaxConstraint,
    OBJECT_PT_bClampToConstraint,
    OBJECT_PT_bTransformConstraint,
    OBJECT_PT_bTransformConstraint_source,
    OBJECT_PT_bTransformConstraint_destination,
    OBJECT_PT_bShrinkwrapConstraint,
    OBJECT_PT_bDampTrackConstraint,
    OBJECT_PT_bPivotConstraint,
    OBJECT_PT_bFollowTrackConstraint,
    OBJECT_PT_bCameraSolverConstraint,
    OBJECT_PT_bObjectSolverConstraint,
    OBJECT_PT_bTransformCacheConstraint,
    OBJECT_PT_bTransformCacheConstraint_time,
    OBJECT_PT_bTransformCacheConstraint_procedural,
    OBJECT_PT_bTransformCacheConstraint_velocity,
    OBJECT_PT_bTransformCacheConstraint_layers,
    OBJECT_PT_bArmatureConstraint,
    OBJECT_PT_bArmatureConstraint_bones,
    # Bone panels
    BONE_PT_bChildOfConstraint,
    BONE_PT_bTrackToConstraint,
    BONE_PT_bKinematicConstraint,
    BONE_PT_bFollowPathConstraint,
    BONE_PT_bRotLimitConstraint,
    BONE_PT_bLocLimitConstraint,
    BONE_PT_bSizeLimitConstraint,
    BONE_PT_bRotateLikeConstraint,
    BONE_PT_bLocateLikeConstraint,
    BONE_PT_bSizeLikeConstraint,
    BONE_PT_bSameVolumeConstraint,
    BONE_PT_bTransLikeConstraint,
    BONE_PT_bActionConstraint,
    BONE_PT_bActionConstraint_target,
    BONE_PT_bActionConstraint_action,
    BONE_PT_bLockTrackConstraint,
    BONE_PT_bDistLimitConstraint,
    BONE_PT_bStretchToConstraint,
    BONE_PT_bMinMaxConstraint,
    BONE_PT_bClampToConstraint,
    BONE_PT_bTransformConstraint,
    BONE_PT_bTransformConstraint_from,
    BONE_PT_bTransformConstraint_to,
    BONE_PT_bShrinkwrapConstraint,
    BONE_PT_bDampTrackConstraint,
    BONE_PT_bSplineIKConstraint,
    BONE_PT_bSplineIKConstraint_fitting,
    BONE_PT_bSplineIKConstraint_chain_scaling,
    BONE_PT_bPivotConstraint,
    BONE_PT_bFollowTrackConstraint,
    BONE_PT_bCameraSolverConstraint,
    BONE_PT_bObjectSolverConstraint,
    BONE_PT_bTransformCacheConstraint,
    BONE_PT_bTransformCacheConstraint_time,
    BONE_PT_bTransformCacheConstraint_procedural,
    BONE_PT_bTransformCacheConstraint_velocity,
    BONE_PT_bTransformCacheConstraint_layers,
    BONE_PT_bArmatureConstraint,
    BONE_PT_bArmatureConstraint_bones,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
