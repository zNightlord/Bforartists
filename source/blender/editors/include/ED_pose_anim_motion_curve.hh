#pragma once

/**
 * ED_pose_anim_motion_curve.hh
 *
 * Place this file in:
 *   source/blender/editors/include/ED_pose_anim_motion_curve.hh
 *
 * That directory is on the include path for EVERY editor module, so both
 * source/blender/editors/armature/  (implements the functions)
 * source/blender/editors/transform/ (calls them)
 * can include it without any CMakeLists changes.
 *
 * In pose_anim_motion_curve.cc  add:
 *   #include "ED_pose_anim_motion_curve.hh"
 *
 * In transform_convert_armature.cc  add:
 *   #include "ED_pose_anim_motion_curve.hh"
 * replacing the old:
 *   #include "pose_anim_motion_curve.hh"
 */
namespace blender {
struct Object;
struct bPoseChannel;

/**
 * Called from pose_grab_with_ik() in transform_convert_armature.cc
 * after the temporary IK constraint has been installed.
 *
 * Signals the motion curve gizmo to switch to its fast-path refresh mode:
 * instead of re-baking all FCurves it just overwrites the current-frame
 * data from the live (IK-solved) scene depsgraph each step.
 *
 * \param ob        Armature object being transformed.
 * \param tip_pchan Chain tip bone (the one being grabbed).
 */
void POSE_motion_curve_auto_ik_begin(Object *ob, bPoseChannel *tip_pchan);

/**
 * Called from special_aftertrans_update__pose() before
 * pose_grab_with_ik_clear(), whether the user confirmed or cancelled.
 *
 * Clears the fast-path mode and forces a full FCurve re-bake so the
 * path reflects the final (committed or rolled-back) animation state.
 *
 * \param cancelled  true when the transform was cancelled (Esc / RMB).
 */
void POSE_motion_curve_auto_ik_end(bool cancelled);
}
