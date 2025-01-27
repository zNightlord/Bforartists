/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BLI_math_quaternion_types.hh"

#include "BKE_attribute.hh"

namespace blender::nodes::xpbd_constraints {

struct ConstraintEvalInputs {
  /* TODO split this by EvaluationTarget, only either (positions + rotations) or (velocities +
   * angular_velocities) are used. */
  VArraySpan<float3> positions;
  VArraySpan<math::Quaternion> rotations;
  Array<float3> positions_buffer;
  Array<math::Quaternion> rotations_buffer;
  VArraySpan<float3> old_positions;
  VArraySpan<math::Quaternion> old_rotations;

  VArraySpan<float3> velocities;
  VArraySpan<float3> angular_velocities;
  Array<float3> velocities_buffer;
  Array<float3> angular_velocities_buffer;
  VArraySpan<float3> orig_velocities;
  VArraySpan<float3> orig_angular_velocities;

  VArraySpan<float> position_weights;
  VArraySpan<float> rotation_weights;

  Span<float4x4> collider_transforms;
  /* Collider transforms of the previous frame for computing velocity constraints. */
  Span<float4x4> old_collider_transforms;

  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<float3> goals;
  } position_goal;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<math::Quaternion> goals;
  } rotation_goal;
  struct {
    VArray<int> solver_group;
    VArraySpan<float> alpha;
    VArraySpan<float> gamma;
    VArraySpan<int> points1;
    VArraySpan<int> points2;
    VArraySpan<float> edge_length;
  } stretch_shear;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points1;
    VArraySpan<int> points2;
    VArraySpan<float3> darboux_vector;
  } bend_twist;
  struct {
    VArray<int> solver_group;
    VArraySpan<int> points;
    VArraySpan<int> collider_index;
    VArraySpan<float3> local_position1;
    VArraySpan<float3> local_position2;
    VArraySpan<float3> normal;

    VArraySpan<float> friction;
    VArraySpan<float> restitution;
  } contact;
};

struct ConstraintEvalOutputs {
  Array<float3> positions;
  Array<math::Quaternion> rotations;
  Array<float3> velocities;
  Array<float3> angular_velocities;

  struct {
    /* Remember active contacts for later velocity update. */
    bke::SpanAttributeWriter<bool> active;
  } contact;
};

}  // namespace blender::nodes::xpbd_constraints
