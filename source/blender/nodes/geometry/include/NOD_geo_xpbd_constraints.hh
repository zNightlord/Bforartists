/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BLI_math_quaternion.hh"
#include "BLI_math_quaternion_types.hh"

#include "BKE_attribute.hh"

namespace blender::nodes::xpbd_constraints {

inline void eval_position_goal(const float3 &position,
                               const float3 &goal,
                               float3 &delta_lambda,
                               float3 &delta_position)
{
  const float3 residual = position - goal;
  /* Gradient is identity transform. */
  delta_lambda = -residual;
  delta_position = delta_lambda;
}

inline void eval_position_stretch_shear(const float3 &position1,
                                        const float3 &position2,
                                        const math::Quaternion &rotation,
                                        const float weight_pos1,
                                        const float weight_pos2,
                                        const float weight_rot,
                                        const float edge_length,
                                        const float alpha,
                                        const float gamma,
                                        float3 &delta_lambda,
                                        float3 &delta_position1,
                                        float3 &delta_position2,
                                        float4 &delta_rotation)
{
  const float inv_edge_length = math::safe_rcp(edge_length);

  const float3 direction = math::transform_point(rotation, float3(0, 0, 1));
  const float3 residual = inv_edge_length * (position2 - position1) - direction;
  const float weight_norm = math::safe_rcp(
      (weight_pos1 + weight_pos2) * inv_edge_length * inv_edge_length + 4.0f * weight_rot);

  delta_lambda = weight_norm * residual;

  delta_position1 = weight_pos1 * inv_edge_length * delta_lambda;
  delta_position2 = weight_pos2 * inv_edge_length * delta_lambda;
  delta_rotation = weight_rot * float4(math::Quaternion(0.0f, delta_lambda) * rotation *
                                       math::Quaternion(0, 0, 0, -1));
}

inline bool eval_position_contact(const float3 &position,
                                  const math::Quaternion &rotation,
                                  const float4x4 &collider_transform,
                                  const float3 &local_position1,
                                  const float3 &local_position2,
                                  const float3 &normal,
                                  float &delta_lambda,
                                  float3 &delta_position,
                                  float4 &delta_rotation)
{
  /* Local positions are relative to moving point and collider respectively.
   * Normal is a fixed shared direction for both participants. */

  const float4x4 point_transform = math::from_loc_rot<float4x4>(position, rotation);
  /* Contact points are computed by applying the transforms to relative local positions. */
  const float3 contact_point1 = math::transform_point(point_transform, local_position1);
  const float3 contact_point2 = math::transform_point(collider_transform, local_position2);

  /* Positional constraint for penetration depth along the normal. */
  const float residual_depth = math::dot(contact_point1 - contact_point2, normal);
  /* Only act on contact. */
  const bool active = residual_depth < 0.0f;
  if (!active) {
    return false;
  }

  /* Gradient is normal, length is 1, no need to compute gradient norm. */
  delta_lambda = -residual_depth;
  const float3 impulse = delta_lambda * normal;
  delta_position = impulse;
  delta_rotation = 0.5f * float4(math::Quaternion(0.0f, math::cross(local_position1, impulse)) *
                                 rotation);
  return true;
}

inline void eval_velocity_contact(const float3 &point_velocity,
                                  const float3 &point_angular_velocity,
                                  const float3 &orig_point_velocity,
                                  const float3 &orig_point_angular_velocity,
                                  const float3 &collider_velocity,
                                  const float3 &collider_angular_velocity,
                                  const float3 &local_position1,
                                  const float3 &local_position2,
                                  const float3 &normal,
                                  const float restitution,
                                  const float friction,
                                  float &delta_lambda_restitution,
                                  float &delta_lambda_friction,
                                  float3 &delta_velocity,
                                  float3 &delta_angular_velocity)
{
  /* Compute velocity of the collider contact point. */
  const float3 point_contact_velocity = point_velocity +
                                        math::cross(point_angular_velocity, local_position1);
  const float3 orig_point_contact_velocity = orig_point_velocity +
                                             math::cross(orig_point_angular_velocity,
                                                         local_position1);
  const float3 collider_contact_velocity = collider_velocity -
                                           math::cross(collider_angular_velocity, local_position2);

  /* Relative contact velocity before and after position corrections. */
  const float3 relative_velocity = point_contact_velocity - collider_contact_velocity;
  const float3 orig_relative_velocity = orig_point_contact_velocity - collider_contact_velocity;

  /* Decompose into normal and tangential velocity. */
  const float normal_velocity = math::dot(relative_velocity, normal);
  const float orig_normal_velocity = math::dot(orig_relative_velocity, normal);
  const float3 surface_velocity = relative_velocity - normal * normal_velocity;

  /* Folded into delta. */
  /* const float residual_restitution = math::dot(relative_velocity, normal); */
  /* const float residual_friction = math::length(surface_velocity); */

  /* Kill normal velocity, then add restitution. */
  delta_lambda_restitution = -std::min(restitution * orig_normal_velocity, 0.0f);
  delta_lambda_friction = -friction * math::length(surface_velocity);
  const float3 impulse_restitution = (-normal_velocity + delta_lambda_restitution) * normal;
  const float3 impulse_friction = -friction * surface_velocity;
  const float3 impulse = impulse_restitution + impulse_friction;

  delta_velocity = impulse;
  delta_angular_velocity = math::cross(local_position1, impulse);
}

struct ConstraintEvalParams {
  /* TODO split this by EvaluationTarget, only either (positions + rotations) or (velocities +
   * angular_velocities) are used. */
  Array<float3> positions;
  Array<math::Quaternion> rotations;
  VArraySpan<float3> old_positions;
  VArraySpan<math::Quaternion> old_rotations;

  Array<float3> velocities;
  Array<float3> angular_velocities;
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

    /* Remember active contacts for later velocity update. */
    bke::SpanAttributeWriter<bool> active;
  } contact;
};

}  // namespace blender::nodes::xpbd_constraints
