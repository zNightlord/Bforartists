/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BLI_math_quaternion.hh"
#include "BLI_math_quaternion_types.hh"

#include "BKE_attribute.hh"

namespace blender::nodes::xpbd_constraints {

enum class ConstraintType {
  /* Set position of a point to a target vector. */
  PositionGoal,
  /* Set orientation of an edge to a target rotation. */
  RotationGoal,
  /* Enforces edge length and aligns forward direction with the edge vector. */
  StretchShear,
  /* Enforces angles between neighboring edges to their relative rest orientation. */
  BendTwist,
  /* Keep contact points from penetrating. */
  Contact,
};
constexpr int NumConstraintTypes = 5;

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
  // TODO
  UNUSED_VARS(alpha, gamma);

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

/**
 * Positional contact constraint based on
 * "Detailed Rigid Body Simulation with Extended Position Based Dynamics", Mueller et al., 2020
 */
inline bool eval_position_contact(const float weight_pos1,
                                  const float weight_pos2,
                                  const float weight_rot1,
                                  const float weight_rot2,
                                  const float3 &local_position1,
                                  const float3 &local_position2,
                                  const float3 &normal,
                                  const float alpha,
                                  float &lambda,
                                  float3 &position1,
                                  float3 &position2,
                                  math::Quaternion &rotation1,
                                  math::Quaternion &rotation2)
{
  /* Local positions are relative to colliders.
   * Normal is a fixed shared direction for both participants. */

  /* Contact points are computed by applying the transforms to relative local positions. */
  const float3 contact_point1 = math::transform_point(rotation1, local_position1) + position1;
  const float3 contact_point2 = math::transform_point(rotation2, local_position2) + position2;

  /* Positional constraint for penetration depth along the normal. */
  const float residual_depth = math::dot(contact_point1 - contact_point2, normal);
  /* Only act on contact. */
  const bool active = residual_depth < 0.0f;
  if (!active) {
    return false;
  }

  /* Effective mass correction to account for the effect of rotation on position displacement.
   * See section 3.3.1 "Positional Constraints" of the paper.
   * We use a simplified rotational weight factor instead of full inverse moment of inertia. */
  const float rot_factor1 = math::length_squared(local_position1) +
                            math::square(math::dot(local_position1, normal));
  const float rot_factor2 = math::length_squared(local_position2) +
                            math::square(math::dot(local_position2, normal));
  const float total_weight = weight_pos1 + weight_rot1 * rot_factor1 + weight_pos2 +
                             weight_rot2 * rot_factor2 + alpha;
  BLI_assert(!math::is_zero(total_weight));

  /* Gradient is normal, length is 1, no need to compute gradient norm. */
  const float delta_lambda = (-residual_depth - alpha * lambda) / total_weight;
  const float3 impulse1 = delta_lambda * normal;
  const float3 impulse2 = -impulse1;
  position1 += impulse1 * weight_pos1;
  position2 += impulse2 * weight_pos2;
  const float3 angular_impulse1 = math::cross(local_position1, impulse1);
  const float3 angular_impulse2 = math::cross(local_position2, impulse2);
  rotation1 = math::normalize(math::Quaternion(1.0f, angular_impulse1) * rotation1);
  rotation2 = math::normalize(math::Quaternion(1.0f, angular_impulse2) * rotation1);
  return true;
}

/**
 * Velocity contact constraint based on
 * "Detailed Rigid Body Simulation with Extended Position Based Dynamics", Mueller et al., 2020
 */
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

}  // namespace blender::nodes::xpbd_constraints
