/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include <atomic>

#include "BLI_math_quaternion.hh"
#include "BLI_math_quaternion_types.hh"

#include "BKE_attribute.hh"

namespace blender::nodes::xpbd_constraints {

namespace error_check {

/* Helper struct for checking that each variable is only written by one constraint.*/
template<bool enable> struct VariableChecker;

template<> struct VariableChecker<false> {
  VariableChecker(const IndexRange /*range*/) {}

  bool claim_variable(const int /*index*/)
  {
    return true;
  }

  bool has_overlap() const
  {
    return false;
  }
};

template<> struct VariableChecker<true> {
  Array<std::atomic_bool> variable_written;
  std::atomic_bool variable_overlap = false;

  VariableChecker(const IndexRange range)
  {
    variable_written.reinitialize(range.size());
    for (const int i : variable_written.index_range()) {
      variable_written[i].store(false, std::memory_order::memory_order_relaxed);
    }
  }

  bool claim_variable(const int index)
  {
    if (variable_written[index].exchange(true, std::memory_order_relaxed)) {
      variable_overlap.store(true, std::memory_order_relaxed);
      return true;
    }
    return false;
  }

  bool has_overlap() const
  {
    return variable_overlap.load(std::memory_order_relaxed);
  }
};

}  // namespace error_check

enum class ConstraintType {
  PositionGoal,
  RotationGoal,
  VelocityGoal,
  AngularVelocityGoal,
  StretchShear,
  BendTwist,
  ContactPosition,
  ContactVelocity,
};

inline void eval_position_goal(const float3 &goal_position,
                               const float alpha,
                               const float gamma,
                               float &lambda,
                               float3 &position)
{
  const float3 old_position = position;
  float residual;
  const float3 gradient = math::normalize_and_get_length(position - goal_position, residual);

  const float delta_lambda = (-residual - alpha * lambda -
                              gamma * math::dot(gradient, position - old_position)) /
                             ((1.0f + gamma) + alpha);
  lambda += delta_lambda;
  position += delta_lambda * gradient;
}

template<bool linearized_quaternion>
inline void eval_rotation_goal(const math::Quaternion &goal_rotation,
                               const float alpha,
                               const float gamma,
                               float &lambda,
                               math::Quaternion &rotation)
{
  const math::Quaternion old_rotation = rotation;
  math::AxisAngle axis_angle = math::to_axis_angle(rotation *
                                                   math::invert_normalized(goal_rotation));
  const float residual = axis_angle.angle().radian();
  const float3 gradient = axis_angle.axis();

  const math::Quaternion diff = math::invert_normalized(old_rotation) * rotation;
  const float delta_lambda = (-residual - alpha * lambda -
                              gamma * 0.5f * math::dot(float4(0.0f, gradient), float4(diff))) /
                             ((1.0f + gamma) + alpha);
  lambda += delta_lambda;

  /* Multiply Quaternion(0, gradient) * rotation. */
  if constexpr (linearized_quaternion) {
    // const float4 q = float4(-math::dot(gradient, rotation.imaginary_part()),
    //                               rotation.w * gradient +
    //                                   math::cross(gradient, rotation.imaginary_part()));
    const float4 q = float4(math::Quaternion(0.0f, gradient) * rotation);
    rotation = math::normalize(math::Quaternion(float4(rotation) + delta_lambda * 0.5f * q));
  }
  else {
    /* Normalize the final result to avoid accumulating errors. */
    constexpr bool normalize_final = true;

    const math::Quaternion q = math::to_quaternion(
        math::AxisAngle(gradient, math::AngleRadian(delta_lambda)));

    rotation = q * rotation;
    if constexpr (normalize_final) {
      rotation = math::normalize(rotation);
    }
  }
}

inline void eval_velocity_goal(const float3 &goal_velocity,
                               const float alpha,
                               float &lambda,
                               float3 &velocity)
{
  float residual;
  const float3 gradient = math::normalize_and_get_length(velocity - goal_velocity, residual);

  const float delta_lambda = (-residual - alpha * lambda) / (1.0f + alpha);
  lambda += delta_lambda;
  velocity += delta_lambda * gradient;
}

inline void eval_angular_velocity_goal(const float3 &goal_angular_velocity,
                                       const float alpha,
                                       float &lambda,
                                       float3 &angular_velocity)
{
  float residual;
  const float3 gradient = math::normalize_and_get_length(angular_velocity - goal_angular_velocity,
                                                         residual);

  const float delta_lambda = (-residual - alpha * lambda) / (1.0f + alpha);
  lambda += delta_lambda;
  angular_velocity += delta_lambda * gradient;
}

template<bool linearized_quaternion>
inline void eval_position_stretch_shear(const float weight_pos1,
                                        const float weight_pos2,
                                        const float weight_rot,
                                        const float edge_length,
                                        const float alpha,
                                        const float gamma,
                                        float3 &lambda,
                                        float3 &position1,
                                        float3 &position2,
                                        math::Quaternion &rotation)
{
  // TODO
  UNUSED_VARS(alpha, gamma);

  const float inv_edge_length = math::safe_rcp(edge_length);

  const float3 direction = math::transform_point(rotation, float3(0, 0, 1));
  const float3 residual = inv_edge_length * (position2 - position1) - direction;
  const float weight_norm = math::safe_rcp(
      (weight_pos1 + weight_pos2) * inv_edge_length * inv_edge_length + 4.0f * weight_rot);

  const float3 delta_lambda = weight_norm * residual;
  lambda = lambda + delta_lambda;

  position1 += weight_pos1 * inv_edge_length * delta_lambda;
  position2 -= weight_pos2 * inv_edge_length * delta_lambda;
  if constexpr (linearized_quaternion) {
    /* XXX is this correct? */
    const float4 delta_rot = weight_rot * float4(math::Quaternion(0.0f, delta_lambda) * rotation *
                                                 math::Quaternion(0, 0, 0, -1));
    rotation = math::normalize(math::Quaternion(float4(rotation) + delta_rot));
  }
  else {
    /* Normalize the final result to avoid accumulating errors. */
    constexpr bool normalize_final = true;

    const math::Quaternion q = math::to_quaternion(
        math::AxisAngle(direction, math::normalize(position2 - position1)));

    rotation = q * rotation;
    if constexpr (normalize_final) {
      rotation = math::normalize(rotation);
    }
  }
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
  const float rot_factor1 = math::length_squared(local_position1) -
                            math::square(math::dot(local_position1, normal));
  const float rot_factor2 = math::length_squared(local_position2) -
                            math::square(math::dot(local_position2, normal));
  const float total_weight = weight_pos1 + weight_rot1 * rot_factor1 + weight_pos2 +
                             weight_rot2 * rot_factor2 + alpha;
  BLI_assert(!math::is_zero(total_weight));

  /* Gradient is normal, length is 1, no need to compute gradient norm. */
  const float delta_lambda = (-residual_depth - alpha * lambda) / total_weight;
  const float3 impulse1 = delta_lambda * normal;
  const float3 impulse2 = -impulse1;

  lambda += delta_lambda;
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
inline void eval_velocity_contact(const float3 &orig_velocity1,
                                  const float3 &orig_velocity2,
                                  const float3 &orig_angular_velocity1,
                                  const float3 &orig_angular_velocity2,
                                  const float3 &local_position1,
                                  const float3 &local_position2,
                                  const float3 &normal,
                                  const float restitution,
                                  const float friction,
                                  float &delta_lambda_restitution,
                                  float &delta_lambda_friction,
                                  float3 &velocity1,
                                  float3 &velocity2,
                                  float3 &angular_velocity1,
                                  float3 &angular_velocity2)
{
  /* Compute velocity of the collider contact point. */
  const float3 contact_velocity1 = velocity1 + math::cross(angular_velocity1, local_position1);
  const float3 contact_velocity2 = velocity2 + math::cross(angular_velocity2, local_position2);
  const float3 orig_contact_velocity1 = orig_velocity1 +
                                        math::cross(orig_angular_velocity1, local_position1);
  const float3 orig_contact_velocity2 = orig_velocity2 -
                                        math::cross(orig_angular_velocity2, local_position2);

  /* Relative contact velocity before and after position corrections. */
  const float3 relative_velocity = contact_velocity1 - contact_velocity2;
  const float3 orig_relative_velocity = orig_contact_velocity1 - orig_contact_velocity2;

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

  velocity1 += impulse;
  velocity2 -= impulse;
  angular_velocity1 += math::cross(local_position1, impulse);
  angular_velocity2 -= math::cross(local_position2, impulse);
}

}  // namespace blender::nodes::xpbd_constraints
