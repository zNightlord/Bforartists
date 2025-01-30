/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_math_rotation.hh"

#include "NOD_xpbd_constraints.hh"

#include "testing/testing.h"

namespace blender::nodes::tests {

TEST(xpbd_constraints, PositionGoal)
{
  float lambda = 0.0f;
  float3 position = {1, 2, 3};

  const float alpha = 0.0f;
  const float gamma = 0.0f;

  xpbd_constraints::eval_position_goal(float3(-2, 0, 2), alpha, gamma, lambda, position);
  EXPECT_NEAR(-math::length(float3(-3, -2, -1)), lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(-2, 0, 2), position, 1e-5f);
}

TEST(xpbd_constraints, InactiveContact)
{
  float lambda = 0.0f;
  float3 position1 = {0, 0, 0};
  float3 position2 = {0, 0, 0};
  math::Quaternion rotation1 = math::Quaternion::identity();
  math::Quaternion rotation2 = math::Quaternion::identity();

  /* One-sided collision. */
  const float weight_pos1 = 1.0f;
  const float weight_pos2 = 0.0f;
  const float weight_rot1 = 1.0f;
  const float weight_rot2 = 0.0f;

  const float3 local_position1 = {0, 0, 1};
  const float3 local_position2 = {0, 0, -1};
  const float3 normal = {0, 0, 1};

  const float alpha = 0.0f;

  bool active = xpbd_constraints::eval_position_contact(weight_pos1,
                                                        weight_pos2,
                                                        weight_rot1,
                                                        weight_rot2,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        alpha,
                                                        lambda,
                                                        position1,
                                                        position2,
                                                        rotation1,
                                                        rotation2);
  EXPECT_FALSE(active);
}

TEST(xpbd_constraints, CentralContact)
{
  float lambda = 0.0f;
  float3 position1 = {0, 0, 0};
  float3 position2 = {0, 0, 0};
  math::Quaternion rotation1 = math::Quaternion::identity();
  math::Quaternion rotation2 = math::Quaternion::identity();

  /* One-sided collision. */
  const float weight_pos1 = 1.0f;
  const float weight_pos2 = 0.0f;
  const float weight_rot1 = 1.0f;
  const float weight_rot2 = 0.0f;

  const float3 local_position1 = {0, 0, -1};
  const float3 local_position2 = {0, 0, 1};
  const float3 normal = {0, 0, 1};

  const float alpha = 0.0f;

  bool active = xpbd_constraints::eval_position_contact(weight_pos1,
                                                        weight_pos2,
                                                        weight_rot1,
                                                        weight_rot2,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        alpha,
                                                        lambda,
                                                        position1,
                                                        position2,
                                                        rotation1,
                                                        rotation2);

  EXPECT_TRUE(active);
  EXPECT_NEAR(2, lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(0, 0, 2), position1, 1e-5f);
  EXPECT_V4_NEAR(float4(1, 0, 0, 0), float4(rotation1), 1e-5f);
}

TEST(xpbd_constraints, OffcenterContact)
{
  float lambda = 0.0f;
  float3 position1 = {0, 0, 0};
  float3 position2 = {0, 0, 0};
  math::Quaternion rotation1 = math::Quaternion::identity();
  math::Quaternion rotation2 = math::Quaternion::identity();

  /* One-sided collision. */
  const float weight_pos1 = 1.0f;
  const float weight_pos2 = 0.0f;
  const float weight_rot1 = 1.0f;
  const float weight_rot2 = 0.0f;

  const float3 local_position1 = {-1, 0, 0};
  const float3 local_position2 = {0, 0, 1};
  const float3 normal = {0, 0, 1};

  const float alpha = 0.0f;

  bool active = xpbd_constraints::eval_position_contact(weight_pos1,
                                                        weight_pos2,
                                                        weight_rot1,
                                                        weight_rot2,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        alpha,
                                                        lambda,
                                                        position1,
                                                        position2,
                                                        rotation1,
                                                        rotation2);

  EXPECT_TRUE(active);
  EXPECT_NEAR(0.5f, lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(0, 0, 0.5f), position1, 1e-5f);
  EXPECT_V4_NEAR(math::normalize(float4(1, 0, 0.5f, 0)), float4(rotation1), 1e-5f);
}

TEST(xpbd_constraints, ContactFriction)
{
  float lambda_restitution = 0.0f;
  float lambda_friction = 0.0f;
  float3 velocity1 = {1, 0, 0};
  float3 velocity2 = {0, 0, 0};
  float3 angular_velocity1 = {0, 0, 0};
  float3 angular_velocity2 = {0, 0, 0};
  const float3 orig_velocity1 = {0, 0, 0};
  const float3 orig_velocity2 = {0, 0, 0};
  const float3 orig_angular_velocity1 = {0, 0, 0};
  const float3 orig_angular_velocity2 = {0, 0, 0};

  const float3 local_position1 = {1, 1, 0};
  const float3 local_position2 = {0, 0, 0};
  const float3 normal = {0, 0, 1};

  const float restitution = 0.0f;
  const float friction = 0.5f;

  xpbd_constraints::eval_velocity_contact(orig_velocity1,
                                          orig_velocity2,
                                          orig_angular_velocity1,
                                          orig_angular_velocity2,
                                          local_position1,
                                          local_position2,
                                          normal,
                                          restitution,
                                          friction,
                                          lambda_restitution,
                                          lambda_friction,
                                          velocity1,
                                          velocity2,
                                          angular_velocity1,
                                          angular_velocity2);
  EXPECT_EQ(0, lambda_restitution);
  EXPECT_NEAR(-0.5f, lambda_friction, 1e-5f);
  EXPECT_V3_NEAR(float3(0.5f, 0, 0), velocity1, 1e-5f);
  /* Friction impulse (-0.5, 0, 0) at point (1, 1, 0) causes torque around the Z axis. */
  EXPECT_V3_NEAR(float3(0, 0, 0.5f), angular_velocity1, 1e-5f);
}

TEST(xpbd_constraints, ContactRestitution)
{
  float lambda_restitution = 0.0f;
  float lambda_friction = 0.0f;
  float3 velocity1 = {1, 0, -1};
  float3 velocity2 = {0, 0, 0};
  float3 angular_velocity1 = {0, 0, 0};
  float3 angular_velocity2 = {0, 0, 0};
  const float3 orig_velocity1 = velocity1;
  const float3 orig_velocity2 = velocity2;
  const float3 orig_angular_velocity1 = angular_velocity1;
  const float3 orig_angular_velocity2 = angular_velocity2;

  const float3 local_position1 = {1, 1, 0};
  const float3 local_position2 = {0, 0, 0};
  const float3 normal = {0, 0, 1};

  const float restitution = 0.5f;
  const float friction = 0.0f;

  xpbd_constraints::eval_velocity_contact(orig_velocity1,
                                          orig_velocity2,
                                          orig_angular_velocity1,
                                          orig_angular_velocity2,
                                          local_position1,
                                          local_position2,
                                          normal,
                                          restitution,
                                          friction,
                                          lambda_restitution,
                                          lambda_friction,
                                          velocity1,
                                          velocity2,
                                          angular_velocity1,
                                          angular_velocity2);
  EXPECT_EQ(0.5f, lambda_restitution);
  EXPECT_NEAR(0, lambda_friction, 1e-5f);
  EXPECT_V3_NEAR(float3(1.0f, 0, 0.5f), velocity1, 1e-5f);
  /* Rebound impulse (0, 0, 0.5) at point (1, 1, 0) causes torque around X and Y. */
  EXPECT_V3_NEAR(float3(1.5f, -1.5f, 0), angular_velocity1, 1e-5f);
}

}  // namespace blender::nodes::tests
