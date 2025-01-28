/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_math_rotation.hh"

#include "NOD_geo_xpbd_constraints.hh"

#include "testing/testing.h"

namespace blender::nodes::tests {

TEST(xpbd_constraints, PositionGoal)
{
  float3 delta_lambda;
  float3 delta_position;

  xpbd_constraints::eval_position_goal(
      float3(1, 2, 3), float3(1, 2, 3), delta_lambda, delta_position);
  EXPECT_V3_NEAR(float3(0.0f), delta_lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(0.0f), delta_position, 1e-5f);

  xpbd_constraints::eval_position_goal(
      float3(1, 2, 3), float3(-2, 0, 2), delta_lambda, delta_position);
  EXPECT_V3_NEAR(float3(-3, -2, -1), delta_lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(-3, -2, -1), delta_position, 1e-5f);
}

TEST(xpbd_constraints, InactiveContact)
{
  const float3 position1 = {0, 0, 0};
  const float3 position2 = {0, 0, 0};
  const math::Quaternion rotation1 = math::Quaternion::identity();
  const math::Quaternion rotation2 = math::Quaternion::identity();
  const float4x4 collider_transform = math::from_loc_rot<float4x4>(position2, rotation2);

  const float3 local_position1 = {0, 0, 1};
  const float3 local_position2 = {0, 0, -1};
  const float3 normal = {0, 0, 1};

  float delta_lambda;
  float3 delta_position;
  float4 delta_rotation;
  bool active = xpbd_constraints::eval_position_contact(position1,
                                                        rotation1,
                                                        collider_transform,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        delta_lambda,
                                                        delta_position,
                                                        delta_rotation);
  EXPECT_FALSE(active);
}

TEST(xpbd_constraints, CentralContact)
{
  const float3 position1 = {0, 0, 0};
  const float3 position2 = {0, 0, 0};
  const math::Quaternion rotation1 = math::Quaternion::identity();
  const math::Quaternion rotation2 = math::Quaternion::identity();
  const float4x4 collider_transform = math::from_loc_rot<float4x4>(position2, rotation2);

  const float3 local_position1 = {0, 0, -1};
  const float3 local_position2 = {0, 0, 1};
  const float3 normal = {0, 0, 1};

  float delta_lambda;
  float3 delta_position;
  float4 delta_rotation;
  bool active = xpbd_constraints::eval_position_contact(position1,
                                                        rotation1,
                                                        collider_transform,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        delta_lambda,
                                                        delta_position,
                                                        delta_rotation);
  EXPECT_NEAR(2, delta_lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(0, 0, 2), delta_position, 1e-5f);
  EXPECT_V4_NEAR(float4(0, 0, 0, 0), delta_rotation, 1e-5f);
}

TEST(xpbd_constraints, OffcenterContact)
{
  const float3 position1 = {0, 0, 0};
  const float3 position2 = {0, 0, 0};
  const math::Quaternion rotation1 = math::Quaternion::identity();
  const math::Quaternion rotation2 = math::Quaternion::identity();
  const float4x4 collider_transform = math::from_loc_rot<float4x4>(position2, rotation2);

  const float3 local_position1 = {-1, 0, 0};
  const float3 local_position2 = {0, 0, 1};
  const float3 normal = {0, 0, 1};

  float delta_lambda;
  float3 delta_position;
  float4 delta_rotation;
  bool active = xpbd_constraints::eval_position_contact(position1,
                                                        rotation1,
                                                        collider_transform,
                                                        local_position1,
                                                        local_position2,
                                                        normal,
                                                        delta_lambda,
                                                        delta_position,
                                                        delta_rotation);
  EXPECT_NEAR(1, delta_lambda, 1e-5f);
  EXPECT_V3_NEAR(float3(0, 0, 1), delta_position, 1e-5f);
  EXPECT_V4_NEAR(float4(0, 0, 0.5f, 0), delta_rotation, 1e-5f);
}

TEST(xpbd_constraints, ContactFriction)
{
  const float3 velocity1 = {1, 0, 0};
  const float3 orig_velocity1 = {0, 0, 0};
  const float3 angular_velocity1 = {0, 0, 0};
  const float3 orig_angular_velocity1 = {0, 0, 0};
  const float3 velocity2 = {0, 0, 0};
  const float3 angular_velocity2 = {0, 0, 0};

  const float3 local_position1 = {0, 0, 0};
  const float3 local_position2 = {0, 0, 0};
  const float3 normal = {0, 0, 1};

  const float restitution = 0.0f;
  const float friction = 0.5f;

  float delta_lambda_restitution, delta_lambda_friction;
  float3 delta_velocity;
  float3 delta_angular_velocity;
  xpbd_constraints::eval_velocity_contact(velocity1,
                                          angular_velocity1,
                                          orig_velocity1,
                                          orig_angular_velocity1,
                                          velocity2,
                                          angular_velocity2,
                                          local_position1,
                                          local_position2,
                                          normal,
                                          restitution,
                                          friction,
                                          delta_lambda_restitution,
                                          delta_lambda_friction,
                                          delta_velocity,
                                          delta_angular_velocity);
  // EXPECT_EQ(0, delta_lambda_restitution);
  // EXPECT_NEAR(, delta_lambda, 1e-5f);
  // EXPECT_V3_NEAR(float3(0, 0, 2), delta_position, 1e-5f);
  // EXPECT_V4_NEAR(float4(0, 0, 0, 0), delta_rotation, 1e-5f);
}

}  // namespace blender::nodes::tests
