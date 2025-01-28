/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

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

}  // namespace blender::nodes::tests
