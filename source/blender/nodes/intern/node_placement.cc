/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include "BLI_listbase_iterator.hh"

#include "DNA_node_types.h"

#include "NOD_node_placement.hh"

namespace blender::nodes {

void place_node_left_of(bNode &node, const bNode &anchor, const float dx, const float dy)
{
  node.location[0] = anchor.location[0] - dx;
  node.location[1] = anchor.location[1] + dy;
}

float estimated_node_height(const bNode &node)
{
  int socket_count = 0;
  for (const bNodeSocket &sock : node.inputs) {
    if (!(sock.flag & (SOCK_UNAVAIL | SOCK_HIDDEN))) {
      socket_count++;
    }
  }
  for (const bNodeSocket &sock : node.outputs) {
    if (!(sock.flag & (SOCK_UNAVAIL | SOCK_HIDDEN))) {
      socket_count++;
    }
  }
  return 60.0f + 25.0f * socket_count;
}

}  // namespace blender::nodes
