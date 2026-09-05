/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * Placement of programmatically created nodes in the node editor, for code
 * that builds node setups on the user's behalf and wants them readable without
 * disturbing the user's own layout.
 */

struct bNode;

namespace blender::nodes {

/** Position #node relative to #anchor, #dx to the left and #dy up.
 * For nodes that belong to #anchor's row (a chain). */
void place_node_left_of(bNode &node, const bNode &anchor, float dx, float dy);

/** Rough on-screen height of #node, for placing nodes before the node editor
 * has drawn them (which is when their real bounds appear). Overestimates, so
 * that placement leaves gaps rather than overlaps. */
float estimated_node_height(const bNode &node);

}  // namespace blender::nodes
