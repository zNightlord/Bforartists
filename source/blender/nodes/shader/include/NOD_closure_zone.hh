/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * Closure zones as used by texture layer adjustments (see
 * #NOD_texture_stack.hh): pairing queries that walk the raw node list so they
 * stay correct while the tree is being edited, and the edits that wrap a
 * node in a closure zone to turn a layer into an adjustment.
 */

struct Main;
struct bNode;
struct bNodeTree;

namespace blender::nodes::closure_zone {

/** Closure Output zone node. */
bool is_output_node(const bNode &node);

/** The paired Closure Input node of #zone_output, or null. */
bNode *input_node(bNodeTree &ntree, const bNode &zone_output);

/** The primary node inside the closure zone: the one feeding #zone_output's
 * first Bundle input (for adjustment layers, the adjustment group). Null when
 * the zone is empty. */
bNode *body_node(const bNodeTree &ntree, bNode &zone_output);

/**
 * Wrap #group_node in a closure zone and wire the zone's Closure output
 * into layer #layer_index on #stack_node, which becomes closure-typed:
 *
 *   Closure Input -> group -> Closure Output --(Closure)--> stack layer
 *
 * The inliner evaluates the closure with the stack accumulated below the
 * layer as the zone's bundle input, per the design's adjustment semantics.
 */
bool wrap_group(
    Main &bmain, bNodeTree &ntree, bNode &group_node, bNode &stack_node, int layer_index);

/**
 * Wrap a nested Texture Layer Stack in a closure zone so it becomes an
 * adjustment group on layer #layer_index of #outer_stack:
 *
 *   Closure Input --(bundle)--> nested stack base
 *   nested stack Result --(bundle)--> Closure Output --(Closure)--> outer layer
 *
 * The closure zone's incoming bundle (everything below the group) feeds the
 * nested stack's base, so the nested layers adjust the accumulated stack. The
 * outer layer becomes closure-typed.
 */
bool wrap_stack(
    Main &bmain, bNodeTree &ntree, bNode &nested_stack, bNode &outer_stack, int layer_index);

}  // namespace blender::nodes::closure_zone
