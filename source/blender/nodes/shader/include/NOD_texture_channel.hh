/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup nodes
 *
 * Texture layer channels: which of a BSDF's inputs are routed through the
 * texture layer stack, and the Combine/Separate Bundle plumbing that carries
 * them. A "BSDF" here is a shader leaf: a node that outputs a shader and takes
 * no shader as input.
 *
 * Which channels are routed is defined by the topology alone: a BSDF input is
 * "textured" when it is fed by the Separate Bundle that the stack's Result
 * drives. The channel universe is the BSDF's fillable inputs; the enabled
 * subset is the Separate Bundle's items linked into the BSDF. These helpers
 * derive and edit that state. Everything walks the raw link lists, so it is
 * safe while the tree is being edited (e.g. from #bNodeType::insert_link).
 */

#include "BLI_string_ref.hh"
#include "BLI_vector.hh"

#include "DNA_node_types.h"

struct NodeCombineBundleItem;

namespace blender::nodes::texture_channel {

/** BSDF (shader leaf): outputs a shader and takes no shader as input. */
bool is_bsdf(const bNode &node);

/** Append every BSDF in #ntree to #r_bsdfs. */
void collect_bsdfs(bNodeTree &ntree, Vector<bNode *> &r_bsdfs);

/** True for inputs a texture-layer channel can drive: color/float/vector and
 * available. */
bool is_fillable_input(const bNodeSocket &socket);

/** The first fillable input on #bsdf, or null. Typically the BSDF's first
 * color input, e.g. Base Color on a Principled BSDF. */
bNodeSocket *first_fillable_input(bNode &bsdf);

/** Separate Bundle node (splits a bundle back into individual sockets). */
bool is_separate_bundle(const bNode &node);

/** How a BSDF input relates to the texture layer stack. */
enum class State {
  /** Plain value, not linked. */
  Free,
  /** Fed by a Separate Bundle that a Texture Layer Stack drives. */
  Textured,
  /** Linked to something else. */
  Linked,
};
State state(const bNodeTree &ntree, const bNodeSocket &bsdf_input);

/** The canonical Separate Bundle for #bsdf: the one feeding its inputs whose
 * bundle input comes from a Texture Layer Stack. When no channel is currently
 * linked, falls back to a stack-driven Separate Bundle that feeds no other
 * BSDF, so re-enabling a channel reuses the existing stack. Null when the
 * material has none. */
bNode *separate_bundle_for_bsdf(const bNodeTree &ntree, const bNode &bsdf);

/** Socket type of the channel named #channel on the BSDF consuming
 * #stack_node. Channels are the BSDF's inputs, matched by name; RGBA when
 * the BSDF or the channel cannot be resolved. */
eNodeSocketDatatype socket_type(const bNodeTree &ntree, bNode &stack_node, StringRef channel);

/** The base (last) layer's Combine Bundle on #stack, or null when the base
 * layer's bundle source is not a Combine Bundle. */
bNode *base_layer_combine_bundle(const bNodeTree &ntree, bNode &stack);

/** The base layer Combine Bundle's input socket carrying #channel, or null. */
bNodeSocket *base_layer_channel_socket(const bNodeTree &ntree, bNode &stack, StringRef channel);

/** Ensure the base layer of #stack carries #channel, seeded with
 * #bsdf_socket's current value, so routing the channel through the stack
 * does not change the render (a bundle key missing from all layers evaluates
 * to the type default, not the BSDF socket's value). Does nothing when the
 * base layer's source is not a Combine Bundle or the channel input is linked. */
void seed_base_value(bNodeTree &ntree,
                     bNode &stack,
                     StringRef channel,
                     const bNodeSocket &bsdf_socket);

/** Route #bsdf_socket through #separate: adds an item named after the
 * socket when missing, links it to the socket and seeds the base layer with
 * the socket's current value. Returns false when #separate is not driven by
 * a Texture Layer Stack. */
bool link(bNodeTree &ntree, bNode &separate, bNode &bsdf, bNodeSocket &bsdf_socket);

/** Undo the routing of a Textured #bsdf_socket: copy the base layer's value
 * back to the socket when it is a plain value, remove the link and drop the
 * Separate Bundle item unless other sockets still use it. The layers keep
 * their channel data, so re-enabling restores the full layered result.
 * Returns false when the socket is not in the Textured state. */
bool unlink(bNodeTree &ntree, bNodeSocket &bsdf_socket);

/** Complete a manually drawn link from #separate to a BSDF input: when the
 * linked item's name matches the target socket (the link enables that
 * channel), seed the stack's base layer with the socket's current value so
 * the render does not change. Name-mismatched links are left alone. Safe to
 * call from #bNodeType::insert_link, before the link is in the tree. */
void handle_manual_link(bNodeTree &ntree, bNode &separate, const bNodeLink &link);

}  // namespace blender::nodes::texture_channel
