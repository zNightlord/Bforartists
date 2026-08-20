/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include "BLI_function_ref.hh"
#include "BLI_listbase_iterator.hh"

#include "DNA_node_types.h"

#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_tree_update.hh"

#include "NOD_closure_zone.hh"
#include "NOD_geo_bundle.hh"
#include "NOD_geo_closure.hh"
#include "NOD_node_placement.hh"
#include "NOD_socket_items.hh"
#include "NOD_texture_stack.hh"

namespace blender::nodes::closure_zone {

bool is_output_node(const bNode &node)
{
  return node.is_type("NodeClosureOutput"_ustr);
}

bNode *input_node(bNodeTree &ntree, const bNode &zone_output)
{
  for (bNode &node : ntree.nodes) {
    if (!node.is_type("NodeClosureInput"_ustr)) {
      continue;
    }
    const auto &storage = *static_cast<const NodeClosureInput *>(node.storage);
    if (storage.output_node_id == zone_output.identifier) {
      return &node;
    }
  }
  return nullptr;
}

bNode *body_node(const bNodeTree &ntree, bNode &zone_output)
{
  bNodeSocket *bundle_in = first_bundle_input(zone_output);
  if (bundle_in == nullptr) {
    return nullptr;
  }
  return bke::node_find_source_node(ntree, *bundle_in);
}

/* Create a paired Closure Input/Output with a Bundle passing through, placed
 * around #body. */
static bool add_zone_around(
    bNodeTree &ntree, bNode &body, const float dx, bNode **r_closure_in, bNode **r_closure_out)
{
  bNode *closure_in = bke::node_add_node(nullptr, ntree, "NodeClosureInput"_ustr);
  bNode *closure_out = bke::node_add_node(nullptr, ntree, "NodeClosureOutput"_ustr);
  if (!closure_in || !closure_out) {
    return false;
  }
  static_cast<NodeClosureInput *>(closure_in->storage)->output_node_id = closure_out->identifier;
  socket_items::add_item_with_socket_type_and_name<ClosureInputItemsAccessor>(
      ntree, *closure_out, SOCK_BUNDLE, "Bundle");
  socket_items::add_item_with_socket_type_and_name<ClosureOutputItemsAccessor>(
      ntree, *closure_out, SOCK_BUNDLE, "Bundle");
  place_node_left_of(*closure_in, body, dx, 0.0f);
  place_node_left_of(*closure_out, body, -dx, 0.0f);
  *r_closure_in = closure_in;
  *r_closure_out = closure_out;
  return true;
}

/* Re-declare the layer's socket as a Closure input (same identifier) and
 * materialize the sockets of all the involved nodes. */
static void make_layer_closure_typed(
    Main &bmain, bNodeTree &ntree, bNode &stack_node, const int layer_index, bNode &in, bNode &out)
{
  layer_stack::storage(stack_node).items[layer_index].item_type = SHADER_LAYER_STACK_ITEM_CLOSURE;
  BKE_ntree_update_tag_node_property(&ntree, &stack_node);
  BKE_ntree_update_tag_node_property(&ntree, &in);
  BKE_ntree_update_tag_node_property(&ntree, &out);
  BKE_main_ensure_invariants(bmain, ntree.id);
}

/* Link #closure_out's Closure output into layer #layer_index. */
static bool connect_zone_to_layer(bNodeTree &ntree,
                                  bNode &closure_out,
                                  bNode &stack_node,
                                  const int layer_index)
{
  bNodeSocket *closure_socket = nullptr;
  for (bNodeSocket &sock : closure_out.outputs) {
    if (sock.type == SOCK_CLOSURE) {
      closure_socket = &sock;
      break;
    }
  }
  bNodeSocket *layer_socket =
      texture_stack::StackLayer{&ntree, &stack_node, layer_index}.bundle_socket();
  if (!closure_socket || !layer_socket) {
    return false;
  }
  bke::node_add_link(ntree, closure_out, *closure_socket, stack_node, *layer_socket);
  return true;
}

/* Put #body inside a closure zone that feeds layer #layer_index of #stack_node:
 * the zone's incoming bundle (the stack accumulated below that layer) goes into
 * the body through #body_bundle_input, the body's #body_bundle_output carries
 * the adjusted bundle back out, and the layer becomes a closure fed by the
 * zone's Closure Output.
 *
 * The body's sockets are looked up through callbacks rather than passed in:
 * #make_layer_closure_typed re-declares the nodes, so socket pointers taken
 * before it may no longer be valid. */
static bool wrap_body_in_zone(Main &bmain,
                              bNodeTree &ntree,
                              bNode &body,
                              const float dx,
                              bNode &stack_node,
                              const int layer_index,
                              const FunctionRef<bNodeSocket *()> body_bundle_input,
                              const FunctionRef<bNodeSocket *()> body_bundle_output)
{
  bNode *closure_in = nullptr;
  bNode *closure_out = nullptr;
  if (!add_zone_around(ntree, body, dx, &closure_in, &closure_out)) {
    return false;
  }
  make_layer_closure_typed(bmain, ntree, stack_node, layer_index, *closure_in, *closure_out);

  bNodeSocket *zone_in_bundle = first_bundle_output(*closure_in);
  bNodeSocket *body_in = body_bundle_input();
  if (zone_in_bundle && body_in) {
    bke::node_add_link(ntree, *closure_in, *zone_in_bundle, body, *body_in);
  }
  bNodeSocket *body_out = body_bundle_output();
  bNodeSocket *zone_out_bundle = first_bundle_input(*closure_out);
  if (body_out && zone_out_bundle) {
    bke::node_add_link(ntree, body, *body_out, *closure_out, *zone_out_bundle);
  }

  return connect_zone_to_layer(ntree, *closure_out, stack_node, layer_index);
}

bool wrap_group(
    Main &bmain, bNodeTree &ntree, bNode &group_node, bNode &stack_node, const int layer_index)
{
  return wrap_body_in_zone(
      bmain,
      ntree,
      group_node,
      200.0f,
      stack_node,
      layer_index,
      [&]() { return first_bundle_input(group_node); },
      [&]() { return first_bundle_output(group_node); });
}

bool wrap_stack(
    Main &bmain, bNodeTree &ntree, bNode &nested_stack, bNode &outer_stack, const int layer_index)
{
  return wrap_body_in_zone(
      bmain,
      ntree,
      nested_stack,
      260.0f,
      outer_stack,
      layer_index,
      [&]() {
        /* The incoming bundle feeds the nested stack's base layer, so the
         * nested layers composite on top of it. */
        const NodeShaderLayerStack &storage = layer_stack::storage(nested_stack);
        return texture_stack::StackLayer{&ntree, &nested_stack, storage.items_num - 1}
            .bundle_socket();
      },
      /* The nested stack's Result is the group's adjusted bundle. */
      [&]() { return bke::node_find_socket(nested_stack, SOCK_OUT, "Result"_ustr); });
}

}  // namespace blender::nodes::closure_zone
