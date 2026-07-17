/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#include <algorithm>

#include "DNA_array_utils.hh"
#include "DNA_node_types.h"

#include "BKE_main_invariants.hh"
#include "BKE_node_tree_update.hh"

#include "NOD_layer_stack.hh"
#include "NOD_mask_stack.hh"
#include "NOD_texture_stack.hh"

namespace blender::nodes::layer_stack {

int active_index_in_range(const bNode &stack_node)
{
  const NodeShaderLayerStack &data = storage(stack_node);
  if (data.active_index < 0 || data.active_index >= data.items_num) {
    return -1;
  }
  return data.active_index;
}

/* Append a new item via the accessor matching #stack_node's kind. */
static void append_layer(bNode &stack_node, const char *name)
{
  if (mask_stack::is_stack(stack_node)) {
    socket_items::add_item_with_name<mask_stack::ItemsAccessor>(stack_node, name);
  }
  else {
    socket_items::add_item_with_name<texture_stack::ItemsAccessor>(stack_node, name);
  }
}

int add_layer_at(
    Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name, const int target)
{
  append_layer(stack_node, name);
  NodeShaderLayerStack &data = storage(stack_node);
  const int count = data.items_num;
  const int index = std::clamp(target, 0, count - 1);
  if (index != count - 1) {
    dna::array::move_index(data.items, count, count - 1, index);
  }
  BKE_ntree_update_tag_node_property(&ntree, &stack_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
  return index;
}

int add_layer_above_active(Main &bmain, bNodeTree &ntree, bNode &stack_node, const char *name)
{
  NodeShaderLayerStack &data = storage(stack_node);
  /* Capture the active index before adding: #append_layer appends the new item
   * and resets active_index to it, so this must be read first to insert the
   * new layer above the previously active one. */
  const int target = std::max(data.active_index, 0);
  append_layer(stack_node, name);
  int new_index = data.items_num - 1;
  if (target < new_index) {
    dna::array::move_index(data.items, data.items_num, new_index, target);
    new_index = target;
  }
  data.active_index = new_index;
  BKE_ntree_update_tag_node_property(&ntree, &stack_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
  return new_index;
}

void remove_layer(Main &bmain, bNodeTree &ntree, bNode &stack_node, const int index)
{
  /* Both stack node kinds share the same storage and the item access and
   * destruction logic in #ItemsAccessorBase, so either accessor removes an
   * item from either node kind. */
  socket_items::SocketItemsRef ref = texture_stack::ItemsAccessor::get_items_from_node(stack_node);
  dna::array::remove_index(ref.items,
                           ref.items_num,
                           ref.active_index,
                           index,
                           texture_stack::ItemsAccessor::destruct_item);
  BKE_ntree_update_tag_node_property(&ntree, &stack_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
}

}  // namespace blender::nodes::layer_stack
