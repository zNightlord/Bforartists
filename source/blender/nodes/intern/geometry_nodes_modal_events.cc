/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include <algorithm>

#include "BLI_map.hh"
#include "BLI_string.hh"

#include "DNA_node_types.h"

#include "BKE_node_runtime.hh"

#include "NOD_geometry_nodes_modal_events.hh"

namespace blender::nodes {

void ModalEvents::merge(const ModalEvents &other)
{
  this->events.extend(other.events);
}

void ModalEvents::deduplicate_and_sort()
{
  std::ranges::sort(this->events, [](const ModalEvent &a, const ModalEvent &b) {
    return BLI_strcasecmp_natural(a.name.c_str(), b.name.c_str()) < 0;
  });
  /* The first node that defines an event also provides its description. */
  const ModalEvent *unique_end = std::unique(
      this->events.begin(), this->events.end(), [](const ModalEvent &a, const ModalEvent &b) {
        return a.name == b.name;
      });
  this->events.resize(unique_end - this->events.begin());
}

static void gather_modal_events(const bNodeTree &ntree,
                                ModalEvents &events,
                                FunctionRef<const ModalEvents *(const bNodeTree &group)> get_group)
{
  ntree.ensure_topology_cache();
  for (const bNode *node : ntree.nodes_by_type("GeometryNodeModalEvent"_ustr)) {
    if (node->is_muted()) {
      continue;
    }
    const auto *storage = static_cast<const GeometryNodeModalEvent *>(node->storage);
    if (!storage->name || storage->name[0] == '\0') {
      continue;
    }
    events.events.append({storage->name, storage->description ? storage->description : ""});
  }
  for (const bNode *node : ntree.group_nodes()) {
    if (!node->id) {
      continue;
    }
    const bNodeTree &group = *reinterpret_cast<const bNodeTree *>(node->id);
    if (const ModalEvents *group_events = get_group(group)) {
      events.merge(*group_events);
    }
  }
}

ModalEvents gather_modal_events_with_cache(const bNodeTree &ntree)
{
  ModalEvents events;
  gather_modal_events(
      ntree, events, [](const bNodeTree &group) { return group.runtime->modal_events.get(); });
  events.deduplicate_and_sort();
  return events;
}

static void gather_modal_events_recursive_impl(const bNodeTree &ntree,
                                               Map<const bNodeTree *, ModalEvents> &events_by_tree)
{
  if (events_by_tree.contains(&ntree)) {
    return;
  }
  /* Add an empty set before recursing, so that node group recursion terminates. */
  events_by_tree.add_new(&ntree, {});
  ModalEvents new_events;
  gather_modal_events(ntree, new_events, [&](const bNodeTree &group) {
    gather_modal_events_recursive_impl(group, events_by_tree);
    return &events_by_tree.lookup(&group);
  });
  new_events.deduplicate_and_sort();
  events_by_tree.add_overwrite(&ntree, std::move(new_events));
}

ModalEvents gather_modal_events_recursive(const bNodeTree &ntree)
{
  Map<const bNodeTree *, ModalEvents> events_by_tree;
  gather_modal_events_recursive_impl(ntree, events_by_tree);
  return events_by_tree.lookup(&ntree);
}

}  // namespace blender::nodes
