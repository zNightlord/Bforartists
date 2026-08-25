/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/**
 * Modal node tools react to a set of named events that is defined by the Modal Event nodes in the
 * node group and all of its nested groups. The set is gathered here so that it can be used before
 * the node group is evaluated, to register a modal keymap for the tool's operator type and to
 * decide which events are worth processing at all.
 */

#include <string>

#include "BLI_vector.hh"

namespace blender {
struct bNodeTree;
}

namespace blender::nodes {

/** One event defined by one or more Modal Event nodes. */
struct ModalEvent {
  std::string name;
  std::string description;

  friend bool operator==(const ModalEvent &a, const ModalEvent &b) = default;
};

struct ModalEvents {
  /**
   * Deduplicated and sorted by name, so that comparing two sets only depends on their
   * content and not on the order the node groups happened to be visited in.
   */
  Vector<ModalEvent> events;

  /** Add all events from the other set that aren't in this one yet. */
  void merge(const ModalEvents &other);
  /** Restore the invariants described above. Called once after all events have been added. */
  void deduplicate_and_sort();

  friend bool operator==(const ModalEvents &a, const ModalEvents &b) = default;
};

/**
 * Gather the events defined by the node group and all of its nested groups. This does not rely on
 * any cached data, so it also works for node groups that are not part of the current file yet.
 */
ModalEvents gather_modal_events_recursive(const bNodeTree &ntree);
/**
 * Same as above, but assumes that the events are already cached on the referenced node groups.
 */
ModalEvents gather_modal_events_with_cache(const bNodeTree &ntree);

}  // namespace blender::nodes
