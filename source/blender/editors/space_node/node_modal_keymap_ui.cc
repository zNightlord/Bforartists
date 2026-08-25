/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spnode
 *
 * Editor for the default modal keymap of a node tool. The events themselves are defined by the
 * Modal Event nodes in the node group and its nested groups, this only assigns a default key
 * binding to each of them. The bindings users actually work with are configured in the keymap
 * editor, this is what a node group ships with.
 */

#include "MEM_guardedalloc.h"

#include "BKE_context.hh"
#include "BKE_node_runtime.hh"
#include "BKE_screen.hh"

#include "BLI_listbase.hh"
#include "BLI_string.hh"
#include "BLI_string_utf8.hh"

#include "BLT_translation.hh"

#include "DNA_screen_types.h"
#include "DNA_space_types.h"

#include "NOD_geometry_nodes_modal_events.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"
#include "RNA_prototypes.hh"

#include "UI_interface_layout.hh"
#include "UI_resources.hh"
#include "UI_tree_view.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "ED_screen.hh"

#include "node_intern.hh"

namespace blender::ed::space_node {

static GeometryNodeModalKeymapItem *find_keymap_item(bNodeTree &tree, const StringRef event_name)
{
  GeometryNodeAssetTraits *traits = tree.geometry_node_asset_traits;
  if (!traits) {
    return nullptr;
  }
  for (const int i : IndexRange(traits->modal_keymap_default_num)) {
    GeometryNodeModalKeymapItem &item = traits->modal_keymap_default[i];
    if (item.event_name && event_name == item.event_name) {
      return &item;
    }
  }
  return nullptr;
}

/* -------------------------------------------------------------------- */
/** \name Add and Remove Operators
 * \{ */

static wmOperatorStatus modal_keymap_item_add_exec(bContext *C, wmOperator *op)
{
  SpaceNode &snode = *CTX_wm_space_node(C);
  bNodeTree &tree = *snode.edittree;
  GeometryNodeAssetTraits &traits = *tree.geometry_node_asset_traits;
  const std::string event_name = RNA_string_get(op->ptr, "event_name");
  if (event_name.empty()) {
    return OPERATOR_CANCELLED;
  }
  if (find_keymap_item(tree, event_name)) {
    return OPERATOR_CANCELLED;
  }

  const int old_num = traits.modal_keymap_default_num;
  auto *new_items = MEM_new_array_zeroed<GeometryNodeModalKeymapItem>(old_num + 1, __func__);
  if (traits.modal_keymap_default) {
    std::copy_n(traits.modal_keymap_default, old_num, new_items);
    MEM_delete(traits.modal_keymap_default);
  }
  new_items[old_num].event_name = BLI_strdupn(event_name.c_str(), event_name.size());
  new_items[old_num].val = KM_PRESS;
  traits.modal_keymap_default = new_items;
  traits.modal_keymap_default_num = old_num + 1;
  traits.modal_keymap_active_index = old_num;

  WM_event_add_notifier(C, NC_NODE | ND_DISPLAY, nullptr);
  return OPERATOR_FINISHED;
}

static bool modal_keymap_edit_poll(bContext *C)
{
  const SpaceNode *snode = CTX_wm_space_node(C);
  if (!snode || !snode->edittree) {
    return false;
  }
  if (snode->edittree->type != NTREE_GEOMETRY) {
    return false;
  }
  return snode->edittree->geometry_node_asset_traits != nullptr;
}

void NODE_OT_modal_keymap_item_add(wmOperatorType *ot)
{
  ot->name = "Add Modal Keymap Item";
  ot->description = "Add a default key binding for one of the node tool's events";
  ot->idname = "NODE_OT_modal_keymap_item_add";

  ot->exec = modal_keymap_item_add_exec;
  ot->poll = modal_keymap_edit_poll;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  RNA_def_string(
      ot->srna, "event_name", nullptr, MAX_NAME, "Name", "Name of the event to bind a key to");
}

static wmOperatorStatus modal_keymap_item_remove_exec(bContext *C, wmOperator *op)
{
  SpaceNode &snode = *CTX_wm_space_node(C);
  bNodeTree &tree = *snode.edittree;
  GeometryNodeAssetTraits &traits = *tree.geometry_node_asset_traits;
  const std::string event_name = RNA_string_get(op->ptr, "event_name");
  GeometryNodeModalKeymapItem *item = find_keymap_item(tree, event_name);
  if (!item) {
    return OPERATOR_CANCELLED;
  }
  const int index = int(item - traits.modal_keymap_default);
  MEM_SAFE_DELETE(item->event_name);

  const int new_num = traits.modal_keymap_default_num - 1;
  GeometryNodeModalKeymapItem *new_items = nullptr;
  if (new_num > 0) {
    new_items = MEM_new_array_zeroed<GeometryNodeModalKeymapItem>(new_num, __func__);
    std::copy_n(traits.modal_keymap_default, index, new_items);
    std::copy_n(traits.modal_keymap_default + index + 1, new_num - index, new_items + index);
  }
  MEM_delete(traits.modal_keymap_default);
  traits.modal_keymap_default = new_items;
  traits.modal_keymap_default_num = new_num;
  traits.modal_keymap_active_index = std::min(traits.modal_keymap_active_index,
                                              std::max(new_num - 1, 0));

  WM_event_add_notifier(C, NC_NODE | ND_DISPLAY, nullptr);
  return OPERATOR_FINISHED;
}

void NODE_OT_modal_keymap_item_remove(wmOperatorType *ot)
{
  ot->name = "Remove Modal Keymap Item";
  ot->description = "Remove the default key binding of one of the node tool's events";
  ot->idname = "NODE_OT_modal_keymap_item_remove";

  ot->exec = modal_keymap_item_remove_exec;
  ot->poll = modal_keymap_edit_poll;

  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  RNA_def_string(ot->srna, "event_name", nullptr, MAX_NAME, "Name", "Name of the event to unbind");
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Tree View
 * \{ */

class ModalEventItem : public ui::AbstractTreeViewItem {
  bNodeTree &tree_;
  /**
   * False for a binding that is left over from a Modal Event node that has been removed or
   * renamed. Those are still listed so that they can be seen and removed, but the node is the
   * source of truth.
   */
  bool event_exists_ = true;

 public:
  ModalEventItem(bNodeTree &tree, const StringRef event_name, const bool event_exists)
      : tree_(tree), event_exists_(event_exists)
  {
    label_ = event_name;
  }

  void build_row(ui::Layout &row) override
  {
    ui::Layout &name_row = row.row(true);
    name_row.active_set(event_exists_);
    name_row.label(label_, ICON_NONE);

    GeometryNodeModalKeymapItem *item = find_keymap_item(tree_, label_);
    ui::Layout &sub = row.row(true);
    sub.alignment_set(ui::LayoutAlign::Right);
    if (!event_exists_) {
      ui::Layout &note = sub.row(true);
      note.active_set(false);
      note.label(IFACE_("No Modal Event node"), ICON_INFO);
    }
    if (item) {
      PointerRNA item_ptr = RNA_pointer_create_discrete(
          &tree_.id, RNA_GeometryNodeModalKeymapItem, item);
      ui::Layout &type_row = sub.row(true);
      type_row.active_set(event_exists_);
      type_row.prop(&item_ptr, "type", UI_ITEM_NONE, "", ICON_NONE);
      PointerRNA op_ptr = sub.op("node.modal_keymap_item_remove", "", ICON_X);
      RNA_string_set(&op_ptr, "event_name", label_.c_str());
    }
    else {
      sub.label(IFACE_("Unbound"), ICON_NONE);
      PointerRNA op_ptr = sub.op("node.modal_keymap_item_add", "", ICON_ADD);
      RNA_string_set(&op_ptr, "event_name", label_.c_str());
    }
  }

  void on_activate(bContext & /*C*/) override
  {
    if (GeometryNodeAssetTraits *traits = tree_.geometry_node_asset_traits) {
      if (const GeometryNodeModalKeymapItem *item = find_keymap_item(tree_, label_)) {
        traits->modal_keymap_active_index = int(item - traits->modal_keymap_default);
      }
    }
  }

  std::optional<bool> should_be_active() const override
  {
    const GeometryNodeAssetTraits &traits = *tree_.geometry_node_asset_traits;
    if (traits.modal_keymap_default_num == 0) {
      return false;
    }
    const GeometryNodeModalKeymapItem *item = find_keymap_item(tree_, label_);
    return item && int(item - traits.modal_keymap_default) == traits.modal_keymap_active_index;
  }
};

class ModalEventsView : public ui::AbstractTreeView {
  bNodeTree &tree_;

 public:
  ModalEventsView(bNodeTree &tree) : tree_(tree) {}

  void build_tree() override
  {
    this->is_flat_ = true;
    Set<StringRef> event_names;
    if (tree_.runtime->modal_events) {
      for (const nodes::ModalEvent &event : tree_.runtime->modal_events->events) {
        event_names.add(event.name);
        this->add_tree_item<ModalEventItem>(tree_, event.name, true);
      }
    }
    /* List bindings whose Modal Event node is gone last, so that they can still be removed. The
     * default keymap being out of date is not an error, those bindings are just ignored. */
    if (const GeometryNodeAssetTraits *traits = tree_.geometry_node_asset_traits) {
      for (const int i : IndexRange(traits->modal_keymap_default_num)) {
        const GeometryNodeModalKeymapItem &item = traits->modal_keymap_default[i];
        if (item.event_name && !event_names.contains(item.event_name)) {
          this->add_tree_item<ModalEventItem>(tree_, item.event_name, false);
        }
      }
    }
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Panel
 * \{ */

static bool modal_keymap_panel_poll(const bContext *C, PanelType * /*pt*/)
{
  const SpaceNode *snode = CTX_wm_space_node(C);
  if (!snode || !snode->edittree) {
    return false;
  }
  const bNodeTree &tree = *snode->edittree;
  if (tree.type != NTREE_GEOMETRY) {
    return false;
  }
  if (!tree.geometry_node_asset_traits) {
    return false;
  }
  if ((tree.geometry_node_asset_traits->flag & GEO_NODE_ASSET_TOOL) == 0) {
    return false;
  }
  /* Show the tab when the node group defines events at any nesting level, but also when it only
   * has bindings left over from removed Modal Event nodes, so that those can still be cleaned up.
   */
  if (tree.runtime->modal_events && !tree.runtime->modal_events->events.is_empty()) {
    return true;
  }
  return tree.geometry_node_asset_traits->modal_keymap_default_num > 0;
}

static void modal_keymap_panel_draw(const bContext *C, Panel *panel)
{
  SpaceNode &snode = *CTX_wm_space_node(C);
  bNodeTree &tree = *snode.edittree;
  ui::Layout &layout = *panel->layout;

  ui::Block *block = layout.block();
  ui::AbstractTreeView *tree_view = block_add_view(
      *block, "Node Tool Modal Events", std::make_unique<ModalEventsView>(tree));
  tree_view->set_default_rows(4);
  ui::TreeViewBuilder::build_tree_view(*C, *tree_view, layout);

  GeometryNodeAssetTraits &traits = *tree.geometry_node_asset_traits;
  if (!traits.modal_keymap_default ||
      !IndexRange(traits.modal_keymap_default_num).contains(traits.modal_keymap_active_index))
  {
    return;
  }
  GeometryNodeModalKeymapItem &item =
      traits.modal_keymap_default[traits.modal_keymap_active_index];
  PointerRNA item_ptr = RNA_pointer_create_discrete(
      &tree.id, RNA_GeometryNodeModalKeymapItem, &item);

  layout.use_property_split_set(true);
  layout.use_property_decorate_set(false);
  layout.prop(&item_ptr, "type", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  layout.prop(&item_ptr, "value", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  layout.prop(&item_ptr, "key_modifier", UI_ITEM_NONE, std::nullopt, ICON_NONE);

  ui::Layout &col = layout.column(true);
  col.prop(&item_ptr, "shift", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  col.prop(&item_ptr, "ctrl", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  col.prop(&item_ptr, "alt", UI_ITEM_NONE, std::nullopt, ICON_NONE);
  col.prop(&item_ptr, "oskey", UI_ITEM_NONE, std::nullopt, ICON_NONE);
}

void node_modal_keymap_panel_register(ARegionType *art)
{
  PanelType *pt = MEM_new_zeroed<PanelType>(__func__);
  STRNCPY_UTF8(pt->idname, "NODE_PT_modal_keymap");
  STRNCPY_UTF8(pt->label, N_("Events"));
  STRNCPY_UTF8(pt->category, "Keymap");
  STRNCPY_UTF8(pt->translation_context, BLT_I18NCONTEXT_DEFAULT_BPYRNA);
  pt->draw = modal_keymap_panel_draw;
  pt->poll = modal_keymap_panel_poll;
  BLI_addtail(&art->paneltypes, pt);
}

/** \} */

}  // namespace blender::ed::space_node
