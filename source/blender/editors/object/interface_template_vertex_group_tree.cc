/* SPDX-FileCopyrightText: 2025 Blender Foundation
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edinterface
 */

#include <fmt/format.h>

#include "BKE_context.hh"
#include "BKE_deform.hh"
#include "BKE_object.hh"
#include "BKE_object_deform.h"

#include "BLI_listbase.h"
#include "BLI_string.h"
#include "BLT_translation.hh"

#include "UI_interface.hh"
#include "UI_interface_layout.hh"
#include "UI_tree_view.hh"

#include "RNA_access.hh"
#include "RNA_prototypes.hh"

#include "DEG_depsgraph.hh"

#include "DNA_object_types.h"

#include "WM_api.hh"
#include "WM_types.hh"

#include "ED_object.hh"
#include "ED_undo.hh"

#include "object_intern.hh"

namespace blender::ed::object::vertexgroup {

/* -------------------------------------------------------------------- */
/** \name Helpers
 * \{ */

/** Returns the first collection that contains dg, or nullptr. */
static bDeformGroupCollection *defgroup_find_collection(const Object *object,
                                                         const bDeformGroup *dg)
{
  for (bDeformGroupCollection *col = static_cast<bDeformGroupCollection *>(
           object->defgroup_collections.first);
       col;
       col = col->next)
  {
    if (collection_find_member(col, dg)) {
      return col;
    }
  }
  return nullptr;
}

/** Remove dg membership from every collection it belongs to.
 * Defined in object_vgroup.cc, declared in object_intern.hh. */

/** \} */

/* -------------------------------------------------------------------- */
/** \name Data carrier
 * \{ */

struct VertexGroupData {
  Object *object;
  bDeformGroup *dg;
  int index; /* 0-based index in defbase */
  /** nullptr if the group is not in any collection. */
  bDeformGroupCollection *collection;
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Tree View
 * \{ */

class VertexGroupTreeView : public ui::AbstractTreeView {
 protected:
  Object &object_;

 public:
  VertexGroupTreeView(Object &ob) : object_(ob) {}

  void build_tree() override;
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Drag Controller
 * \{ */

class VertexGroupDragController : public ui::AbstractViewItemDragController {
 private:
  VertexGroupData drag_vg_;

 public:
  VertexGroupDragController(VertexGroupTreeView &view, VertexGroupData drag_vg)
      : AbstractViewItemDragController(view), drag_vg_(drag_vg)
  {
  }

  std::optional<eWM_DragDataType> get_drag_type() const override
  {
    return WM_DRAG_VERTEX_GROUP;
  }

  void *create_drag_data() const override
  {
    /* Collect all selected groups from defbase. Groups always live in defbase
     * regardless of collection membership, so this is always correct. */
    const ListBaseT<bDeformGroup> *defbase = BKE_object_defgroup_list(drag_vg_.object);

    int selected_count = [&]() -> int {
      int count = 0;
      for (const bDeformGroup *dg = static_cast<const bDeformGroup *>(defbase->first); dg;
           dg = dg->next)
      {
        count += (dg->flag & DG_SEL) != 0;
      }
      return count;
    }();

    bDeformGroup **selected_groups = MEM_new_array_zeroed<bDeformGroup *>(selected_count + 1,
                                                                          "Selected Def Groups");
    selected_count = 0;
    for (bDeformGroup *dg = static_cast<bDeformGroup *>(defbase->first); dg; dg = dg->next) {
      if (dg->flag & DG_SEL) {
        selected_groups[selected_count] = dg;
        selected_count++;
      }
    }
    BLI_assert_msg(selected_groups[selected_count] == nullptr,
                   "Expected last element to be null (null-delimiter)");
    return selected_groups;
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Drop Target (reorder in defbase)
 * \{ */

class VertexGroupDropTarget : public ui::TreeViewItemDropTarget {
 private:
  bDeformGroup &drop_dg_;
  int drop_index_;
  /** Collection context of the drop target item, nullptr if in defbase view. */
  bDeformGroupCollection *drop_collection_;

 public:
  VertexGroupDropTarget(ui::AbstractTreeViewItem &item,
                        ui::DropBehavior behavior,
                        bDeformGroup &drop_dg,
                        int index,
                        bDeformGroupCollection *collection)
      : TreeViewItemDropTarget(item, behavior),
        drop_dg_(drop_dg),
        drop_index_(index),
        drop_collection_(collection)
  {
  }

  bool can_drop(const wmDrag &drag, const char ** /*r_disabled_hint*/) const override
  {
    if (drag.type != WM_DRAG_VERTEX_GROUP) {
      return false;
    }
    const bDeformGroup **drag_groups = static_cast<const bDeformGroup **>(drag.poin);
    return drag_groups && drag_groups[0];
  }

  std::string drop_tooltip(const ui::DragInfo &drag_info) const override
  {
    const StringRef drag_name = TIP_("Selected Groups");
    const StringRef drop_name = drop_dg_.name;

    switch (drag_info.drop_location) {
      case ui::DropLocation::Into:
        BLI_assert_unreachable();
        break;
      case ui::DropLocation::Before:
        return fmt::format(fmt::runtime(TIP_("Move {} above {}")), drag_name, drop_name);
      case ui::DropLocation::After:
        return fmt::format(fmt::runtime(TIP_("Move {} below {}")), drag_name, drop_name);
      default:
        BLI_assert_unreachable();
        break;
    }
    return "";
  }

  bool on_drop(bContext *C, const ui::DragInfo &drag_info) const override
  {
    const bDeformGroup **drag_groups = static_cast<const bDeformGroup **>(
        drag_info.drag_data.poin);

    Object *object = CTX_data_active_object(C);
    const ListBaseT<bDeformGroup> *defbase = BKE_object_defgroup_list(object);

    if (drop_collection_ != nullptr) {
      /* Dropping onto an item inside a collection —
       * each group belongs to at most one collection: evict from any previous
       * collection before adding to the new one. */
      for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
        bDeformGroup *dg = const_cast<bDeformGroup *>(drag_groups[i]);
        if (!collection_find_member(drop_collection_, dg)) {
          /* Remove from whichever collection currently owns this group. */
          defgroup_remove_from_all_collections(object, dg);
          bDeformGroupMember *m = MEM_new<bDeformGroupMember>("bDeformGroupMember");
          m->dg = dg;
          BLI_addtail(&drop_collection_->members, m);
        }
      }

      /* Reorder members relative to drop target. */
      bDeformGroupMember *drop_member = collection_find_member(drop_collection_, &drop_dg_);
      for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
        bDeformGroupMember *m = collection_find_member(drop_collection_, drag_groups[i]);
        if (!m) {
          continue;
        }
        BLI_remlink(&drop_collection_->members, m);
        switch (drag_info.drop_location) {
          case ui::DropLocation::Before:
            BLI_insertlinkbefore(&drop_collection_->members, drop_member, m);
            break;
          case ui::DropLocation::After:
            BLI_insertlinkafter(&drop_collection_->members, drop_member, m);
            break;
          default:
            BLI_addtail(&drop_collection_->members, m);
            break;
        }
      }
    }
    else {
      /* Dropping onto a defbase item — remove from any collection first
       * so the group becomes a plain ungrouped item, then reorder. */
      for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
        defgroup_remove_from_all_collections(object,
                                             const_cast<bDeformGroup *>(drag_groups[i]));
      }

      const int first_drag_index = BLI_findindex(defbase, drag_groups[0]);
      int drop_index = BLI_findindex(defbase, &drop_dg_);

      switch (drag_info.drop_location) {
        case ui::DropLocation::Into:
          BLI_assert_unreachable();
          break;
        case ui::DropLocation::Before:
          drop_index -= int(first_drag_index < drop_index);
          break;
        case ui::DropLocation::After:
          drop_index += int(first_drag_index > drop_index);
          break;
      }

      for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
        const int drag_index = BLI_findindex(defbase, drag_groups[i]);
        if (drag_index == -1) {
          continue;
        }
        if (i > 0) {
          drop_index += int(drag_index > drop_index);
        }
        object_defgroup_move(object, drag_index, drop_index);
      }
    }

    DEG_id_tag_update(&object->id, ID_RECALC_GEOMETRY);
    WM_event_add_notifier(C, NC_GEOM | ND_VERTEX_GROUP, object);
    ED_undo_push(C, "Drop Vertex Group");

    return true;
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Collection Drop Target
 * \{ */

class VertexGroupCollectionDropTarget : public ui::TreeViewItemDropTarget {
 private:
  bDeformGroupCollection *collection_;

 public:
  VertexGroupCollectionDropTarget(ui::AbstractTreeViewItem &item,
                                  ui::DropBehavior behavior,
                                  bDeformGroupCollection *collection)
      : TreeViewItemDropTarget(item, behavior), collection_(collection)
  {
  }

  bool can_drop(const wmDrag &drag, const char ** /*r_disabled_hint*/) const override
  {
    if (drag.type == WM_DRAG_VERTEX_GROUP_COLLECTION) {
      const bDeformGroupCollection *drag_col = *static_cast<bDeformGroupCollection **>(
          drag.poin);
      return drag_col != collection_;
    }
    if (drag.type == WM_DRAG_VERTEX_GROUP) {
      const bDeformGroup **drag_groups = static_cast<const bDeformGroup **>(drag.poin);
      return drag_groups && drag_groups[0];
    }
    return false;
  }

  std::string drop_tooltip(const ui::DragInfo &drag_info) const override
  {
    const StringRef drop_name = collection_->name;

    if (drag_info.drag_data.type == WM_DRAG_VERTEX_GROUP_COLLECTION) {
      const StringRef drag_name =
          (*static_cast<bDeformGroupCollection **>(drag_info.drag_data.poin))->name;
      switch (drag_info.drop_location) {
        case ui::DropLocation::Before:
          return fmt::format(fmt::runtime(TIP_("Move {} above {}")), drag_name, drop_name);
        case ui::DropLocation::After:
          return fmt::format(fmt::runtime(TIP_("Move {} below {}")), drag_name, drop_name);
        default:
          break;
      }
      return "";
    }

    const StringRef drag_name = TIP_("Selected Groups");
    switch (drag_info.drop_location) {
      case ui::DropLocation::Into:
        return fmt::format(fmt::runtime(TIP_("Add {} to {}")), drag_name, drop_name);
      case ui::DropLocation::Before:
        return fmt::format(fmt::runtime(TIP_("Move {} above {}")), drag_name, drop_name);
      case ui::DropLocation::After:
        return fmt::format(fmt::runtime(TIP_("Move {} below {}")), drag_name, drop_name);
      default:
        break;
    }
    return "";
  }

  bool on_drop(bContext *C, const ui::DragInfo &drag_info) const override
  {
    Object *object = CTX_data_active_object(C);

    if (drag_info.drag_data.type == WM_DRAG_VERTEX_GROUP_COLLECTION) {
      /* Reorder collections. */
      bDeformGroupCollection *drag_col = *static_cast<bDeformGroupCollection **>(
          drag_info.drag_data.poin);
      BLI_remlink(&object->defgroup_collections, drag_col);
      switch (drag_info.drop_location) {
        case ui::DropLocation::Before:
          BLI_insertlinkbefore(&object->defgroup_collections, collection_, drag_col);
          break;
        case ui::DropLocation::After:
          BLI_insertlinkafter(&object->defgroup_collections, collection_, drag_col);
          break;
        default:
          BLI_addtail(&object->defgroup_collections, drag_col);
          break;
      }
      DEG_id_tag_update(&object->id, ID_RECALC_GEOMETRY);
      WM_event_add_notifier(C, NC_GEOM | ND_VERTEX_GROUP, object);
      ED_undo_push(C, "Reorder Vertex Group Collection");
      return true;
    }

    /* Vertex group drag — add membership, never move out of defbase. */
    bDeformGroup **drag_groups = static_cast<bDeformGroup **>(drag_info.drag_data.poin);

    switch (drag_info.drop_location) {
      case ui::DropLocation::Into: {
        /* Add groups to this collection. Each group belongs to at most one
         * collection — evict from any previous collection first. */
        for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
          bDeformGroup *dg = drag_groups[i];
          if (!collection_find_member(collection_, dg)) {
            defgroup_remove_from_all_collections(object, dg);
            bDeformGroupMember *m = MEM_new<bDeformGroupMember>("bDeformGroupMember");
            m->dg = dg;
            BLI_addtail(&collection_->members, m);
          }
        }
        break;
      }
      case ui::DropLocation::Before:
      case ui::DropLocation::After: {
        /* Remove membership from this collection — group stays in defbase. */
        for (int8_t i = 0; drag_groups[i] != nullptr; i++) {
          bDeformGroupMember *m = collection_find_member(collection_, drag_groups[i]);
          if (m) {
            BLI_remlink(&collection_->members, m);
            MEM_delete(m);
          }
        }
        break;
      }
      default:
        break;
    }

    DEG_id_tag_update(&object->id, ID_RECALC_GEOMETRY);
    WM_event_add_notifier(C, NC_GEOM | ND_VERTEX_GROUP, object);
    ED_undo_push(C, "Drop Vertex Group into Collection");

    return true;
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Tree View Item
 * \{ */

class VertexGroupItem : public ui::AbstractTreeViewItem {
 private:
  VertexGroupData vg_data_;

 public:
  VertexGroupItem(Object *object,
                  bDeformGroup *dg,
                  int index,
                  bDeformGroupCollection *collection)
  {
    label_ = dg->name;
    vg_data_.object = object;
    vg_data_.dg = dg;
    vg_data_.index = index;
    vg_data_.collection = collection;
  }

  void build_row(ui::Layout &row) override
  {
    uiItemL_ex(&row, this->label_, ICON_GROUP_VERTEX, false, false);

    ui::Layout &sub = row.row(true);
    sub.use_property_decorate_set(false);

    PointerRNA vg_ptr = RNA_pointer_create_discrete(
        &vg_data_.object->id, RNA_VertexGroup, vg_data_.dg);

    const BIFIconID icon = (vg_data_.dg->flag & DG_LOCK_WEIGHT) ? ICON_LOCKED : ICON_UNLOCKED;
    sub.prop(&vg_ptr, "lock_weight", ui::ITEM_R_ICON_ONLY, std::nullopt, icon);
  }

  std::optional<bool> should_be_active() const override
  {
    return BKE_object_defgroup_active_index_get(vg_data_.object) == vg_data_.index + 1;
  }

  void on_activate(bContext &C) override
  {
    BKE_object_defgroup_active_index_set(vg_data_.object, vg_data_.index + 1);
    DEG_id_tag_update(&vg_data_.object->id, ID_RECALC_GEOMETRY);
    WM_event_add_notifier(&C, NC_GEOM | ND_VERTEX_GROUP, vg_data_.object);
    ED_undo_push(&C, "Set Active Vertex Group");
  }

  std::optional<bool> should_be_selected() const override
  {
    return vg_data_.dg->flag & DG_SEL;
  }

  void set_selected(const bool select) override
  {
    AbstractViewItem::set_selected(select);
    SET_FLAG_FROM_TEST(vg_data_.dg->flag, select, DG_SEL);
  }

  bool supports_renaming() const override
  {
    return true;
  }

  bool rename(const bContext &C, StringRefNull new_name) override
  {
    PointerRNA vg_ptr = RNA_pointer_create_discrete(
        &vg_data_.object->id, RNA_VertexGroup, vg_data_.dg);
    RNA_string_set(&vg_ptr, "name", new_name.c_str());
    label_ = vg_data_.dg->name;
    ED_undo_push(const_cast<bContext *>(&C), "Rename Vertex Group");
    return true;
  }

  StringRef get_rename_string() const override
  {
    return label_;
  }

  void delete_item(bContext *C) override
  {
    /* Delegate to the operator so OPTYPE_UNDO handles undo registration
     * cleanly — pushing undo from inside the tree view redraw cycle crashes. */
    WM_operator_name_call(
        C, "OBJECT_OT_vertex_group_remove", wm::OpCallContext::InvokeDefault, nullptr, nullptr);
  }

  void build_context_menu(bContext &C, ui::Layout &layout) const override
  {
    MenuType *mt = WM_menutype_find("OBJECT_MT_vertex_group_tree_context_menu", true);
    if (!mt) {
      return;
    }
    ui::menutype_draw(&C, mt, &layout);
  }

  std::unique_ptr<ui::AbstractViewItemDragController> create_drag_controller() const override
  {
    return std::make_unique<VertexGroupDragController>(
        static_cast<VertexGroupTreeView &>(get_tree_view()), vg_data_);
  }

  std::unique_ptr<ui::TreeViewItemDropTarget> create_drop_target() override
  {
    return std::make_unique<VertexGroupDropTarget>(
        *this, ui::DropBehavior::Reorder, *vg_data_.dg, vg_data_.index, vg_data_.collection);
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Collection Drag Controller
 * \{ */

class VertexGroupCollectionDragController : public ui::AbstractViewItemDragController {
 private:
  bDeformGroupCollection *collection_;

 public:
  VertexGroupCollectionDragController(VertexGroupTreeView &view,
                                      bDeformGroupCollection *collection)
      : AbstractViewItemDragController(view), collection_(collection)
  {
  }

  std::optional<eWM_DragDataType> get_drag_type() const override
  {
    return WM_DRAG_VERTEX_GROUP_COLLECTION;
  }

  void *create_drag_data() const override
  {
    /* Allocate wrapper so framework can free poin without corrupting DNA pointer. */
    bDeformGroupCollection **wrapper = MEM_new_array_zeroed<bDeformGroupCollection *>(
        1, "Drag Collection");
    wrapper[0] = collection_;
    return wrapper;
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Collection Item
 * \{ */

class VertexGroupCollectionItem : public ui::AbstractTreeViewItem {
 private:
  Object *object_;
  bDeformGroupCollection *collection_;

 public:
  VertexGroupCollectionItem(Object *object, bDeformGroupCollection *collection)
  {
    /* Prefix to avoid label_ matching conflicts with VertexGroupItem of same name. */
    label_ = std::string("collection:") + collection->name;
    object_ = object;
    collection_ = collection;
  }

  std::optional<bool> should_be_collapsed() const override
  {
    return (collection_->flag & DG_COLLECTION_EXPANDED) == 0;
  }

  bool set_collapsed(const bool collapsed) override
  {
    if (!AbstractTreeViewItem::set_collapsed(collapsed)) {
      return false;
    }
    SET_FLAG_FROM_TEST(collection_->flag, !collapsed, DG_COLLECTION_EXPANDED);
    return true;
  }

  void build_row(ui::Layout &row) override
  {
    const int label_icon = (collection_->flag & DG_COLLECTION_LAYER) ?
                               ICON_MOD_MASK :
                               ICON_GROUP;
    uiItemL_ex(&row, collection_->name, label_icon, false, false);

    ui::Layout &sub = row.row(true);
    sub.use_property_decorate_set(false);

    /* Layer toggle — calls operator on the active collection (DG_COLLECTION_ACTIVE).
     * TODO: replace with sub.prop(&col_ptr, "is_layer", ...) once
     *       RNA_VertexGroupCollection is registered. */
    {
      const int icon = (collection_->flag & DG_COLLECTION_LAYER) ?
                           ICON_MOD_MASK :
                           ICON_GROUP_VERTEX;
      ui::Layout &layer_sub = sub.row(true);
      layer_sub.active_set(bool(collection_->flag & DG_COLLECTION_ACTIVE));
      layer_sub.op("OBJECT_OT_vertex_group_collection_layer_toggle", "", icon);
    }



    /* Lock button — mirrors bone collection is_visible pattern.
     * Grey out when collection has no members (nothing to lock),
     * analogous to BONE_COLLECTION_ANCESTORS_VISIBLE greying out
     * when a parent collection is hidden.
     * TODO: replace representative_dg prop with sub.prop(&col_ptr, "lock_weight", ...)
     *       once RNA_VertexGroupCollection is registered. The RNA setter should
     *       iterate col->members and write DG_LOCK_WEIGHT to each m->dg. */
    if (collection_->representative_dg != nullptr) {
      ui::Layout &lock_sub = sub.row(true);
      lock_sub.active_set(!BLI_listbase_is_empty(&collection_->members));
      const int icon = (collection_->representative_dg->flag & DG_LOCK_WEIGHT) ?
                           ICON_LOCKED :
                           ICON_UNLOCKED;
      PointerRNA dg_ptr = RNA_pointer_create_discrete(
          &object_->id, RNA_VertexGroup, collection_->representative_dg);
      lock_sub.prop(&dg_ptr, "lock_weight", ui::ITEM_R_ICON_ONLY, "", icon);
    }
  }

  std::optional<bool> should_be_active() const override
  {
    return collection_->flag & DG_COLLECTION_ACTIVE;
  }

  void on_activate(bContext &C) override
  {
    /* Clear active from all other collections. */
    for (bDeformGroupCollection *col = static_cast<bDeformGroupCollection *>(
             object_->defgroup_collections.first);
         col;
         col = col->next)
    {
      col->flag &= ~DG_COLLECTION_ACTIVE;
    }
    collection_->flag |= DG_COLLECTION_ACTIVE;

    /* Activate the collection's representative group for weight painting.
     * It lives in defbase so actdef, MDeformVert indices and all operators
     * work without any special casing. */
    if (collection_->representative_dg != nullptr) {
      const ListBaseT<bDeformGroup> *defbase = BKE_object_defgroup_list(object_);
      const int index = BLI_findindex(defbase, collection_->representative_dg);
      if (index != -1) {
        BKE_object_defgroup_active_index_set(object_, index + 1);
      }
    }

    DEG_id_tag_update(&object_->id, ID_RECALC_GEOMETRY);
    WM_event_add_notifier(&C, NC_GEOM | ND_VERTEX_GROUP, object_);
    ED_undo_push(&C, "Set Active Vertex Group Collection");
  }

  std::optional<bool> should_be_selected() const override
  {
    return collection_->flag & DG_COLLECTION_SEL;
  }

  void set_selected(const bool select) override
  {
    AbstractViewItem::set_selected(select);
    SET_FLAG_FROM_TEST(collection_->flag, select, DG_COLLECTION_SEL);
  }

  bool supports_renaming() const override
  {
    return true;
  }

  bool rename(const bContext &C, StringRefNull new_name) override
  {
    /* TODO: replace with RNA_string_set once RNA_VertexGroupCollection is registered. */
    BLI_strncpy(collection_->name, new_name.c_str(), sizeof(collection_->name));
    label_ = std::string("collection:") + collection_->name;

    /* Keep the representative group name in sync so weight paint shows the
     * updated collection name. Use RNA so uniqueness is enforced. */
    if (collection_->representative_dg != nullptr) {
      char rep_name[MAX_NAME];
      SNPRINTF(rep_name, "col:%s", collection_->name);
      PointerRNA dg_ptr = RNA_pointer_create_discrete(
          &object_->id, RNA_VertexGroup, collection_->representative_dg);
      RNA_string_set(&dg_ptr, "name", rep_name);
    }

    ED_undo_push(const_cast<bContext *>(&C), "Rename Vertex Group Collection");
    return true;
  }

  StringRef get_rename_string() const override
  {
    return collection_->name;
  }

  void delete_item(bContext *C) override
  {
    /* Remove all member links — member groups stay in defbase untouched. */
    bDeformGroupMember *m = static_cast<bDeformGroupMember *>(collection_->members.first);
    while (m) {
      bDeformGroupMember *next = m->next;
      BLI_remlink(&collection_->members, m);
      MEM_delete(m);
      m = next;
    }

    /* Remove the representative group from defbase — this also removes its
     * weight data from all mesh vertices via BKE_object_defgroup_remove. */
    if (collection_->representative_dg != nullptr) {
      BKE_object_defgroup_remove(object_, collection_->representative_dg);
      collection_->representative_dg = nullptr;
    }

    BLI_remlink(&object_->defgroup_collections, collection_);
    MEM_delete(collection_);
    collection_ = nullptr;

    DEG_id_tag_update(&object_->id, ID_RECALC_GEOMETRY);
    WM_event_add_notifier(C, NC_GEOM | ND_VERTEX_GROUP, object_);
    ED_undo_grouped_push(C, "Delete Vertex Group Collection");
  }

  std::unique_ptr<ui::AbstractViewItemDragController> create_drag_controller() const override
  {
    return std::make_unique<VertexGroupCollectionDragController>(
        static_cast<VertexGroupTreeView &>(get_tree_view()), collection_);
  }

  std::unique_ptr<ui::TreeViewItemDropTarget> create_drop_target() override
  {
    return std::make_unique<VertexGroupCollectionDropTarget>(
        *this, ui::DropBehavior::ReorderAndInsert, collection_);
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name build_tree
 * \{ */

void VertexGroupTreeView::build_tree()
{
  const ListBaseT<bDeformGroup> *defbase = BKE_object_defgroup_list(&object_);

  /* Pass 1: collections — members are references into defbase. */
  for (bDeformGroupCollection *col = static_cast<bDeformGroupCollection *>(
           object_.defgroup_collections.first);
       col;
       col = col->next)
  {
    auto &col_item = this->add_tree_item<VertexGroupCollectionItem>(&object_, col);

    for (bDeformGroupMember *m = static_cast<bDeformGroupMember *>(col->members.first); m;
         m = m->next)
    {
      if (!m->dg) {
        continue; /* Stale reference — group was deleted. */
      }
      /* Verify the group still exists in defbase. */
      if (BLI_findindex(defbase, m->dg) == -1) {
        continue;
      }
      const int index = BLI_findindex(defbase, m->dg);
      col_item.add_tree_item<VertexGroupItem>(&object_, m->dg, index, col);
    }
  }

  /* Pass 2: all groups in defbase that are NOT already shown as collection members,
   * and excluding representative groups (identified by "col:" prefix).
   * Groups that belong to a collection appear only under their collection in pass 1. */
  int index = 0;
  for (bDeformGroup *dg = static_cast<bDeformGroup *>(defbase->first); dg;
       dg = dg->next, index++)
  {
    if (STRPREFIX(dg->name, "col:")) {
      continue;
    }
    if (defgroup_find_collection(&object_, dg) != nullptr) {
      /* Already shown as a member under its collection. */
      continue;
    }
    this->add_tree_item<VertexGroupItem>(&object_, dg, index, nullptr);
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Public Template Entry Point
 * \{ */

void template_tree(ui::Layout *layout, bContext *C)
{
  Object *ob = CTX_data_active_object(C);
  if (ob == nullptr) {
    return;
  }

  if (!OB_TYPE_SUPPORT_VGROUP(ob->type)) {
    return;
  }

  ui::Block *block = layout->block();

  ui::AbstractTreeView *tree_view = block_add_view(
      *block,
      "Vertex Group Tree View",
      std::make_unique<VertexGroupTreeView>(*ob));

  tree_view->set_context_menu_title("Vertex Group");
  tree_view->set_default_rows(4);
  tree_view->allow_multiselect_items();

  ui::TreeViewBuilder::build_tree_view(*C, *tree_view, *layout);
}

/** \} */

}  // namespace blender::ed::object::vertexgroup