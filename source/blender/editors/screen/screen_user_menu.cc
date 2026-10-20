/* SPDX-FileCopyrightText: 2009 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spview3d
 */

#include <cfloat>
#include <cstring>

#include "MEM_guardedalloc.h"

#include "AS_asset_library.hh"
#include "AS_asset_representation.hh"

#include "BLT_translation.hh"

#include "BKE_blender_user_menu.hh"
#include "BKE_context.hh"
#include "BKE_idprop.hh"
#include "BKE_screen.hh"

#include "BLI_listbase.hh"
#include "BLI_string_utf8.hh"
#include "BLI_vector.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "ED_asset_list.hh"
#include "ED_asset_menu_utils.hh"
#include "ED_screen.hh"

#include "UI_interface.hh"
#include "UI_interface_c.hh"
#include "UI_interface_layout.hh"
#include "UI_resources.hh"
#include "UI_tree_view.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"
#include "RNA_path.hh"
#include "RNA_prototypes.hh"

namespace blender {

using namespace ui;

/* -------------------------------------------------------------------- */
/** \name Internal Utilities
 * \{ */

static const char *screen_menu_context_string(const bContext *C, const SpaceLink *sl)
{
  if (sl->spacetype == SPACE_NODE) {
    const SpaceNode *snode = reinterpret_cast<const SpaceNode *>(sl);
    return snode->tree_idname;
  }
  return CTX_data_mode_string(C);
}

/* BFA: Quick Favorites - shared visibility filter (space type + mode). */
static bool user_menu_item_is_visible(const bContext *C, const bUserMenuItem &umi)
{
  const SpaceLink *sl = CTX_wm_space_data(C);
  if (sl == nullptr) {
    return false;
  }
  const int space_type = sl->spacetype;
  const int mode = int(CTX_data_mode_enum(C));

  const bool mode_ok = (umi.mode == 0) || (uchar(umi.mode) == 255) || (int(umi.mode) == mode);
  if (!mode_ok) {
    return false;
  }
  if (umi.space_type == space_type) {
    return true;
  }
  /* Topbar and Toolbar items are also shown in the 3D View. */
  return (space_type == SPACE_VIEW3D) &&
         ELEM(umi.space_type, SPACE_TOPBAR, SPACE_TOOLBAR);
}

/* BFA: Quick Favorites - resolve display name + icon for the editor. */
static StringRefNull user_menu_item_label(const bContext *C,
                                          bUserMenuItem &umi,
                                          char (&r_buf)[256],
                                          int &r_icon)
{
  r_icon = ICON_NONE;
  const StringRefNull custom_name = umi.ui_name;

  if (umi.type == USER_MENU_TYPE_OPERATOR) {
    bUserMenuItem_Op *umi_op = reinterpret_cast<bUserMenuItem_Op *>(&umi);
    wmOperatorType *ot = WM_operatortype_find(umi_op->op_idname, false);
    if (ot == nullptr) {
      SNPRINTF_UTF8(r_buf, RPT_("Missing: %s"), umi_op->op_idname);
      return r_buf;
    }
    r_icon = umi_op->icon;
    if (!custom_name.is_empty()) {
      return custom_name;
    }
    /* `WM_operatortype_name` returns a temporary std::string, so copy it into the buffer. */
    SNPRINTF_UTF8(r_buf, "%s", WM_operatortype_name(ot, nullptr).c_str());
    return r_buf;
  }

  if (umi.type == USER_MENU_TYPE_MENU) {
    bUserMenuItem_Menu *umi_mt = reinterpret_cast<bUserMenuItem_Menu *>(&umi);
    MenuType *mt = WM_menutype_find(umi_mt->mt_idname, false);
    if (mt == nullptr) {
      SNPRINTF_UTF8(r_buf, RPT_("Missing: %s"), umi_mt->mt_idname);
      return r_buf;
    }
    r_icon = umi_mt->icon;
    if (!custom_name.is_empty()) {
      return custom_name;
    }
    return CTX_IFACE_(mt->translation_context, mt->label);
  }

  if (umi.type == USER_MENU_TYPE_PROP) {
    bUserMenuItem_Prop *umi_pr = reinterpret_cast<bUserMenuItem_Prop *>(&umi);
    r_icon = umi_pr->icon;
    if (!custom_name.is_empty()) {
      return custom_name;
    }

    /* Work on a local copy of the path, never modify the stored one. */
    char path[sizeof(umi_pr->context_data_path)];
    STRNCPY_UTF8(path, umi_pr->context_data_path);
    const char *tail = nullptr;
    if (char *dot = strchr(path, '.')) {
      *dot = '\0';
      tail = dot + 1;
    }

    PointerRNA ptr = CTX_data_pointer_get(C, path);
    if (ptr.type == nullptr) {
      PointerRNA ctx_ptr = RNA_pointer_create_discrete(nullptr, RNA_Context, (void *)C);
      if (!RNA_path_resolve_full(&ctx_ptr, path, &ptr, nullptr, nullptr)) {
        ptr.type = nullptr;
      }
    }
    if (ptr.type != nullptr) {
      PointerRNA prop_ptr = ptr;
      if (tail == nullptr || RNA_path_resolve_full(&ptr, tail, &prop_ptr, nullptr, nullptr)) {
        if (PropertyRNA *prop = RNA_struct_find_property(&prop_ptr, umi_pr->prop_id)) {
          const char *prop_name = RNA_property_ui_name(prop);
          if (prop_name && prop_name[0]) {
            return prop_name; /* Static RNA string, safe to point at. */
          }
        }
      }
    }
    SNPRINTF_UTF8(r_buf, "%s.%s", umi_pr->context_data_path, umi_pr->prop_id);
    return r_buf;
  }

  if (umi.type == USER_MENU_TYPE_SEP) {
    return "--- Separator ---";
  }
  return custom_name;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Menu Type
 * \{ */

bUserMenu **ED_screen_user_menus_find(const bContext *C, uint *r_len)
{
  SpaceLink *sl = CTX_wm_space_data(C);

  if (sl == nullptr) {
    *r_len = 0;
    return nullptr;
  }

  /* Return the global menu for all operations - filtering happens during display. */
  uint array_len = 1;
  bUserMenu **um_array = MEM_new_array_zeroed<bUserMenu *>(array_len, __func__);
  um_array[0] = BKE_blender_user_menu_find(&U.user_menus, SPACE_VIEW3D, "");

  *r_len = array_len;
  return um_array;
}

bUserMenu *ED_screen_user_menu_ensure(bContext *C)
{
  SpaceLink *sl = CTX_wm_space_data(C);
  const char *context = screen_menu_context_string(C, sl);
  
  /* Store all items in a global menu but track their original context. */
  bUserMenu *global_menu = BKE_blender_user_menu_ensure(&U.user_menus, SPACE_VIEW3D, "");
  
  /* Find existing menu to migrate items from context-specific storage. */
  bUserMenu *old_menu = BKE_blender_user_menu_find(&U.user_menus, sl->spacetype, context);
  if (old_menu && old_menu != global_menu) {
    /* Move items from old menu to global menu. */
    while (bUserMenuItem *umi = reinterpret_cast<bUserMenuItem *>(BLI_pophead(&old_menu->items))) {
      /* Store original context in the item. */
      umi->space_type = sl->spacetype;
      STRNCPY_UTF8(umi->context, context);
      BLI_addtail(&global_menu->items, umi);
    }
  }
  
  return global_menu;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Quick Favorites Editor Dialog (tree view)
 * \{ */

/* Currently selected item in the editor. Must be declared before the classes below.
 * Only compared by pointer, and reset on remove/invoke. */
static bUserMenuItem *g_user_menu_editor_selected = nullptr;

static bUserMenu *user_menu_global_find()
{
  return BKE_blender_user_menu_find(&U.user_menus, SPACE_VIEW3D, "");
}

/* BFA: Quick Favorites - every item of the global menu. The editors list everything;
 * context/mode filtering only applies to the menu itself. */
static blender::Vector<bUserMenuItem *> user_menu_all_items()
{
  blender::Vector<bUserMenuItem *> result;
  if (bUserMenu *um = user_menu_global_find()) {
    for (bUserMenuItem &umi : um->items) {
      result.append(&umi);
    }
  }
  return result;
}

static blender::Vector<bUserMenuItem *> user_menu_siblings(const bUserMenuItem *umi)
{
  blender::Vector<bUserMenuItem *> result;
  if (umi == nullptr) {
    return result;
  }
  for (bUserMenuItem *other : user_menu_all_items()) {
    if (other->space_type == umi->space_type) {
      result.append(other);
    }
  }
  return result;
}

static const char *user_menu_space_name(const int space_type)
{
  switch (space_type) {
    case SPACE_TOPBAR:
      return IFACE_("Topbar");
    case SPACE_VIEW3D:
      return IFACE_("3D View");
    case SPACE_PROPERTIES:
      return IFACE_("Properties");
    case SPACE_TOOLBAR:
      return IFACE_("Toolbar");
    default:
      return IFACE_("Other");
  }
}

/* Hover text: where the item was added, plus what it runs. */
static std::string user_menu_item_tooltip(const bUserMenuItem &umi)
{
  std::string tip = std::string(IFACE_("Editor")) + ": " + user_menu_space_name(umi.space_type);
  if (umi.context[0] != '\0') {
    tip += "\n";
    tip += IFACE_("Context");
    tip += ": ";
    tip += umi.context;
  }

  if (umi.type == USER_MENU_TYPE_OPERATOR) {
    const bUserMenuItem_Op *umi_op = reinterpret_cast<const bUserMenuItem_Op *>(&umi);
    tip += "\n";
    tip += IFACE_("Operator");
    tip += ": ";
    tip += umi_op->op_idname;
  }
  else if (umi.type == USER_MENU_TYPE_MENU) {
    const bUserMenuItem_Menu *umi_mt = reinterpret_cast<const bUserMenuItem_Menu *>(&umi);
    tip += "\n";
    tip += IFACE_("Menu");
    tip += ": ";
    tip += umi_mt->mt_idname;
  }
  else if (umi.type == USER_MENU_TYPE_PROP) {
    const bUserMenuItem_Prop *umi_pr = reinterpret_cast<const bUserMenuItem_Prop *>(&umi);
    tip += "\n";
    tip += IFACE_("Property");
    tip += ": ";
    tip += umi_pr->context_data_path;
    tip += ".";
    tip += umi_pr->prop_id;
  }
  return tip;
}

/* Tooltip callback. `arg` is the bUserMenuItem, resolved at hover time (no allocation to free).
 * Check the signature against an existing button_func_tooltip_set user in your fork. */
static std::string user_menu_tooltip_fn(bContext * /*C*/, void *arg, StringRef /*tip*/)
{
  return user_menu_item_tooltip(*static_cast<const bUserMenuItem *>(arg));
}

/* Make sure the selection points at an item in the list (defaults to the first). */
static void user_menu_editor_validate_selection(const blender::Vector<bUserMenuItem *> &items)
{
  if (items.is_empty()) {
    g_user_menu_editor_selected = nullptr;
    return;
  }
  if (!items.contains(g_user_menu_editor_selected)) {
    g_user_menu_editor_selected = items.first();
  }
}

static void user_menu_editor_refresh(bContext *C)
{
  if (ARegion *region = CTX_wm_region(C)) {
    ED_region_tag_redraw(region);
  }
}

struct UserMenuDragData {
  bUserMenuItem *umi;
};

/* Take the return type from the base class so the override always matches. */
using UserMenuDragType = decltype(std::declval<const ui::AbstractViewItemDragController &>()
                                      .get_drag_type());

class QuickFavoritesDragController : public ui::AbstractViewItemDragController {
  bUserMenuItem *umi_;

 public:
  QuickFavoritesDragController(ui::AbstractView &view, bUserMenuItem *umi)
      : ui::AbstractViewItemDragController(view), umi_(umi)
  {
  }

  UserMenuDragType get_drag_type() const override
  {
    return WM_DRAG_USER_MENU_ITEM;
  }

  void *create_drag_data() const override
  {
    /* Same ownership as the asset catalog controller: WM frees this when the drag ends. */
    UserMenuDragData *data = MEM_new_zeroed<UserMenuDragData>(__func__);
    data->umi = umi_;
    return data;
  }
};

class QuickFavoritesDropTarget : public ui::TreeViewItemDropTarget {
  bUserMenuItem *target_;

 public:
  QuickFavoritesDropTarget(ui::AbstractTreeViewItem &item, bUserMenuItem *target)
      : ui::TreeViewItemDropTarget(item, ui::DropBehavior::Reorder), target_(target)
  {
  }

  static bUserMenuItem *dragged_item(const wmDrag &drag)
  {
    return static_cast<const UserMenuDragData *>(drag.poin)->umi;
  }

  /* Reordering is only allowed inside a group (same editor space). */
  bool can_drop(const wmDrag &drag, const char ** /*r_disabled_hint*/) const override
  {
    if (drag.type != WM_DRAG_USER_MENU_ITEM) {
      return false;
    }
    const bUserMenuItem *dragged = dragged_item(drag);
    return dragged != target_;
  }

  std::string drop_tooltip(const ui::DragInfo & /*drag_info*/) const override
  {
    return TIP_("Move Quick Favorites item here");
  }

  bool on_drop(bContext * /*C*/, const ui::DragInfo &drag_info) const override
  {
    bUserMenu *um = user_menu_global_find();
    bUserMenuItem *dragged = dragged_item(drag_info.drag_data);
    if (um == nullptr || dragged == nullptr || dragged == target_) {
      return false;
    }
    BLI_remlink(&um->items, dragged);
    if (drag_info.drop_location == ui::DropLocation::Before) {
      BLI_insertlinkbefore(&um->items, target_, dragged);
    }
    else {
      BLI_insertlinkafter(&um->items, target_, dragged);
    }
    g_user_menu_editor_selected = dragged;
    U.runtime.is_dirty = true;
    return true;
  }
};

class QuickFavoritesTreeItem : public ui::BasicTreeViewItem {
  bUserMenuItem *umi_;

 public:
  QuickFavoritesTreeItem(StringRef label, int icon, bUserMenuItem *umi)
      : ui::BasicTreeViewItem(label, BIFIconID(icon)), umi_(umi)
  {
  }

  void on_activate(bContext &C) override
  {
    g_user_menu_editor_selected = umi_;
    ED_region_tag_redraw(CTX_wm_region(&C));
  }

  std::optional<bool> should_be_active() const override
  {
    return umi_ == g_user_menu_editor_selected;
  }

  bool supports_collapsing() const override
  {
    return false;
  }

  bool matches(const ui::AbstractViewItem &other) const override
  {
    const QuickFavoritesTreeItem *other_item = dynamic_cast<const QuickFavoritesTreeItem *>(
        &other);
    return other_item && (other_item->umi_ == umi_);
  }

  std::unique_ptr<ui::AbstractViewItemDragController> create_drag_controller() const override
  {
    return std::make_unique<QuickFavoritesDragController>(get_tree_view(), umi_);
  }

  std::unique_ptr<ui::TreeViewItemDropTarget> create_drop_target() override
  {
    return std::make_unique<QuickFavoritesDropTarget>(*this, umi_);
  }
};

/* Display-only header row for one editor space. Not backed by a bUserMenuItem. */
class QuickFavoritesGroupItem : public ui::BasicTreeViewItem {
  int space_type_;

 public:
  QuickFavoritesGroupItem(StringRef label, const int space_type)
      : ui::BasicTreeViewItem(label, BIFIconID(ICON_NONE)), space_type_(space_type)
  {
  }

  /* Headers are never highlighted or selected. */
  std::optional<bool> should_be_active() const override
  {
    return false;
  }

  bool matches(const ui::AbstractViewItem &other) const override
  {
    const QuickFavoritesGroupItem *other_group = dynamic_cast<const QuickFavoritesGroupItem *>(
        &other);
    return other_group && (other_group->space_type_ == space_type_);
  }
};

class QuickFavoritesTreeView : public ui::AbstractTreeView {
  const bContext *C_;

 public:
  explicit QuickFavoritesTreeView(const bContext *C) : C_(C) {}

  void build_tree() override
  {
    const blender::Vector<bUserMenuItem *> items = user_menu_all_items();

    /* Groups in order of first appearance. */
    blender::Vector<int> space_order;
    for (const bUserMenuItem *umi : items) {
      if (!space_order.contains(int(umi->space_type))) {
        space_order.append(int(umi->space_type));
      }
    }

    for (const int space_type : space_order) {
      QuickFavoritesGroupItem &group = add_tree_item<QuickFavoritesGroupItem>(
          user_menu_space_name(space_type), space_type);
      group.uncollapse_by_default(); /* Start expanded. */

      for (bUserMenuItem *umi : items) {
        if (umi->space_type != space_type) {
          continue;
        }
        char name_buf[256];
        int icon = ICON_NONE;
        const StringRefNull name = user_menu_item_label(C_, *umi, name_buf, icon);
        group.add_tree_item<QuickFavoritesTreeItem>(name, icon, umi);
      }
    }
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Menu Item
 * \{ */

bUserMenuItem_Op *ED_screen_user_menu_item_find_operator(ListBaseT<bUserMenuItem> *lb,
                                                         const wmOperatorType *ot,
                                                         IDProperty *prop,
                                                         const char *op_prop_enum,
                                                         wm::OpCallContext opcontext)
{
  for (bUserMenuItem &umi : *lb) {
    if (umi.type == USER_MENU_TYPE_OPERATOR) {
      bUserMenuItem_Op *umi_op = reinterpret_cast<bUserMenuItem_Op *>(&umi);
      const bool is_strict = prop && umi_op->prop;
      const bool ok_idprop = IDP_EqualsProperties_ex(prop, umi_op->prop, is_strict);
      const bool ok_prop_enum = (umi_op->op_prop_enum[0] != '\0') ?
                                    STREQ(umi_op->op_prop_enum, op_prop_enum) :
                                    true;
      if (STREQ(ot->idname, umi_op->op_idname) &&
          (opcontext == wm::OpCallContext(umi_op->opcontext)) && ok_idprop && ok_prop_enum)
      {
        return umi_op;
      }
    }
  }
  return nullptr;
}

bUserMenuItem_Menu *ED_screen_user_menu_item_find_menu(ListBaseT<bUserMenuItem> *lb,
                                                       const MenuType *mt)
{
  for (bUserMenuItem &umi : *lb) {
    if (umi.type == USER_MENU_TYPE_MENU) {
      bUserMenuItem_Menu *umi_mt = reinterpret_cast<bUserMenuItem_Menu *>(&umi);
      if (STREQ(mt->idname, umi_mt->mt_idname)) {
        return umi_mt;
      }
    }
  }
  return nullptr;
}

bUserMenuItem_Prop *ED_screen_user_menu_item_find_prop(ListBaseT<bUserMenuItem> *lb,
                                                       const char *context_data_path,
                                                       const char *prop_id,
                                                       int prop_index)
{
  for (bUserMenuItem &umi : *lb) {
    if (umi.type == USER_MENU_TYPE_PROP) {
      bUserMenuItem_Prop *umi_pr = reinterpret_cast<bUserMenuItem_Prop *>(&umi);
      if (STREQ(context_data_path, umi_pr->context_data_path) && STREQ(prop_id, umi_pr->prop_id) &&
          (prop_index == umi_pr->prop_index))
      {
        return umi_pr;
      }
    }
  }
  return nullptr;
}

/* BFA: Quick Favorites - Add operator item with context-aware filtering */
void ED_screen_user_menu_item_add_operator(const bContext *C,
                                           ListBaseT<bUserMenuItem> *lb,
                                           const char *ui_name,
                                           const wmOperatorType *ot,
                                           const IDProperty *prop,
                                           const char *op_prop_enum,
                                           wm::OpCallContext opcontext,
                                           int icon)
{
  bUserMenuItem_Op *umi_op = reinterpret_cast<bUserMenuItem_Op *>(
      BKE_blender_user_menu_item_add(lb, USER_MENU_TYPE_OPERATOR));
  umi_op->opcontext = int8_t(opcontext);
  if (!STREQ(ui_name, ot->name)) {
    STRNCPY_UTF8(umi_op->item.ui_name, ui_name);
  }
  STRNCPY_UTF8(umi_op->op_idname, ot->idname);
  STRNCPY_UTF8(umi_op->op_prop_enum, op_prop_enum);
  umi_op->prop = prop ? IDP_CopyProperty(prop) : nullptr;
  umi_op->icon = icon;
  
  /* Set context information from the current context */
  SpaceLink *sl = CTX_wm_space_data(C);
  if (sl) {
    /* Always try to get the context from the active area, not from space_data */
    ScrArea *area = CTX_wm_area(C);
    if (area && area->spacetype != SPACE_TOPBAR) {
      /* Get the actual active editor context */
      sl = area->spacedata.first_as<SpaceLink>();
    }
    
    umi_op->item.space_type = sl->spacetype;
    const char *context = screen_menu_context_string(C, sl);
    STRNCPY_UTF8(umi_op->item.context, context);
    
    /* Store the current mode for mode-specific filtering */
    umi_op->item.mode = CTX_data_mode_enum(C);
    
    } else {
    /* No space data available, keep default values */
  }
}

/* BFA: Quick Favorites - Add menu item with context-aware filtering */
void ED_screen_user_menu_item_add_menu(const bContext *C,
                                       ListBaseT<bUserMenuItem> *lb,
                                       const char *ui_name,
                                       const MenuType *mt,
                                       int icon)
{
  bUserMenuItem_Menu *umi_mt = reinterpret_cast<bUserMenuItem_Menu *>(
      BKE_blender_user_menu_item_add(lb, USER_MENU_TYPE_MENU));
  if (!STREQ(ui_name, mt->label)) {
    STRNCPY_UTF8(umi_mt->item.ui_name, ui_name);
  }
  STRNCPY_UTF8(umi_mt->mt_idname, mt->idname);
  umi_mt->icon = icon;
  
  /* Set context information from the current context */
  SpaceLink *sl = CTX_wm_space_data(C);
  if (sl) {
    /* Always try to get the context from the active area, not from space_data */
    ScrArea *area = CTX_wm_area(C);
    if (area && area->spacetype != SPACE_TOPBAR) {
      /* Get the actual active editor context */
      sl = area->spacedata.first_as<SpaceLink>();
    }
    
    umi_mt->item.space_type = sl->spacetype;
    const char *context = screen_menu_context_string(C, sl);
    STRNCPY_UTF8(umi_mt->item.context, context);
    
    /* Store the current mode for mode-specific filtering */
    umi_mt->item.mode = CTX_data_mode_enum(C);
  }
}

/* BFA: Quick Favorites - Add property item with context-aware filtering */
void ED_screen_user_menu_item_add_prop(const bContext *C,
                                       ListBaseT<bUserMenuItem> *lb,
                                       const char *ui_name,
                                       const char *context_data_path,
                                       const char *prop_id,
                                       int prop_index,
                                       int icon)
{
  bUserMenuItem_Prop *umi_pr = reinterpret_cast<bUserMenuItem_Prop *>(
      BKE_blender_user_menu_item_add(lb, USER_MENU_TYPE_PROP));
  STRNCPY_UTF8(umi_pr->item.ui_name, ui_name);
  STRNCPY_UTF8(umi_pr->context_data_path, context_data_path);
  STRNCPY_UTF8(umi_pr->prop_id, prop_id);
  umi_pr->prop_index = prop_index;
  umi_pr->icon = icon;
  
  /* Set context information from the current context */
  SpaceLink *sl = CTX_wm_space_data(C);
  if (sl) {
    /* Always try to get the context from the active area, not from space_data */
    ScrArea *area = CTX_wm_area(C);
    if (area && area->spacetype != SPACE_TOPBAR) {
      /* Get the actual active editor context */
      sl = area->spacedata.first_as<SpaceLink>();
    }
    
    umi_pr->item.space_type = sl->spacetype;
    const char *context = screen_menu_context_string(C, sl);
    STRNCPY_UTF8(umi_pr->item.context, context);
    
    /* Store the current mode for mode-specific filtering */
    umi_pr->item.mode = CTX_data_mode_enum(C);
  }
}

void ED_screen_user_menu_item_remove(ListBaseT<bUserMenuItem> *lb, bUserMenuItem *umi)
{
  BLI_remlink(lb, umi);
  BKE_blender_user_menu_item_free(umi);
}

void ED_screen_user_menu_item_move_up(ListBaseT<bUserMenuItem> *lb, bUserMenuItem *umi)
{
  if (umi->prev) {
    BLI_listbase_swaplinks(lb, umi, umi->prev);
  }
}

void ED_screen_user_menu_item_move_down(ListBaseT<bUserMenuItem> *lb, bUserMenuItem *umi)
{
  if (umi->next) {
    BLI_listbase_swaplinks(lb, umi, umi->next);
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Quick Favorites Editor Dialog
 * \{
 */

/* BFA*/
static void user_menu_editor_draw(const bContext *C, ui::Layout &layout)
{
  const blender::Vector<bUserMenuItem *> items = user_menu_all_items();
  if (items.is_empty()) {
    ui::Layout *col = &layout.column(false);
    col->label(RPT_("No Quick Favorites items found."), ICON_INFO);
    col->label(RPT_("Right-click on buttons in menus to add them here."), ICON_NONE);
    return;
  }

  user_menu_editor_validate_selection(items);
  const int selected_index = int(items.first_index_of_try(g_user_menu_editor_selected));
  bUserMenuItem *selected = selected_index >= 0 ? items[selected_index] : nullptr;

  /* Flat-list move enablement (any item, ignoring groups). */
  const bool can_move_up = selected_index > 0;
  const bool can_move_down = selected_index >= 0 && selected_index < int(items.size()) - 1;
  const bool can_remove = selected != nullptr;

  ui::Layout &main_row = layout.row(false);

  /* Left: tree view (groups preserved for display). */
  ui::Layout &list_col = main_row.column(false);
  ui::Block *block = list_col.block();
  ui::AbstractTreeView *tree_view = block_add_view(
      *block, "Quick Favorites Tree View", std::make_unique<QuickFavoritesTreeView>(C));
  ui::TreeViewBuilder::build_tree_view(*C, *tree_view, list_col, true);

  /* Right: control buttons. */
  ui::Layout &control_col = main_row.column(true);
  control_col.alignment_set(ui::LayoutAlign::Right);

  ui::Layout &up_row = control_col.row(true);
  up_row.enabled_set(can_move_up);
  PointerRNA ptr_up = up_row.op("SCREEN_OT_user_menu_item_move",
                                "",
                                ICON_TRIA_UP,
                                wm::OpCallContext::InvokeDefault,
                                ITEM_R_ICON_ONLY);
  RNA_enum_set(&ptr_up, "direction", 0);

  ui::Layout &down_row = control_col.row(true);
  down_row.enabled_set(can_move_down);
  PointerRNA ptr_down = down_row.op("SCREEN_OT_user_menu_item_move",
                                    "",
                                    ICON_TRIA_DOWN,
                                    wm::OpCallContext::InvokeDefault,
                                    ITEM_R_ICON_ONLY);
  RNA_enum_set(&ptr_down, "direction", 1);

  control_col.separator();

  ui::Layout &remove_row = control_col.row(true);
  remove_row.enabled_set(can_remove);
  remove_row.op("SCREEN_OT_user_menu_item_remove",
                "",
                ICON_X,
                wm::OpCallContext::InvokeDefault,
                ITEM_R_ICON_ONLY);
}

/* BFA: Quick Favorites - Context-aware edit menu with mode-specific filtering */
static void user_menu_edit_draw(bContext *C, wmOperator *op)
{
  user_menu_editor_draw(C, *op->layout);
}

static wmOperatorStatus user_menu_edit_exec(bContext * /*C*/, wmOperator * /*op*/)
{
  return OPERATOR_FINISHED;
}

static wmOperatorStatus user_menu_edit_invoke(bContext *C, wmOperator *op, const wmEvent * /*event*/)
{
  const blender::Vector<bUserMenuItem *> items = user_menu_all_items();
  g_user_menu_editor_selected = items.is_empty() ? nullptr : items.first();
  return WM_operator_props_dialog_popup(C, op, 300, std::nullopt);
}

static void SCREEN_OT_user_menu_edit(wmOperatorType *ot)
{
  ot->name = "Edit Quick Favorites";
  ot->idname = "SCREEN_OT_user_menu_edit";
  ot->description = "Open editor to reorder and manage Quick Favorites items";

  ot->invoke = user_menu_edit_invoke;
  ot->exec = user_menu_edit_exec;
  ot->poll = ED_operator_screenactive;

  ot->flag = OPTYPE_UNDO;

  RNA_def_int(ot->srna,
              "selected_index",
              0,
              0,
              INT_MAX,
              "Selected Index",
              "Index of the currently selected item",
              0,
              INT_MAX);

  ot->ui = user_menu_edit_draw;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Remove Item Operator
 * \{
 */

/* BFA: Quick Favorites - Remove item with context-aware filtering */
static wmOperatorStatus user_menu_item_remove_exec(bContext *C, wmOperator * /*op*/)
{
  bUserMenu *um = user_menu_global_find();
  if (um == nullptr) {
    return OPERATOR_CANCELLED;
  }

  const blender::Vector<bUserMenuItem *> items = user_menu_all_items();
  const int index = int(items.first_index_of_try(g_user_menu_editor_selected));
  if (index < 0) {
    return OPERATOR_CANCELLED;
  }

  bUserMenuItem *target = items[index];

  /* Prefer next item, else previous (flat list). */
  bUserMenuItem *new_selection = nullptr;
  if (index + 1 < int(items.size())) {
    new_selection = items[index + 1];
  }
  else if (index > 0) {
    new_selection = items[index - 1];
  }

  U.runtime.is_dirty = true;
  ED_screen_user_menu_item_remove(&um->items, target);
  g_user_menu_editor_selected = new_selection;
  user_menu_editor_refresh(C);
  return OPERATOR_FINISHED;
}

static void SCREEN_OT_user_menu_item_remove(wmOperatorType *ot)
{
  ot->name = "Remove Quick Favorites Item";
  ot->idname = "SCREEN_OT_user_menu_item_remove";
  ot->description = "Remove an item from Quick Favorites";

  ot->exec = user_menu_item_remove_exec;
  ot->poll = ED_operator_screenactive;

  ot->flag = OPTYPE_UNDO;

  RNA_def_int(ot->srna,
              "index",
              0,
              0,
              INT_MAX,
              "Index",
              "Index of the item to remove",
              0,
              INT_MAX);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Menu Definition
 * \{ */

static bool all_loading_finished()
{
  AssetLibraryReference all_library_ref = asset_system::all_library_reference();
  return ed::asset::list::is_loaded(&all_library_ref);
}

/**
 * When adding operators that reference an asset (see
 * `ed::asset::operator_asset_reference_props_is_set`) make sure the asset libraries are loaded.
 */
static void handle_operator_asset_reference_props(const bContext &C,
                                                  bUserMenuItem_Op &umi_op,
                                                  ui::Layout &row,
                                                  wmOperatorType *ot,
                                                  int &r_icon,
                                                  bool &r_add_operator)
{
  if (!umi_op.prop) {
    return;
  }
  PointerRNA opptr = WM_operator_properties_create_ptr(ot);
  opptr.data = bke::idprop::create_group("wmOperatorProperties").release();
  IDP_CopyPropertyContent(opptr.data_as<IDProperty>(), umi_op.prop);
  if (ed::asset::operator_asset_reference_props_is_set(opptr)) {
    const bool loading_finished = all_loading_finished();
    if (!loading_finished) {
      row.label(IFACE_("Loading Asset Libraries"), ICON_STATUS_INFO);
      r_add_operator = false;
    }
    else {
      /* Set `check_context_asset` to false because we're setting the context pointer after getting
       * the asset from the operator properties. */
      const asset_system::AssetRepresentation *asset =
          ed::asset::operator_asset_reference_props_get_asset_from_all_library(
              C, opptr, CTX_wm_reports(&C));
      if (asset) {
        if (asset->is_online_only()) {
          r_icon = ICON_INTERNET;
        }
        PointerRNA asset_ptr = RNA_pointer_create_discrete(
            nullptr,
            RNA_AssetRepresentation,
            const_cast<asset_system::AssetRepresentation *>(asset));
        row.context_ptr_set("asset", &asset_ptr);
      }
      else {
        r_add_operator = false;
      }
    }
  }
  WM_operator_properties_free(&opptr);
}

static void screen_user_menu_draw(const bContext *C, Menu *menu)
{
  /* Enable when we have the ability to edit menus. */
  const bool show_missing = false;
  char label[512];

  uint um_array_len;
  bUserMenu **um_array = ED_screen_user_menus_find(C, &um_array_len);
  bool is_empty = true;
  
  /* Get current context for filtering */
  SpaceLink *sl = CTX_wm_space_data(C);
  char current_space_type = sl->spacetype;
  
  for (int um_index = 0; um_index < um_array_len; um_index++) {
    bUserMenu *um = um_array[um_index];
    if (um == nullptr) {
      continue;
    }
    for (bUserMenuItem &umi : um->items) {
      /* Filter items based on current context - show only items that work in this context */
      if (!user_menu_item_is_visible(C, umi)) {
        continue;
      }
      bool show_item = false;
      
      /* Show items from current space type AND mode */
      if (umi.space_type == current_space_type) {
        /* Check if modes match, or if item has no specific mode (0 or uninitialized) */
        if (umi.mode == 0 || umi.mode == CTX_data_mode_enum(C) || umi.mode == 255) {
          show_item = true;
        }
      }
      /* Allow Topbar and Toolbar items to show in 3D View context */
      else if (current_space_type == SPACE_VIEW3D && 
               (umi.space_type == SPACE_TOPBAR || umi.space_type == SPACE_TOOLBAR)) {
        /* Check if modes match, or if item has no specific mode (0 or uninitialized) */
        if (umi.mode == 0 || umi.mode == CTX_data_mode_enum(C) || umi.mode == 255) {
          show_item = true;
        }
      }
      
      if (!show_item) {
        continue;  /* Skip this item */
      }
      std::optional<StringRefNull> ui_name = umi.ui_name[0] ?
                                                 std::make_optional<StringRefNull>(umi.ui_name) :
                                                 std::nullopt;
      if (umi.type == USER_MENU_TYPE_OPERATOR) {
        bUserMenuItem_Op *umi_op = reinterpret_cast<bUserMenuItem_Op *>(&umi);
        if (wmOperatorType *ot = WM_operatortype_find(umi_op->op_idname, false)) {
          if (ui_name) {
            ui_name = CTX_IFACE_(ot->translation_context, ui_name->c_str());
          }
          if (umi_op->op_prop_enum[0] == '\0') {
            ui::Layout &row = menu->layout->row(true);
            int icon = umi_op->icon; /*BFA - start with user-set icon*/
            bool add_operator = true;
            handle_operator_asset_reference_props(*C, *umi_op, row, ot, icon, add_operator);
            if (add_operator) {
              PointerRNA ptr = row.op(
                  ot, ui_name, icon, wm::OpCallContext(umi_op->opcontext), UI_ITEM_NONE);
              if (umi_op->prop) {
                IDP_CopyPropertyContent(ptr.data_as<IDProperty>(), umi_op->prop);
              }
            }
          }
          else {
            /* umi_op->prop could be used to set other properties but it's currently unsupported.
             */
            menu->layout->op_menu_enum(C, ot, umi_op->op_prop_enum, ui_name, umi_op->icon);
          }
          is_empty = false;
        }
        else {
          if (show_missing) {
            SNPRINTF_UTF8(label, RPT_("Missing: %s"), umi_op->op_idname);
            menu->layout->label(label, ICON_NONE);
          }
        }
      }
      else if (umi.type == USER_MENU_TYPE_MENU) {
        bUserMenuItem_Menu *umi_mt = reinterpret_cast<bUserMenuItem_Menu *>(&umi);
        MenuType *mt = WM_menutype_find(umi_mt->mt_idname, false);
        if (mt != nullptr) {
          menu->layout->menu(mt, ui_name, umi_mt->icon);
          is_empty = false;
        }
        else {
          if (show_missing) {
            SNPRINTF_UTF8(label, RPT_("Missing: %s"), umi_mt->mt_idname);
            menu->layout->label(label, ICON_NONE);
          }
        }
      }
      else if (umi.type == USER_MENU_TYPE_PROP) {
        bUserMenuItem_Prop *umi_pr = reinterpret_cast<bUserMenuItem_Prop *>(&umi);

        char *data_path = strchr(umi_pr->context_data_path, '.');
        if (data_path) {
          *data_path = '\0';
        }
        PointerRNA ptr = CTX_data_pointer_get(C, umi_pr->context_data_path);
        if (!ptr.has_type()) {
          PointerRNA ctx_ptr = RNA_pointer_create_discrete(nullptr, RNA_Context, (void *)C);
          if (!RNA_path_resolve_full(&ctx_ptr, umi_pr->context_data_path, &ptr, nullptr, nullptr))
          {
            ptr.type = nullptr;
          }
        }
        if (data_path) {
          *data_path = '.';
          data_path += 1;
        }

        bool ok = false;
        if (ptr.has_type()) {
          PropertyRNA *prop = nullptr;
          PointerRNA prop_ptr = ptr;
          if ((data_path == nullptr) ||
              RNA_path_resolve_full(&ptr, data_path, &prop_ptr, nullptr, nullptr))
          {
            prop = RNA_struct_find_property(&prop_ptr, umi_pr->prop_id);
            if (prop) {
              ok = true;
              menu->layout->prop(
                  &prop_ptr, prop, umi_pr->prop_index, 0, UI_ITEM_NONE, ui_name, umi_pr->icon);
              is_empty = false;
            }
          }
        }
        if (!ok) {
          if (show_missing) {
            SNPRINTF_UTF8(
                label, RPT_("Missing: %s.%s"), umi_pr->context_data_path, umi_pr->prop_id);
            menu->layout->label(label, ICON_NONE);
          }
        }
      }
      else if (umi.type == USER_MENU_TYPE_SEP) {
        menu->layout->separator();
      }
    }
  }
  if (um_array) {
    MEM_delete(um_array);
  }

  if (is_empty) {
    menu->layout->label(RPT_("No menu items found"), ICON_NONE);
    menu->layout->label(RPT_("Right click on buttons to add them to this menu"), ICON_NONE);
  }
  else {
    /* Add Edit button at the bottom of the menu. */
    menu->layout->separator();
    menu->layout->op("SCREEN_OT_user_menu_edit",
                     IFACE_("Edit Quick Favorites..."),
                     ICON_LASTOPERATOR,
                     wm::OpCallContext::InvokeDefault,
                     UI_ITEM_NONE);
  }
}

/* -------------------------------------------------------------------- */
/** \name Move Item Operator
 * \{
 */

/* BFA: Quick Favorites - Move item with context-aware filtering */
static wmOperatorStatus user_menu_item_move_exec(bContext *C, wmOperator *op)
{
  const int direction = RNA_enum_get(op->ptr, "direction");
  bUserMenu *um = user_menu_global_find();
  if (um == nullptr) {
    return OPERATOR_CANCELLED;
  }

  const blender::Vector<bUserMenuItem *> items = user_menu_all_items();
  const int index = int(items.first_index_of_try(g_user_menu_editor_selected));
  if (index < 0) {
    return OPERATOR_CANCELLED;
  }

  bUserMenuItem *target = items[index];
  const int neighbor = (direction == 0) ? index - 1 : index + 1;
  if (neighbor < 0 || neighbor >= int(items.size())) {
    return OPERATOR_CANCELLED;
  }

  /* Flat-list swap – groups are only visual. */
  BLI_listbase_swaplinks(&um->items, target, items[neighbor]);
  U.runtime.is_dirty = true;
  user_menu_editor_refresh(C);
  return OPERATOR_FINISHED;
}

static void SCREEN_OT_user_menu_item_move(wmOperatorType *ot)
{
  static const EnumPropertyItem direction_items[] = {
      {0, "UP", 0, "Up", "Move item up"},
      {1, "DOWN", 0, "Down", "Move item down"},
      {0, nullptr, 0, nullptr, nullptr},
  };

  ot->name = "Move Quick Favorites Item";
  ot->idname = "SCREEN_OT_user_menu_item_move";
  ot->description = "Move an item in the Quick Favorites menu up or down";

  ot->exec = user_menu_item_move_exec;
  ot->poll = ED_operator_screenactive;

  ot->flag = OPTYPE_UNDO;

  RNA_def_enum(ot->srna, "direction", direction_items, 0, "Direction", "Direction to move");
  RNA_def_int(ot->srna,
              "index",
              0,
              0,
              INT_MAX,
              "Index",
              "Index of the item to move",
              0,
              INT_MAX);
}

/** \} */

namespace ed::screen {

/* Preferences panel entry point, referenced directly by RNA (`layout.template_quick_favorites()`).
 * The parameter order must be (layout, C): that is the order RNA calls it in when
 * FUNC_USE_CONTEXT is set. */
void template_quick_favorites(ui::Layout *layout, bContext *C)
{
  user_menu_editor_draw(C, *layout);
}

}  // namespace ed::screen

void ED_screen_user_menu_register()
{
  MenuType *mt = MEM_new_zeroed<MenuType>(__func__);
  STRNCPY_UTF8(mt->idname, "SCREEN_MT_user_menu");
  STRNCPY_UTF8(mt->label, N_("Quick Favorites"));
  STRNCPY_UTF8(mt->translation_context, BLT_I18NCONTEXT_DEFAULT_BPYRNA);
  mt->draw = screen_user_menu_draw;
  mt->listener = ed::asset::list::asset_reading_region_listen_fn;
  WM_menutype_add(mt);

  WM_operatortype_append(SCREEN_OT_user_menu_item_move);
  WM_operatortype_append(SCREEN_OT_user_menu_item_remove);
  WM_operatortype_append(SCREEN_OT_user_menu_edit);
}

/** \} */

}  // namespace blender
