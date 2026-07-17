/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup shdnodes
 */

#include "NOD_mask_stack.hh"
#include "NOD_socket_items_ops.hh"
#include "NOD_socket_items_ui.hh"

#include "RNA_prototypes.hh"

#include "node_shader_util.hh"
#include "node_util.hh"

namespace blender {

namespace nodes::node_shader_mask_stack_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  const bNodeTree *tree = b.tree_or_null();
  const bNode *node = b.node_or_null();
  if (node) {
    const NodeShaderLayerStack &storage = layer_stack::storage(*node);
    const int last_index = storage.items_num - 1;
    for (const int i : IndexRange(storage.items_num)) {
      const NodeShaderLayerStackItem &item = storage.items[i];
      const UString name = item.name ? UString(item.name) : ""_ustr;
      const UString mask_id(mask_stack::ItemsAccessor::socket_identifier_for_item(item));
      auto &mask_decl = b.add_input<decl::Float>(name, mask_id)
                            .default_value(0.0f)
                            .min(0.0f)
                            .max(1.0f)
                            .subtype(PROP_FACTOR)
                            .structure_type(StructureType::Single);
      if (tree) {
        mask_decl.socket_name_ptr(&tree->id, *mask_stack::ItemsAccessor::item_srna, &item, "name");
      }
      /* Every mask layer exposes Opacity. For non-base layers it weights the
       * layer over the layers below; for the base (last) layer it blends
       * between no mask (fully unmasked, 1.0) and this layer's mask. */
      const UString opacity_id(
          mask_stack::ItemsAccessor::opacity_socket_identifier_for_item(item));
      b.add_input<decl::Float>("Opacity"_ustr, opacity_id)
          .default_value(1.0f)
          .min(0.0f)
          .max(1.0f)
          .subtype(PROP_FACTOR)
          .description(i != last_index ?
                           "Weight of this mask layer over the layers below" :
                           "Blend between no mask (fully unmasked) and this mask layer");
    }
  }
  /* Trailing extend-socket: dragging any link onto this hollow socket adds a
   * new mask layer at the end. The new item becomes the base. */
  b.add_input<decl::Extend>(""_ustr, "__extend__"_ustr)
      .custom_draw(socket_items::ui::draw_extend_socket_fn<mask_stack::ItemsAccessor>());
  b.add_output<decl::Float>("Result"_ustr).structure_type(StructureType::Single);
}

static void node_operators()
{
  socket_items::ops::make_common_operators<mask_stack::ItemsAccessor>();
}

/* Sidebar layout: the mask items list with add/remove/move, plus blend mode
 * and mute of the active mask. The base (last) mask has no blend mode. */
static void node_layout_ex(ui::Layout &layout, bContext *C, PointerRNA *node_ptr)
{
  bNodeTree &ntree = *reinterpret_cast<bNodeTree *>(node_ptr->owner_id);
  const bNode &node = *static_cast<const bNode *>(node_ptr->data);
  socket_items::ui::draw_items_list_with_operators<mask_stack::ItemsAccessor>(
      C, &layout, ntree, node);
  const NodeShaderLayerStack &storage = layer_stack::storage(node);
  socket_items::ui::draw_active_item_props<mask_stack::ItemsAccessor>(
      ntree, node, [&](PointerRNA *item_ptr) {
        layout.use_property_split_set(true);
        layout.use_property_decorate_set(false);
        if (storage.active_index != storage.items_num - 1) {
          layout.prop(item_ptr, "blend_type", UI_ITEM_NONE, std::nullopt, ICON_NONE);
        }
        layout.prop(item_ptr, "mute", UI_ITEM_NONE, std::nullopt, ICON_NONE);
      });
}

}  // namespace nodes::node_shader_mask_stack_cc

void register_node_type_sh_mask_stack()
{
  namespace file_ns = nodes::node_shader_mask_stack_cc;
  using Accessor = nodes::mask_stack::ItemsAccessor;

  static bke::bNodeType ntype;
  common_node_type_base(&ntype, "ShaderNodeMaskStack"_ustr, SH_NODE_MASK_STACK);
  ntype.ui_name = "Mask Stack";
  ntype.ui_description =
      "Composite a stack of mask layers (floats) with blend modes and opacity. "
      "Used as the source of a texture layer's Mask socket so masks are themselves "
      "stackable";
  ntype.enum_name_legacy = "MASK_STACK";
  ntype.nclass = NODE_CLASS_OP_COLOR;
  ntype.declare = file_ns::node_declare;
  ntype.initfunc = nodes::layer_stack::node_init;
  ntype.insert_link = nodes::layer_stack::node_insert_link<Accessor>;
  ntype.draw_buttons_ex = file_ns::node_layout_ex;
  ntype.register_operators = file_ns::node_operators;
  ntype.blend_write_storage_content = nodes::layer_stack::node_blend_write<Accessor>;
  ntype.blend_data_read_storage_content = nodes::layer_stack::node_blend_read<Accessor>;
  bke::node_type_storage(ntype,
                         "NodeShaderLayerStack",
                         nodes::layer_stack::node_free_storage<Accessor>,
                         nodes::layer_stack::node_copy_storage<Accessor>);

  bke::node_register_type(ntype);
}

namespace nodes {

StructRNA **mask_stack::ItemsAccessor::item_srna = &RNA_ShaderMaskStackItem;

}  // namespace nodes

}  // namespace blender
