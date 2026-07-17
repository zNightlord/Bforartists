/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup shdnodes
 */

#include "NOD_sh_layer_stack.hh"
#include "NOD_socket_items_ops.hh"
#include "NOD_socket_items_ui.hh"

#include "RNA_prototypes.hh"

#include "node_shader_util.hh"
#include "node_util.hh"

namespace blender {

namespace nodes::node_shader_texture_layer_stack_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  const bNodeTree *tree = b.tree_or_null();
  const bNode *node = b.node_or_null();
  if (node) {
    const NodeShaderLayerStack &storage = layer_stack_storage(*node);
    const int last_index = storage.items_num - 1;
    for (const int i : IndexRange(storage.items_num)) {
      const NodeShaderLayerStackItem &item = storage.items[i];
      const UString name = item.name ? UString(item.name) : ""_ustr;
      const UString layer_id(TextureLayerStackItemsAccessor::socket_identifier_for_item(item));
      /* Closure-typed layers (adjustments) take a Closure that the inliner
       * evaluates with the accumulated stack below as its bundle input;
       * Bundle-typed layers (generators) take the bundle directly. Untyped
       * layers expose a hollow extend socket that takes its type from the
       * first link dropped onto it. */
      BaseSocketDeclarationBuilder *layer_decl = nullptr;
      switch (eShaderLayerStackItemType(item.item_type)) {
        case SHADER_LAYER_STACK_ITEM_CLOSURE:
          layer_decl = &b.add_input<decl::Closure>(name, layer_id);
          break;
        case SHADER_LAYER_STACK_ITEM_BUNDLE:
          layer_decl =
              &b.add_input<decl::Bundle>(name, layer_id).structure_type(StructureType::Single);
          break;
        case SHADER_LAYER_STACK_ITEM_UNSET:
        default:
          layer_decl = &b.add_input<decl::Extend>(name, layer_id);
          break;
      }
      if (tree) {
        layer_decl->socket_name_ptr(
            &tree->id, *TextureLayerStackItemsAccessor::item_srna, &item, "name");
      }
      /* Non-base layers expose Opacity and Mask. The base (last) layer has no
       * Opacity (nothing below it to blend over), but still exposes a Mask so
       * a masked layer keeps its mask when it becomes the base: there the mask
       * blends between the plain pre-stack value and this layer. */
      if (i != last_index) {
        const UString opacity_id(
            TextureLayerStackItemsAccessor::opacity_socket_identifier_for_item(item));
        b.add_input<decl::Float>("Opacity"_ustr, opacity_id)
            .default_value(1.0f)
            .min(0.0f)
            .max(1.0f)
            .subtype(PROP_FACTOR)
            .description("Opacity of this layer over the layers below");
      }
      const UString mask_id(TextureLayerStackItemsAccessor::mask_socket_identifier_for_item(item));
      b.add_input<decl::Float>("Mask"_ustr, mask_id)
          .default_value(1.0f)
          .min(0.0f)
          .max(1.0f)
          .subtype(PROP_FACTOR)
          .hide_value()
          .description(
              i != last_index ?
                  "Mask multiplied with Opacity to form the final blend weight. Defaults to 1.0 "
                  "when nothing is connected so it has no effect" :
                  "Mask blending between the plain pre-stack value and this base layer. Defaults "
                  "to 1.0 when nothing is connected so it has no effect");
    }
  }
  /* Trailing extend-socket: dragging any link onto this hollow socket adds a
   * new layer (Combine-Bundle style). #node_insert_link appends the new item at
   * the end, so it becomes the new base and the previously-last item gains its
   * blend-weight sockets. */
  b.add_input<decl::Extend>(""_ustr, "__extend__"_ustr)
      .custom_draw(socket_items::ui::draw_extend_socket_fn<TextureLayerStackItemsAccessor>());
  b.add_output<decl::Bundle>("Result"_ustr).structure_type(StructureType::Single);
}

/* Keep each item's type in sync with its link state:
 * - A typed item whose layer socket lost its links resets to UNSET, so the
 *   socket becomes a hollow extend socket that can be re-typed by the next
 *   dropped link.
 * - An UNSET item whose socket is linked infers its type from the link. This
 *   covers links created without the insert-link callback (Python's
 *   links.new) and upgrades files saved before the type existed (their
 *   zeroed #item_type reads as UNSET while the sockets are linked bundles).
 *
 * Links survive the socket re-declaration: the socket keeps its identifier
 * and the refresh migrates links onto the re-built socket. */
static void node_update(bNodeTree *ntree, bNode *node)
{
  NodeShaderLayerStack &storage = layer_stack_storage(*node);
  bool changed = false;
  for (const int i : IndexRange(storage.items_num)) {
    NodeShaderLayerStackItem &item = storage.items[i];
    const std::string identifier = TextureLayerStackItemsAccessor::socket_identifier_for_item(
        item);
    /* Walk the links directly: the topology cache may be dirty here. */
    const bNodeSocket *from_socket = nullptr;
    for (const bNodeLink &link : ntree->links) {
      if (link.tonode == node && link.tosock && identifier == link.tosock->identifier) {
        from_socket = link.fromsock;
        break;
      }
    }
    if (item.item_type == SHADER_LAYER_STACK_ITEM_UNSET && from_socket) {
      if (ELEM(from_socket->type, SOCK_BUNDLE, SOCK_CLOSURE)) {
        item.item_type = (from_socket->type == SOCK_CLOSURE) ? SHADER_LAYER_STACK_ITEM_CLOSURE :
                                                               SHADER_LAYER_STACK_ITEM_BUNDLE;
        changed = true;
      }
    }
    else if (item.item_type != SHADER_LAYER_STACK_ITEM_UNSET && !from_socket) {
      item.item_type = SHADER_LAYER_STACK_ITEM_UNSET;
      /* The extend declaration accepts any socket with a matching identifier,
       * so the typed socket must be removed explicitly for the refresh below
       * to build the hollow one in its place. */
      if (bNodeSocket *sock = bke::node_find_socket(*node, SOCK_IN, UString(identifier))) {
        bke::node_remove_socket(*ntree, *node, *sock);
      }
      changed = true;
    }
  }
  if (changed) {
    /* Refresh the sockets right away: this runs after the updater's own
     * declaration pass, so a deferred tag alone would leave the sockets
     * stale until some later update. Links keep their sockets by identifier. */
    update_node_declaration_and_sockets(*ntree, *node);
    BKE_ntree_update_tag_node_property(ntree, node);
  }
}

static void node_operators()
{
  socket_items::ops::make_common_operators<TextureLayerStackItemsAccessor>();
}

/* Sidebar layout: the layer items list with add/remove/move, plus blend mode
 * and mute of the active layer. The base (last) layer has no blend mode. */
static void node_layout_ex(ui::Layout &layout, bContext *C, PointerRNA *node_ptr)
{
  bNodeTree &ntree = *reinterpret_cast<bNodeTree *>(node_ptr->owner_id);
  const bNode &node = *static_cast<const bNode *>(node_ptr->data);
  socket_items::ui::draw_items_list_with_operators<TextureLayerStackItemsAccessor>(
      C, &layout, ntree, node);
  const NodeShaderLayerStack &storage = layer_stack_storage(node);
  socket_items::ui::draw_active_item_props<TextureLayerStackItemsAccessor>(
      ntree, node, [&](PointerRNA *item_ptr) {
        layout.use_property_split_set(true);
        layout.use_property_decorate_set(false);
        if (storage.active_index != storage.items_num - 1) {
          layout.prop(item_ptr, "blend_type", UI_ITEM_NONE, std::nullopt, ICON_NONE);
        }
        layout.prop(item_ptr, "mute", UI_ITEM_NONE, std::nullopt, ICON_NONE);
      });
}

}  // namespace nodes::node_shader_texture_layer_stack_cc

void register_node_type_sh_texture_layer_stack()
{
  namespace file_ns = nodes::node_shader_texture_layer_stack_cc;
  using Accessor = nodes::TextureLayerStackItemsAccessor;

  static bke::bNodeType ntype;
  common_node_type_base(&ntype, "ShaderNodeTextureLayerStack"_ustr, SH_NODE_TEXTURE_LAYER_STACK);
  ntype.ui_name = "Texture Layer Stack";
  ntype.ui_description =
      "Composite a stack of texture layer bundles, channel by channel. "
      "Expands into per-channel Mix Color / Mix Vector / Mix Float nodes during shader node "
      "inlining";
  ntype.enum_name_legacy = "TEXTURE_LAYER_STACK";
  ntype.nclass = NODE_CLASS_OP_COLOR;
  ntype.declare = file_ns::node_declare;
  ntype.initfunc = nodes::layer_stack::node_init;
  ntype.updatefunc = file_ns::node_update;
  ntype.insert_link = nodes::layer_stack::node_insert_link<Accessor>;
  ntype.draw_buttons_ex = file_ns::node_layout_ex;
  ntype.register_operators = file_ns::node_operators;
  ntype.blend_write_storage_content = nodes::layer_stack::node_blend_write<Accessor>;
  ntype.blend_data_read_storage_content = nodes::layer_stack::node_blend_read<Accessor>;
  bke::node_type_storage(ntype,
                         "NodeShaderLayerStack",
                         nodes::layer_stack::node_free_storage<Accessor>,
                         nodes::layer_stack::node_copy_storage<Accessor>);
  ntype.is_bundle_join = true;

  bke::node_register_type(ntype);
}

namespace nodes {

StructRNA **TextureLayerStackItemsAccessor::item_srna = &RNA_ShaderTextureLayerStackItem;

}  // namespace nodes

Span<NodeShaderLayerStackItem> NodeShaderLayerStack::items_span() const
{
  return Span<NodeShaderLayerStackItem>(items, items_num);
}

MutableSpan<NodeShaderLayerStackItem> NodeShaderLayerStack::items_span()
{
  return MutableSpan<NodeShaderLayerStackItem>(items, items_num);
}

}  // namespace blender
