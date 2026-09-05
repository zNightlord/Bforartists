/* SPDX-FileCopyrightText: 2026 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edrend
 *
 * Operators that edit the Texture Layer Stack from the Material Surface panel's
 * Shader Layers tree view: add/remove layers, Fill and mask layers, grouping
 * and drag & drop reordering (reparent). This is a thin layer: context
 * resolution, polls, properties and notifiers live here, while the graph
 * interpretation and manipulation live in #NOD_texture_stack.hh and its
 * sibling shader-layer headers, so the operators, the tree view and the node
 * editor agree on how the node graph maps to layers.
 */

#include <algorithm>
#include <climits>

#include "BLI_listbase_iterator.hh"
#include "BLI_string_ref.hh"
#include "BLI_string_utf8.hh"
#include "BLI_vector.hh"
#include "BLI_vector_set.hh"

#include "DNA_array_utils.hh"
#include "DNA_image_types.h"
#include "DNA_material_types.h"
#include "DNA_node_types.h"
#include "DNA_userdef_types.h"

#include "BKE_context.hh"
#include "BKE_image.hh"
#include "BKE_lib_id.hh"
#include "BKE_library.hh"
#include "BKE_main.hh"
#include "BKE_main_invariants.hh"
#include "BKE_node.hh"
#include "BKE_node_runtime.hh"
#include "BKE_node_tree_update.hh"
#include "BKE_report.hh"

#include "RNA_access.hh"
#include "RNA_define.hh"
#include "RNA_prototypes.hh"

#include "NOD_closure_zone.hh"
#include "NOD_geo_bundle.hh"
#include "NOD_layer_stack.hh"
#include "NOD_mask_stack.hh"
#include "NOD_node_placement.hh"
#include "NOD_socket_items.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"

#include "WM_api.hh"
#include "WM_types.hh"

#include "render_intern.hh"

namespace blender::ed::render {

namespace layer_stack = nodes::layer_stack;
namespace texture_stack = nodes::texture_stack;
namespace mask_stack = nodes::mask_stack;
namespace texture_channel = nodes::texture_channel;
namespace closure_zone = nodes::closure_zone;

using nodes::CombineBundleItemsAccessor;

/* -------------------------------------------------------------------- */
/** \name Shared Helpers
 * \{ */

Material *active_material(const bContext &C)
{
  PointerRNA mat_ptr = CTX_data_pointer_get_type(&C, "material", RNA_Material);
  return static_cast<Material *>(mat_ptr.data);
}

static bNodeTree *active_material_ntree(const bContext &C)
{
  Material *mat = active_material(C);
  return (mat && mat->nodetree) ? mat->nodetree : nullptr;
}

/* True when the active material's embedded node tree may be structurally edited. Linked or
 * library-override materials must not be edited: changes are lost on reload or corrupt overrides.
 * Gates the layer operators and tree-view edit callbacks (the panels still draw for viewing). */
bool active_material_editable(const bContext &C)
{
  const Material *mat = active_material(C);
  return mat && ID_IS_EDITABLE(mat) && !ID_IS_OVERRIDE_LIBRARY(mat);
}

/* Common operator gate: the Texture Layers feature is enabled (it is experimental, and its
 * panels/menus are only drawn then) and the active material can be edited. Keeps the operators
 * off operator search when the feature is disabled. */
static bool texture_layers_poll(const bContext &C)
{
  return U.experimental.use_texture_layers && active_material_editable(C);
}

std::optional<ActiveStackContext> resolve_active_stack(const bContext &C)
{
  Material *mat = active_material(C);
  if (mat == nullptr || mat->nodetree == nullptr) {
    return std::nullopt;
  }
  bNode *stack = texture_stack::active(*mat->nodetree);
  if (stack == nullptr) {
    return std::nullopt;
  }
  return ActiveStackContext{
      mat, mat->nodetree, stack, texture_stack::active_layer_index(*mat->nodetree, *stack)};
}

void tree_changed(bContext &C, bNodeTree &ntree, Material *mat)
{
  BKE_main_ensure_invariants(*CTX_data_main(&C), ntree.id);
  WM_event_add_notifier(&C, NC_NODE | NA_EDITED, &ntree.id);
  if (mat) {
    WM_event_add_notifier(&C, NC_MATERIAL | ND_SHADING, mat);
  }
}

bNode *layer_active_node(bNodeTree &ntree, bNode &stack, const int index)
{
  if (index < 0 || index >= layer_stack::storage(stack).items_num) {
    return &stack;
  }
  if (mask_stack::is_stack(stack)) {
    bNode *source = mask_stack::layer_source(ntree, stack, index);
    return source ? source : &stack;
  }
  const texture_stack::StackLayer layer{&ntree, &stack, index};
  if (layer.nested_stack()) {
    /* A group layer: its nested stack stands for the group's children, so the
     * group itself is identified by the stack node holding it. */
    return &stack;
  }
  bNode *source = layer.display_source();
  return source ? source : &stack;
}

void select_layer_nodes(bNodeTree &ntree, bNode &stack, const int index)
{
  VectorSet<bNode *> nodes;
  if (index >= 0 && index < layer_stack::storage(stack).items_num) {
    texture_stack::nodes_for_layer({&ntree, &stack, index}, nodes);
  }
  nodes.add(layer_active_node(ntree, stack, index));
  for (bNode &node : ntree.nodes) {
    bke::node_set_selected(node, nodes.contains(&node));
  }
}

void set_active_layer(bNodeTree &ntree, bNode &stack, const int index)
{
  layer_stack::storage(stack).active_index = index;
  bke::node_set_active(ntree, *layer_active_node(ntree, stack, index));
  select_layer_nodes(ntree, stack, index);
}

/* The Mask Stack the active node belongs to, with a valid active item; else
 * none (a selected mask row in the tree view). */
struct ActiveMask {
  bNode *stack;
  int index;
};
static std::optional<ActiveMask> active_mask(bNodeTree &ntree)
{
  bNode *active = bke::node_get_active(ntree);
  if (active == nullptr) {
    return std::nullopt;
  }
  const std::optional<texture_stack::StackLayer> layer = texture_stack::layer_for_node(ntree,
                                                                                       *active);
  if (!layer || !mask_stack::is_stack(*layer->stack)) {
    return std::nullopt;
  }
  return ActiveMask{layer->stack, layer->index};
}

/* Remove mask item #index of #mask_stack_node, its source nodes, and the
 * Mask Stack node itself when it becomes empty (so the layer's Mask input
 * returns to its default, fully opaque). */
static void remove_mask(bContext &C, bNodeTree &ntree, bNode &mask_stack_node, const int index)
{
  Main &bmain = *CTX_data_main(&C);
  if (index < 0 || index >= layer_stack::storage(mask_stack_node).items_num) {
    return;
  }
  bNode *source = mask_stack::layer_source(ntree, mask_stack_node, index);
  /* The layer stack the mask belongs to, resolved while the mask's nodes (which
   * the active node is one of) are still around. */
  bNode *layer_stack_node = texture_stack::active(ntree);
  layer_stack::remove_layer(bmain, ntree, mask_stack_node, index);

  VectorSet<bNode *> deleted;
  texture_stack::delete_layer_source_recursive(bmain, ntree, source, deleted);

  if (layer_stack::storage(mask_stack_node).items_num > 0) {
    set_active_layer(ntree, mask_stack_node, layer_stack::active_index_in_range(mask_stack_node));
    return;
  }
  bke::node_remove_node(&bmain, ntree, mask_stack_node, true);
  /* The last mask is gone: fall back to the layer that owned the mask stack. */
  if (layer_stack_node) {
    set_active_layer(
        ntree, *layer_stack_node, layer_stack::active_index_in_range(*layer_stack_node));
  }
}

struct FillChannel {
  eNodeSocketDatatype type;
  std::string name;
};

/* The channels a bootstrapped Fill layer covers when the material has no textured channels yet:
 * common Principled BSDF inputs, so they map cleanly to channels. The first is the base channel.
 * Listed in the BSDF's own input order, so the links into it do not cross.
 */
static Vector<FillChannel> default_fill_channels()
{
  return {{SOCK_RGBA, "Base Color"},
          {SOCK_FLOAT, "Metallic"},
          {SOCK_FLOAT, "Roughness"},
          {SOCK_VECTOR, "Normal"}};
}

/* Mirror the BSDF input #channel's UI (subtype, soft range, default) onto a Combine Bundle #item,
 * so a Fill layer's channel input shows the BSDF's slider instead of a bare value. */
static void inherit_channel_ui(bNode *bsdf, NodeCombineBundleItem &item, const StringRef channel)
{
  if (bsdf == nullptr) {
    return;
  }
  for (bNodeSocket &sock : bsdf->inputs) {
    if (StringRef(sock.name) == channel && texture_channel::is_fillable_input(sock)) {
      nodes::combine_bundle_item_copy_socket_data(item, sock);
      break;
    }
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Polls
 * \{ */

static bool active_stack_poll(bContext *C)
{
  return texture_layers_poll(*C) && resolve_active_stack(*C).has_value();
}

/* Unlike #active_stack_poll, this also accepts a material with no Texture
 * Layer Stack yet: #texture_layer_add_exec creates the root stack on demand. */
static bool active_material_ntree_poll(bContext *C)
{
  return texture_layers_poll(*C) && active_material_ntree(*C) != nullptr;
}

/* There is a layer to operate on: the selected mask row, or the active layer of
 * the active stack. Also gates adding masks, which every layer accepts
 * (including the base, which has a Mask socket too). */
static bool active_layer_poll(bContext *C)
{
  if (!texture_layers_poll(*C)) {
    return false;
  }
  bNodeTree *ntree = active_material_ntree(*C);
  if (ntree && active_mask(*ntree)) {
    return true;
  }
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  return ctx && ctx->layer_index != -1;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add / Remove Layer
 * \{ */

static wmOperatorStatus texture_layer_add_exec(bContext *C, wmOperator *op)
{
  Main &bmain = *CTX_data_main(C);
  Material *mat = active_material(*C);
  if (mat == nullptr || mat->nodetree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *mat->nodetree;
  bNode *stack = texture_stack::active(ntree);
  const bool bootstrapping = stack == nullptr;
  Vector<bNode *> bsdfs;
  if (bootstrapping) {
    texture_channel::collect_bsdfs(ntree, bsdfs);
    if (bsdfs.is_empty()) {
      /* Nothing to wire a stack into; don't leave an unwired orphan stack node behind. */
      BKE_report(op->reports, RPT_WARNING, "No shader output to attach a texture layer stack to");
      return OPERATOR_CANCELLED;
    }
    stack = texture_stack::create_root(ntree, bsdfs[0]);
    if (stack == nullptr) {
      return OPERATOR_CANCELLED;
    }
  }
  char name[MAX_NAME];
  RNA_string_get(op->ptr, "name", name);
  const int index = texture_stack::add_layer(bmain, ntree, *stack, name);
  if (bootstrapping) {
    texture_stack::wire_root_into_bsdf(bmain, ntree, *stack, index, *bsdfs[0]);
  }
  set_active_layer(ntree, *stack, index);
  tree_changed(*C, ntree, mat);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_add(wmOperatorType *ot)
{
  ot->name = "Add Texture Layer";
  ot->description =
      "Add a new layer to the active Texture Layer Stack, creating one first if needed";
  ot->idname = "MATERIAL_OT_texture_layer_add";

  ot->exec = texture_layer_add_exec;
  ot->poll = active_material_ntree_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  RNA_def_string(ot->srna, "name", "Layer", MAX_NAME, "Name", "Name of the new layer");
}

/* TODO: hardcoded default channels for a bootstrapped stack; make configurable. */
static wmOperatorStatus texture_layer_add_default_exec(bContext *C, wmOperator *op)
{
  Main &bmain = *CTX_data_main(C);
  Material *mat = active_material(*C);
  if (mat == nullptr || mat->nodetree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *mat->nodetree;

  /* This operator bootstraps a stack from scratch. Refuse when the material already has one
   * (even a stack disconnected from any BSDF): creating a second would orphan the
   * existing stack and all its layers. */
  for (const bNode &node : ntree.nodes) {
    if (texture_stack::is_stack(node)) {
      BKE_report(op->reports, RPT_WARNING, "Material already has a texture layer stack");
      return OPERATOR_CANCELLED;
    }
  }

  Vector<bNode *> bsdfs;
  texture_channel::collect_bsdfs(ntree, bsdfs);
  if (bsdfs.is_empty()) {
    return OPERATOR_CANCELLED;
  }
  bNode &bsdf = *bsdfs[0];

  bNode *stack = texture_stack::create_root(ntree, &bsdf);
  if (stack == nullptr) {
    return OPERATOR_CANCELLED;
  }
  const int index = texture_stack::add_layer(bmain, ntree, *stack, "Fill");

  auto free_channel_socket = [&](const StringRef name) -> bNodeSocket * {
    for (bNodeSocket &sock : bsdf.inputs) {
      if (StringRef(sock.name) == name && texture_channel::is_fillable_input(sock) &&
          texture_channel::state(ntree, sock) == texture_channel::State::Free)
      {
        return &sock;
      }
    }
    return nullptr;
  };

  /* Bootstrap the base Fill with the first default channel, then route the rest. */
  const Vector<FillChannel> defaults = default_fill_channels();
  texture_stack::wire_root_into_bsdf(
      bmain, ntree, *stack, index, bsdf, free_channel_socket(defaults[0].name));

  for (const FillChannel &channel : defaults.as_span().drop_front(1)) {
    if (bNodeSocket *sock = free_channel_socket(channel.name)) {
      if (bNode *separate = texture_channel::separate_bundle_for_bsdf(ntree, bsdf)) {
        texture_channel::link(ntree, *separate, bsdf, *sock);
      }
    }
  }

  set_active_layer(ntree, *stack, index);
  tree_changed(*C, ntree, mat);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_add_default(wmOperatorType *ot)
{
  ot->name = "Add Texture Layers";
  ot->description =
      "Create a Texture Layer Stack with a Fill layer covering the material's main channels";
  ot->idname = "MATERIAL_OT_texture_layer_add_default";

  ot->exec = texture_layer_add_default_exec;
  ot->poll = active_material_ntree_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static wmOperatorStatus texture_layer_remove_exec(bContext *C, wmOperator * /*op*/)
{
  Main &bmain = *CTX_data_main(C);

  /* A selected mask row removes that mask, so the panel's single remove
   * button works on whatever the tree view has active. */
  if (bNodeTree *mask_ntree = active_material_ntree(*C)) {
    if (const std::optional<ActiveMask> mask = active_mask(*mask_ntree)) {
      remove_mask(*C, *mask_ntree, *mask->stack, mask->index);
      tree_changed(*C, *mask_ntree, active_material(*C));
      return OPERATOR_FINISHED;
    }
  }

  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;
  const int index = ctx->layer_index;
  if (index == -1) {
    return OPERATOR_CANCELLED;
  }

  /* Trace the layer's source nodes before removing the item — once the item is
   * gone its sockets disappear and the links can't be followed. */
  const texture_stack::StackLayer layer{&ntree, &stack, index};
  bNode *bundle_source = layer.source();
  bNodeSocket *mask_sock = layer.mask_socket();
  bNode *mask_source = mask_sock ? bke::node_find_source_node(ntree, *mask_sock) : nullptr;

  /* Removing the last layer of a root stack tears down the whole stack: the
   * stack node and its Separate Bundle go too, restoring the BSDF to plain
   * inputs. A nested stack (a group's contents) is left empty instead, so
   * removing a group's last layer never touches the outer stack's Separate
   * Bundle and disconnects the whole material. Resolve this while the topology
   * is still intact. */
  bool remove_stack = false;
  bNode *separate = nullptr;
  if (layer_stack::storage(stack).items_num <= 1) {
    if (bNode *bsdf = texture_stack::bsdf_for(ntree, stack)) {
      Vector<bNode *> roots;
      texture_stack::roots_for_bsdf(ntree, *bsdf, roots);
      if (roots.contains(&stack)) {
        remove_stack = true;
        separate = texture_channel::separate_bundle_for_bsdf(ntree, *bsdf);
      }
    }
  }

  layer_stack::remove_layer(bmain, ntree, stack, index);

  VectorSet<bNode *> deleted;
  texture_stack::delete_layer_source_recursive(bmain, ntree, bundle_source, deleted);
  texture_stack::delete_layer_source_recursive(bmain, ntree, mask_source, deleted);

  if (remove_stack) {
    if (separate) {
      bke::node_remove_node(&bmain, ntree, *separate, true);
    }
    bke::node_remove_node(&bmain, ntree, stack, true);
  }
  else {
    /* The removed layer's nodes are gone with it, so point the active node at
     * the layer that took its place. */
    set_active_layer(ntree, stack, layer_stack::active_index_in_range(stack));
  }

  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_remove(wmOperatorType *ot)
{
  ot->name = "Remove Texture Layer";
  ot->description = "Remove the active layer from the Texture Layer Stack";
  ot->idname = "MATERIAL_OT_texture_layer_remove";

  ot->exec = texture_layer_remove_exec;
  ot->poll = active_layer_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Channel Toggle
 * \{ */

static wmOperatorStatus texture_channel_toggle_exec(bContext *C, wmOperator *op)
{
  Main &bmain = *CTX_data_main(C);
  Material *mat = active_material(*C);
  if (mat == nullptr || mat->nodetree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *mat->nodetree;

  char node_name[MAX_NAME];
  char socket_id[MAX_NAME];
  RNA_string_get(op->ptr, "node", node_name);
  RNA_string_get(op->ptr, "socket", socket_id);
  bNode *bsdf = bke::node_find_node_by_name(ntree, node_name);
  if (bsdf == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNodeSocket *socket = bke::node_find_socket(*bsdf, SOCK_IN, UString(socket_id));
  if (socket == nullptr || !texture_channel::is_fillable_input(*socket)) {
    return OPERATOR_CANCELLED;
  }

  switch (texture_channel::state(ntree, *socket)) {
    case texture_channel::State::Textured:
      if (!texture_channel::unlink(ntree, *socket)) {
        return OPERATOR_CANCELLED;
      }
      break;
    case texture_channel::State::Free:
      if (bNode *separate = texture_channel::separate_bundle_for_bsdf(ntree, *bsdf)) {
        if (!texture_channel::link(ntree, *separate, *bsdf, *socket)) {
          return OPERATOR_CANCELLED;
        }
      }
      else {
        /* No stack routes into this BSDF yet: bootstrap one whose base fill
         * layer carries this channel. */
        bNode *stack = texture_stack::create_root(ntree, bsdf);
        if (stack == nullptr) {
          return OPERATOR_CANCELLED;
        }
        const int index = texture_stack::add_layer(bmain, ntree, *stack, "Fill");
        texture_stack::wire_root_into_bsdf(bmain, ntree, *stack, index, *bsdf, socket);
        set_active_layer(ntree, *stack, index);
      }
      break;
    case texture_channel::State::Linked:
      BKE_report(op->reports, RPT_WARNING, "Input is linked to another node");
      return OPERATOR_CANCELLED;
  }

  tree_changed(*C, ntree, mat);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_channel_toggle(wmOperatorType *ot)
{
  ot->name = "Toggle Texture Channel";
  ot->description =
      "Route the shader input through the material's texture layer stack, or restore it "
      "to a plain value. The layers keep their channel data when disabled";
  ot->idname = "MATERIAL_OT_texture_channel_toggle";

  ot->exec = texture_channel_toggle_exec;
  ot->poll = active_material_ntree_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  RNA_def_string(ot->srna, "node", nullptr, MAX_NAME, "Node", "Name of the BSDF node");
  RNA_def_string(
      ot->srna, "socket", nullptr, MAX_NAME, "Socket", "Identifier of the input socket");
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add Fill Layer
 * \{ */

/* The channels a new Fill layer carries: the ones the material textures
 * (BSDF inputs routed through the stack). When none are, fall back to the
 * default channels so they map cleanly to the BSDF. */
static Vector<FillChannel> fill_layer_channels(const bNodeTree &ntree, bNode &stack)
{
  Vector<FillChannel> channels;
  if (bNode *bsdf = texture_stack::bsdf_for(ntree, stack)) {
    for (bNodeSocket &sock : bsdf->inputs) {
      if (!texture_channel::is_fillable_input(sock)) {
        continue;
      }
      if (texture_channel::state(ntree, sock) == texture_channel::State::Textured) {
        channels.append({eNodeSocketDatatype(sock.type), sock.name});
      }
    }
  }
  if (channels.is_empty()) {
    channels = default_fill_channels();
  }
  return channels;
}

static wmOperatorStatus texture_layer_add_fill_exec(bContext *C, wmOperator * /*op*/)
{
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;

  const Vector<FillChannel> channels = fill_layer_channels(ntree, stack);
  const int index = texture_stack::add_layer(*CTX_data_main(C), ntree, stack, "Fill");

  bNode *combine = bke::node_add_node(C, ntree, "NodeCombineBundle"_ustr);
  if (combine == nullptr) {
    /* Roll back the layer item added above so no empty, source-less layer is left behind. */
    layer_stack::remove_layer(*CTX_data_main(C), ntree, stack, index);
    return OPERATOR_CANCELLED;
  }
  STRNCPY_UTF8(combine->label, "Fill");
  bNode *bsdf = texture_stack::bsdf_for(ntree, stack);
  for (const FillChannel &channel : channels) {
    NodeCombineBundleItem *item =
        nodes::socket_items::add_item_with_socket_type_and_name<CombineBundleItemsAccessor>(
            ntree, *combine, channel.type, channel.name.c_str());
    if (item) {
      inherit_channel_ui(bsdf, *item, channel.name);
    }
  }
  BKE_ntree_update_tag_node_property(&ntree, combine);
  BKE_main_ensure_invariants(*CTX_data_main(C), ntree.id);
  /* Placed once its sockets exist, so its height is known. */
  texture_stack::place_layer_source(*combine, {&ntree, &stack, index}, 280.0f);

  texture_stack::connect_bundle({&ntree, &stack, index}, *combine, "Bundle");
  set_active_layer(ntree, stack, index);
  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_add_fill(wmOperatorType *ot)
{
  ot->name = "Add Fill Layer";
  ot->description = "Add a new Fill texture layer (Combine Bundle with default channels)";
  ot->idname = "MATERIAL_OT_texture_layer_add_fill";

  ot->exec = texture_layer_add_fill_exec;
  ot->poll = active_stack_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add Paint Layer
 * \{ */

/* Channel layout and initial fill of the image pass backing #channel. The pass
 * names are the channel names, so the Image Texture node's Passes bundle
 * output feeds the stack's channels directly. The fill mirrors a Fill layer:
 * the BSDF input's default value (in linear space, matching the float image's
 * colorspace), so a fresh Paint layer carries the channel defaults until
 * painted. */
static ImageGeneratedPass image_pass_for_channel(bNode *bsdf, const FillChannel &channel)
{
  ImageGeneratedPass pass;
  pass.name = channel.name;
  switch (channel.type) {
    case SOCK_FLOAT:
      pass.chan_id = "X";
      pass.channels_num = 1;
      break;
    case SOCK_VECTOR:
      pass.chan_id = "XYZ";
      pass.channels_num = 3;
      break;
    default:
      pass.chan_id = "RGBA";
      pass.channels_num = 4;
      break;
  }

  if (bsdf == nullptr) {
    return pass;
  }
  for (bNodeSocket &sock : bsdf->inputs) {
    if (StringRef(sock.name) != channel.name || !texture_channel::is_fillable_input(sock)) {
      continue;
    }
    switch (sock.type) {
      case SOCK_FLOAT: {
        const float value = sock.default_value_typed<bNodeSocketValueFloat>()->value;
        pass.color[0] = pass.color[1] = pass.color[2] = value;
        pass.color[3] = 1.0f;
        break;
      }
      case SOCK_VECTOR: {
        const float *value = sock.default_value_typed<bNodeSocketValueVector>()->value;
        pass.color[0] = value[0];
        pass.color[1] = value[1];
        pass.color[2] = value[2];
        pass.color[3] = 1.0f;
        break;
      }
      case SOCK_RGBA: {
        const float *value = sock.default_value_typed<bNodeSocketValueRGBA>()->value;
        pass.color[0] = value[0];
        pass.color[1] = value[1];
        pass.color[2] = value[2];
        pass.color[3] = value[3];
        break;
      }
      default:
        break;
    }
    break;
  }
  return pass;
}

static wmOperatorStatus texture_layer_add_paint_exec(bContext *C, wmOperator *op)
{
  Main &bmain = *CTX_data_main(C);
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;

  const Vector<FillChannel> channels = fill_layer_channels(ntree, stack);
  bNode *bsdf = texture_stack::bsdf_for(ntree, stack);
  Vector<ImageGeneratedPass> passes;
  for (const FillChannel &channel : channels) {
    passes.append(image_pass_for_channel(bsdf, channel));
  }

  char name[MAX_NAME];
  RNA_string_get(op->ptr, "name", name);

  /* One multi-layer image whose passes are the layer's channels: each pass is
   * its own paintable buffer, created on first use. */
  Image *image = BKE_image_add_generated_multilayer(
      &bmain, RNA_int_get(op->ptr, "width"), RNA_int_get(op->ptr, "height"), name, passes);
  if (image == nullptr) {
    return OPERATOR_CANCELLED;
  }

  const int index = texture_stack::add_layer(bmain, ntree, stack, name);

  bNode *image_node = bke::node_add_node(C, ntree, "ShaderNodeTexImage"_ustr);
  if (image_node == nullptr) {
    layer_stack::remove_layer(bmain, ntree, stack, index);
    BKE_id_delete(&bmain, image);
    return OPERATOR_CANCELLED;
  }
  /* The image's initial user from creation becomes the node's user. */
  image_node->id = &image->id;
  STRNCPY_UTF8(image_node->label, name);

  /* Rebuild the node's declaration for the newly assigned multi-layer image,
   * so the per-pass outputs and the Passes bundle output exist to connect. */
  bke::node_tag_update_id(*image_node);
  BKE_ntree_update_tag_node_property(&ntree, image_node);
  BKE_main_ensure_invariants(bmain, ntree.id);
  texture_stack::place_layer_source(*image_node, {&ntree, &stack, index}, 320.0f);

  texture_stack::connect_bundle({&ntree, &stack, index}, *image_node, "Passes");
  set_active_layer(ntree, stack, index);
  /* Make the fresh image the paint target right away (the tree view keeps it
   * in sync when other layers are selected). */
  bke::node_set_active_texture(ntree, *image_node);
  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_add_paint(wmOperatorType *ot)
{
  ot->name = "Add Paint Layer";
  ot->description =
      "Add a texture layer backed by a new multi-layer image with a paintable image pass "
      "for each of the layer's channels";
  ot->idname = "MATERIAL_OT_texture_layer_add_paint";

  ot->exec = texture_layer_add_paint_exec;
  ot->poll = active_stack_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

  RNA_def_string(ot->srna, "name", "Paint", MAX_NAME, "Name", "Name of the new layer and image");
  RNA_def_int(ot->srna, "width", 1024, 1, INT_MAX, "Width", "Image width", 1, 16384);
  RNA_def_int(ot->srna, "height", 1024, 1, INT_MAX, "Height", "Image height", 1, 16384);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Add Mask Layers (White / Black / Paint)
 * \{ */

bNode *ensure_active_layer_mask_stack(bContext &C, bNodeTree &ntree)
{
  /* A selected mask row targets its own Mask Stack: the new mask is added
   * above the active one. */
  if (const std::optional<ActiveMask> mask = active_mask(ntree)) {
    return mask->stack;
  }
  bNode *stack = texture_stack::active(ntree);
  if (stack == nullptr) {
    return nullptr;
  }
  /* Every layer, including the base, has a Mask socket. */
  const int index = texture_stack::active_layer_index(ntree, *stack);
  if (index == -1) {
    return nullptr;
  }
  const NodeShaderLayerStack &storage = layer_stack::storage(*stack);
  return mask_stack::ensure_for_layer(*CTX_data_main(&C), ntree, *stack, storage.items[index]);
}

static wmOperatorStatus texture_layer_add_constant_mask(bContext *C,
                                                        const char *name,
                                                        const float value)
{
  bNodeTree *ntree = active_material_ntree(*C);
  if (ntree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNode *mask_stack_node = ensure_active_layer_mask_stack(*C, *ntree);
  if (mask_stack_node == nullptr) {
    return OPERATOR_CANCELLED;
  }
  const int index = mask_stack::add_layer(
      *CTX_data_main(C), *ntree, *mask_stack_node, name, value);
  set_active_layer(*ntree, *mask_stack_node, index);
  tree_changed(*C, *ntree, active_material(*C));
  return OPERATOR_FINISHED;
}

static wmOperatorStatus texture_layer_add_white_mask_exec(bContext *C, wmOperator * /*op*/)
{
  return texture_layer_add_constant_mask(C, "White", 1.0f);
}

static wmOperatorStatus texture_layer_add_black_mask_exec(bContext *C, wmOperator * /*op*/)
{
  return texture_layer_add_constant_mask(C, "Black", 0.0f);
}

static wmOperatorStatus texture_layer_add_paint_mask_exec(bContext *C, wmOperator * /*op*/)
{
  bNodeTree *ntree = active_material_ntree(*C);
  if (ntree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNode *mask_stack_node = ensure_active_layer_mask_stack(*C, *ntree);
  if (mask_stack_node == nullptr) {
    return OPERATOR_CANCELLED;
  }

  /* Source the mask from a fresh Image Texture (paintable once assigned). */
  bNode *image_node = bke::node_add_node(C, *ntree, "ShaderNodeTexImage"_ustr);
  if (image_node == nullptr) {
    return OPERATOR_CANCELLED;
  }
  /* The painted mask value is the image's Color, as for an Image Texture added
   * from the mask menu. */
  bNode *mask_source = nullptr;
  bNodeSocket *mask_out = mask_stack::find_source_output(
      *CTX_data_main(C), *ntree, *image_node, &mask_source);
  if (mask_out == nullptr) {
    /* Should not happen (an Image Texture always has a Color output); clean up the freshly
     * added node rather than leave it orphaned in the tree. */
    bke::node_remove_node(CTX_data_main(C), *ntree, *image_node, true);
    return OPERATOR_CANCELLED;
  }
  const int index = mask_stack::add_layer_from_source(
      *CTX_data_main(C), *ntree, *mask_stack_node, "Paint", *mask_source, *mask_out);
  /* Placed once the mask layer exists, so it lands in that layer's slot. */
  texture_stack::place_layer_source(*image_node, {ntree, mask_stack_node, index}, 320.0f);
  set_active_layer(*ntree, *mask_stack_node, index);
  tree_changed(*C, *ntree, active_material(*C));
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_add_white_mask(wmOperatorType *ot)
{
  ot->name = "Add White Mask";
  ot->description = "Add an opaque (1.0) mask layer to the active layer's mask stack";
  ot->idname = "MATERIAL_OT_texture_layer_add_white_mask";
  ot->exec = texture_layer_add_white_mask_exec;
  ot->poll = active_layer_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static void MATERIAL_OT_texture_layer_add_black_mask(wmOperatorType *ot)
{
  ot->name = "Add Black Mask";
  ot->description = "Add a fully transparent (0.0) mask layer to the active layer's mask stack";
  ot->idname = "MATERIAL_OT_texture_layer_add_black_mask";
  ot->exec = texture_layer_add_black_mask_exec;
  ot->poll = active_layer_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static void MATERIAL_OT_texture_layer_add_paint_mask(wmOperatorType *ot)
{
  ot->name = "Add Paint Mask";
  ot->description = "Add a paintable image-texture mask layer to the active layer's mask stack";
  ot->idname = "MATERIAL_OT_texture_layer_add_paint_mask";
  ot->exec = texture_layer_add_paint_mask_exec;
  ot->poll = active_layer_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Remove Mask Layer
 * \{ */

static bool active_mask_item_poll(bContext *C)
{
  if (!texture_layers_poll(*C)) {
    return false;
  }
  bNodeTree *ntree = active_material_ntree(*C);
  return ntree && active_mask(*ntree).has_value();
}

static wmOperatorStatus texture_layer_mask_remove_exec(bContext *C, wmOperator * /*op*/)
{
  bNodeTree *ntree = active_material_ntree(*C);
  if (ntree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  const std::optional<ActiveMask> mask = active_mask(*ntree);
  if (!mask) {
    return OPERATOR_CANCELLED;
  }
  remove_mask(*C, *ntree, *mask->stack, mask->index);
  tree_changed(*C, *ntree, active_material(*C));
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_mask_remove(wmOperatorType *ot)
{
  ot->name = "Remove Mask Layer";
  ot->description = "Remove the active mask layer from its Mask Stack";
  ot->idname = "MATERIAL_OT_texture_layer_mask_remove";

  ot->exec = texture_layer_mask_remove_exec;
  ot->poll = active_mask_item_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Group / Convert / Ungroup
 * \{ */

/* What a group layer does with the stack nested inside it: a generator group
 * contributes its result as a layer of its own, an adjustment group applies its
 * layers to the bundle accumulated below it (through a closure zone). */
enum class GroupKind {
  Generator,
  Adjustment,
};

/* The name of layer #index, used to label the nodes backing it. */
static std::string layer_item_name(const bNode &stack,
                                   const int index,
                                   const StringRefNull fallback)
{
  const char *name = layer_stack::storage(stack).items[index].name;
  return name ? name : std::string(fallback);
}

/* Create the nested Texture Layer Stack node for group layer #index of
 * #stack_node: labeled #label, placed in that layer's slot, and with one item
 * per entry of #item_names (top to bottom). The node's sockets are
 * materialized, so the caller can link to the new items right away. */
static bNode *create_nested_stack_node(bContext &C,
                                       bNodeTree &ntree,
                                       bNode &stack_node,
                                       const int index,
                                       const StringRefNull label,
                                       const Span<StringRefNull> item_names)
{
  bNode *nested = bke::node_add_node(&C, ntree, "ShaderNodeTextureLayerStack"_ustr);
  if (nested == nullptr) {
    return nullptr;
  }
  STRNCPY_UTF8(nested->label, label.c_str());
  texture_stack::place_layer_source(*nested, {&ntree, &stack_node, index}, 320.0f);
  for (const StringRefNull item_name : item_names) {
    nodes::socket_items::add_item_with_name<texture_stack::ItemsAccessor>(*nested,
                                                                          item_name.c_str());
  }
  BKE_ntree_update_tag_node_property(&ntree, nested);
  BKE_main_ensure_invariants(*CTX_data_main(&C), ntree.id);
  return nested;
}

/* Add an empty group layer of #kind, named #name, above the active layer. */
static wmOperatorStatus texture_layer_add_group(bContext *C,
                                                const StringRefNull name,
                                                const GroupKind kind)
{
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  Main &bmain = *CTX_data_main(C);
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;

  const int index = texture_stack::add_layer(bmain, ntree, stack, name.c_str());
  const std::string item_name = layer_item_name(stack, index, name);

  /* An adjustment group's nested stack starts with the base item that receives
   * the incoming (below) bundle through the closure zone; the user then adds
   * adjustment/generator layers above it. A generator group starts empty. */
  Vector<StringRefNull> item_names;
  if (kind == GroupKind::Adjustment) {
    item_names.append("Input");
  }
  bNode *nested = create_nested_stack_node(*C, ntree, stack, index, item_name, item_names);
  if (nested == nullptr) {
    return OPERATOR_CANCELLED;
  }

  if (kind == GroupKind::Adjustment) {
    if (!closure_zone::wrap_stack(bmain, ntree, *nested, stack, index)) {
      return OPERATOR_CANCELLED;
    }
  }
  else {
    texture_stack::connect_bundle({&ntree, &stack, index}, *nested, "Result");
  }
  set_active_layer(ntree, stack, index);
  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static wmOperatorStatus texture_layer_add_group_exec(bContext *C, wmOperator * /*op*/)
{
  return texture_layer_add_group(C, "Group", GroupKind::Generator);
}

static void MATERIAL_OT_texture_layer_add_group(wmOperatorType *ot)
{
  ot->name = "Add Texture Layer Group";
  ot->description = "Add a generator group layer that contains its own stack of layers";
  ot->idname = "MATERIAL_OT_texture_layer_add_group";
  ot->exec = texture_layer_add_group_exec;
  ot->poll = active_stack_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static wmOperatorStatus texture_layer_add_adjustment_group_exec(bContext *C, wmOperator * /*op*/)
{
  return texture_layer_add_group(C, "Adjustment Group", GroupKind::Adjustment);
}

static void MATERIAL_OT_texture_layer_add_adjustment_group(wmOperatorType *ot)
{
  ot->name = "Add Adjustment Group";
  ot->description =
      "Add an adjustment group that applies its own stack of adjustments to the layers below";
  ot->idname = "MATERIAL_OT_texture_layer_add_adjustment_group";
  ot->exec = texture_layer_add_adjustment_group_exec;
  ot->poll = active_stack_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static bool convert_to_group_poll(bContext *C)
{
  if (!texture_layers_poll(*C)) {
    return false;
  }
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return false;
  }
  const int index = ctx->layer_index;
  if (index == -1) {
    return false;
  }
  /* The layer must have a source (a bundle generator or a closure adjustment)
   * and not already be a group. A bundle layer becomes a generator group; a
   * closure (adjustment) layer becomes an adjustment group. */
  const texture_stack::StackLayer layer{ctx->ntree, ctx->stack, index};
  if (layer.nested_stack() != nullptr) {
    return false;
  }
  return layer.source() != nullptr;
}

static wmOperatorStatus texture_layer_convert_to_group_exec(bContext *C, wmOperator * /*op*/)
{
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;
  const NodeShaderLayerStack &storage = layer_stack::storage(stack);
  const int index = ctx->layer_index;
  if (index == -1) {
    return OPERATOR_CANCELLED;
  }
  const texture_stack::StackLayer layer{&ntree, &stack, index};
  if (layer.nested_stack() != nullptr) {
    return OPERATOR_CANCELLED;
  }

  /* Capture the current source edge feeding the layer's Bundle/Closure input,
   * before the tree updates below can rebuild the stack's sockets. */
  bNodeSocket *layer_socket = layer.bundle_socket();
  const bNodeLink *outer_link = layer_socket ? bke::node_find_incoming_link(ntree, *layer_socket) :
                                               nullptr;
  if (outer_link == nullptr) {
    return OPERATOR_CANCELLED;
  }
  bNode *source_node = outer_link->fromnode;
  bNodeSocket *source_socket = outer_link->fromsock;

  Main &bmain = *CTX_data_main(C);
  const std::string item_name = layer_item_name(stack, index, "Group");

  /* A closure (adjustment) layer becomes an adjustment group, a bundle layer a
   * generator group. */
  const GroupKind kind = (storage.items[index].item_type == SHADER_LAYER_STACK_ITEM_CLOSURE) ?
                             GroupKind::Adjustment :
                             GroupKind::Generator;

  /* The layer's own content becomes the nested stack's top layer. An adjustment
   * group also gets the base item that receives the group's incoming bundle
   * (the outer stack accumulated below) through its closure zone. */
  Vector<StringRefNull> item_names = {StringRefNull(item_name)};
  if (kind == GroupKind::Adjustment) {
    item_names.append("Input");
  }
  bNode *nested = create_nested_stack_node(*C, ntree, stack, index, item_name, item_names);
  if (nested == nullptr) {
    return OPERATOR_CANCELLED;
  }

  /* Move the layer's source onto the nested stack's top layer. For an
   * adjustment its type is inferred as a closure from the link. */
  bNodeSocket *outer_socket = layer.bundle_socket();
  if (const bNodeLink *link = outer_socket ? bke::node_find_incoming_link(ntree, *outer_socket) :
                                             nullptr)
  {
    bke::node_remove_link(&ntree, const_cast<bNodeLink &>(*link));
  }
  bNodeSocket *nested_item_socket = texture_stack::StackLayer{&ntree, nested, 0}.bundle_socket();
  if (nested_item_socket && source_node && source_socket) {
    bke::node_add_link(ntree, *source_node, *source_socket, *nested, *nested_item_socket);
  }

  if (kind == GroupKind::Adjustment) {
    /* Wrap the nested stack in the group's closure zone; the outer layer keeps
     * its blend mode, opacity and mask and is re-fed from Closure Output. */
    if (!closure_zone::wrap_stack(bmain, ntree, *nested, stack, index)) {
      return OPERATOR_CANCELLED;
    }
  }
  else {
    /* Feed the layer from the nested stack's Result instead. */
    bNodeSocket *result_out = bke::node_find_socket(*nested, SOCK_OUT, "Result"_ustr);
    if (outer_socket && result_out) {
      bke::node_add_link(ntree, *nested, *result_out, stack, *outer_socket);
    }
  }

  set_active_layer(ntree, stack, index);
  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_convert_to_group(wmOperatorType *ot)
{
  ot->name = "Make Group";
  ot->description = "Wrap the active layer in a Texture Layer Group";
  ot->idname = "MATERIAL_OT_texture_layer_convert_to_group";
  ot->exec = texture_layer_convert_to_group_exec;
  ot->poll = convert_to_group_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static bool ungroup_poll(bContext *C)
{
  if (!texture_layers_poll(*C)) {
    return false;
  }
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return false;
  }
  const int index = ctx->layer_index;
  if (index == -1) {
    return false;
  }
  return texture_stack::StackLayer{ctx->ntree, ctx->stack, index}.nested_stack() != nullptr;
}

static wmOperatorStatus texture_layer_ungroup_exec(bContext *C, wmOperator * /*op*/)
{
  Main &bmain = *CTX_data_main(C);
  const std::optional<ActiveStackContext> ctx = resolve_active_stack(*C);
  if (!ctx) {
    return OPERATOR_CANCELLED;
  }
  bNodeTree &ntree = *ctx->ntree;
  bNode &stack = *ctx->stack;
  const int insert_at = ctx->layer_index;
  if (insert_at == -1) {
    return OPERATOR_CANCELLED;
  }
  const texture_stack::StackLayer group_layer{&ntree, &stack, insert_at};
  bNode *nested = group_layer.nested_stack();
  if (nested == nullptr) {
    return OPERATOR_CANCELLED;
  }

  /* An adjustment group's nested stack has a hidden base "input" item, fed by
   * the group's Closure Input, that is dropped (not flattened). Capture the
   * group's closure zone nodes so they can be deleted once the layers move out;
   * the nested adjustment layers keep their own zones and are re-linked. */
  const bool adjustment_group = texture_stack::is_adjustment_group(ntree, *nested);
  bNode *group_closure_out = nullptr;
  bNode *group_closure_in = nullptr;
  if (adjustment_group) {
    group_closure_out = group_layer.source();
    if (group_closure_out && closure_zone::is_output_node(*group_closure_out)) {
      group_closure_in = closure_zone::input_node(ntree, *group_closure_out);
    }
  }

  /* Snapshot each nested item (top-to-bottom) so it can be recreated outside,
   * excluding an adjustment group's hidden base input. */
  Vector<texture_stack::LayerState> nested_states;
  {
    const NodeShaderLayerStack &nested_storage = layer_stack::storage(*nested);
    const int snapshot_num = adjustment_group ? nested_storage.items_num - 1 :
                                                nested_storage.items_num;
    for (const int i : IndexRange(snapshot_num)) {
      nested_states.append(texture_stack::capture_layer_state({&ntree, nested, i}));
    }
  }

  /* The group's own mask, captured so it can be re-attached to the top
   * flattened layer instead of being lost with the group. */
  int32_t outer_mask_node_id = 0;
  std::string outer_mask_socket;
  bNodeSocket *outer_mask_sock = group_layer.mask_socket();
  if (const bNodeLink *mask_link = outer_mask_sock ?
                                       bke::node_find_incoming_link(ntree, *outer_mask_sock) :
                                       nullptr)
  {
    if (mask_link->fromnode && mask_link->fromsock) {
      outer_mask_node_id = mask_link->fromnode->identifier;
      outer_mask_socket = mask_link->fromsock->identifier;
    }
  }

  /* Remove the group layer, the nested stack node, and (for an adjustment
   * group) the now-orphaned closure zone nodes. The group's mask source is left
   * in place for now (re-attached or cleaned up below). */
  layer_stack::remove_layer(bmain, ntree, stack, insert_at);
  bke::node_remove_node(&bmain, ntree, *nested, true);
  if (group_closure_out) {
    bke::node_remove_node(&bmain, ntree, *group_closure_out, true);
  }
  if (group_closure_in) {
    bke::node_remove_node(&bmain, ntree, *group_closure_in, true);
  }
  BKE_main_ensure_invariants(bmain, ntree.id);

  /* Recreate the nested items in the outer stack at the group's position,
   * preserving order (first nested item ends up at insert_at). */
  int first_new_index = -1;
  for (const int offset : nested_states.index_range()) {
    const texture_stack::LayerState &state = nested_states[offset];
    const int target = layer_stack::add_layer_at(
        bmain, ntree, stack, state.name.c_str(), insert_at + offset);
    texture_stack::restore_layer_state({&ntree, &stack, target}, state);
    if (first_new_index == -1) {
      first_new_index = target;
    }
  }
  if (first_new_index != -1) {
    set_active_layer(ntree, stack, first_new_index);
  }

  /* Preserve the group's mask by re-attaching it to the top flattened layer,
   * when that layer has no mask of its own and can take one. Otherwise the
   * now-orphaned mask source is removed so it doesn't dangle. */
  if (outer_mask_node_id != 0) {
    bNode *mask_from = bke::node_find_node_by_identifier(ntree, outer_mask_node_id);
    bNodeSocket *mask_in =
        (first_new_index != -1) ?
            texture_stack::StackLayer{&ntree, &stack, first_new_index}.mask_socket() :
            nullptr;
    bool reattached = false;
    if (mask_from && mask_in && bke::node_find_incoming_link(ntree, *mask_in) == nullptr) {
      if (bNodeSocket *mask_out = bke::node_find_socket(
              *mask_from, SOCK_OUT, UString(outer_mask_socket)))
      {
        bke::node_add_link(ntree, *mask_from, *mask_out, stack, *mask_in);
        reattached = true;
      }
    }
    if (!reattached && mask_from) {
      VectorSet<bNode *> orphan_deleted;
      texture_stack::delete_layer_source_recursive(bmain, ntree, mask_from, orphan_deleted);
    }
  }

  tree_changed(*C, ntree, ctx->material);
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_ungroup(wmOperatorType *ot)
{
  ot->name = "Ungroup Texture Layer";
  ot->description = "Flatten the active group's nested layers into the parent stack";
  ot->idname = "MATERIAL_OT_texture_layer_ungroup";
  ot->exec = texture_layer_ungroup_exec;
  ot->poll = ungroup_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Reparent (drag & drop reorder / move between stacks)
 * \{ */

static wmOperatorStatus texture_layer_reparent_exec(bContext *C, wmOperator *op)
{
  bNodeTree *ntree = active_material_ntree(*C);
  if (ntree == nullptr) {
    return OPERATOR_CANCELLED;
  }
  char from_name[MAX_NAME], to_name[MAX_NAME];
  RNA_string_get(op->ptr, "from_stack", from_name);
  RNA_string_get(op->ptr, "to_stack", to_name);
  const int from_index = RNA_int_get(op->ptr, "from_index");
  const int to_index_prop = RNA_int_get(op->ptr, "to_index");

  bNode *src = bke::node_find_node_by_name(*ntree, from_name);
  bNode *dst = bke::node_find_node_by_name(*ntree, to_name);
  if (src == nullptr || dst == nullptr) {
    return OPERATOR_CANCELLED;
  }
  /* Masks move between Mask Stacks, layers between Texture Layer Stacks —
   * never across the two kinds. */
  const bool is_mask_move = mask_stack::is_stack(*src);
  if (is_mask_move ? !mask_stack::is_stack(*dst) :
                     (!texture_stack::is_stack(*src) || !texture_stack::is_stack(*dst)))
  {
    return OPERATOR_CANCELLED;
  }
  NodeShaderLayerStack &src_storage = layer_stack::storage(*src);
  if (from_index < 0 || from_index >= src_storage.items_num) {
    return OPERATOR_CANCELLED;
  }

  if (src == dst) {
    const int count = src_storage.items_num;
    const int target = std::clamp(to_index_prop, 0, count - 1);
    if (target == from_index) {
      return OPERATOR_CANCELLED;
    }
    dna::array::move_index(src_storage.items, count, from_index, target);
    /* Highlight the moved row, matching the cross-stack path (dragging a row of a non-active
     * stack should activate it). */
    set_active_layer(*ntree, *src, target);
    /* The per-item sockets follow the item order and the base (last) layer
     * has no Opacity socket: rebuild the declaration. A mask on a layer moved
     * to the base is preserved: the base keeps a Mask socket. */
    BKE_ntree_update_tag_node_property(ntree, src);
    tree_changed(*C, *ntree, active_material(*C));
    return OPERATOR_FINISHED;
  }

  /* Cross-stack move: refuse dropping a group into its own nested tree. */
  bNode *src_nested = is_mask_move ?
                          nullptr :
                          texture_stack::StackLayer{ntree, src, from_index}.nested_stack();
  if (src_nested && texture_stack::is_stack(*src_nested) &&
      texture_stack::contains(*ntree, *src_nested, *dst))
  {
    return OPERATOR_CANCELLED;
  }

  const texture_stack::LayerState state = texture_stack::capture_layer_state(
      {ntree, src, from_index});

  const int target = layer_stack::add_layer_at(
      *CTX_data_main(C), *ntree, *dst, state.name.c_str(), to_index_prop);
  layer_stack::storage(*dst).active_index = target;
  texture_stack::restore_layer_state({ntree, dst, target}, state);

  /* Remove the source item (this also drops its now-dangling links). */
  layer_stack::remove_layer(*CTX_data_main(C), *ntree, *src, from_index);

  set_active_layer(*ntree, *dst, target);
  tree_changed(*C, *ntree, active_material(*C));
  return OPERATOR_FINISHED;
}

static void MATERIAL_OT_texture_layer_reparent(wmOperatorType *ot)
{
  ot->name = "Move Texture Layer";
  ot->description =
      "Move a texture layer between or within Texture Layer Stack nodes, or a mask "
      "layer between or within Mask Stack nodes";
  ot->idname = "MATERIAL_OT_texture_layer_reparent";
  ot->exec = texture_layer_reparent_exec;
  /* The source and target stacks come from the operator properties, so only a
   * material node tree is required (the active node may be anything). */
  ot->poll = active_material_ntree_poll;
  ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO | OPTYPE_INTERNAL;

  RNA_def_string(ot->srna, "from_stack", nullptr, MAX_NAME, "From Stack", "");
  RNA_def_int(ot->srna, "from_index", -1, -1, INT32_MAX, "From Index", "", -1, INT32_MAX);
  RNA_def_string(ot->srna, "to_stack", nullptr, MAX_NAME, "To Stack", "");
  RNA_def_int(ot->srna, "to_index", -1, -1, INT32_MAX, "To Index", "", -1, INT32_MAX);
}

/** \} */

void material_texture_layers_register()
{
  WM_operatortype_append(MATERIAL_OT_texture_layer_add);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_default);
  WM_operatortype_append(MATERIAL_OT_texture_layer_remove);
  WM_operatortype_append(MATERIAL_OT_texture_channel_toggle);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_fill);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_paint);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_white_mask);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_black_mask);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_paint_mask);
  WM_operatortype_append(MATERIAL_OT_texture_layer_mask_remove);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_group);
  WM_operatortype_append(MATERIAL_OT_texture_layer_add_adjustment_group);
  WM_operatortype_append(MATERIAL_OT_texture_layer_convert_to_group);
  WM_operatortype_append(MATERIAL_OT_texture_layer_ungroup);
  WM_operatortype_append(MATERIAL_OT_texture_layer_reparent);
}

}  // namespace blender::ed::render
