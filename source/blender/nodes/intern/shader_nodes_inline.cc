/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include <fmt/format.h>
#include <variant>

#include "BKE_compute_context_cache.hh"
#include "BKE_image.hh"
#include "BKE_lib_id.hh"
#include "BKE_node.hh"
#include "BKE_node_tree_zones.hh"
#include "BKE_type_conversions.hh"

#include "BLI_listbase.hh"
#include "BLI_math_vector_c.hh"
#include "BLI_stack.hh"
#include "BLI_string.hh"
#include "BLI_string_utf8.hh"

#include "DNA_image_types.h"
#include "DNA_material_types.h"
#include "DNA_node_types.h"

#include "NOD_layer_stack.hh"
#include "NOD_menu_value.hh"
#include "NOD_multi_function.hh"
#include "NOD_node_declaration.hh"
#include "NOD_node_in_compute_context.hh"
#include "NOD_shader_nodes_inline.hh"
#include "NOD_socket.hh"
#include "NOD_texture_channel.hh"
#include "NOD_texture_stack.hh"
#include "NOD_shader_nodes_multi_function.hh"

namespace blender::nodes {
namespace {

struct SocketValue;
struct BundleSocketValue;
using BundleSocketValuePtr = std::shared_ptr<BundleSocketValue>;

struct FallbackValue {};

/** This indicates that the value should be ignored when it is linked to an input socket. */
struct DanglingValue {};

struct NodeAndSocket {
  bNode *node = nullptr;
  bNodeSocket *socket = nullptr;
};

struct PrimitiveSocketValue {
  std::variant<int, float, bool, ColorGeometry4f, float3, MenuValue, std::string> value;

  const void *buffer() const
  {
    return std::visit([](auto &&value) -> const void * { return &value; }, value);
  }

  void *buffer()
  {
    return const_cast<void *>(const_cast<const PrimitiveSocketValue *>(this)->buffer());
  }

  static PrimitiveSocketValue from_value(const GPointer value)
  {
    const CPPType &type = *value.type();
    if (type.is<int>()) {
      return {*static_cast<const int *>(value.get())};
    }
    if (type.is<float>()) {
      return {*static_cast<const float *>(value.get())};
    }
    if (type.is<bool>()) {
      return {*static_cast<const bool *>(value.get())};
    }
    if (type.is<ColorGeometry4f>()) {
      return {*static_cast<const ColorGeometry4f *>(value.get())};
    }
    if (type.is<float3>()) {
      return {*static_cast<const float3 *>(value.get())};
    }
    if (type.is<MenuValue>()) {
      return {*static_cast<const MenuValue *>(value.get())};
    }
    if (type.is<std::string>()) {
      return {*static_cast<const std::string *>(value.get())};
    }
    BLI_assert_unreachable();
    return {};
  }
};

/** References an output socket in the generated node tree. */
struct LinkedSocketValue {
  bNode *node = nullptr;
  bNodeSocket *socket = nullptr;
};

/** References an input socket in the source node tree. */
struct InputSocketValue {
  const bNodeSocket *socket = nullptr;
};

struct ClosureZoneValue {
  const bke::bNodeTreeZone *zone = nullptr;
  const ComputeContext *closure_creation_context = nullptr;
};

struct MultiInputValue {
  Vector<SocketValue, 0> values;
};

struct SocketValue {
  /**
   * The value of an arbitrary socket value can have one of many different types. At a high level
   * it can either have a specific constant-folded value, or it references a socket that can't be
   * constant-folded.
   */
  std::variant<FallbackValue,
               DanglingValue,
               LinkedSocketValue,
               InputSocketValue,
               PrimitiveSocketValue,
               ClosureZoneValue,
               BundleSocketValuePtr,
               MultiInputValue>
      value;

  /** Try to get the value as a primitive value. */
  std::optional<PrimitiveSocketValue> to_primitive(const bke::bNodeSocketType &type) const
  {
    if (const auto *primitive_value = std::get_if<PrimitiveSocketValue>(&this->value)) {
      return *primitive_value;
    }
    if (const auto *input_socket_value = std::get_if<InputSocketValue>(&this->value)) {
      const bNodeSocket &socket = *input_socket_value->socket;
      BLI_assert(socket.type == type.type);
      if (!socket.runtime->declaration) {
        return std::nullopt;
      }
      if (socket.runtime->declaration->default_input_type != NODE_DEFAULT_INPUT_VALUE) {
        return std::nullopt;
      }
      switch (socket.typeinfo->type) {
        case SOCK_FLOAT:
          return PrimitiveSocketValue{socket.default_value_typed<bNodeSocketValueFloat>()->value};
        case SOCK_INT:
          return PrimitiveSocketValue{socket.default_value_typed<bNodeSocketValueInt>()->value};
        case SOCK_BOOLEAN:
          return PrimitiveSocketValue{
              bool(socket.default_value_typed<bNodeSocketValueBoolean>()->value)};
        case SOCK_VECTOR:
          return PrimitiveSocketValue{
              float3(socket.default_value_typed<bNodeSocketValueVector>()->value)};
        case SOCK_RGBA:
          return PrimitiveSocketValue{
              ColorGeometry4f(socket.default_value_typed<bNodeSocketValueRGBA>()->value)};
        case SOCK_MENU:
          return PrimitiveSocketValue{
              MenuValue(socket.default_value_typed<bNodeSocketValueMenu>()->value)};
        case SOCK_STRING:
          return PrimitiveSocketValue{
              std::string(socket.default_value_typed<bNodeSocketValueString>()->value)};
        default:
          return std::nullopt;
      }
    }
    if (std::get_if<FallbackValue>(&this->value)) {
      switch (type.type) {
        case SOCK_INT:
        case SOCK_BOOLEAN:
        case SOCK_VECTOR:
        case SOCK_RGBA:
        case SOCK_FLOAT:
        case SOCK_STRING:
        case SOCK_MENU:
          return PrimitiveSocketValue::from_value(
              {type.base_cpp_type, type.base_cpp_type->default_value()});
        default:
          return std::nullopt;
      }
    }
    return std::nullopt;
  }
};

struct BundleSocketValue {
  struct Item {
    std::string key;
    SocketValue value;
    const bke::bNodeSocketType *socket_type = nullptr;
  };

  Vector<Item> items;
};

struct PreservedZone {
  bNode *input_node = nullptr;
  bNode *output_node = nullptr;
};

/**
 * When a repeat zone is preserved, the corresponding #RepeatZoneComputeContext will use this
 * iteration. It should be a value that is not a valid iteration index. This allows optionally
 * inlining the repeat zone later without reusing the same context for different values.
 */
constexpr int preserved_repeat_zone_iteration = -1;

class ShaderNodesInliner {
 private:
  /** Cache for intermediate values used during the inline process. */
  ResourceScope scope_;
  /** The original tree the has to be inlined. */
  const bNodeTree &src_tree_;
  /** The tree where the inlined nodes will be added. */
  bNodeTree &dst_tree_;
  /** Parameters passed in by the caller. */
  InlineShaderNodeTreeParams &params_;
  /** Simplifies building the all the compute contexts for nodes in zones and groups. */
  bke::ComputeContextCache compute_context_cache_;
  /**
   * Stores compute context of the direct parent of each zone. In most cases, this is just the
   * parent compute context directly, except for closures.
   */
  Map<const ComputeContext *, const ComputeContext *> parent_zone_contexts_;
  /** Stores the computed value for each socket. The final value for each socket may be constant */
  Map<SocketInContext, SocketValue> value_by_socket_;
  /**
   * Remember zone nodes that have been copied to the destination so that they can be connected
   * again in the end.
   */
  Map<NodeInContext, PreservedZone> copied_zone_by_zone_output_node_;
  /** Sockets that still have to be evaluated. */
  Stack<SocketInContext> scheduled_sockets_stack_;
  /* Running composition state per #ShaderNodeTextureLayerStack invocation.
   * Closure-typed layers evaluate their body asynchronously, so composition
   * can span multiple handler invocations; tracking progress here ensures
   * each layer's Mix nodes are emitted exactly once. */
  struct StackCompositionState {
    /* Number of layers already composed, counting from the base upwards. */
    int layers_done = 0;
    BundleSocketValuePtr accumulator;
  };
  Map<NodeInContext, StackCompositionState> stack_states_;
  /* Closure-layer contexts whose closure-input values have been seeded. */
  Set<const ComputeContext *> seeded_closure_layer_contexts_;
  /**
   * Single-pass Image Texture node copies created when inlining a multi-pass
   * Image Texture node, deduplicated per (source node, pass name) so that a pass
   * socket and its synthetic Alpha share one copy (design phase 16).
   */
  Map<NodeInContext, Map<std::string, bNode *>> single_pass_image_copies_;
  /** Knows how to compute between different data types. */
  const bke::DataTypeConversions &data_type_conversions_;
  /** This is used to generate unique names and ids. */
  int dst_node_counter_ = 0;

  struct FailedToInlineRepeatZone {
    const ComputeContext *parent_ctx = nullptr;
    int output_node_id = 0;

    friend bool operator==(const FailedToInlineRepeatZone &a,
                           const FailedToInlineRepeatZone &b) = default;

    uint64_t hash() const
    {
      return get_default_hash(this->parent_ctx, this->output_node_id);
    }
  };

  /**
   * Sometimes the inliner first attempts to preserve a repeat zone without knowing whether it
   * will succeed. In the general case (with passed in bundles and closures) that cannot be
   * determined in advance with 100% certainty. Since preserved repeat zones are an important
   * optimization, it's not enough to only have them when we can "prove" beforehand that it will
   * work. This proof might not be possible or expensive in some cases where preserving the repeat
   * zone would be possible.
   *
   * To solve this, the inliner can "backtrack" an earlier decision to preserve a repeat zone.
   * Basically, when it notices later that there is a specific error in the preserved repeat zone,
   * it can switch to inlining it instead.
   *
   * The backtracking is not perfect currently, it could leave behind some nodes in the generated
   * tree, but those are not connected to the final output. Those could still be removed in a
   * separate pass.
   */
  Set<FailedToInlineRepeatZone> repeat_zones_to_force_inline_;

 public:
  ShaderNodesInliner(const bNodeTree &src_tree,
                     bNodeTree &dst_tree,
                     InlineShaderNodeTreeParams &params)
      : src_tree_(src_tree),
        dst_tree_(dst_tree),
        params_(params),
        data_type_conversions_(bke::get_implicit_type_conversions())
  {
  }

  bool do_inline()
  {
    src_tree_.ensure_topology_cache();
    if (src_tree_.has_available_link_cycle()) {
      return false;
    }

    const Vector<SocketInContext> final_output_sockets = this->find_final_output_sockets();

    /* Evaluation starts at the final output sockets which will request the evaluation of whether
     * sockets are linked to them. */
    for (const SocketInContext &socket : final_output_sockets) {
      this->schedule_socket(socket);
    }

    /* Evaluate until all scheduled sockets have a value. While evaluating a single socket, it may
     * either end up having a value, or request more other sockets that need to be evaluated first.
     *
     * This uses an explicit stack instead of recursion to avoid stack overflows which can easily
     * happen when there are long chains of nodes (or e.g. repeat zones with many iterations). */
    while (!scheduled_sockets_stack_.is_empty()) {
      const SocketInContext socket = scheduled_sockets_stack_.peek();
      const int old_stack_size = scheduled_sockets_stack_.size();

      this->handle_socket(socket);

      if (scheduled_sockets_stack_.size() == old_stack_size) {
        /* No additional dependencies were pushed, so this socket is fully handled and can be
         * popped from the stack. */
        BLI_assert(socket == scheduled_sockets_stack_.peek());
        scheduled_sockets_stack_.pop();
      }
    }

    /* Create actual output nodes. */
    Map<NodeInContext, bNode *> final_output_nodes;
    for (const SocketInContext &socket : final_output_sockets) {
      const NodeInContext src_node = socket.owner_node();
      bNode *copied_node = final_output_nodes.lookup_or_add_cb(src_node, [&]() {
        Map<const bNodeSocket *, bNodeSocket *> socket_map;
        bNode *copied_node = bke::node_copy_with_mapping(&dst_tree_,
                                                         *src_node.node,
                                                         this->node_copy_flag(),
                                                         std::nullopt,
                                                         this->get_next_node_identifier(),
                                                         socket_map);
        copied_node->parent = nullptr;
        return copied_node;
      });
      bNodeSocket *copied_socket = static_cast<bNodeSocket *>(
          BLI_findlink(&copied_node->inputs, socket.socket->index()));
      this->set_input_socket_value(
          src_node, *copied_node, *copied_socket, value_by_socket_.lookup(socket));
    }

    this->restore_zones_in_output_tree();
    this->position_nodes_in_output_tree();
    return true;
  }

  Vector<SocketInContext> find_final_output_sockets()
  {
    Vector<TreeInContext> trees;
    this->find_trees_potentially_containing_shader_outputs_recursive(nullptr, src_tree_, trees);

    auto get_engine_target = [](const bNode *output_node) {
      if (STR_ELEM(output_node->idname,
                   "ShaderNodeOutputMaterial",
                   "ShaderNodeOutputLight",
                   "ShaderNodeOutputWorld"))
      {
        return NodeShaderOutputTarget(output_node->custom1);
      }
      return SHD_OUTPUT_ALL;
    };

    Vector<SocketInContext> output_sockets;
    auto add_output_type = [&](const UString output_type) {
      for (const TreeInContext &tree : trees) {
        const bke::bNodeTreeZones &zones = *tree->zones();
        for (const bNode *node : tree->nodes_by_type(output_type)) {
          if (!ELEM(get_engine_target(node), SHD_OUTPUT_ALL, params_.target_engine_)) {
            continue;
          }
          const bke::bNodeTreeZone *zone = zones.get_zone_by_node(node->identifier);
          if (zone) {
            this->report_error({tree.context, node},
                               TIP_("Output node must not be in zone"),
                               NodeWarningType::Error);
            continue;
          }
          for (const bNodeSocket *socket : node->input_sockets()) {
            output_sockets.append({tree.context, socket});
          }
        }
      }
    };

    /* owner_id can be null for DefaultSurfaceNodeTree. */
    ID_Type tree_type = src_tree_.owner_id ? src_tree_.owner_id->id_type() : ID_MA;

    switch (tree_type) {
      case ID_MA:
        add_output_type("ShaderNodeOutputMaterial"_ustr);
        add_output_type("ShaderNodeOutputLight"_ustr);
        add_output_type("ShaderNodeOutputAOV"_ustr);
        break;
      case ID_WO:
        add_output_type("ShaderNodeOutputWorld"_ustr);
        add_output_type("ShaderNodeOutputAOV"_ustr);
        break;
      case ID_LA:
        add_output_type("ShaderNodeOutputLight"_ustr);
        break;
      case ID_LS:
        add_output_type("ShaderNodeOutputLineStyle"_ustr);
        break;
      default:
        BLI_assert_unreachable();
    }

    return output_sockets;
  }

  void find_trees_potentially_containing_shader_outputs_recursive(const ComputeContext *context,
                                                                  const bNodeTree &tree,
                                                                  Vector<TreeInContext> &r_trees)
  {
    const bke::bNodeTreeZones *zones = src_tree_.zones();
    if (!zones) {
      return;
    }
    if (tree.has_available_link_cycle()) {
      return;
    }
    r_trees.append({context, &tree});
    for (const bNode *group_node : tree.group_nodes()) {
      if (group_node->is_muted()) {
        continue;
      }
      const bNodeTree *group = id_cast<const bNodeTree *>(group_node->id);
      if (!group || ID_MISSING(&group->id)) {
        continue;
      }
      group->ensure_topology_cache();
      const bke::bNodeTreeZone *zone = zones->get_zone_by_node(group_node->identifier);
      if (zone) {
        /* Node groups in zones are ignored. */
        continue;
      }
      const ComputeContext &group_context = compute_context_cache_.for_group_node(
          context, group_node->identifier, &tree);
      this->find_trees_potentially_containing_shader_outputs_recursive(
          &group_context, *group, r_trees);
    }
  }

  void handle_socket(const SocketInContext &socket)
  {
    if (!socket->is_available()) {
      return;
    }
    if (value_by_socket_.contains(socket)) {
      /* The socket already has a value, so there is nothing to do. */
      return;
    }
    if (socket->is_input()) {
      this->handle_input_socket(socket);
    }
    else {
      this->handle_output_socket(socket);
    }
  }

  void handle_input_socket(const SocketInContext &socket)
  {
    if (socket->is_multi_input()) {
      this->handle_multi_input_socket(socket);
      return;
    }

    const bNodeLink *used_link = nullptr;
    for (const bNodeLink *link : socket->directly_linked_links()) {
      if (!link->is_used() || link->fromnode == nullptr || link->fromnode->is_undefined()) {
        continue;
      }
      used_link = link;
    }
    if (!used_link) {
      /* If there is no link on the input, use the value of the socket directly. */
      if (this->input_socket_may_have_dangling_value(socket)) {
        this->store_socket_value(socket, {DanglingValue{}});
      }
      else {
        this->store_socket_value(socket, {InputSocketValue{socket.socket}});
      }
      return;
    }

    const ComputeContext *from_context = this->get_link_source_context(*used_link, socket);
    const SocketInContext origin_socket = {from_context, used_link->fromsock};
    if (const auto *value = value_by_socket_.lookup_ptr(origin_socket)) {
      if (std::holds_alternative<DanglingValue>(value->value)) {
        if (this->input_socket_may_have_dangling_value(socket)) {
          this->store_socket_value(socket, {DanglingValue{}});
        }
        else {
          /* If the input value is dangling, use the value of the socket itself. */
          this->store_socket_value(socket, {InputSocketValue{socket.socket}});
        }
        return;
      }
      /* If the socket linked to the input has a value already, copy that value to the current
       * socket, potentially with an implicit conversion. */
      this->store_socket_value(socket,
                               this->handle_implicit_conversion(*value,
                                                                *used_link->fromsock->typeinfo,
                                                                *used_link->tosock->typeinfo));
      return;
    }
    /* If the origin socket does not have a value yet, only schedule it for evaluation for now. */
    this->schedule_socket(origin_socket);
  }

  void handle_multi_input_socket(const SocketInContext &socket)
  {
    bool all_links_ready = true;
    Vector<SocketValue, 0> values;
    for (const bNodeLink *link : socket->directly_linked_links()) {
      if (!link->is_used()) {
        continue;
      }
      const ComputeContext *from_context = this->get_link_source_context(*link, socket);
      const SocketInContext origin_socket = {from_context, link->fromsock};
      const SocketValue *value = value_by_socket_.lookup_ptr(origin_socket);
      if (!value) {
        this->schedule_socket(origin_socket);
        all_links_ready = false;
        continue;
      }
      values.append(*value);
    }
    if (!all_links_ready) {
      return;
    }
    this->store_socket_value(socket, {MultiInputValue{std::move(values)}});
  }

  /**
   * Generally, input values of a node should never be dangling because otherwise the node can't be
   * evaluated. However, if a node is never evaluated anyway, then its inputs can be dangling. This
   * allows the dangling-state to be properly forwarded through the node.
   */
  bool input_socket_may_have_dangling_value(const SocketInContext &socket)
  {
    BLI_assert(socket->is_input());
    const NodeInContext node = socket.owner_node();
    return node->is_reroute() || node->is_muted();
  }

  const ComputeContext *get_link_source_context(const bNodeLink &link,
                                                const SocketInContext &to_socket)
  {
    const bNodeTree &tree = to_socket->owner_tree();
    const bke::bNodeTreeZones *zones = tree.zones();
    if (!zones) {
      return nullptr;
    }
    const bke::bNodeTreeZone *to_zone = zones->get_zone_by_socket(*to_socket);
    const bke::bNodeTreeZone *from_zone = zones->get_zone_by_socket(*link.fromsock);
    const ComputeContext *context = to_socket.context;
    for (const bke::bNodeTreeZone *zone = to_zone; zone != from_zone; zone = zone->parent_zone) {
      context = parent_zone_contexts_.lookup(context);
    }
    return context;
  }

  void handle_output_socket(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    if (node->is_reroute()) {
      this->handle_output_socket__reroute(socket);
      return;
    }
    if (node->is_muted()) {
      if (!this->handle_output_socket__internal_links(socket)) {
        /* The output socket does not have a corresponding input, so the value is ignored. */
        this->store_socket_value_dangling(socket);
      }
      return;
    }
    if (node->is_group()) {
      this->handle_output_socket__group(socket);
      return;
    }
    if (node->is_group_input()) {
      this->handle_output_socket__group_input(socket);
      return;
    }
    if (node->is_type("GeometryNodeRepeatOutput"_ustr)) {
      if (this->should_preserve_repeat_zone_node(node.context, *node)) {
        this->handle_output_socket__preserved_repeat_output(socket);
        return;
      }
      this->handle_output_socket__repeat_output(socket);
      return;
    }
    if (node->is_type("GeometryNodeRepeatInput"_ustr)) {
      const auto *context = dynamic_cast<const bke::RepeatZoneComputeContext *>(socket.context);
      if (!context) {
        /* This socket is expected to be in the repeat zone, so it should have a context. */
        this->store_socket_value_fallback(socket);
        return;
      }
      if (context->iteration() == preserved_repeat_zone_iteration) {
        this->handle_output_socket__preserved_repeat_input(socket);
        return;
      }
      this->handle_output_socket__repeat_input(socket);
      return;
    }
    if (node->is_type("NodeClosureOutput"_ustr)) {
      this->handle_output_socket__closure_output(socket);
      return;
    }
    if (node->is_type("NodeClosureInput"_ustr)) {
      this->handle_output_socket__closure_input(socket);
      return;
    }
    if (node->is_type("NodeEvaluateClosure"_ustr)) {
      this->handle_output_socket__evaluate_closure(socket);
      return;
    }
    if (node->is_type("NodeCombineBundle"_ustr)) {
      this->handle_output_socket__combine_bundle(socket);
      return;
    }
    if (node->is_type("NodeSeparateBundle"_ustr)) {
      this->handle_output_socket__separate_bundle(socket);
      return;
    }
    if (node->is_type("GeometryNodeMenuSwitch"_ustr)) {
      this->handle_output_socket__menu_switch(socket);
      return;
    }
    if (node->is_type("GeometryNodeIndexSwitch"_ustr)) {
      this->handle_output_socket__index_switch(socket);
      return;
    }
    if (node->is_type("GeometryNodeSwitch"_ustr)) {
      this->handle_output_socket__switch(socket);
      return;
    }
    if (node->is_type("FunctionNodeInputMenu"_ustr)) {
      this->handle_output_socket__input_menu(socket);
      return;
    }
    if (node->is_type("NodeJoinBundle"_ustr)) {
      this->handle_output_socket__join_bundle(socket);
      return;
    }
    if (node->is_type("ShaderNodeTextureLayerStack"_ustr)) {
      this->handle_output_socket__texture_layer_stack(socket);
      return;
    }
    if (node->is_type("ShaderNodeMaskStack"_ustr)) {
      this->handle_output_socket__mask_stack(socket);
      return;
    }
    if (node->is_type("NodeImplicitConversion"_ustr)) {
      this->handle_output_socket__implicit_conversion(socket);
      return;
    }
    this->handle_output_socket__eval(socket);
  }

  void handle_output_socket__reroute(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const SocketInContext input_socket = node.input_socket(0);
    this->forward_value_or_schedule(socket, input_socket);
  }

  /* Returns whether the socket was handled. */
  [[nodiscard]] bool handle_output_socket__internal_links(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    for (const bNodeLink &internal_link : node->internal_links()) {
      if (internal_link.tosock == socket.socket) {
        const SocketInContext src_socket = {socket.context, internal_link.fromsock};
        if (src_socket->is_multi_input()) {
          const bNodeLink *src_link = nullptr;
          for (const bNodeLink *link : src_socket->directly_linked_links()) {
            if (link->is_used()) {
              src_link = link;
              break;
            }
          }
          if (!src_link) {
            return false;
          }
          const ComputeContext *from_context = this->get_link_source_context(*src_link,
                                                                             src_socket);
          const SocketInContext origin_socket = {from_context, src_link->fromsock};
          this->forward_value_or_schedule(socket, origin_socket);
          return true;
        }
        if (const SocketValue *value = value_by_socket_.lookup_ptr(src_socket)) {
          /* Pass the value of the internally linked input socket, with an implicit conversion if
           * necessary. */
          this->store_socket_value(
              socket,
              this->handle_implicit_conversion(
                  *value, *internal_link.fromsock->typeinfo, *internal_link.tosock->typeinfo));
          return true;
        }
        this->schedule_socket(src_socket);
        return true;
      }
    }
    return false;
  }

  void handle_output_socket__group(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const bNodeTree *group = reinterpret_cast<const bNodeTree *>(node->id);
    if (!group || ID_MISSING(&group->id)) {
      this->store_socket_value_fallback(socket);
      return;
    }
    if (socket.context && socket.context->parents_num() >= U.nodes_stack_limit) {
      this->store_socket_value_fallback(socket);
      this->report_error(node,
                         TIP_("Nodes stack limit reached (too many levels of nested nodes)"),
                         NodeWarningType::Error);
      return;
    }
    group->ensure_interface_cache();
    group->ensure_topology_cache();
    const bNode *group_output_node = group->group_output_node();
    if (!group_output_node) {
      this->store_socket_value_fallback(socket);
      return;
    }
    /* Get the value of an output of a group node by evaluating the corresponding output of the
     * node group. Since this socket is in a different tree, the compute context is different. */
    const ComputeContext &group_compute_context = compute_context_cache_.for_group_node(
        socket.context, node->identifier, &node->owner_tree());
    const SocketInContext group_output_socket_ctx = {
        &group_compute_context, &group_output_node->input_socket(socket->index())};
    this->forward_value_or_schedule(socket, group_output_socket_ctx);
  }

  void handle_output_socket__group_input(const SocketInContext &socket)
  {
    if (const auto *group_node_compute_context =
            dynamic_cast<const bke::GroupNodeComputeContext *>(socket.context))
    {
      /* Get the value of a group input from the corresponding input socket of the parent group
       * node. */
      const ComputeContext *parent_compute_context = group_node_compute_context->parent();
      const bNode *group_node = group_node_compute_context->node();
      BLI_assert(group_node);
      const bNodeSocket &group_node_input = group_node->input_socket(socket->index());
      const SocketInContext group_input_socket_ctx = {parent_compute_context, &group_node_input};
      this->forward_value_or_schedule(socket, group_input_socket_ctx);
      return;
    }
    this->store_socket_value_fallback(socket);
  }

  /**
   * Returns whether the given repeat zone should be preserved, which is an important optimization
   * for render engines that support it (currently only EEVEE).
   *
   * \note A repeat zone may first be preserved and can later change to be inlined anyway if it
   * failed. Also see #repeat_zones_to_force_inline_.
   */
  bool should_preserve_repeat_zone_node(const ComputeContext *outer_ctx,
                                        const bNode &repeat_zone_node) const
  {
    BLI_assert(repeat_zone_node.is_type("GeometryNodeRepeatOutput"_ustr) ||
               repeat_zone_node.is_type("GeometryNodeRepeatInput"_ustr));
    if (!params_.allow_preserving_repeat_zones) {
      return false;
    }
    const bNodeTree &tree = repeat_zone_node.owner_tree();
    const bke::bNodeTreeZones *zones = tree.zones();
    if (!zones) {
      return false;
    }
    const bke::bNodeTreeZone *zone = zones->get_zone_by_node(repeat_zone_node.identifier);
    if (!zone) {
      return false;
    }
    const bNode *repeat_zone_input_node = zone->input_node();
    const bNode *repeat_zone_output_node = zone->output_node();
    if (!repeat_zone_input_node || !repeat_zone_output_node) {
      return false;
    }
    if (repeat_zones_to_force_inline_.contains({outer_ctx, *zone->output_node_id})) {
      return false;
    }
    const auto &storage = *static_cast<const NodeGeometryRepeatOutput *>(
        repeat_zone_output_node->storage);
    for (const int i : IndexRange(storage.items_num)) {
      const NodeRepeatItem &item = storage.items[i];
      if (!ELEM(item.socket_type, SOCK_INT, SOCK_FLOAT, SOCK_BOOLEAN, SOCK_RGBA, SOCK_VECTOR)) {
        /* Repeat zones with more special types have to be inlined. */
        return false;
      }
    }
    return true;
  }

  bool is_inliner_evaluation_node(const bNode &node) const
  {
    const UString idname = node.typeinfo->idname;
    return ELEM(idname,
                "FunctionNodeFindInString"_ustr,
                "FunctionNodeFormatString"_ustr,
                "FunctionNodeInputSpecialCharacters"_ustr,
                "FunctionNodeInputString"_ustr,
                "FunctionNodeMatchString"_ustr,
                "FunctionNodeReplaceString"_ustr,
                "FunctionNodeReverseString"_ustr,
                "FunctionNodeSetStringCase"_ustr,
                "FunctionNodeSliceString"_ustr,
                "FunctionNodeStringLength"_ustr,
                "FunctionNodeStringToValue"_ustr,
                "FunctionNodeTrimString"_ustr,
                "FunctionNodeValueToString"_ustr);
  }

  void handle_output_socket__repeat_output(const SocketInContext &socket)
  {
    const bNode &repeat_output_node = socket->owner_node();
    const bNodeTree &tree = socket->owner_tree();

    const bke::bNodeTreeZones *zones = tree.zones();
    if (!zones) {
      this->store_socket_value_fallback(socket);
      return;
    }
    const bke::bNodeTreeZone *zone = zones->get_zone_by_node(repeat_output_node.identifier);
    if (!zone) {
      this->store_socket_value_fallback(socket);
      return;
    }
    const NodeInContext repeat_input_node = {socket.context, zone->input_node()};
    const SocketInContext iterations_input = repeat_input_node.input_socket(0);
    const SocketValue *iterations_socket_value = value_by_socket_.lookup_ptr(iterations_input);
    if (!iterations_socket_value) {
      /* The number of iterations is not known yet, so only schedule that socket for now. */
      this->schedule_socket(iterations_input);
      return;
    }
    const std::optional<PrimitiveSocketValue> iterations_value_opt =
        iterations_socket_value->to_primitive(*iterations_input->typeinfo);
    if (!iterations_value_opt) {
      this->add_dynamic_repeat_zone_iterations_error(repeat_input_node);
    }
    const int iterations = iterations_value_opt.has_value() ?
                               std::get<int>(iterations_value_opt->value) :
                               0;
    if (iterations <= 0) {
      /* If the number of iterations is zero, the values are copied directly from the repeat input
       * node. */
      const SocketInContext origin_socket = repeat_input_node.input_socket(1 + socket->index());
      this->forward_value_or_schedule(socket, origin_socket);
      return;
    }
    /* Otherwise, the value is copied from the output of the last iteration. */
    const ComputeContext &last_iteration_context = compute_context_cache_.for_repeat_zone(
        socket.context, repeat_output_node, iterations - 1);
    parent_zone_contexts_.add(&last_iteration_context, socket.context);
    const SocketInContext origin_socket = {&last_iteration_context,
                                           &repeat_output_node.input_socket(socket->index())};
    this->forward_value_or_schedule(socket, origin_socket);
  }

  void handle_output_socket__preserved_repeat_output(const SocketInContext &socket)
  {
    const bNodeTree &tree = socket->owner_tree();
    const NodeInContext repeat_output_node = socket.owner_node();
    const bke::bNodeTreeZones &zones = *tree.zones();
    const bke::bNodeTreeZone &zone = *zones.get_zone_by_node(repeat_output_node->identifier);
    const bNode &repeat_input_node = *zone.input_node();

    const ComputeContext *out_context = socket.context;
    const ComputeContext &in_context = compute_context_cache_.for_repeat_zone(
        socket.context, repeat_output_node->identifier, preserved_repeat_zone_iteration);
    parent_zone_contexts_.add(&in_context, out_context);

    const EnsureInputsResult ensured_inputs = this->ensure_node_inputs(
        {&in_context, &socket->owner_node()});
    if (ensured_inputs.has_missing_inputs) {
      /* The node can only be evaluated if all inputs values are known. */
      return;
    }
    const NodeInContext node = socket.owner_node();
    bNode &copied_node = this->handle_output_socket__eval_copy_node(
        *node, &in_context, out_context);
    PreservedZone &preserved_zone = copied_zone_by_zone_output_node_.lookup_or_add_default(
        repeat_output_node);
    preserved_zone.output_node = &copied_node;
    /* Ensure that the repeat input node is created as well. */
    this->schedule_socket({&in_context, &repeat_input_node.output_socket(0)});
  }

  void handle_output_socket__preserved_repeat_input(const SocketInContext &socket)
  {
    const ComputeContext &out_context = *socket.context;
    const ComputeContext *in_context = out_context.parent();
    BLI_assert(dynamic_cast<const bke::RepeatZoneComputeContext &>(out_context).iteration() ==
               preserved_repeat_zone_iteration);

    const EnsureInputsResult ensured_inputs = this->ensure_node_inputs(
        {in_context, &socket->owner_node()});
    if (ensured_inputs.has_missing_inputs) {
      /* The node can only be evaluated if all inputs values are known. */
      return;
    }
    const bNodeTree &tree = socket->owner_tree();
    const NodeInContext node = socket.owner_node();
    bNode &copied_node = this->handle_output_socket__eval_copy_node(
        *node, in_context, &out_context);
    const auto &storage = *static_cast<const NodeGeometryRepeatInput *>(node->storage);
    const NodeInContext repeat_output_node{in_context, tree.node_by_id(storage.output_node_id)};
    PreservedZone &preserved_zone = copied_zone_by_zone_output_node_.lookup_or_add_default(
        repeat_output_node);
    preserved_zone.input_node = &copied_node;
  }

  void add_dynamic_repeat_zone_iterations_error(const NodeInContext &repeat_input_node)
  {
    this->report_error(repeat_input_node,
                       TIP_("Iterations input has to be a constant value"),
                       NodeWarningType::Error);
  }

  void handle_output_socket__repeat_input(const SocketInContext &socket)
  {
    const bNode &repeat_input_node = socket->owner_node();
    const auto *repeat_zone_context = dynamic_cast<const bke::RepeatZoneComputeContext *>(
        socket.context);
    if (!repeat_zone_context) {
      this->store_socket_value_fallback(socket);
      return;
    }
    /* The index of the current iteration comes from the context. */
    const int iteration = repeat_zone_context->iteration();

    if (socket->index() == 0) {
      /* The first output is the current iteration index. */
      this->store_socket_value(socket, {PrimitiveSocketValue{iteration}});
      return;
    }

    if (iteration == 0) {
      /* In the first iteration, the values are copied from the corresponding input socket. */
      const SocketInContext origin_socket = {repeat_zone_context->parent(),
                                             &repeat_input_node.input_socket(socket->index())};
      this->forward_value_or_schedule(socket, origin_socket);
      return;
    }
    /* For later iterations, the values are copied from the corresponding output of the previous
     * iteration. */
    const bNode &repeat_output_node = *repeat_input_node.owner_tree().node_by_id(
        repeat_zone_context->output_node_id());
    const int previous_iteration = iteration - 1;
    const ComputeContext &previous_iteration_context = compute_context_cache_.for_repeat_zone(
        repeat_zone_context->parent(), repeat_output_node, previous_iteration);
    parent_zone_contexts_.add(&previous_iteration_context, repeat_zone_context->parent());
    const SocketInContext origin_socket = {&previous_iteration_context,
                                           &repeat_output_node.input_socket(socket->index() - 1)};
    this->forward_value_or_schedule(socket, origin_socket);
  }

  void handle_output_socket__closure_output(const SocketInContext &socket)
  {
    const bNode &node = socket->owner_node();
    const bke::bNodeTreeZones *zones = node.owner_tree().zones();
    if (!zones) {
      this->store_socket_value_fallback(socket);
      return;
    }
    const bke::bNodeTreeZone *zone = zones->get_zone_by_node(node.identifier);
    if (!zone) {
      this->store_socket_value_fallback(socket);
      return;
    }
    /* Just store a reference to the closure. */
    this->store_socket_value(socket, {ClosureZoneValue{zone, socket.context}});
  }

  void handle_output_socket__evaluate_closure(const SocketInContext &socket)
  {
    const NodeInContext evaluate_closure_node = socket.owner_node();
    const SocketInContext closure_input_socket = evaluate_closure_node.input_socket(0);
    const SocketValue *closure_input_value = value_by_socket_.lookup_ptr(closure_input_socket);
    if (!closure_input_value) {
      /* The closure to evaluate is not known yet, so schedule the closure input before it can be
       * evaluated. */
      this->schedule_socket(closure_input_socket);
      return;
    }
    const ClosureZoneValue *closure_zone_value = std::get_if<ClosureZoneValue>(
        &closure_input_value->value);
    if (!closure_zone_value) {
      /* If the closure is null, the node behaves as if it is muted. */
      if (!this->handle_output_socket__internal_links(socket)) {
        this->store_socket_value_fallback(socket);
      }
      return;
    }
    if (socket.context && socket.context->parents_num() >= U.nodes_stack_limit) {
      this->store_socket_value_fallback(socket);
      this->report_error(evaluate_closure_node,
                         TIP_("Nodes stack limit reached (too many levels of nested nodes)"),
                         NodeWarningType::Error);
      return;
    }

    const auto *evaluate_closure_storage = static_cast<const NodeEvaluateClosure *>(
        evaluate_closure_node->storage);
    const bNode &closure_output_node = *closure_zone_value->zone->output_node();
    const auto &closure_storage = *static_cast<const NodeClosureOutput *>(
        closure_output_node.storage);
    const StringRef key = evaluate_closure_storage->output_items.items[socket->index()].name;

    const bNodeTree &closure_tree = closure_output_node.owner_tree();
    const ClosureSourceLocation closure_source_location{
        &closure_tree,
        closure_output_node.identifier,
        closure_zone_value->closure_creation_context ?
            closure_zone_value->closure_creation_context->hash() :
            ComputeContextHash{},
        closure_zone_value->closure_creation_context};
    const bke::EvaluateClosureComputeContext &closure_eval_context =
        compute_context_cache_.for_evaluate_closure(socket.context,
                                                    evaluate_closure_node->identifier,
                                                    &socket->owner_tree(),
                                                    closure_source_location);
    parent_zone_contexts_.add(&closure_eval_context, closure_zone_value->closure_creation_context);

    for (const int i : IndexRange(closure_storage.output_items.items_num)) {
      const NodeClosureOutputItem &item = closure_storage.output_items.items[i];
      if (key != item.name) {
        continue;
      }
      /* Get the value of the output by evaluating the corresponding output in the closure zone. */
      const SocketInContext origin_socket = {&closure_eval_context,
                                             &closure_output_node.input_socket(i)};
      this->forward_value_or_schedule(socket, origin_socket);
      return;
    }
    this->store_socket_value_fallback(socket);
  }

  void handle_output_socket__closure_input(const SocketInContext &socket)
  {
    const bNode &closure_input_node = socket->owner_node();
    const auto *closure_eval_context = dynamic_cast<const bke::EvaluateClosureComputeContext *>(
        socket.context);
    if (!closure_eval_context) {
      this->store_socket_value_fallback(socket);
      return;
    }
    const bNode &closure_output_node = *closure_input_node.owner_tree().node_by_id(
        closure_eval_context->closure_source_location()->closure_output_node_id);
    const NodeInContext closure_eval_node = {closure_eval_context->parent(),
                                             closure_eval_context->node()};

    const auto &closure_storage = *static_cast<const NodeClosureOutput *>(
        closure_output_node.storage);
    const auto &eval_closure_storage = *static_cast<const NodeEvaluateClosure *>(
        closure_eval_node->storage);

    const StringRef key = closure_storage.input_items.items[socket->index()].name;
    for (const int i : IndexRange(eval_closure_storage.input_items.items_num)) {
      const NodeEvaluateClosureInputItem &item = eval_closure_storage.input_items.items[i];
      if (key != item.name) {
        continue;
      }
      /* The input of a closure zone gets its value from the corresponding input of the Evaluate
       * Closure node that evaluates it. */
      const SocketInContext origin_socket = closure_eval_node.input_socket(i + 1);
      this->forward_value_or_schedule(socket, origin_socket);
      return;
    }
    this->store_socket_value_fallback(socket);
  }

  void handle_output_socket__combine_bundle(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeCombineBundle *>(node->storage);

    bool all_inputs_available = true;
    for (const bNodeSocket *input_socket : node->input_sockets()) {
      const SocketInContext input_socket_ctx = {socket.context, input_socket};
      if (!value_by_socket_.lookup_ptr(input_socket_ctx)) {
        this->schedule_socket(input_socket_ctx);
        all_inputs_available = false;
      }
    }
    if (!all_inputs_available) {
      /* Can't create the bundle yet. Wait until all inputs are available. */
      return;
    }
    /* Build the actual bundle socket value from the input values. */
    auto bundle_value = std::make_shared<BundleSocketValue>();
    for (const int i : IndexRange(storage.items_num)) {
      const SocketInContext input_socket = node.input_socket(i);
      const NodeCombineBundleItem &item = storage.items[i];
      const StringRef key = item.name;
      const auto &socket_value = value_by_socket_.lookup(input_socket);
      bundle_value->items.append({key, socket_value, input_socket->typeinfo});
    }
    this->store_socket_value(socket, {bundle_value});
  }

  void handle_output_socket__separate_bundle(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeSeparateBundle *>(node->storage);

    const SocketInContext input_socket = node.input_socket(0);
    const SocketValue *socket_value = value_by_socket_.lookup_ptr(input_socket);
    if (!socket_value) {
      /* The input bundle is not known yet, so schedule it for now. */
      this->schedule_socket(input_socket);
      return;
    }
    const auto *bundle_value_ptr = std::get_if<BundleSocketValuePtr>(&socket_value->value);
    if (!bundle_value_ptr) {
      /* The bundle is empty. Use the fallback value. */
      this->store_socket_value_fallback(socket);
      return;
    }
    const BundleSocketValue &bundle_value = **bundle_value_ptr;

    const StringRef key = storage.items[socket->index()].name;
    for (const BundleSocketValue::Item &item : bundle_value.items) {
      if (key != item.key) {
        continue;
      }
      /* Extract the value from the bundle. */
      const SocketValue converted_value = this->handle_implicit_conversion(
          item.value, *item.socket_type, *socket->typeinfo);
      this->store_socket_value(socket, converted_value);
      return;
    }
    /* The bundle does not contain the requested key, so use the fallback value. */
    this->store_socket_value_fallback(socket);
  }

  void handle_output_socket__join_bundle(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const SocketInContext input_socket = node.input_socket(0);
    const SocketValue *socket_value = value_by_socket_.lookup_ptr(input_socket);
    if (!socket_value) {
      /* The input bundles are not known yet, so schedule them for now. */
      this->schedule_socket(input_socket);
      return;
    }
    const auto &multi_input_value = *std::get_if<MultiInputValue>(&socket_value->value);
    if (multi_input_value.values.is_empty()) {
      /* The input is empty, so use the fallback value. */
      this->store_socket_value_fallback(socket);
      return;
    }

    Set<StringRef> existing_keys;
    auto joined_bundle = std::make_shared<BundleSocketValue>();
    for (const SocketValue &value : multi_input_value.values) {
      const auto *bundle_value = std::get_if<BundleSocketValuePtr>(&value.value);
      if (!bundle_value || !*bundle_value) {
        /* Ignore invalid values. */
        continue;
      }
      for (const BundleSocketValue::Item &item : (*bundle_value)->items) {
        if (existing_keys.add(item.key)) {
          joined_bundle->items.append(item);
        }
      }
    }
    this->store_socket_value(socket, {BundleSocketValuePtr{joined_bundle}});
  }

  /* Emit a uniform-factor Mix node of the given data type mixing #b over
   * #a and return its Result output. #blend_type applies to color mixes.
   * When #a_type / #b_type are given, the inputs are implicitly converted
   * to the Mix node's socket type first. */
  SocketValue emit_mix(const NodeInContext &owner_node,
                       const eNodeSocketDatatype type,
                       const int blend_type,
                       const SocketValue &a,
                       const SocketValue &b,
                       const SocketValue &factor,
                       const bke::bNodeSocketType *a_type = nullptr,
                       const bke::bNodeSocketType *b_type = nullptr)
  {
    bNode *mix_node = this->add_node("ShaderNodeMix"_ustr);
    auto &mix_storage = *static_cast<NodeShaderMix *>(mix_node->storage);
    mix_storage.data_type = type;
    mix_storage.factor_mode = NODE_MIX_MODE_UNIFORM;
    mix_storage.clamp_factor = 0;
    mix_storage.clamp_result = 0;
    mix_storage.blend_type = (type == SOCK_RGBA) ? blend_type : MA_RAMP_BLEND;
    if (mix_node->typeinfo->updatefunc) {
      mix_node->typeinfo->updatefunc(&dst_tree_, mix_node);
    }

    UString a_identifier;
    UString b_identifier;
    UString result_identifier;
    switch (type) {
      case SOCK_FLOAT:
        a_identifier = "A_Float"_ustr;
        b_identifier = "B_Float"_ustr;
        result_identifier = "Result_Float"_ustr;
        break;
      case SOCK_VECTOR:
        a_identifier = "A_Vector"_ustr;
        b_identifier = "B_Vector"_ustr;
        result_identifier = "Result_Vector"_ustr;
        break;
      case SOCK_RGBA:
        a_identifier = "A_Color"_ustr;
        b_identifier = "B_Color"_ustr;
        result_identifier = "Result_Color"_ustr;
        break;
      default:
        BLI_assert_unreachable();
    }

    bNodeSocket *mix_factor = bke::node_find_socket(*mix_node, SOCK_IN, "Factor_Float"_ustr);
    bNodeSocket *mix_a = bke::node_find_socket(*mix_node, SOCK_IN, a_identifier);
    bNodeSocket *mix_b = bke::node_find_socket(*mix_node, SOCK_IN, b_identifier);
    bNodeSocket *mix_result = bke::node_find_socket(*mix_node, SOCK_OUT, result_identifier);

    const SocketValue a_converted = a_type ? this->handle_implicit_conversion(
                                                 a, *a_type, *mix_a->typeinfo) :
                                             a;
    const SocketValue b_converted = b_type ? this->handle_implicit_conversion(
                                                 b, *b_type, *mix_b->typeinfo) :
                                             b;

    this->set_input_socket_value(owner_node, *mix_node, *mix_a, a_converted);
    this->set_input_socket_value(owner_node, *mix_node, *mix_b, b_converted);
    this->set_input_socket_value(owner_node, *mix_node, *mix_factor, factor);

    return SocketValue{LinkedSocketValue{mix_node, mix_result}};
  }

  /* Blend two bundles channel-by-channel, with the #top bundle composited
   * over the #base bundle using #blend_factor as the per-channel weight and
   * #blend_type as the color blend mode. Returns the resulting bundle as a
   * BundleSocketValue. Channels present on only one side pass through. */
  BundleSocketValuePtr blend_layer_bundles(const NodeInContext &owner_node,
                                           const BundleSocketValue &base,
                                           const BundleSocketValue &top,
                                           const SocketValue &blend_factor,
                                           const int blend_type)
  {
    Map<StringRef, const BundleSocketValue::Item *> top_items_by_key;
    for (const BundleSocketValue::Item &item : top.items) {
      top_items_by_key.add(item.key, &item);
    }

    auto result_bundle = std::make_shared<BundleSocketValue>();
    Set<StringRef> handled_keys;

    for (const BundleSocketValue::Item &base_item : base.items) {
      handled_keys.add(base_item.key);
      const BundleSocketValue::Item *const *top_lookup = top_items_by_key.lookup_ptr(
          base_item.key);
      if (!top_lookup) {
        result_bundle->items.append(base_item);
        continue;
      }
      const BundleSocketValue::Item &top_item = **top_lookup;

      const eNodeSocketDatatype type = eNodeSocketDatatype(base_item.socket_type->type);
      if (!ELEM(type, SOCK_FLOAT, SOCK_VECTOR, SOCK_RGBA)) {
        result_bundle->items.append(base_item);
        continue;
      }

      const SocketValue mixed = this->emit_mix(owner_node,
                                               type,
                                               blend_type,
                                               base_item.value,
                                               top_item.value,
                                               blend_factor,
                                               base_item.socket_type,
                                               top_item.socket_type);
      const LinkedSocketValue &linked = std::get<LinkedSocketValue>(mixed.value);
      result_bundle->items.append({base_item.key, mixed, linked.socket->typeinfo});
    }

    /* Pass through items only in the top bundle. */
    for (const BundleSocketValue::Item &top_item : top.items) {
      if (!handled_keys.contains(top_item.key)) {
        result_bundle->items.append(top_item);
      }
    }

    return BundleSocketValuePtr{result_bundle};
  }

  /* Ensure all input socket values of #node are computed. Returns true if
   * any were missing and had to be scheduled — the caller should return and
   * wait to be re-run. */
  bool schedule_missing_inputs(const NodeInContext &node)
  {
    bool missing_input = false;
    for (const bNodeSocket *bsock : node->input_sockets()) {
      const SocketInContext s = {node.context, bsock};
      if (!value_by_socket_.lookup_ptr(s)) {
        this->schedule_socket(s);
        missing_input = true;
      }
    }
    return missing_input;
  }

  /* Evaluate a closure-typed layer: run the closure body with #accumulator
   * as its bundle input, in a compute context unique to this stack/layer pair
   * so multiple closure layers do not collide. Returns false when the body
   * had to be scheduled — the caller must return and resume when the handler
   * re-runs. On success #r_result holds the closure's output bundle, or
   * null when the closure is unusable and the layer should pass through. */
  bool evaluate_closure_layer(const NodeInContext &stack_node,
                              const NodeShaderLayerStackItem &item,
                              const ClosureZoneValue &closure_value,
                              const BundleSocketValuePtr &accumulator,
                              BundleSocketValuePtr *r_result)
  {
    *r_result = {};
    const bNode *closure_output_node = closure_value.zone->output_node();
    const bNode *closure_input_node = closure_value.zone->input_node();
    if (!closure_output_node || !closure_input_node) {
      return true;
    }
    const auto &closure_storage = *static_cast<const NodeClosureOutput *>(
        closure_output_node->storage);

    /* The closure is expected to take and produce a bundle; further items are
     * ignored and receive fallback values. */
    int bundle_in_index = -1;
    for (const int i : IndexRange(closure_storage.input_items.items_num)) {
      if (closure_storage.input_items.items[i].socket_type == SOCK_BUNDLE) {
        bundle_in_index = i;
        break;
      }
    }
    int bundle_out_index = -1;
    for (const int i : IndexRange(closure_storage.output_items.items_num)) {
      if (closure_storage.output_items.items[i].socket_type == SOCK_BUNDLE) {
        bundle_out_index = i;
        break;
      }
    }
    if (bundle_out_index == -1) {
      params_.r_error_messages.append(
          {&*stack_node, TIP_("Adjustment closure has no bundle output")});
      return true;
    }

    const ClosureSourceLocation source_location{
        &closure_output_node->owner_tree(),
        closure_output_node->identifier,
        closure_value.closure_creation_context ? closure_value.closure_creation_context->hash() :
                                                 ComputeContextHash{},
        closure_value.closure_creation_context};
    const bke::ClosureLayerComputeContext &ctx = compute_context_cache_.for_closure_layer(
        stack_node.context,
        stack_node->identifier,
        item.identifier,
        &stack_node->owner_tree(),
        source_location);

    if (ctx.is_recursive()) {
      params_.r_error_messages.append(
          {&*stack_node, TIP_("Recursive closures are not supported")});
      return true;
    }

    if (seeded_closure_layer_contexts_.add(&ctx)) {
      parent_zone_contexts_.add(&ctx, closure_value.closure_creation_context);
      /* Seed every output of the zone's input node before anything inside the
       * closure body is scheduled: the bundle input receives the running
       * accumulator, everything else its fallback. This must cover all
       * outputs so #handle_output_socket__closure_input never runs in this
       * context — that handler expects an Evaluate Closure node as the
       * context node, but here it is the stack node. */
      for (const int i : closure_input_node->output_sockets().index_range()) {
        const SocketInContext out_socket{&ctx, &closure_input_node->output_socket(i)};
        if (i == bundle_in_index) {
          BundleSocketValuePtr input_bundle = accumulator;
          if (!input_bundle) {
            input_bundle = BundleSocketValuePtr{std::make_shared<BundleSocketValue>()};
          }
          this->store_socket_value(out_socket, {input_bundle});
        }
        else {
          this->store_socket_value_fallback(out_socket);
        }
      }
    }

    const SocketInContext body_result{&ctx, &closure_output_node->input_socket(bundle_out_index)};
    const SocketValue *value = value_by_socket_.lookup_ptr(body_result);
    if (!value) {
      this->schedule_socket(body_result);
      return false;
    }
    if (const auto *bundle = std::get_if<BundleSocketValuePtr>(&value->value)) {
      *r_result = *bundle;
    }
    return true;
  }

  /* Build the plain (pre-stack) value bundle for a Texture Layer Stack: the
   * BSDF's input default for each channel the base layer carries. Null when
   * the stack is not a root (its Result does not feed the BSDF's Separate
   * Bundle), reaches no BSDF, or none of the base channels map to a fillable
   * BSDF input. Only a root stack has a meaningful plain value: nested group
   * stacks contribute a bundle, not the material's pre-stack channel values. */
  BundleSocketValuePtr build_stack_plain_bundle(const NodeInContext &node,
                                                const BundleSocketValue *base_bundle)
  {
    if (base_bundle == nullptr) {
      return nullptr;
    }
    bNode &stack_bnode = const_cast<bNode &>(*node.node);
    /* Only for a root stack, whose Result feeds a Separate Bundle. */
    bool is_root = false;
    if (bNodeSocket *result = bke::node_find_socket(stack_bnode, SOCK_OUT, "Result"_ustr)) {
      for (const bNodeSocket *target : result->directly_linked_sockets()) {
        if (texture_channel::is_separate_bundle(target->owner_node())) {
          is_root = true;
          break;
        }
      }
    }
    if (!is_root) {
      return nullptr;
    }
    bNode *bsdf = texture_stack::bsdf_for(stack_bnode.owner_tree(), stack_bnode);
    if (bsdf == nullptr) {
      return nullptr;
    }
    auto plain = std::make_shared<BundleSocketValue>();
    for (const BundleSocketValue::Item &base_item : base_bundle->items) {
      for (bNodeSocket &sock : bsdf->inputs) {
        if (base_item.key != sock.name || !texture_channel::is_fillable_input(sock)) {
          continue;
        }
        std::optional<PrimitiveSocketValue> prim;
        switch (sock.type) {
          case SOCK_FLOAT:
            prim = PrimitiveSocketValue{sock.default_value_typed<bNodeSocketValueFloat>()->value};
            break;
          case SOCK_VECTOR:
            prim = PrimitiveSocketValue{
                float3(sock.default_value_typed<bNodeSocketValueVector>()->value)};
            break;
          case SOCK_RGBA:
            prim = PrimitiveSocketValue{
                ColorGeometry4f(sock.default_value_typed<bNodeSocketValueRGBA>()->value)};
            break;
          default:
            break;
        }
        if (prim) {
          plain->items.append({base_item.key, SocketValue{*prim}, base_item.socket_type});
        }
        break;
      }
    }
    if (plain->items.is_empty()) {
      return nullptr;
    }
    return BundleSocketValuePtr{plain};
  }

  /* Whether #value is a constant 1.0 (a primitive or an unlinked socket's
   * default), used to skip the base layer's plain-value blend when its mask has
   * no effect. */
  static bool socket_value_is_one(const SocketValue &value)
  {
    if (const auto *prim = std::get_if<PrimitiveSocketValue>(&value.value)) {
      if (const float *f = std::get_if<float>(&prim->value)) {
        return *f == 1.0f;
      }
    }
    if (const auto *isv = std::get_if<InputSocketValue>(&value.value)) {
      if (isv->socket->type == SOCK_FLOAT) {
        return isv->socket->default_value_typed<bNodeSocketValueFloat>()->value == 1.0f;
      }
    }
    return false;
  }

  void handle_output_socket__texture_layer_stack(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const NodeShaderLayerStack &storage = nodes::layer_stack::storage(*node);
    const int items_num = storage.items_num;

    if (items_num == 0) {
      this->store_socket_value_fallback(socket);
      return;
    }

    /* Sockets are declared per item: for items 0..N-2 the order is
     * Layer, Opacity, Mask. The last item (N-1, base layer) only declares
     * Layer. */
    if (this->schedule_missing_inputs(node)) {
      return;
    }

    /* Gather per-item socket values from the declaration order. Muted items
     * have their bundle dropped so they pass the accumulator through unchanged. */
    struct Layer {
      const BundleSocketValue *bundle = nullptr;
      std::optional<ClosureZoneValue> closure;
      SocketValue opacity{};
      SocketValue mask{};
      int blend_type = MA_RAMP_BLEND;
    };
    Vector<Layer> layers(items_num);
    int sock_index = 0;
    for (const int i : IndexRange(items_num)) {
      const SocketInContext layer_sock = node.input_socket(sock_index++);
      const SocketValue &layer_value = value_by_socket_.lookup(layer_sock);
      const bool is_muted = (storage.items[i].flag & SHADER_LAYER_STACK_ITEM_MUTED) != 0;
      if (!is_muted) {
        if (const auto *bp = std::get_if<BundleSocketValuePtr>(&layer_value.value)) {
          if (*bp) {
            layers[i].bundle = bp->get();
          }
        }
        else if (const auto *cz = std::get_if<ClosureZoneValue>(&layer_value.value)) {
          if (cz->zone) {
            layers[i].closure = *cz;
          }
        }
      }
      layers[i].blend_type = storage.items[i].blend_type;
      const bool is_base = (i == items_num - 1);
      if (!is_base) {
        const SocketInContext opacity_sock = node.input_socket(sock_index++);
        const SocketInContext mask_sock = node.input_socket(sock_index++);
        layers[i].opacity = value_by_socket_.lookup(opacity_sock);
        layers[i].mask = value_by_socket_.lookup(mask_sock);
      }
      else {
        /* The base layer has a Mask but no Opacity socket. */
        const SocketInContext mask_sock = node.input_socket(sock_index++);
        layers[i].mask = value_by_socket_.lookup(mask_sock);
      }
    }

    /* Plain (pre-stack) values for the base channels, read from the BSDF's
     * input defaults. The base layer blends over these by its Mask, so a
     * masked-out or disabled base channel reveals the plain value instead of
     * the socket type's default (black). */
    const BundleSocketValuePtr plain_bundle = this->build_stack_plain_bundle(node,
                                                                             layers.last().bundle);

    /* Compose the final result. Start from the base layer at the bottom,
     * then composite each layer above it using opacity * mask as the
     * per-channel weight.
     * Closure-typed layers evaluate their body with the running accumulator
     * as its bundle input, which may require scheduling the body and
     * resuming later — the state records progress across handler runs. */
    StackCompositionState &state = stack_states_.lookup_or_add_default(node);
    while (state.layers_done < items_num) {
      const int i = items_num - 1 - state.layers_done;
      const Layer &layer = layers[i];

      const BundleSocketValue *layer_bundle = layer.bundle;
      BundleSocketValuePtr closure_result;
      if (layer.closure) {
        if (!this->evaluate_closure_layer(
                node, storage.items[i], *layer.closure, state.accumulator, &closure_result))
        {
          /* The closure body was scheduled; resume on the next handler run. */
          return;
        }
        layer_bundle = closure_result.get();
      }

      state.layers_done++;
      if (!layer_bundle) {
        /* This layer is empty: pass the accumulator through. */
        continue;
      }
      /* Drop the channels the user disabled on this layer, so it only
       * contributes to the enabled ones. */
      BundleSocketValuePtr filtered_holder;
      if (storage.items[i].disabled_channels_num > 0) {
        auto filtered = std::make_shared<BundleSocketValue>();
        for (const BundleSocketValue::Item &bundle_item : layer_bundle->items) {
          if (!nodes::layer_stack::channel_disabled(storage.items[i], bundle_item.key)) {
            filtered->items.append(bundle_item);
          }
        }
        filtered_holder = BundleSocketValuePtr{filtered};
        layer_bundle = filtered_holder.get();
      }
      if (!state.accumulator) {
        /* This layer becomes the start of the stack. When it is the base layer
         * and the plain bundle is known, reveal the plain pre-stack value where
         * the base does not fully cover it. With a constant 1.0 mask the base
         * fully covers its channels, so only the channels it drops (disabled or
         * not carried) fall back to the plain value — no Mix nodes needed.
         * Otherwise blend the base over the plain values by its mask. */
        if (i == items_num - 1 && plain_bundle) {
          if (socket_value_is_one(layer.mask)) {
            auto merged = std::make_shared<BundleSocketValue>(*layer_bundle);
            Set<StringRef> keys;
            for (const BundleSocketValue::Item &it : merged->items) {
              keys.add(it.key);
            }
            for (const BundleSocketValue::Item &plain_item : plain_bundle->items) {
              if (!keys.contains(plain_item.key)) {
                merged->items.append(plain_item);
              }
            }
            state.accumulator = BundleSocketValuePtr{merged};
          }
          else {
            state.accumulator = this->blend_layer_bundles(
                node, *plain_bundle, *layer_bundle, layer.mask, MA_RAMP_BLEND);
          }
        }
        else {
          state.accumulator = BundleSocketValuePtr{
              std::make_shared<BundleSocketValue>(*layer_bundle)};
        }
        continue;
      }
      const SocketValue blend_factor_value = this->multiply_float_socket_values(
          node, layer.opacity, layer.mask);
      state.accumulator = this->blend_layer_bundles(
          node, *state.accumulator, *layer_bundle, blend_factor_value, layer.blend_type);
    }

    const BundleSocketValuePtr accumulator = state.accumulator;
    stack_states_.remove(node);
    if (!accumulator) {
      this->store_socket_value_fallback(socket);
      return;
    }
    this->store_socket_value(socket, {accumulator});
  }

  void handle_output_socket__mask_stack(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const NodeShaderLayerStack &storage = nodes::layer_stack::storage(*node);
    const int items_num = storage.items_num;

    if (items_num == 0) {
      /* An empty Mask Stack contributes no mask. Like the all-muted case below, it must act as
       * fully unmasked (1.0), not the socket default of 0.0 which would hide the consuming
       * layer entirely. */
      this->store_socket_value(socket, {PrimitiveSocketValue{1.0f}});
      return;
    }

    /* Sockets are declared per item, each item declaring Mask then Opacity. */
    if (this->schedule_missing_inputs(node)) {
      return;
    }

    struct Layer {
      SocketValue value{};
      bool has_value = false;
      SocketValue opacity{};
      int blend_type = MA_RAMP_BLEND;
    };
    Vector<Layer> layers(items_num);
    int sock_index = 0;
    for (const int i : IndexRange(items_num)) {
      const SocketInContext mask_sock = node.input_socket(sock_index++);
      const bool is_muted = (storage.items[i].flag & SHADER_LAYER_STACK_ITEM_MUTED) != 0;
      if (!is_muted) {
        layers[i].value = value_by_socket_.lookup(mask_sock);
        layers[i].has_value = true;
      }
      layers[i].blend_type = storage.items[i].blend_type;
      const SocketInContext opacity_sock = node.input_socket(sock_index++);
      layers[i].opacity = value_by_socket_.lookup(opacity_sock);
    }

    /* Start at the base mask, then composite each layer above using opacity
     * as the per-mask weight and the layer's blend_type as the operation.
     * The base's own opacity blends it between fully unmasked (1.0) and its
     * mask value, so the last mask can be faded out like the ones above it. */
    SocketValue accumulator{};
    bool have_acc = false;
    if (layers.last().has_value) {
      const SocketValue unmasked{PrimitiveSocketValue{1.0f}};
      accumulator = this->blend_mask_floats(
          node, unmasked, layers.last().value, layers.last().opacity, MA_RAMP_BLEND);
      have_acc = true;
    }

    for (int i = items_num - 2; i >= 0; --i) {
      const Layer &layer = layers[i];
      if (!layer.has_value) {
        continue;
      }
      if (!have_acc) {
        accumulator = layer.value;
        have_acc = true;
        continue;
      }
      accumulator = this->blend_mask_floats(
          node, accumulator, layer.value, layer.opacity, layer.blend_type);
    }

    if (!have_acc) {
      /* Every mask layer is muted, so the stack contributes no mask. A muted
       * mask must act as if it were not there, which for the consuming layer
       * means fully unmasked (1.0), not the socket default of 0.0 (which would
       * hide the layer entirely). */
      this->store_socket_value(socket, {PrimitiveSocketValue{1.0f}});
      return;
    }
    this->store_socket_value(socket, accumulator);
  }

  /* Composite #top onto #base using #factor and #blend_type. Generic
   * formula:
   *     blended = blend_op(base, top)   (blend_op is identity for MIX)
   *     result  = mix(base, blended, factor)
   * which expands the four blend modes the design calls out (Mix, Add,
   * Multiply, Subtract) into 0–1 Math nodes plus 1 Mix(Float) node. Other
   * MA_RAMP_* modes fall back to plain Mix — they don't have a clean float
   * meaning and the user can model them via dedicated math nodes when needed. */
  SocketValue blend_mask_floats(const NodeInContext &owner_node,
                                const SocketValue &base,
                                const SocketValue &top,
                                const SocketValue &factor,
                                const int blend_type)
  {
    SocketValue blended = top;
    if (blend_type == MA_RAMP_ADD || blend_type == MA_RAMP_SUB || blend_type == MA_RAMP_MULT ||
        blend_type == MA_RAMP_DIV)
    {
      bNode *math_node = this->add_node("ShaderNodeMath"_ustr);
      switch (blend_type) {
        case MA_RAMP_ADD:
          math_node->custom1 = NODE_MATH_ADD;
          break;
        case MA_RAMP_SUB:
          math_node->custom1 = NODE_MATH_SUBTRACT;
          break;
        case MA_RAMP_MULT:
          math_node->custom1 = NODE_MATH_MULTIPLY;
          break;
        case MA_RAMP_DIV:
          math_node->custom1 = NODE_MATH_DIVIDE;
          break;
      }
      bNodeSocket *in_a = static_cast<bNodeSocket *>(math_node->inputs.first);
      bNodeSocket *in_b = in_a->next;
      bNodeSocket *out = static_cast<bNodeSocket *>(math_node->outputs.first);
      this->set_input_socket_value(owner_node, *math_node, *in_a, base);
      this->set_input_socket_value(owner_node, *math_node, *in_b, top);
      blended = SocketValue{LinkedSocketValue{math_node, out}};
    }

    return this->emit_mix(owner_node, SOCK_FLOAT, MA_RAMP_BLEND, base, blended, factor);
  }

  /* Multiply two float-valued socket values. Collapses to a literal when both
   * sides are known constants (or defaults on unlinked sockets), skips the
   * multiplication entirely on identity (×1) or zero (×0) shortcuts, and
   * otherwise inserts a single Math Multiply node. */
  SocketValue multiply_float_socket_values(const NodeInContext &owner_node,
                                           const SocketValue &a,
                                           const SocketValue &b)
  {
    auto as_float = [](const SocketValue &v) -> std::optional<float> {
      if (const auto *prim = std::get_if<PrimitiveSocketValue>(&v.value)) {
        if (const float *f = std::get_if<float>(&prim->value)) {
          return *f;
        }
      }
      if (const auto *isv = std::get_if<InputSocketValue>(&v.value)) {
        if (isv->socket->type == SOCK_FLOAT) {
          return isv->socket->default_value_typed<bNodeSocketValueFloat>()->value;
        }
      }
      return std::nullopt;
    };

    const std::optional<float> a_f = as_float(a);
    const std::optional<float> b_f = as_float(b);
    if (a_f && b_f) {
      return SocketValue{PrimitiveSocketValue{*a_f * *b_f}};
    }
    if ((a_f && *a_f == 0.0f) || (b_f && *b_f == 0.0f)) {
      return SocketValue{PrimitiveSocketValue{0.0f}};
    }
    if (a_f && *a_f == 1.0f) {
      return b;
    }
    if (b_f && *b_f == 1.0f) {
      return a;
    }

    bNode *math_node = this->add_node("ShaderNodeMath"_ustr);
    math_node->custom1 = NODE_MATH_MULTIPLY;
    bNodeSocket *in_a = static_cast<bNodeSocket *>(math_node->inputs.first);
    bNodeSocket *in_b = in_a->next;
    bNodeSocket *out = static_cast<bNodeSocket *>(math_node->outputs.first);
    this->set_input_socket_value(owner_node, *math_node, *in_a, a);
    this->set_input_socket_value(owner_node, *math_node, *in_b, b);
    return SocketValue{LinkedSocketValue{math_node, out}};
  }

  void handle_output_socket__menu_switch(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeMenuSwitch *>(node->storage);

    const SocketInContext menu_input = node.input_socket(0);
    const SocketValue *menu_socket_value = value_by_socket_.lookup_ptr(menu_input);
    if (!menu_socket_value) {
      /* The menu value is not known yet, so schedule it for now. */
      this->schedule_socket(menu_input);
      return;
    }

    const std::optional<PrimitiveSocketValue> menu_value_opt = menu_socket_value->to_primitive(
        *menu_input->typeinfo);
    if (!menu_value_opt) {
      /* This limitation may be lifted in the future. Menu Switch nodes could be supported natively
       * by render engines or we convert them to a bunch of mix nodes. */
      this->store_socket_value_fallback(socket);
      this->report_required_constant_input_or_backtrack(
          node, TIP_("Menu value has to be a constant value"));
      return;
    }
    const MenuValue menu_value = std::get<MenuValue>(menu_value_opt->value);
    /* Find the selected item index. */
    std::optional<int> selected_index;
    for (const int item_i : IndexRange(storage.enum_definition.items_num)) {
      const NodeEnumItem &item = storage.enum_definition.items_array[item_i];
      if (MenuValue(item.identifier) == menu_value) {
        selected_index = item_i;
        break;
      }
    }
    if (!selected_index.has_value()) {
      /* The input value does not exist in the menu. */
      this->store_socket_value_fallback(socket);
      return;
    }
    if (socket->index() == 0) {
      /* Handle forwarding the selected value. */
      this->forward_value_or_schedule(socket, node.input_socket(*selected_index + 1));
      return;
    }
    /* Set the value of the mask output. */
    const bool is_selected = selected_index == socket->index() - 1;
    this->store_socket_value(socket, {PrimitiveSocketValue{is_selected}});
  }

  struct MixNodeInfo {
    bNode *node = nullptr;
    bNodeSocket *factor_in = nullptr;
    bNodeSocket *a_in = nullptr;
    bNodeSocket *b_in = nullptr;
    bNodeSocket *result_out = nullptr;
  };

  MixNodeInfo create_mix_node(const eNodeSocketDatatype socket_type)
  {
    MixNodeInfo mix;
    switch (socket_type) {
      case SOCK_FLOAT:
      case SOCK_VECTOR:
      case SOCK_RGBA: {
        /* The ShaderNodeMix node uses different socket identifiers based on the data type, so find
         * the correct socket based on the name instead. */
        auto find_socket_by_name =
            [](bNode &node, const eNodeSocketInOut in_out, const StringRef name) -> bNodeSocket * {
          for (bNodeSocket &socket : in_out == SOCK_IN ? node.inputs : node.outputs) {
            if (socket.is_available()) {
              if (socket.name == name) {
                return &socket;
              }
            }
          }
          BLI_assert_unreachable();
          return nullptr;
        };

        mix.node = this->add_node("ShaderNodeMix"_ustr);
        auto &mix_storage = *static_cast<NodeShaderMix *>(mix.node->storage);
        mix_storage.data_type = socket_type;
        mix_storage.clamp_result = false;
        /* Use clamping of the mix node to avoid the need for a separate clamp node. */
        mix_storage.clamp_factor = true;
        /* This makes the right sockets for the given data type available. */
        mix.node->typeinfo->updatefunc(&dst_tree_, mix.node);
        mix.factor_in = find_socket_by_name(*mix.node, SOCK_IN, "Factor");
        mix.a_in = find_socket_by_name(*mix.node, SOCK_IN, "A");
        mix.b_in = find_socket_by_name(*mix.node, SOCK_IN, "B");
        mix.result_out = find_socket_by_name(*mix.node, SOCK_OUT, "Result");
        break;
      }
      case SOCK_SHADER: {
        mix.node = this->add_node("ShaderNodeMixShader"_ustr);
        mix.factor_in = static_cast<bNodeSocket *>(mix.node->inputs.first);
        mix.a_in = mix.factor_in->next;
        mix.b_in = mix.a_in->next;
        mix.result_out = static_cast<bNodeSocket *>(mix.node->outputs.first);
        break;
      }
      default: {
        BLI_assert_unreachable();
      }
    }
    return mix;
  }

  /**
   * Returns the socket type that's used internally to emulate switch nodes. Not all types are
   * natively supported by render engines.
   */
  std::optional<eNodeSocketDatatype> get_internal_mix_socket_type(
      const eNodeSocketDatatype socket_type) const
  {
    switch (socket_type) {
      case SOCK_FLOAT:
      case SOCK_INT:
      case SOCK_BOOLEAN:
        return SOCK_FLOAT;
      case SOCK_VECTOR:
        return SOCK_VECTOR;
      case SOCK_RGBA:
        return SOCK_RGBA;
      case SOCK_SHADER:
        return SOCK_SHADER;
      default:
        return std::nullopt;
    }
  }

  void handle_output_socket__index_switch(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeIndexSwitch *>(node->storage);

    if (storage.items_num == 0) {
      this->store_socket_value_fallback(socket);
      return;
    }

    const SocketInContext index_input = node.input_socket(0);
    const SocketValue *index_input_value = value_by_socket_.lookup_ptr(index_input);
    if (!index_input_value) {
      /* The index is not known yet, so schedule it for now. */
      this->schedule_socket(index_input);
      return;
    }

    if (const std::optional<PrimitiveSocketValue> primitive_index_value_opt =
            index_input_value->to_primitive(*index_input->typeinfo))
    {
      const int index = std::get<int>(primitive_index_value_opt->value);
      if (index < 0 || index >= storage.items_num) {
        this->store_socket_value_fallback(socket);
        this->report_error(node,
                           fmt::format("{}: {}", TIP_("Index out of range"), index),
                           NodeWarningType::Error);
        return;
      }
      this->forward_value_or_schedule(socket, node.input_socket(index + 1));
      return;
    }

    /* If the index is not a constant value, we immitate the index switch node using a chain of
     * mix nodes. This allows renderers using the Index Switch node with rendering backends which
     * don't support it natively. */
    const std::optional<eNodeSocketDatatype> internal_mix_type =
        this->get_internal_mix_socket_type(storage.data_type);
    if (!internal_mix_type) {
      this->report_required_constant_input_or_backtrack(
          node, TIP_("Index must be a constant value for this data type"));
      this->store_socket_value_fallback(socket);
      return;
    }

    const EnsureInputsResult ensured_inputs = this->ensure_node_inputs(node);
    if (ensured_inputs.has_missing_inputs) {
      /* Wait until all inputs values are available. */
      return;
    }

    /* Use a truncate node to turn the input value into an int. */
    bNode &truncate_math_node = *this->add_node("ShaderNodeMath"_ustr);
    truncate_math_node.custom1 = NODE_MATH_TRUNC;
    bNodeSocket &truncate_input = *static_cast<bNodeSocket *>(truncate_math_node.inputs.first);
    bNodeSocket &truncate_output = *static_cast<bNodeSocket *>(truncate_math_node.outputs.first);
    this->set_input_socket_value(node, truncate_math_node, truncate_input, *index_input_value);

    bNode *prev_mix = nullptr;
    bNodeSocket *prev_mix_result = nullptr;
    for (const int i : IndexRange(storage.items_num + 1)) {
      /* Use a math node to turn the index into a factor for the mix node.*/
      const int index_to_factor_offset = 1 - i;

      bNode *factor_node = nullptr;
      bNodeSocket *factor_out = nullptr;
      if (index_to_factor_offset == 0) {
        /* No need to add an extra math node here that just computes +0. */
        factor_node = &truncate_math_node;
        factor_out = &truncate_output;
      }
      else {
        bNode &add_math_node = *this->add_node("ShaderNodeMath"_ustr);
        add_math_node.custom1 = NODE_MATH_ADD;
        bNodeSocket &add_in_1 = *static_cast<bNodeSocket *>(add_math_node.inputs.first);
        bNodeSocket &add_in_2 = *add_in_1.next;
        bke::node_add_link(
            dst_tree_, truncate_math_node, truncate_output, add_math_node, add_in_1);
        static_cast<bNodeSocketValueFloat *>(add_in_2.default_value)->value =
            index_to_factor_offset;
        factor_node = &add_math_node;
        factor_out = static_cast<bNodeSocket *>(add_math_node.outputs.first);
      }

      const MixNodeInfo mix = this->create_mix_node(*internal_mix_type);
      bke::node_add_link(dst_tree_, *factor_node, *factor_out, *mix.node, *mix.factor_in);
      if (i == 0) {
        this->set_input_socket_value(node, *mix.node, *mix.a_in, {FallbackValue{}});
      }
      else {
        bke::node_add_link(dst_tree_, *prev_mix, *prev_mix_result, *mix.node, *mix.a_in);
      }
      if (i < storage.items_num) {
        const SocketInContext input_socket = node.input_socket(i + 1);
        const SocketValue &input_value = value_by_socket_.lookup(input_socket);
        const SocketValue converted_value = this->handle_implicit_conversion(
            input_value, *input_socket->typeinfo, *mix.b_in->typeinfo);
        this->set_input_socket_value(node, *mix.node, *mix.b_in, converted_value);
      }
      else {
        this->set_input_socket_value(node, *mix.node, *mix.b_in, {FallbackValue{}});
      }

      prev_mix = mix.node;
      prev_mix_result = mix.result_out;
    }

    this->store_socket_value(socket, {LinkedSocketValue{prev_mix, prev_mix_result}});
  }

  void handle_output_socket__switch(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeSwitch *>(node->storage);

    const SocketInContext switch_input = node.input_socket(0);
    const SocketValue *switch_input_value = value_by_socket_.lookup_ptr(switch_input);
    if (!switch_input_value) {
      /* The switch condition is not known yet, so schedule it for now. */
      this->schedule_socket(switch_input);
      return;
    }

    if (const std::optional<PrimitiveSocketValue> primitive_switch_value_opt =
            switch_input_value->to_primitive(*switch_input->typeinfo))
    {
      const bool condition = std::get<bool>(primitive_switch_value_opt->value);
      this->forward_value_or_schedule(socket, node.input_socket(condition ? 2 : 1));
      return;
    }

    /* If the switch condition is not a constant value, imitate the Switch node using a mix node.
     * This allows using the Switch node with rendering backends which don't support it natively.
     */
    const std::optional<eNodeSocketDatatype> internal_mix_type =
        this->get_internal_mix_socket_type(storage.input_type);
    if (!internal_mix_type) {
      this->report_required_constant_input_or_backtrack(
          node, TIP_("Switch must be a constant value for this data type"));
      this->store_socket_value_fallback(socket);
      return;
    }

    const EnsureInputsResult ensured_inputs = this->ensure_node_inputs(node);
    if (ensured_inputs.has_missing_inputs) {
      /* Wait until all inputs values are available. */
      return;
    }

    /* Convert the condition to either 0 or 1 so that it can be used as factor input without
     * accidentally mixing between the two input values. */
    bNode &to_bool_math_node = *this->add_node("ShaderNodeMath"_ustr);
    to_bool_math_node.custom1 = NODE_MATH_GREATER_THAN;
    bNodeSocket &to_bool_in_1 = *static_cast<bNodeSocket *>(to_bool_math_node.inputs.first);
    bNodeSocket &to_bool_in_2 = *to_bool_in_1.next;
    bNodeSocket &to_bool_out = *static_cast<bNodeSocket *>(to_bool_math_node.outputs.first);
    this->set_input_socket_value(node, to_bool_math_node, to_bool_in_1, *switch_input_value);
    static_cast<bNodeSocketValueFloat *>(to_bool_in_2.default_value)->value = 0.0f;

    const MixNodeInfo mix = this->create_mix_node(*internal_mix_type);
    bke::node_add_link(dst_tree_, to_bool_math_node, to_bool_out, *mix.node, *mix.factor_in);

    const SocketInContext false_input = node.input_socket(1);
    const SocketInContext true_input = node.input_socket(2);
    const SocketValue &false_value = value_by_socket_.lookup(false_input);
    const SocketValue &true_value = value_by_socket_.lookup(true_input);
    this->set_input_socket_value(node,
                                 *mix.node,
                                 *mix.a_in,
                                 this->handle_implicit_conversion(
                                     false_value, *false_input->typeinfo, *mix.a_in->typeinfo));
    this->set_input_socket_value(
        node,
        *mix.node,
        *mix.b_in,
        this->handle_implicit_conversion(true_value, *true_input->typeinfo, *mix.b_in->typeinfo));

    this->store_socket_value(socket, {LinkedSocketValue{mix.node, mix.result_out}});
  }

  void handle_output_socket__input_menu(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const auto &storage = *static_cast<const NodeInputMenu *>(node->storage);
    SocketInContext output_socket = node.output_socket(0);
    this->store_socket_value(output_socket,
                             {PrimitiveSocketValue::from_value(
                                 {output_socket->typeinfo->base_cpp_type, &storage.value})});
  }

  void handle_output_socket__implicit_conversion(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();

    const SocketInContext input_socket = node.input_socket(0);
    const SocketValue *socket_value = value_by_socket_.lookup_ptr(input_socket);
    if (!socket_value) {
      /* The input bundle is not known yet, so schedule it for now. */
      this->schedule_socket(input_socket);
      return;
    }
    const SocketValue converted_value = this->handle_implicit_conversion(
        *socket_value, *socket->typeinfo, *socket->typeinfo);
    this->store_socket_value(socket, converted_value);
  }

  /**
   * Evaluate a node to compute the value of the given output socket. This may also compute all the
   * other outputs of the node.
   */
  void handle_output_socket__eval(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const EnsureInputsResult ensured_inputs = this->ensure_node_inputs(node);
    if (ensured_inputs.has_missing_inputs) {
      /* The node can only be evaluated if all inputs values are known. */
      return;
    }
    const bke::bNodeType &node_type = *node->typeinfo;
    if (node_type.build_multi_function && ensured_inputs.all_inputs_primitive) {
      /* Do constant folding. */
      this->handle_output_socket__eval_multi_function(node);
      return;
    }
    if (is_multi_pass_image_node(node)) {
      /* Expand a multi-pass Image Texture node into single-pass copies, one per
       * used pass, so the render engine only ever sees ordinary single-pass
       * Image Texture nodes (design phase 16). */
      this->handle_output_socket__multi_pass_image(socket);
      return;
    }
    /* Check if this node is supported by the renderer. If not, either try inlining it harder to
     * eliminate the node, or report an error and skip it. */
    if (this->is_inliner_evaluation_node(*node)) {
      this->report_required_constant_input_or_backtrack(
          node, TIP_("Node requires constant input values"));
      this->store_socket_value_fallback(socket);
      return;
    }
    /* The node can't be constant-folded. So copy it to the destination tree instead. */
    this->handle_output_socket__eval_copy_node(*node, node.context, node.context);
  }

  /**
   * Can be called when an input socket that only supports constant values got a non-constant one.
   * This may either be caused by an issue in the original tree, or because a repeat-zone has been
   * attempted to preserve that can't be preserved. In the latter case, no error is reported yet.
   * Instead, the repeat zone is force-inlined which generally either fixes the issue or makes it
   * obvious that the issue is in fact in the original tree (and then this function will be called
   * again and the error is actually reported).
   */
  void report_required_constant_input_or_backtrack(const NodeInContext node,
                                                   const StringRef message)
  {
    /* Before reporting an error, attempt to inline a surrounding repeat zone. */
    for (const ComputeContext *ctx = node.context; ctx; ctx = ctx->parent()) {
      if (const auto *repeat_zone_ctx = dynamic_cast<const bke::RepeatZoneComputeContext *>(ctx)) {
        if (repeat_zone_ctx->iteration() == preserved_repeat_zone_iteration) {
          repeat_zones_to_force_inline_.add(
              {repeat_zone_ctx->parent(), repeat_zone_ctx->output_node_id()});
          return;
        }
      }
    }
    this->report_error(node, message, NodeWarningType::Error);
  }

  /**
   * The multi-layer catalog (#Image.layers) is built lazily on first image
   * acquire and is not persisted in the blend file, so a freshly-loaded image
   * has an empty catalog until something acquires it. In a background render
   * nothing does before inlining runs, and the pass lookup below would then fail
   * and wrongly fall back to a single-pass node. Force the catalog to build, the
   * same way the node declaration's prepare_image() does.
   */
  static void ensure_multilayer_catalog(Image &image, const ImageUser &iuser)
  {
    if (!BLI_listbase_is_empty(&image.layers)) {
      return;
    }
    ImageUser local_iuser = iuser;
    local_iuser.framenr = BKE_image_sequence_guess_offset(&image);
    ImBuf *ibuf = BKE_image_acquire_ibuf(&image, &local_iuser, nullptr);
    BKE_image_release_ibuf(&image, ibuf, nullptr);
  }

  /**
   * True for an Image Texture node assigned a multi-layer image, which declares
   * one output socket per pass instead of the fixed Color/Alpha pair. Mirrors
   * the gate used by the declaration (see node_shader_tex_image.cc): both must
   * agree on which nodes are multi-pass.
   */
  static bool is_multi_pass_image_node(const NodeInContext &node)
  {
    if (!node->is_type("ShaderNodeTexImage"_ustr)) {
      return false;
    }
    Image *image = reinterpret_cast<Image *>(node->id);
    if (image == nullptr) {
      return false;
    }
    const NodeTexImage *tex = static_cast<const NodeTexImage *>(node->storage);
    if (tex == nullptr || (tex->flag & SHD_TEX_IMAGE_SINGLE_PASS)) {
      return false;
    }
    if (!BKE_image_is_multilayer(image)) {
      return false;
    }
    ensure_multilayer_catalog(*image, tex->iuser);
    /* The declaration falls back to single-layer Color/Alpha when the layer
     * index doesn't resolve; match that fallback here too. */
    return BLI_findlink(&image->layers, tex->iuser.layer) != nullptr;
  }

  void handle_output_socket__multi_pass_image(const SocketInContext &socket)
  {
    /* The Bundle output collects every pass; a regular output is a single pass. */
    if (socket->type == SOCK_BUNDLE) {
      this->handle_output_socket__multi_pass_image_bundle(socket);
      return;
    }
    this->store_socket_value(socket, this->resolve_multi_pass_output(socket));
  }

  /**
   * Resolves one per-pass output socket of a multi-pass Image Texture node to a
   * single-pass copy's Color/Alpha output (design phase 16).
   */
  SocketValue resolve_multi_pass_output(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    const Image *image = reinterpret_cast<const Image *>(node->id);
    const NodeTexImage *tex = static_cast<const NodeTexImage *>(node->storage);

    /* Determine which pass this output socket represents. A socket named "Alpha"
     * is the synthetic alpha generated from the Combined pass unless the layer
     * has a real pass of that name. */
    const StringRefNull socket_name = socket->name;
    const ImageLayer *layer = static_cast<const ImageLayer *>(
        BLI_findlink(&image->layers, tex->iuser.layer));
    bool layer_has_socket_pass = false;
    if (layer != nullptr) {
      for (const ImagePass &pass : layer->passes) {
        if (socket_name == pass.name) {
          layer_has_socket_pass = true;
          break;
        }
      }
    }
    const bool is_synthetic_alpha = socket_name == "Alpha" && !layer_has_socket_pass;
    const std::string pass_name = is_synthetic_alpha ? std::string("Combined") :
                                                       std::string(socket_name);

    bNode &copied_node = this->get_or_create_single_pass_image_copy(node, pass_name);
    /* The synthetic alpha uses the copy's Alpha output; every real pass uses its
     * Color output, with the inliner's implicit conversion bridging the type. */
    const UString output_identifier = is_synthetic_alpha ? "Alpha"_ustr : "Color"_ustr;
    bNodeSocket *copied_output = bke::node_find_socket(copied_node, SOCK_OUT, output_identifier);
    BLI_assert(copied_output != nullptr);
    return {LinkedSocketValue{&copied_node, copied_output}};
  }

  /**
   * Builds a bundle of every pass of a multi-layer Image Texture node, one named
   * item per pass output, for use with a Separate Bundle node (design phase 15b).
   */
  void handle_output_socket__multi_pass_image_bundle(const SocketInContext &socket)
  {
    const NodeInContext node = socket.owner_node();
    auto bundle_value = std::make_shared<BundleSocketValue>();
    for (const bNodeSocket *output : node->output_sockets()) {
      if (output->type == SOCK_BUNDLE || !output->is_available()) {
        continue;
      }
      const SocketInContext pass_socket = {node.context, output};
      bundle_value->items.append({std::string(output->name),
                                  this->resolve_multi_pass_output(pass_socket),
                                  output->typeinfo});
    }
    this->store_socket_value(socket, {std::move(bundle_value)});
  }

  bNode &get_or_create_single_pass_image_copy(const NodeInContext &node,
                                              const std::string &pass_name)
  {
    Map<std::string, bNode *> &copies_by_pass = single_pass_image_copies_.lookup_or_add_default(
        node);
    if (bNode *const *existing = copies_by_pass.lookup_ptr(pass_name)) {
      return **existing;
    }

    Map<const bNodeSocket *, bNodeSocket *> socket_map;
    bNode &copied_node = *bke::node_copy_with_mapping(&dst_tree_,
                                                      *node.node,
                                                      this->node_copy_flag(),
                                                      std::nullopt,
                                                      this->get_next_node_identifier(),
                                                      socket_map);
    copied_node.parent = nullptr;

    /* Pin the copy to the single pass and force the standard Color/Alpha
     * declaration. */
    NodeTexImage &tex = *static_cast<NodeTexImage *>(copied_node.storage);
    tex.flag |= SHD_TEX_IMAGE_SINGLE_PASS;
    this->pin_image_user_to_pass(copied_node, pass_name);
    /* node_copy_with_mapping copies the source (multi-pass) sockets verbatim, so
     * rebuild the declaration now that the single-pass flag is set. */
    update_node_declaration_and_sockets(dst_tree_, copied_node);

    /* Wire the input sockets, mirroring handle_output_socket__eval_copy_node. */
    for (const bNodeSocket *src_input_socket : node->input_sockets()) {
      if (!src_input_socket->is_available()) {
        continue;
      }
      bNodeSocket *dst_input_socket = bke::node_find_socket(
          copied_node, SOCK_IN, UString(src_input_socket->identifier));
      if (dst_input_socket == nullptr) {
        continue;
      }
      const SocketInContext input_socket_ctx = {node.context, src_input_socket};
      const SocketValue &value = value_by_socket_.lookup(input_socket_ctx);
      this->set_input_socket_value(node, copied_node, *dst_input_socket, value);
    }

    copies_by_pass.add_new(pass_name, &copied_node);
    return copied_node;
  }

  /** Sets the copied node's #ImageUser to select the named pass of its layer.
   * The integer indices and the name fields are written together, so the iuser
   * never ends up with a name that doesn't match its index. */
  void pin_image_user_to_pass(bNode &copied_node, const std::string &pass_name)
  {
    const Image *image = reinterpret_cast<const Image *>(copied_node.id);
    NodeTexImage &tex = *static_cast<NodeTexImage *>(copied_node.storage);
    ImageUser &iuser = tex.iuser;

    const ImageLayer *layer = static_cast<const ImageLayer *>(
        BLI_findlink(&image->layers, iuser.layer));
    if (layer == nullptr) {
      return;
    }
    int pass_index = 0;
    for (const ImagePass &pass : layer->passes) {
      if (pass_name == pass.name) {
        iuser.pass = short(pass_index);
        STRNCPY(iuser.pass_name, pass_name.c_str());
        STRNCPY(iuser.layer_name, layer->name);
        return;
      }
      pass_index++;
    }
  }

  struct EnsureInputsResult {
    bool has_missing_inputs = false;
    bool all_inputs_primitive = false;
  };

  EnsureInputsResult ensure_node_inputs(const NodeInContext &node)
  {
    EnsureInputsResult result;
    result.has_missing_inputs = false;
    result.all_inputs_primitive = true;
    for (const bNodeSocket *input_socket : node->input_sockets()) {
      if (this->socket_is_ignored(*input_socket)) {
        continue;
      }
      const SocketInContext input_socket_ctx = {node.context, input_socket};
      const SocketValue *value = value_by_socket_.lookup_ptr(input_socket_ctx);
      if (!value) {
        this->schedule_socket(input_socket_ctx);
        result.has_missing_inputs = true;
        continue;
      }
      if (!value->to_primitive(*input_socket->typeinfo)) {
        result.all_inputs_primitive = false;
      }
    }
    return result;
  }

  void handle_output_socket__eval_multi_function(const NodeInContext &node)
  {
    NodeMultiFunctionBuilder builder{*node.node, node->owner_tree()};
    node->typeinfo->build_multi_function(builder);
    const mf::MultiFunction &fn = builder.function();
    mf::ContextBuilder context;
    ShaderNodesMultiFunctionUserData user_data;
    context.user_data(&user_data);
    IndexMask mask(1);
    mf::ParamsBuilder params{fn, &mask};

    /* Prepare inputs to the multi-function evaluation. */
    for (const bNodeSocket *input_socket : node->input_sockets()) {
      if (this->socket_is_ignored(*input_socket)) {
        continue;
      }
      const SocketInContext input_socket_ctx = {node.context, input_socket};
      const PrimitiveSocketValue value =
          *value_by_socket_.lookup(input_socket_ctx).to_primitive(*input_socket->typeinfo);
      params.add_readonly_single_input(
          GVArray::from_single(*input_socket->typeinfo->base_cpp_type, 1, value.buffer()));
    }

    /* Prepare output buffers. */
    Vector<void *> output_values;
    for (const bNodeSocket *output_socket : node->output_sockets()) {
      if (this->socket_is_ignored(*output_socket)) {
        continue;
      }
      void *value = scope_.allocate_owned(*output_socket->typeinfo->base_cpp_type);
      output_values.append(value);
      params.add_uninitialized_single_output(
          GMutableSpan(output_socket->typeinfo->base_cpp_type, value, 1));
    }

    fn.call(mask, params, context);

    /* Store constant-folded values for the output sockets. */
    int current_output_i = 0;
    for (const bNodeSocket *output_socket : node->output_sockets()) {
      if (this->socket_is_ignored(*output_socket)) {
        continue;
      }
      const void *value = output_values[current_output_i++];
      this->store_socket_value(
          {node.context, output_socket},
          {PrimitiveSocketValue::from_value({output_socket->typeinfo->base_cpp_type, value})});
    }

    for (const NodeWarning &warning : user_data.warnings) {
      this->report_error(node, warning.message, warning.type);
    }
  }

  bNode &handle_output_socket__eval_copy_node(const bNode &node,
                                              const ComputeContext *in_context,
                                              const ComputeContext *out_context)
  {
    Map<const bNodeSocket *, bNodeSocket *> socket_map;
    /* We generate our own identifier and name here to get unique values without having to scan all
     * already existing nodes. */
    const int identifier = this->get_next_node_identifier();
    const std::string unique_name = fmt::format("{}_{}", identifier, node.name);
    bNode &copied_node = *bke::node_copy_with_mapping(
        &dst_tree_,
        node,
        this->node_copy_flag(),
        unique_name.size() < sizeof(bNode::name) ? std::make_optional<StringRefNull>(unique_name) :
                                                   std::nullopt,
        identifier,
        socket_map);

    /* Clear the parent frame pointer, because it does not exist in the destination tree. */
    copied_node.parent = nullptr;

    /* Setup input sockets for the copied node. */
    for (const bNodeSocket *src_input_socket : node.input_sockets()) {
      if (this->socket_is_ignored(*src_input_socket)) {
        continue;
      }
      bNodeSocket &dst_input_socket = *socket_map.lookup(src_input_socket);
      const SocketInContext input_socket_ctx = {in_context, src_input_socket};
      const SocketValue &value = value_by_socket_.lookup(input_socket_ctx);
      this->set_input_socket_value({in_context, &node}, copied_node, dst_input_socket, value);
    }
    for (const bNodeSocket *src_output_socket : node.output_sockets()) {
      if (this->socket_is_ignored(*src_output_socket)) {
        continue;
      }
      bNodeSocket &dst_output_socket = *socket_map.lookup(src_output_socket);
      const SocketInContext output_socket_ctx = {out_context, src_output_socket};
      this->store_socket_value(output_socket_ctx,
                               {LinkedSocketValue{&copied_node, &dst_output_socket}});
    }
    return copied_node;
  }

  /** Converts the given socket value if necessary. */
  SocketValue handle_implicit_conversion(const SocketValue &src_value,
                                         const bke::bNodeSocketType &from_socket_type,
                                         const bke::bNodeSocketType &to_socket_type)
  {
    if (from_socket_type.type == to_socket_type.type) {
      return src_value;
    }
    if (std::get_if<LinkedSocketValue>(&src_value.value)) {
      return src_value;
    }
    if (std::get_if<DanglingValue>(&src_value.value)) {
      return src_value;
    }
    const std::optional<PrimitiveSocketValue> src_primitive_value = src_value.to_primitive(
        from_socket_type);
    if (src_primitive_value && to_socket_type.base_cpp_type) {
      if (data_type_conversions_.is_convertible(*from_socket_type.base_cpp_type,
                                                *to_socket_type.base_cpp_type))
      {
        const void *src_buffer = src_primitive_value->buffer();
        BUFFER_FOR_CPP_TYPE_VALUE(*to_socket_type.base_cpp_type, dst_buffer);
        data_type_conversions_.convert_to_uninitialized(*from_socket_type.base_cpp_type,
                                                        *to_socket_type.base_cpp_type,
                                                        src_buffer,
                                                        dst_buffer);
        return {
            PrimitiveSocketValue::from_value(GPointer{to_socket_type.base_cpp_type, dst_buffer})};
      }
    }
    if (src_primitive_value && to_socket_type.type == SOCK_SHADER) {
      const CPPType &from_cpp_type = *from_socket_type.base_cpp_type;
      const CPPType &to_cpp_type = CPPType::get<ColorGeometry4f>();
      if (from_cpp_type == to_cpp_type ||
          data_type_conversions_.is_convertible(from_cpp_type, to_cpp_type))
      {
        /* Insert a Color node when converting a primitive value to a shader. */
        bNode *color_node = this->add_node("ShaderNodeRGB"_ustr);
        const void *src_buffer = src_primitive_value->buffer();
        ColorGeometry4f color;
        data_type_conversions_.convert_to_uninitialized(
            from_cpp_type, to_cpp_type, src_buffer, &color);
        bNodeSocket *output_socket = static_cast<bNodeSocket *>(color_node->outputs.first);
        auto *socket_storage = static_cast<bNodeSocketValueRGBA *>(output_socket->default_value);
        copy_v3_v3(socket_storage->value, color);
        socket_storage->value[3] = 1.0f;
        return {LinkedSocketValue{color_node, output_socket}};
      }
    }

    return SocketValue{FallbackValue{}};
  }

  void set_input_socket_value(const NodeInContext &original_node,
                              bNode &dst_node,
                              bNodeSocket &dst_socket,
                              const SocketValue &value)
  {
    BLI_assert(dst_socket.is_input());
    if (dst_socket.flag & SOCK_HIDE_VALUE) {
      if (const auto *input_socket_value = std::get_if<InputSocketValue>(&value.value)) {
        if (input_socket_value->socket->flag & SOCK_HIDE_VALUE) {
          /* Don't add a value or link of the source and destination sockets don't have a value. */
          return;
        }
      }
    }
    if (const std::optional<PrimitiveSocketValue> primitive_value = value.to_primitive(
            *dst_socket.typeinfo))
    {
      if (dst_socket.flag & SOCK_HIDE_VALUE) {
        /* Can't store the primitive value directly on the socket. So create a new input node and
         * link it instead. */
        const NodeAndSocket node_and_socket = this->primitive_value_to_output_socket(
            *primitive_value);
        if (dst_tree_.typeinfo->validate_link(node_and_socket.socket->typeinfo->type,
                                              dst_socket.typeinfo->type))
        {
          bke::node_add_link(
              dst_tree_, *node_and_socket.node, *node_and_socket.socket, dst_node, dst_socket);
        }
      }
      else {
        this->set_primitive_value_on_socket(dst_socket, *primitive_value);
      }
      return;
    }
    if (!params_.allow_preserving_repeat_zones) {
      const bool is_iterations_input = dst_node.inputs.first == &dst_socket &&
                                       dst_node.is_type("GeometryNodeRepeatInput"_ustr);
      if (is_iterations_input) {
        this->add_dynamic_repeat_zone_iterations_error(original_node);
        this->set_primitive_value_on_socket(dst_socket, PrimitiveSocketValue{0});
        return;
      }
    }
    if (std::get_if<InputSocketValue>(&value.value)) {
      /* Cases were the input has a primitive value are handled above. */
      return;
    }
    if (std::get_if<FallbackValue>(&value.value)) {
      /* Cases were the input has a primitive fallback value are handled above. */
      return;
    }
    if (std::get_if<DanglingValue>(&value.value)) {
      /* Input sockets should never have a dangling value, because they are replaced by the socket
       * value in #handle_input_socket. */
      BLI_assert_unreachable();
      return;
    }
    if (std::get_if<BundleSocketValuePtr>(&value.value)) {
      /* This type can't be assigned to a socket. The bundle has to be separated first. */
      BLI_assert_unreachable();
      return;
    }
    if (std::get_if<ClosureZoneValue>(&value.value)) {
      /* This type can't be assigned to a socket. One has to evaluate a closure. */
      BLI_assert_unreachable();
      return;
    }
    if (std::get_if<MultiInputValue>(&value.value)) {
      /* This type can't be assigned to a socket. */
      BLI_assert_unreachable();
      return;
    }
    if (const auto *src_socket_value = std::get_if<LinkedSocketValue>(&value.value)) {
      if (dst_tree_.typeinfo->validate_link(src_socket_value->socket->typeinfo->type,
                                            dst_socket.typeinfo->type))
      {
        bke::node_add_link(
            dst_tree_, *src_socket_value->node, *src_socket_value->socket, dst_node, dst_socket);
      }
      return;
    }
    BLI_assert_unreachable();
  }

  NodeAndSocket primitive_value_to_output_socket(const PrimitiveSocketValue &value)
  {
    if (const float *value_float = std::get_if<float>(&value.value)) {
      bNode *node = this->add_node("ShaderNodeValue"_ustr);
      bNodeSocket *socket = static_cast<bNodeSocket *>(node->outputs.first);
      socket->default_value_typed<bNodeSocketValueFloat>()->value = *value_float;
      return {node, socket};
    }
    if (const int *value_int = std::get_if<int>(&value.value)) {
      bNode *node = this->add_node("FunctionNodeInputInt"_ustr);
      auto &storage = *static_cast<NodeInputInt *>(node->storage);
      storage.integer = *value_int;
      return {node, static_cast<bNodeSocket *>(node->outputs.first)};
    }
    if (const bool *value_bool = std::get_if<bool>(&value.value)) {
      bNode *node = this->add_node("FunctionNodeInputBool"_ustr);
      auto &storage = *static_cast<NodeInputBool *>(node->storage);
      storage.boolean = int(*value_bool);
      return {node, static_cast<bNodeSocket *>(node->outputs.first)};
    }
    if (const float3 *value_float3 = std::get_if<float3>(&value.value)) {
      bNode *node = this->add_node("FunctionNodeInputVector"_ustr);
      auto &storage = *static_cast<NodeInputVector *>(node->storage);
      copy_v3_v3(storage.vector, *value_float3);
      storage.dimensions = 3;
      return {node, static_cast<bNodeSocket *>(node->outputs.first)};
    }
    if (const ColorGeometry4f *value_color = std::get_if<ColorGeometry4f>(&value.value)) {
      bNode *node = this->add_node("ShaderNodeRGB"_ustr);
      bNodeSocket *output_socket = static_cast<bNodeSocket *>(node->outputs.first);
      auto *socket_storage = static_cast<bNodeSocketValueRGBA *>(output_socket->default_value);
      copy_v3_v3(socket_storage->value, *value_color);
      socket_storage->value[3] = 1.0f;
      return {node, output_socket};
    }
    BLI_assert_unreachable();
    return {};
  }

  bNode *add_node(const UString idname)
  {
    return bke::node_add_node(nullptr, dst_tree_, idname, this->get_next_node_identifier());
  }

  int get_next_node_identifier()
  {
    return ++dst_node_counter_;
  }

  void set_primitive_value_on_socket(bNodeSocket &socket, const PrimitiveSocketValue &value)
  {
    switch (socket.type) {
      case SOCK_FLOAT: {
        socket.default_value_typed<bNodeSocketValueFloat>()->value = std::get<float>(value.value);
        break;
      }
      case SOCK_INT: {
        socket.default_value_typed<bNodeSocketValueInt>()->value = std::get<int>(value.value);
        break;
      }
      case SOCK_BOOLEAN: {
        socket.default_value_typed<bNodeSocketValueBoolean>()->value = std::get<bool>(value.value);
        break;
      }
      case SOCK_VECTOR: {
        copy_v3_v3(socket.default_value_typed<bNodeSocketValueVector>()->value,
                   std::get<float3>(value.value));
        break;
      }
      case SOCK_RGBA: {
        copy_v4_v4(socket.default_value_typed<bNodeSocketValueRGBA>()->value,
                   std::get<ColorGeometry4f>(value.value));
        break;
      }
      case SOCK_MENU: {
        socket.default_value_typed<bNodeSocketValueMenu>()->value =
            std::get<nodes::MenuValue>(value.value).value;
        break;
      }
      case SOCK_STRING: {
        STRNCPY_UTF8(socket.default_value_typed<bNodeSocketValueString>()->value,
                     std::get<std::string>(value.value).c_str());
        break;
      }
      default: {
        BLI_assert_unreachable();
        break;
      }
    }
  }

  void restore_zones_in_output_tree()
  {
    for (const PreservedZone &copied_zone : copied_zone_by_zone_output_node_.values()) {
      if (!copied_zone.input_node || !copied_zone.output_node) {
        continue;
      }
      const bke::bNodeZoneType *zone_type = bke::zone_type_by_node_type(
          copied_zone.input_node->type_legacy);
      if (!zone_type) {
        continue;
      }
      int &output_id = zone_type->get_corresponding_output_id(*copied_zone.input_node);
      output_id = copied_zone.output_node->identifier;
    }
  }

  void position_nodes_in_output_tree()
  {
    bNodeTree &tree = dst_tree_;
    tree.ensure_topology_cache();

    Map<int, int> num_by_depth;
    Map<bNode *, int> depth_by_node;

    /* Simple algorithm that does a very rough layout of the generated tree. This does not produce
     * great results generally, but is usually good enough when debugging smaller node trees. */
    for (bNode *node : tree.toposort_right_to_left()) {
      int depth = 0;
      for (bNodeSocket *socket : node->output_sockets()) {
        for (bNodeSocket *target : socket->directly_linked_sockets()) {
          depth = std::max(depth, depth_by_node.lookup(&target->owner_node()) + 1);
        }
      }
      depth_by_node.add_new(node, depth);
      const int index_at_depth = num_by_depth.lookup_or_add(depth, 0)++;
      node->location[0] = 200 - depth * 200;
      node->location[1] = -index_at_depth * 300;
    }
  }

  /**
   * Utility to that copies the value of the origin socket to the current socket. If the origin
   * value does not exist yet, the origin socket is only scheduled.
   */
  void forward_value_or_schedule(const SocketInContext &socket, const SocketInContext &origin)
  {
    if (const SocketValue *value = value_by_socket_.lookup_ptr(origin)) {
      if (socket->type == origin->type) {
        this->store_socket_value(socket, *value);
        return;
      }
      this->store_socket_value(
          socket, this->handle_implicit_conversion(*value, *origin->typeinfo, *socket->typeinfo));
      return;
    }
    this->schedule_socket(origin);
  }

  void store_socket_value(const SocketInContext &socket, SocketValue value)
  {
    value_by_socket_.add_new(socket, std::move(value));
  }

  void store_socket_value_fallback(const SocketInContext &socket)
  {
    value_by_socket_.add_new(socket, {FallbackValue{}});
  }

  void store_socket_value_dangling(const SocketInContext &socket)
  {
    value_by_socket_.add_new(socket, {DanglingValue{}});
  }

  void schedule_socket(const SocketInContext &socket)
  {
    scheduled_sockets_stack_.push(socket);
  }

  int node_copy_flag() const
  {
    const bool use_refcounting = !(dst_tree_.id.tag & ID_TAG_NO_MAIN);
    return use_refcounting ? 0 : LIB_ID_CREATE_NO_USER_REFCOUNT;
  }

  bool socket_is_ignored(const bNodeSocket &socket) const
  {
    return !socket.is_available() || socket.typeinfo->idname == "NodeSocketVirtual"_ustr;
  }

  void report_error(const NodeInContext &node, const StringRef message, const NodeWarningType type)
  {
    Vector<NodeInContext> nodes;
    nodes.append(node);
    for (const ComputeContext *context = node.context; context; context = context->parent()) {
      if (const auto *group_context = dynamic_cast<const bke::GroupNodeComputeContext *>(context))
      {
        nodes.append({context->parent(), group_context->node()});
      }
    }
    for (const NodeInContext &node : nodes) {
      params_.r_error_messages.append({&*node, message, type});
    }
  }
};

}  // namespace

bool inline_shader_node_tree(const bNodeTree &src_tree,
                             bNodeTree &dst_tree,
                             InlineShaderNodeTreeParams &params)
{
  ShaderNodesInliner inliner(src_tree, dst_tree, params);

  if (inliner.do_inline()) {
    /* Update deprecated bNodeSocket.link pointers because some code still depends on it. */
    for (bNode &node : dst_tree.nodes) {
      for (bNodeSocket &sock : node.inputs) {
        sock.link = nullptr;
      }
      for (bNodeSocket &sock : node.outputs) {
        sock.link = nullptr;
      }
    }
    for (bNodeLink &link : dst_tree.links) {
      link.tosock->link = &link;
      BLI_assert(dst_tree.typeinfo->validate_link(link.fromsock->typeinfo->type,
                                                  link.tosock->typeinfo->type));
      link.flag |= NODE_LINK_VALID;
    }
    return true;
  }

  return false;
}

}  // namespace blender::nodes
