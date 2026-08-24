/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "node_geometry_util.hh"

#include <cstring>

#include "BLI_string.hh"
#include "BLI_string_utf8.hh"

#include "BLT_translation.hh"

#include "BLO_read_write.hh"

#include "NOD_geometry_nodes_lazy_function.hh"
#include "NOD_rna_define.hh"

#include "RNA_enum_types.hh"

#include "UI_interface_layout.hh"
#include "UI_resources.hh"

namespace blender {

namespace nodes::node_geo_modal_event_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Bool>("Enable"_ustr)
      .default_value(true)
      .structure_type(StructureType::Single)
      .description(
          "Keep the tool running so that it can react to this event instead of ending execution");
  b.add_output<decl::Bool>("Is Active"_ustr)
      .structure_type(StructureType::Single)
      .description("True when this execution was started by this event");
}

static void node_layout(ui::Layout &layout, bContext * /*C*/, PointerRNA *ptr)
{
  layout.prop(ptr, "event_name", UI_ITEM_NONE, "", ICON_NONE);
  layout.prop(ptr, "description", UI_ITEM_NONE, "", ICON_NONE);
}

static void node_label(const bNodeTree * /*tree*/,
                       const bNode *node,
                       char *label,
                       const int label_maxncpy)
{
  const auto &storage = *static_cast<const GeometryNodeModalEvent *>(node->storage);
  if (storage.name && storage.name[0] != '\0') {
    BLI_strncpy_utf8(label, storage.name, label_maxncpy);
  }
  else {
    BLI_strncpy_utf8(label, CTX_N_(BLT_I18NCONTEXT_ID_NODETREE, "Modal Event"), label_maxncpy);
  }
}

class LazyFunctionForModalEventNode : public LazyFunction {
  const bNode &node_;

 public:
  LazyFunctionForModalEventNode(const bNode &node, MutableSpan<int> r_lf_index_by_bsocket)
      : node_(node)
  {
    debug_name_ = node.name;
    r_lf_index_by_bsocket[node.input_socket(0).index_in_tree()] = inputs_.append_and_get_index_as(
        "Enable", CPPType::get<SocketValueVariant>());
    r_lf_index_by_bsocket[node.output_socket(0).index_in_tree()] =
        outputs_.append_and_get_index_as("Is Active", CPPType::get<SocketValueVariant>());
  }

  void execute_impl(lf::Params &params, const lf::Context &context) const override
  {
    const GeoNodesUserData &user_data = *static_cast<GeoNodesUserData *>(context.user_data);
    const GeoNodesOperatorData *operator_data = user_data.call_data->operator_data;
    if (!operator_data || !operator_data->modal_requested) {
      GeoNodesLocalUserData &local_user_data = *static_cast<GeoNodesLocalUserData *>(
          context.local_user_data);
      if (eval_log::NodeTreeLogger *tree_logger = local_user_data.try_get_tree_logger(user_data)) {
        tree_logger->node_warnings.append(
            *tree_logger->allocator,
            {node_.identifier, {NodeWarningType::Error, TIP_("Node must be run as tool")}});
      }
      set_default_remaining_node_outputs(params, node_);
      return;
    }

    /* The node is evaluated as a side-effect node, so the "Enable" input is always available. */
    const SocketValueVariant &enable_variant = params.get_input<SocketValueVariant>(0);
    const bool enabled = enable_variant.is_single() && enable_variant.get<bool>();
    if (enabled) {
      /* Other modal nodes may set this at the same time, so it must never be reset here. */
      operator_data->modal_requested->store(true, std::memory_order_relaxed);
    }

    const auto &storage = *static_cast<const GeometryNodeModalEvent *>(node_.storage);
    const StringRef name = storage.name ? StringRef(storage.name) : StringRef();
    const bool is_active = enabled && !name.is_empty() && name == operator_data->active_event;
    params.set_output(0, SocketValueVariant(is_active));
  }
};

static void node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->storage = MEM_new<GeometryNodeModalEvent>(__func__);
}

static void node_storage_free(bNode *node)
{
  auto *storage = static_cast<GeometryNodeModalEvent *>(node->storage);
  if (storage == nullptr) {
    return;
  }
  MEM_delete(storage->name);
  MEM_delete(storage->description);
  MEM_delete(storage);
}

static void node_storage_copy(bNodeTree * /*dst_tree*/, bNode *dst_node, const bNode *src_node)
{
  const auto *src_storage = static_cast<const GeometryNodeModalEvent *>(src_node->storage);
  auto *dst_storage = MEM_dupalloc(src_storage);
  dst_storage->name = BLI_strdup_null(src_storage->name);
  dst_storage->description = BLI_strdup_null(src_storage->description);
  dst_node->storage = dst_storage;
}

static void node_blend_write(const bNodeTree & /*tree*/, const bNode &node, BlendWriter &writer)
{
  const auto &storage = *static_cast<const GeometryNodeModalEvent *>(node.storage);
  writer.write_string(storage.name);
  writer.write_string(storage.description);
}

static void node_blend_read(bNodeTree & /*tree*/, bNode &node, BlendDataReader &reader)
{
  auto &storage = *static_cast<GeometryNodeModalEvent *>(node.storage);
  BLO_read_string(&reader, &storage.name);
  BLO_read_string(&reader, &storage.description);
}

/**
 * The strings are stored as pointers in the node storage, so they need runtime accessors instead
 * of the DNA based ones.
 */
template<char *GeometryNodeModalEvent::*member>
static std::string storage_string_get(PointerRNA *ptr, PropertyRNA * /*prop*/)
{
  const bNode &node = *static_cast<const bNode *>(ptr->data);
  const auto &storage = *static_cast<const GeometryNodeModalEvent *>(node.storage);
  const char *value = storage.*member;
  return value ? value : "";
}

template<char *GeometryNodeModalEvent::*member>
static int storage_string_length(PointerRNA *ptr, PropertyRNA * /*prop*/)
{
  const bNode &node = *static_cast<const bNode *>(ptr->data);
  const auto &storage = *static_cast<const GeometryNodeModalEvent *>(node.storage);
  const char *value = storage.*member;
  return value ? int(std::strlen(value)) : 0;
}

template<char *GeometryNodeModalEvent::*member>
static void storage_string_set(PointerRNA *ptr, PropertyRNA * /*prop*/, const std::string &value)
{
  bNode &node = *static_cast<bNode *>(ptr->data);
  auto &storage = *static_cast<GeometryNodeModalEvent *>(node.storage);
  MEM_SAFE_DELETE(storage.*member);
  storage.*member = BLI_strdupn(value.c_str(), value.size());
}

// TODO: REPLACE WITH OPERATOR MODAL KEYMAP
static void event_name_search(const bContext * /*C*/,
                              PointerRNA * /*ptr*/,
                              PropertyRNA * /*prop*/,
                              const char * /*edit_text*/,
                              FunctionRef<void(StringPropertySearchVisitParams)> visit_fn)
{
  for (const EnumPropertyItem *item = rna_enum_event_type_items; item->identifier; item++) {
    if (item->identifier[0] == '\0') {
      /* Skip the separators between the event type categories. */
      continue;
    }
    StringPropertySearchVisitParams params{};
    params.text = item->identifier;
    params.info = item->name;
    visit_fn(params);
  }
}

static void node_rna(StructRNA *srna)
{
  PropertyRNA *prop;

  prop = RNA_def_property(srna, "event_name", PROP_STRING, PROP_NONE);
  RNA_def_property_string_funcs_runtime(prop,
                                        storage_string_get<&GeometryNodeModalEvent::name>,
                                        storage_string_length<&GeometryNodeModalEvent::name>,
                                        storage_string_set<&GeometryNodeModalEvent::name>,
                                        nullptr,
                                        nullptr);
  RNA_def_property_string_search_func_runtime(
      prop, event_name_search, PROP_STRING_SEARCH_SORT | PROP_STRING_SEARCH_SUGGESTION);
  RNA_def_property_ui_text(
      prop, "Name", "Name of the event that the node tool reacts to while it keeps running");
  RNA_def_property_update_runtime(prop, rna_Node_update);
  RNA_def_property_update_notifier(prop, NC_NODE | NA_EDITED);

  prop = RNA_def_property(srna, "description", PROP_STRING, PROP_NONE);
  RNA_def_property_string_funcs_runtime(
      prop,
      storage_string_get<&GeometryNodeModalEvent::description>,
      storage_string_length<&GeometryNodeModalEvent::description>,
      storage_string_set<&GeometryNodeModalEvent::description>,
      nullptr,
      nullptr);
  RNA_def_property_ui_text(
      prop, "Description", "Explanation of what the node tool does when the event is processed");
  RNA_def_property_update_runtime(prop, rna_Node_update);
  RNA_def_property_update_notifier(prop, NC_NODE | NA_EDITED);
}

static void node_register()
{
  static bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeModalEvent"_ustr);
  ntype.ui_name = "Modal Event";
  ntype.ui_description =
      "Define an event that the tool reacts to while it keeps running, and detect when that event "
      "is being processed";
  ntype.nclass = NODE_CLASS_INPUT;
  ntype.declare = node_declare;
  ntype.initfunc = node_init;
  ntype.labelfunc = node_label;
  ntype.draw_buttons = node_layout;
  bke::node_type_storage(ntype, "GeometryNodeModalEvent", node_storage_free, node_storage_copy);
  ntype.blend_write_storage_content = node_blend_write;
  ntype.blend_data_read_storage_content = node_blend_read;
  ntype.gather_link_search_ops = search_link_ops_for_tool_node;
  bke::node_register_type(ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace nodes::node_geo_modal_event_cc

namespace nodes {

std::unique_ptr<LazyFunction> get_modal_event_node_lazy_function(
    const bNode &node, GeometryNodesLazyFunctionGraphInfo &own_lf_graph_info)
{
  using namespace node_geo_modal_event_cc;
  return std::make_unique<LazyFunctionForModalEventNode>(
      node, own_lf_graph_info.mapping.lf_index_by_bsocket);
}

}  // namespace nodes
}  // namespace blender
