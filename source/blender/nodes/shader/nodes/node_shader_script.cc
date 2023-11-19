/* SPDX-FileCopyrightText: 2005 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup shdnodes
 */

#include "node_shader_util.hh"

#include "RNA_access.hh"

#include "UI_interface.hh"
#include "UI_resources.hh"

namespace blender::nodes::node_shader_script_cc {

static void node_shader_buts_script(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiLayout *row;

  row = uiLayoutRow(layout, false);
  uiItemR(row, ptr, "mode", UI_ITEM_R_SPLIT_EMPTY_NAME | UI_ITEM_R_EXPAND, nullptr, ICON_NONE);

  row = uiLayoutRow(layout, true);

  if (RNA_enum_get(ptr, "mode") == NODE_SCRIPT_INTERNAL) {
    uiItemR(row, ptr, "script", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
  }
  else {
    uiItemR(row, ptr, "filepath", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
  }

  uiItemO(row, "", ICON_FILE_REFRESH, "node.shader_script_update");

  row = uiLayoutRow(layout, true);
  uiItemR(row, ptr, "label", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
}

static void node_shader_buts_script_ex(uiLayout *layout, bContext *C, PointerRNA *ptr)
{
  uiItemS(layout);

  node_shader_buts_script(layout, C, ptr);

#if 0 /* not implemented yet */
  if (RNA_enum_get(ptr, "mode") == NODE_SCRIPT_EXTERNAL) {
    uiItemR(layout, ptr, "use_auto_update", UI_ITEM_R_SPLIT_EMPTY_NAME, nullptr, ICON_NONE);
  }
#endif
}

static void init(bNodeTree * /*ntree*/, bNode *node)
{
  NodeShaderScript *nss = MEM_cnew<NodeShaderScript>("shader script node");
  node->storage = nss;
}

static void node_free_script(bNode *node)
{
  NodeShaderScript *nss = static_cast<NodeShaderScript *>(node->storage);

  if (nss) {
    if (nss->bytecode) {
      MEM_freeN(nss->bytecode);
    }

    MEM_freeN(nss);
  }
}

static void node_copy_script(bNodeTree * /*dst_ntree*/, bNode *dest_node, const bNode *src_node)
{
  NodeShaderScript *src_nss = static_cast<NodeShaderScript *>(src_node->storage);
  NodeShaderScript *dest_nss = static_cast<NodeShaderScript *>(MEM_dupallocN(src_nss));

  if (src_nss->bytecode) {
    dest_nss->bytecode = static_cast<char *>(MEM_dupallocN(src_nss->bytecode));
  }

  dest_node->storage = dest_nss;
}

static int node_shader_gpu_script(GPUMaterial *mat,
                                  bNode *node,
                                  bNodeExecData * /*execdata*/,
                                  GPUNodeStack *in,
                                  GPUNodeStack *out)
{
  // TODO External script loading
  Text* text = reinterpret_cast<Text *>(node->id);
  if (text == nullptr) {
    return 0;
  }

  // Register text as runtime library, but not if it's already registered
  GPU_material_library_text_add(text, false);

  return GPU_stack_link(mat, node, node->label, in, out);
}

void node_declare_dynamic(const bNodeTree &node_tree,
                          const bNode &node,
                          NodeDeclaration &r_declaration)
{
  const NodeShaderScript *nss = static_cast<NodeShaderScript *>(node.storage);
  Text* text = reinterpret_cast<Text *>(node.id);
  if (text == nullptr) {
    return;
  }

  // Register text as runtime library, but only if it has not already been registered.
  // If the text has changed since, then this will *not* reflect the updated version
  GPU_material_library_text_add(text, false);
  GPU_material_node_decl_from_function(&r_declaration, node.label);
}

}  // namespace blender::nodes::node_shader_script_cc

void register_node_type_sh_script()
{
  namespace file_ns = blender::nodes::node_shader_script_cc;

  static bNodeType ntype;

  sh_node_type_base(&ntype, SH_NODE_SCRIPT, "Script", NODE_CLASS_SCRIPT);
  ntype.draw_buttons = file_ns::node_shader_buts_script;
  ntype.draw_buttons_ex = file_ns::node_shader_buts_script_ex;
  ntype.initfunc = file_ns::init;
  ntype.declare_dynamic = file_ns::node_declare_dynamic;
  ntype.gpu_fn = file_ns::node_shader_gpu_script;

  node_type_storage(
      &ntype, "NodeShaderScript", file_ns::node_free_script, file_ns::node_copy_script);

  nodeRegisterType(&ntype);
}
