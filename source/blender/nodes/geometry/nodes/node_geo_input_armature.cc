/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_listbase.h"

#include "BKE_pointcloud.hh"

#include "UI_interface.hh"
#include "UI_resources.hh"

#include "node_geometry_util.hh"

namespace blender::nodes::node_geo_input_armature_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Object>("Armature").custom_draw([](CustomSocketDrawParams &params) {
    params.layout.alignment_set(blender::ui::LayoutAlign::Expand);
    params.layout.prop(&params.node_ptr, "armature", UI_ITEM_NONE, "", ICON_NONE);
  });
  ;
  b.add_output<decl::Geometry>("Armature", "Geometry");
}

static void node_geo_exec(GeoNodeExecParams params)
{
  const Object *armature_obj = reinterpret_cast<Object *>(params.node().id);
  if (!armature_obj || armature_obj->type != OB_ARMATURE) {
    params.error_message_add(NodeWarningType::Error, TIP_("Object is not an Armature"));
    params.set_default_remaining_outputs();
    return;
  }

  const bPose *pose = armature_obj->pose;
  if (!pose) {
    params.set_default_remaining_outputs();
    return;
  }
  const int count = BLI_listbase_count(&pose->chanbase);
  const Span<bPoseChannel *> pose_channels{pose->chan_array, count};

  PointCloud *armature_points = BKE_pointcloud_new_nomain(count);
  MutableAttributeAccessor attributes = armature_points->attributes_for_write();
  MutableSpan<float3> positions = armature_points->positions_for_write();
  SpanAttributeWriter<float4x4> rest_transforms =
      attributes.lookup_or_add_for_write_only_span<float4x4>("rest_transform",
                                                             bke::AttrDomain::Point);
  SpanAttributeWriter<float4x4> pose_transforms =
      attributes.lookup_or_add_for_write_only_span<float4x4>("pose_transform",
                                                             bke::AttrDomain::Point);
  SpanAttributeWriter<int> parents = attributes.lookup_or_add_for_write_only_span<int>(
      "parent", bke::AttrDomain::Point);
  SpanAttributeWriter<bool> deform = attributes.lookup_or_add_for_write_only_span<bool>(
      "deform", bke::AttrDomain::Point);

  for (const int bone_i : pose_channels.index_range()) {
    bPoseChannel *pchan = pose_channels[bone_i];
    // const int parent_bone = pose_channels.first_index_try(pchan->parent);
    positions[bone_i] = float3(pchan->pose_head);
    rest_transforms.span[bone_i] = float4x4(pchan->bone->arm_mat);
    pose_transforms.span[bone_i] = float4x4(pchan->pose_mat);
    parents.span[bone_i] = pose_channels.first_index_try(pchan->parent);
    deform.span[bone_i] = (pchan->bone->flag & BONE_NO_DEFORM) == 0;
  }
  rest_transforms.finish();
  pose_transforms.finish();
  parents.finish();
  deform.finish();

  params.set_output("Geometry", GeometrySet::from_pointcloud(armature_points));
}

static void node_register()
{
  static blender::bke::bNodeType ntype;

  geo_node_type_base(&ntype, "GeometryNodeInputArmature");
  ntype.ui_name = "Armature";
  ntype.ui_description = "Get armature as point cloud geometry";
  ntype.nclass = NODE_CLASS_INPUT;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.declare = node_declare;
  blender::bke::node_register_type(ntype);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_input_armature_cc
