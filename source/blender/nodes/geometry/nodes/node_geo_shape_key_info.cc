/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_math_matrix.hh"
#include "BLI_listbase_iterator.hh"

#include "BKE_geometry_set.hh"
#include "BKE_mesh.hh"
#include "BKE_curve.hh"
#include "BKE_curves.hh"
#include "BKE_curve_legacy_convert.hh"
#include "BKE_mesh_runtime.hh"
#include "BKE_geometry_fields.hh"
#include "BKE_key.hh"
#include "BKE_object.hh"
#include "BKE_idtype.hh"
#include "BKE_deform.hh"

#include "GEO_transform.hh"

#include "NOD_geometry_nodes_list.hh"
#include "DEG_depsgraph_query.hh"

#include "NOD_rna_define.hh"
#include "NOD_socket_search_link.hh"

#include "UI_interface_layout.hh"
#include "UI_resources.hh"

#include "node_geometry_util.hh"

#include "DNA_key_types.h"
#include "DNA_mesh_types.h"
#include "DNA_curve_types.h"
#include "DNA_object_types.h"

#include "FN_field.hh"
#include "FN_multi_function.hh"

enum GEO_SHAPE_KEY_INFO_MODE : int {
  GEO_NODE_SHAPE_KEY_INFO_MODE_MESH = 0,
  GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE = 1,
};

namespace blender::nodes::node_geo_shape_key_info_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Vector>("Position"_ustr)
    .structure_type(StructureType::Field)
    .description("Position of each element in the geometry for the selected Shape Key.");

  const bNode * node = b.node_or_null();
  if (node) {
    if (node->custom1 == GEO_NODE_SHAPE_KEY_INFO_MODE_MESH) {
      b.add_output<decl::Float>("Vertex Group"_ustr)
        .structure_type(StructureType::Field)
        .description("Vertex group weights associated with the selected Shape Key.");
    } 
    else {
      b.add_output<decl::Vector>("Handle Left"_ustr)
          .structure_type(StructureType::Field)
          .description("Position of the left handle of each control point for the selected Shape Key.");
      b.add_output<decl::Vector>("Handle Right"_ustr)
          .structure_type(StructureType::Field)
          .description("Position of the right handle of each control point for the selected Shape Key.");
      b.add_output<decl::Float>("Radius"_ustr)
          .structure_type(StructureType::Field)
          .description("Radius of each control point for the selected Shape Key.");
      b.add_output<decl::Float>("Tilt"_ustr)
          .structure_type(StructureType::Field)
          .description("Tilt angle of each control point for the selected Shape Key.");
    }
  }

  b.add_output<decl::Float>("Value"_ustr)
    .structure_type(StructureType::Single)
    .description("Current value of the selected Shape Key.");
  
  b.add_output<decl::Float>("Range Min"_ustr)
    .structure_type(StructureType::Single)
    .description("Current range min of the selected Shape Key.");

  b.add_output<decl::Float>("Range Max"_ustr)
    .structure_type(StructureType::Single)
    .description("Current range max of the selected Shape Key.");

  b.add_output<decl::Bool>("Exists"_ustr)
    .description("True if the specified Shape Key exists in the object.");

  b.add_input<decl::Object>("Object"_ustr)
    .optional_label()
    .description("Object containing the Shape Keys. Uses the current object if not specified.");

  b.add_input<decl::String>("Shape Key Name"_ustr)
    .optional_label()
    .description("Name of the Shape Key to retrieve.");

}

static void node_layout(ui::Layout &layout, bContext * /*C*/, PointerRNA *ptr)
{
  layout.prop(ptr, "mode", UI_ITEM_NONE, "", ICON_NONE);
}

static void node_gather_link_search_ops(GatherLinkSearchOpParams &params)
{
  const eNodeSocketDatatype other_type = params.other_socket().type;

  if (params.in_out() == SOCK_OUT) {
    if (other_type == SOCK_VECTOR) {
      params.add_item(IFACE_("Position"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Position"_ustr);
      });
      params.add_item(IFACE_("Handle Left"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        node.custom1 = GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE;
        params.update_and_connect_available_socket(node, "Handle Left"_ustr);
      });
      params.add_item(IFACE_("Handle Right"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        node.custom1 = GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE;
        params.update_and_connect_available_socket(node, "Handle Right"_ustr);
      });
    }
    if (other_type == SOCK_FLOAT) {
      params.add_item(IFACE_("Vertex Group"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Vertex Group"_ustr);
      });
      params.add_item(IFACE_("Radius"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        node.custom1 = GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE;
        params.update_and_connect_available_socket(node, "Radius"_ustr);
      });
      params.add_item(IFACE_("Tilt"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        node.custom1 = GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE;
        params.update_and_connect_available_socket(node, "Tilt"_ustr);
      });
      params.add_item(IFACE_("Value"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Value"_ustr);
      });
    }
    if (other_type == SOCK_BOOLEAN) {
      params.add_item(IFACE_("Exists"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Exists"_ustr);
      }); 
    }
  } 
  else {
    if (other_type == SOCK_OBJECT) {
      params.add_item(IFACE_("Object"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Object"_ustr);
      }); 
    }
    if (other_type == SOCK_STRING) {
      params.add_item(IFACE_("Shape Key Name"), [](LinkSearchOpParams &params) {
        bNode &node = params.add_node("GeometryNodeShapeKeyInfo"_ustr);
        params.update_and_connect_available_socket(node, "Shape Key Name"_ustr);
      }); 
    }
  }

}

static void node_node_init(bNodeTree * /*tree*/, bNode *node)
{
  node->custom1 = GEO_NODE_SHAPE_KEY_INFO_MODE_MESH;
}

struct ShapeKeyInfo {
  KeyBlock *key_block;
  int point_count;
  float value;
  float range_min;
  float range_max;
  bool exists;
};

struct ShapeKeyMeshInfo{
  bool is_match;
  Field<float3> position;
  Field<float> vertex_group;
};

struct ShapeKeyCurveInfo {
  bool is_match = false;
  Field<float3> position;
  Field<float3> handle_left;
  Field<float3> handle_right;
  Field<float> radius;
  Field<float> tilt;
};

class ShapeKeyVectorsFunction : public mf::MultiFunction {
  private:
    std::shared_ptr<std::vector<float3>> data_;
    mf::Signature signature_;

  public:
    ShapeKeyVectorsFunction(std::shared_ptr<std::vector<float3>> data, const char *function_name)
        : data_(std::move(data))
    {
      mf::SignatureBuilder builder{function_name, signature_};
      
      builder.single_input<int>("Index");
      builder.single_output<float3>("Position");
      this->set_signature(&signature_);
    }

    void call(const IndexMask &mask, mf::Params params, mf::Context /*context*/) const override
    {
      const VArray<int> &indices = params.readonly_single_input<int>(0);
      MutableSpan<float3> vetors = params.uninitialized_single_output<float3>(1);

      mask.foreach_index([&](const int64_t i) {
        const int index = indices[i];
        if (index < 0 || index >= int(data_->size())) {
          new (&vetors[i]) float3(0.0f);
          return;
        }
        new (&vetors[i]) float3((*data_)[index]);
      });
    }
};

class ShapeKeyValuesFunctions : public mf::MultiFunction {
 private:
  std::shared_ptr<std::vector<float>> data_;
  mf::Signature signature_;

 public:
  ShapeKeyValuesFunctions(std::shared_ptr<std::vector<float>> data,
                          const char *function_name)
      : data_(std::move(data))
  {
    mf::SignatureBuilder builder{function_name, signature_};

    builder.single_input<int>("Index");
    builder.single_output<float>("Value");

    this->set_signature(&signature_);
  }
  void call(const IndexMask &mask,
            mf::Params params,
            mf::Context /*context*/) const override
  {
    const VArray<int> &indices = params.readonly_single_input<int>(0);
    MutableSpan<float> values = params.uninitialized_single_output<float>(1);

    mask.foreach_index([&](const int64_t i) {
      const int index = indices[i];

      if (index < 0 || index >= int(data_->size())) {
        new (&values[i]) float(0.0f);
        return;
      }

      new (&values[i]) float((*data_)[index]);
    });
  }
};

/* get shape keys. */
static ShapeKeyInfo compute_shape_key_info(const Object &object,
                                           const std::string shape_key_name) 
{

  Key *key = BKE_key_from_object(const_cast<Object *>(&object));

  if (key == nullptr) {
    return {};
  }

  KeyBlock *key_block =
      BKE_keyblock_find_name(key, shape_key_name.c_str());

  if (key_block == nullptr) {
    return {};
  }

  return {
    key_block,            /* key block    */
    key_block->totelem,   /* point count  */
    key_block->curval,    /* value        */
    key_block->slidermin, /* range min    */
    key_block->slidermax, /* range max    */
    true                  /* exists       */
  };
}

/* get shape key mesh info */
static ShapeKeyMeshInfo compute_shape_key_mesh_info(const Object &object,
                                                    const ShapeKeyInfo *skape_key_info
                                                    ) {
  
  if (object.type != OB_MESH) {
    return {};
  }

  const Mesh *mesh = id_cast<const Mesh *>(object.data);   
  if (!mesh) {
    return {};
  }
  const int mesh_point_count = mesh->verts_num;

  if (mesh_point_count == 0) {
    return {};
  }

  const Span<MDeformVert> mesh_deforms = mesh->deform_verts(); 

  KeyBlock *keyblock = skape_key_info->key_block;
  if (!keyblock->data) {
    return {};
  }

  const float3 *data = static_cast<const float3 *>(keyblock->data);

  std::shared_ptr<std::vector<float3>> positions = 
    std::make_shared<std::vector<float3>>(
      data, data + skape_key_info->point_count);

  Field<float3> field_position(
      fn::FieldOperation::from(
          std::make_shared<ShapeKeyVectorsFunction>(
              positions, "Get Shape Key Positions"),
          {fn::IndexFieldInput::get_field()}));
      
  /* get the vertex group weights. */
  std::vector<float> vertex_group_weights(mesh_point_count, 1.0f);

  if (keyblock->vgroup[0] != '\0') {
    const int group_index = BKE_object_defgroup_name_index(
        const_cast<Object *>(&object), keyblock->vgroup);

    if (group_index >= 0) {
      std::fill(vertex_group_weights.begin(), vertex_group_weights.end(), 0.0f);

      for (int i = 0; i < mesh_point_count && i < mesh_deforms.size(); i++) {
        const MDeformVert &deform_vert = mesh_deforms[i];

        for (const MDeformWeight &weight :
            Span<MDeformWeight>(deform_vert.dw, deform_vert.totweight))
        {
          if (weight.def_nr == group_index) {
            vertex_group_weights[i] = weight.weight;
            break;
          }
        }
      }
    }
  }

  /* create the Vertex group field. */
  auto weights = std::make_shared<std::vector<float>>(
      std::move(vertex_group_weights));

  Field<float> field_vertex_group(
      fn::FieldOperation::from(
          std::make_shared<ShapeKeyValuesFunctions>(
              weights, "Get Shape Key Weights"),
          {fn::IndexFieldInput::get_field()}));

  return {
    skape_key_info->point_count == mesh_point_count,
    std::move(field_position),
    std::move(field_vertex_group), 
  };
}

/* get shape key curve info */
static ShapeKeyCurveInfo compute_shape_key_curve_info(
    const Object &object,
    const ShapeKeyInfo *shape_key_info)
{
  if (object.type != OB_CURVES_LEGACY || !shape_key_info ||
      !shape_key_info->key_block)
  {
    return {};
  }

  const Curve *curve = id_cast<const Curve *>(object.data);

  if (!curve) {
    return {};
  }

  const KeyBlock *key_block = shape_key_info->key_block;

  if (!key_block->data) {
    return {};
  }

  const int curve_point_count = BKE_nurbList_verts_count(&curve->nurb);

  if (curve_point_count == 0) {
    return {};
  }

  const float *data = static_cast<const float *>(key_block->data);

  auto positions = std::make_shared<std::vector<float3>>();
  auto handles_left = std::make_shared<std::vector<float3>>();
  auto handles_right = std::make_shared<std::vector<float3>>();
  auto radii = std::make_shared<std::vector<float>>();
  auto tilts = std::make_shared<std::vector<float>>();

  positions->reserve(curve_point_count);
  handles_left->reserve(curve_point_count);
  handles_right->reserve(curve_point_count);
  radii->reserve(curve_point_count);
  tilts->reserve(curve_point_count);

  int offset = 0;

  for (const Nurb *nu = static_cast<const Nurb *>(curve->nurb.first());
       nu;
       nu = nu->next)
  {
    if (nu->bezt) {
      for (int i = 0; i < nu->pntsu; i++) {

        // Left handle.
        handles_left->emplace_back(
            data[offset],
            data[offset + 1],
            data[offset + 2]);

        // Control point position.
        positions->emplace_back(
            data[offset + 3],
            data[offset + 4],
            data[offset + 5]);

        // Right handle.
        handles_right->emplace_back(
            data[offset + 6],
            data[offset + 7],
            data[offset + 8]);

        tilts->emplace_back(data[offset + 9]);
        radii->emplace_back(data[offset + 10]);

        offset += KEYELEM_FLOAT_LEN_BEZTRIPLE;
      }
    }
    else if (nu->bp) {
      const int count = nu->pntsu * nu->pntsv;

      for (int i = 0; i < count; i++) {
        positions->emplace_back(
            data[offset],
            data[offset + 1],
            data[offset + 2]);

        handles_left->emplace_back(
            data[offset],
            data[offset + 1],
            data[offset + 2]);

        handles_right->emplace_back(
            data[offset],
            data[offset + 1],
            data[offset + 2]);

        tilts->emplace_back(data[offset + 3]);
        radii->emplace_back(data[offset + 4]);

        offset += KEYELEM_FLOAT_LEN_BPOINT;
      }
    }
  }

  auto make_vector_field = [](
      const std::shared_ptr<std::vector<float3>> &values,
      const char *name)
  {
    return Field<float3>(
        fn::FieldOperation::from(
            std::make_shared<ShapeKeyVectorsFunction>(values, name),
            {fn::IndexFieldInput::get_field()}));
  };

  auto make_float_field = [](
      const std::shared_ptr<std::vector<float>> &values,
      const char *name)
  {
    return Field<float>(
        fn::FieldOperation::from(
            std::make_shared<ShapeKeyValuesFunctions>(values, name),
            {fn::IndexFieldInput::get_field()}));
  };

  return {
      true,
      make_vector_field(positions, "Shape Key Position"),
      make_vector_field(handles_left, "Shape Key Handle Left"),
      make_vector_field(handles_right, "Shape Key Handle Right"),
      make_float_field(radii, "Shape Key Radius"),
      make_float_field(tilts, "Shape Key Tilt"),
  };
}

static void node_geo_exec(GeoNodeExecParams params)
{
  const Object *object = params.extract_input<Object *>("Object"_ustr);
  const std::string shape_key_name =
      params.extract_input<std::string>("Shape Key Name"_ustr);

  auto set_failure = [&]() {
    params.set_default_remaining_outputs();
    params.set_output("Exists"_ustr, false);
  };

  /* get self object is socket input is null */
  if (object == nullptr) {
    object = params.self_object();
  }

  if (object == nullptr) {
    set_failure();
    return;
  }

  if (!ELEM(object->type, OB_MESH, OB_CURVES_LEGACY)) {
    set_failure();
    params.error_message_add(
      NodeWarningType::Error, 
      TIP_("Object type does not support editable Shape Keys.")
    );
    return;
  }


  const ShapeKeyInfo skape_key_info = compute_shape_key_info(*object, shape_key_name);

  if (skape_key_info.exists == false) {
    set_failure();
    return;
  }

  const int mode = params.node().custom1;

  switch (mode) {
    case GEO_NODE_SHAPE_KEY_INFO_MODE_MESH: {
      const ShapeKeyMeshInfo shape_key_mesh_info = 
        compute_shape_key_mesh_info(*object, &skape_key_info);

      if (object->type == OB_CURVES_LEGACY) {
        set_failure();
        params.error_message_add(
          NodeWarningType::Warning, 
          TIP_("The object and the mode do not match.")
        );
        return;
      }
      
      if (shape_key_mesh_info.is_match == false) {
        params.error_message_add(
          NodeWarningType::Info, 
          TIP_("The Shape Key data does not match the geometry element count.")
        );
      }
      params.set_output("Position"_ustr, 
        std::move(shape_key_mesh_info.position));
      params.set_output("Vertex Group"_ustr, 
        std::move(shape_key_mesh_info.vertex_group));
      
      break;
    }

    case GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE: {
      const ShapeKeyCurveInfo shape_key_curve_info = 
        compute_shape_key_curve_info(*object, &skape_key_info);

      if (object->type == OB_MESH) {
        set_failure();
        params.error_message_add(
          NodeWarningType::Warning, 
          TIP_("The object and the mode do not match.")
        );
        return;
      }
      
      if (shape_key_curve_info.is_match == false) {
        params.error_message_add(
          NodeWarningType::Info, 
          TIP_("The Shape Key data does not match the geometry element count.")
        );
      }
      params.set_output("Position"_ustr, 
        std::move(shape_key_curve_info.position));
      params.set_output("Handle Left"_ustr, 
        std::move(shape_key_curve_info.handle_left));
      params.set_output("Handle Right"_ustr, 
        std::move(shape_key_curve_info.handle_right));
      params.set_output("Radius"_ustr, 
        std::move(shape_key_curve_info.radius));
      params.set_output("Tilt"_ustr, 
        std::move(shape_key_curve_info.tilt));
      break;
    }

    default:
      break;
  }

  params.set_output("Value"_ustr, 
    skape_key_info.value);
  params.set_output("Range Min"_ustr, 
    skape_key_info.range_min);
  params.set_output("Range Max"_ustr, 
    skape_key_info.range_min);
  params.set_output("Exists"_ustr, 
    true);
}

static void node_rna(StructRNA *srna)
{
  static const EnumPropertyItem rna_node_geometry_shape_key_info_mode[] = {
    {GEO_NODE_SHAPE_KEY_INFO_MODE_MESH, "MESH", ICON_MESH_DATA, "Mesh", "..."}, 
    {GEO_NODE_SHAPE_KEY_INFO_MODE_CURVE, "CURVE", ICON_CURVE_DATA, "Curve", "..."}, 
     {0, nullptr, 0, nullptr, nullptr},
  } ;

  PropertyRNA *prop_mode = RNA_def_node_enum(srna, 
                                            "mode",
                                            "Mode",
                                            "...",
                                            rna_node_geometry_shape_key_info_mode, 
                                            NOD_inline_enum_accessors(custom1), 
                                            GEO_NODE_SHAPE_KEY_INFO_MODE_MESH);
  RNA_def_property_update_runtime(prop_mode, rna_Node_update_relations);
}

static void node_register()
{
  static bke::bNodeType ntype;
  geo_node_type_base(&ntype, "GeometryNodeShapeKeyInfo"_ustr);
  ntype.ui_name = "Shape Key Info";
  ntype.ui_description =
      "Retrieve information about a Shape Key, including its positions, vertex group weights, and current value";
  ntype.nclass = NODE_CLASS_INPUT;
  ntype.declare = node_declare;
  ntype.initfunc = node_node_init;
  ntype.draw_buttons = node_layout;
  ntype.geometry_node_execute = node_geo_exec;
  ntype.gather_link_search_ops = node_gather_link_search_ops;
  bke::node_register_type(ntype);

  node_rna(ntype.rna_ext.srna);
}
NOD_REGISTER_NODE(node_register)

}  // namespace blender::nodes::node_geo_shape_key_info_cc
