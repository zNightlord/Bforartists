/* SPDX-FileCopyrightText: 2005 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup shdnodes
 */

#include "node_shader_util.hh"
#include "node_util.hh"

#include "BLI_listbase.hh"
#include "BLI_string.hh"
#include "BLI_string_ref.hh"

#include "BKE_image.hh"
#include "BKE_node_runtime.hh"
#include "BKE_texture.h"

#include "DNA_image_types.h"
#include "DNA_scene_enums.h"

#include "IMB_colormanagement.hh"

#include "DEG_depsgraph_query.hh"

namespace blender {

namespace nodes::node_shader_tex_image_cc {

/* Declares the fixed Color/Alpha outputs used for single-layer images and when
 * no multi-layer image is assigned. */
static void declare_single_layer(NodeDeclarationBuilder &b)
{
  b.add_output<decl::Color>("Color"_ustr).no_muted_links();
  b.add_output<decl::Float>("Alpha"_ustr).no_muted_links();
}

/* Declares an output that matches the type of the given pass. Mirrors the
 * compositor Image node's per-pass declaration (design §15). */
static void declare_pass(NodeDeclarationBuilder &b, const ImagePass &pass)
{
  const UString name(pass.name);
  switch (pass.channels_num) {
    case 1:
      b.add_output<decl::Float>(name).no_muted_links();
      return;
    case 2:
      b.add_output<decl::Vector>(name).dimensions(2).no_muted_links();
      return;
    case 3:
      if (STR_ELEM(pass.chan_id, "RGB", "rgb")) {
        b.add_output<decl::Color>(name).no_muted_links();
        return;
      }
      b.add_output<decl::Vector>(name).dimensions(3).no_muted_links();
      return;
    case 4:
      if (STR_ELEM(pass.chan_id, "RGBA", "rgba")) {
        b.add_output<decl::Color>(name).no_muted_links();
        return;
      }
      b.add_output<decl::Vector>(name).dimensions(4).no_muted_links();
      return;
  }
  BLI_assert_unreachable();
}

/* Declares one output socket per pass of the selected layer, plus a synthetic
 * Alpha output when the layer has none of its own, and a Bundle output carrying
 * every pass (design §15a, §15b). Returns false when the layer is unknown, in
 * which case the caller falls back to the single-layer declaration. */
static bool declare_multi_layer(NodeDeclarationBuilder &b, Image *image, const ImageUser *iuser)
{
  ImageLayer *layer = static_cast<ImageLayer *>(BLI_findlink(&image->layers, iuser->layer));
  if (!layer) {
    return false;
  }

  bool has_alpha_pass = false;
  for (const ImagePass &pass : layer->passes) {
    if (StringRef(pass.name) == "Alpha") {
      has_alpha_pass = true;
      break;
    }
  }

  for (const ImagePass &pass : layer->passes) {
    declare_pass(b, pass);

    /* Generate an Alpha output from the Combined pass when no dedicated Alpha
     * pass exists, matching the compositor Image node. */
    if (!has_alpha_pass && StringRef(pass.name) == RE_PASSNAME_COMBINED &&
        pass.channels_num == 4 && StringRef(pass.chan_id) == "RGBA")
    {
      b.add_output<decl::Float>("Alpha"_ustr).no_muted_links();
    }
  }

  /* Bundle of every pass, extracted downstream with a Separate Bundle node. */
  b.add_output<decl::Bundle>("Passes"_ustr);
  return true;
}

/* Re-declares the outputs that already exist on the node. Used when the
 * declaration is rebuilt for a reason unrelated to the image data, so the
 * existing sockets and their links are retained unchanged. */
static void declare_existing(NodeDeclarationBuilder &b, const bNode &node)
{
  for (const bNodeSocket &output : node.outputs) {
    if (output.type == SOCK_VECTOR) {
      const int dimensions = output.default_value_typed<bNodeSocketValueVector>()->dimensions;
      b.add_output<decl::Vector>(UString(output.name)).dimensions(dimensions).no_muted_links();
    }
    else {
      b.add_output(eNodeSocketDatatype(output.type), UString(output.name)).no_muted_links();
    }
  }
}

/* The image structure of a multi-layer image is only known once an image buffer
 * has been acquired. Acquire and immediately release a dummy buffer so the
 * Image.layers catalog is populated as a side effect. */
static void prepare_image(Image *image, const ImageUser *iuser)
{
  const int start_frame_offset = BKE_image_sequence_guess_offset(image);
  ImageUser initial_frame_iuser = *iuser;
  initial_frame_iuser.framenr = start_frame_offset;

  ImBuf *ibuf = BKE_image_acquire_ibuf(image, &initial_frame_iuser, nullptr);
  BKE_image_release_ibuf(image, ibuf, nullptr);
}

/* A multi-pass node carries a Bundle output, which is not valid on a function
 * node (whose sockets are all forced to a dynamic structure type). So a node is
 * only declared as a function node when it exposes the plain Color/Alpha pair. */
static void declare_single_layer_function(NodeDeclarationBuilder &b)
{
  b.is_function_node();
  declare_single_layer(b);
}

static void sh_node_tex_image_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Vector>("Vector"_ustr).default_input_type(NODE_DEFAULT_INPUT_POSITION_FIELD);

  const bNode *node = b.node_or_null();
  if (!node) {
    declare_single_layer_function(b);
    return;
  }

  Image *image = reinterpret_cast<Image *>(node->id);
  const NodeTexImage *tex = static_cast<const NodeTexImage *>(node->storage);
  if (!image || !tex) {
    declare_single_layer_function(b);
    return;
  }

  /* A single-pass node (produced by shader node inlining) always exposes the
   * standard Color/Alpha outputs; its pass is pinned in iuser (design §16). */
  if (tex->flag & SHD_TEX_IMAGE_SINGLE_PASS) {
    declare_single_layer_function(b);
    return;
  }

  /* Only the Image/Image User data affects the declared sockets. */
  if (!(node->runtime->update & NODE_UPDATE_ID)) {
    /* Keep the function-node flag consistent with the existing sockets: a
     * multi-pass node must not be a function. Use the same structural gate as
     * the multi-pass declaration path below, so both agree. */
    const bool is_multi_pass = BKE_image_is_multilayer(image) &&
                               BLI_findlink(&image->layers, tex->iuser.layer) != nullptr;
    if (!is_multi_pass) {
      b.is_function_node();
    }
    declare_existing(b, *node);
    return;
  }

  prepare_image(image, &tex->iuser);

  if (!BKE_image_is_multilayer(image) || !declare_multi_layer(b, image, &tex->iuser)) {
    declare_single_layer_function(b);
    return;
  }
}

static void node_shader_init_tex_image(bNodeTree * /*ntree*/, bNode *node)
{
  NodeTexImage *tex = MEM_new<NodeTexImage>(__func__);
  BKE_texture_mapping_default(&tex->base.tex_mapping, TEXMAP_TYPE_POINT);
  BKE_texture_colormapping_default(&tex->base.color_mapping);
  BKE_imageuser_default(&tex->iuser);

  node->storage = tex;
}

static int node_shader_gpu_tex_image(GPUMaterial *mat,
                                     bNode *node,
                                     bNodeExecData * /*execdata*/,
                                     GPUNodeStack *in,
                                     GPUNodeStack *out)
{
  Image *ima = id_cast<Image *>(node->id);
  NodeTexImage *tex = static_cast<NodeTexImage *>(node->storage);

  /* We get the image user from the original node, since GPU image keeps
   * a pointer to it and the dependency refreshes the original. */
  bNode *node_original = node->runtime->original ? node->runtime->original : node;
  NodeTexImage *tex_original = static_cast<NodeTexImage *>(node_original->storage);
  ImageUser *iuser = &tex_original->iuser;

  if (!ima) {
    return GPU_stack_link(mat, node, "node_tex_image_empty", in, out);
  }

  GPUNodeLink **texco = &in[0].link;
  if (!*texco) {
    *texco = GPU_attribute(mat, CD_AUTO_FROM_NAME, "");
    node_shader_gpu_bump_tex_coord(mat, node, texco);
  }

  node_shader_gpu_tex_mapping(mat, node, in, out);

  GPUSamplerState sampler_state = GPUSamplerState::default_sampler();

  switch (tex->extension) {
    case SHD_IMAGE_EXTENSION_EXTEND:
      sampler_state.extend_x = GPU_SAMPLER_EXTEND_MODE_EXTEND;
      sampler_state.extend_yz = GPU_SAMPLER_EXTEND_MODE_EXTEND;
      break;
    case SHD_IMAGE_EXTENSION_REPEAT:
      sampler_state.extend_x = GPU_SAMPLER_EXTEND_MODE_REPEAT;
      sampler_state.extend_yz = GPU_SAMPLER_EXTEND_MODE_REPEAT;
      break;
    case SHD_IMAGE_EXTENSION_CLIP:
      sampler_state.extend_x = GPU_SAMPLER_EXTEND_MODE_CLAMP_TO_BORDER;
      sampler_state.extend_yz = GPU_SAMPLER_EXTEND_MODE_CLAMP_TO_BORDER;
      break;
    case SHD_IMAGE_EXTENSION_MIRROR:
      sampler_state.extend_x = GPU_SAMPLER_EXTEND_MODE_MIRRORED_REPEAT;
      sampler_state.extend_yz = GPU_SAMPLER_EXTEND_MODE_MIRRORED_REPEAT;
      break;
    default:
      break;
  }

  if (tex->interpolation != SHD_INTERP_CLOSEST) {
    /* TODO(fclem): For now assume mipmap is always enabled. */
    /* Setting the GPU_SAMPLER_FILTERING_ANISOTROPIC_ENABLE enables anisotropic filtering. The
     * exact number of samples are being determined at bind time by the engine.
     * See #blender::draw::PassBase<T>::material_set */
    sampler_state.filtering = GPU_SAMPLER_FILTERING_ANISOTROPIC_ENABLE |
                              GPU_SAMPLER_FILTERING_LINEAR | GPU_SAMPLER_FILTERING_MIPMAP;
  }
  const bool use_cubic = ELEM(tex->interpolation, SHD_INTERP_CUBIC, SHD_INTERP_SMART);

  /* Only use UDIM tiles if projection is flat.
   * Otherwise treat the first tile as a single image. (See #141776). */
  const bool use_udim = ima->source == IMA_SRC_TILED && tex->projection == SHD_PROJ_FLAT;
  if (use_udim) {
    const char *gpu_node_name = use_cubic ? "node_tex_tile_cubic" : "node_tex_tile_linear";
    GPUNodeLink *gpu_image, *gpu_image_tile_mapping;
    GPU_image_tiled(mat, ima, iuser, sampler_state, &gpu_image, &gpu_image_tile_mapping);
    /* UDIM tiles needs a `sampler2DArray` and `sampler1DArray` for tile mapping. */
    GPU_stack_link(mat, node, gpu_node_name, in, out, gpu_image, gpu_image_tile_mapping);
  }
  else {
    const char *gpu_node_name = use_cubic ? "node_tex_image_cubic" : "node_tex_image_linear";

    switch (tex->projection) {
      case SHD_PROJ_FLAT: {
        GPUNodeLink *gpu_image = GPU_image(mat, ima, iuser, sampler_state);
        GPU_stack_link(mat, node, gpu_node_name, in, out, gpu_image);
        break;
      }
      case SHD_PROJ_BOX: {
        gpu_node_name = use_cubic ? "tex_box_sample_cubic" : "tex_box_sample_linear";
        GPUNodeLink *vnor, *wnor, *col1, *col2, *col3;
        GPUNodeLink *blend = GPU_uniform(&tex->projection_blend);
        GPUNodeLink *gpu_image = GPU_image(mat, ima, iuser, sampler_state);
        GPU_link(mat, "world_normals_get", &vnor);
        GPU_link(mat, "normal_transform_world_to_object", vnor, &wnor);
        GPU_link(mat, gpu_node_name, in[0].link, wnor, gpu_image, &col1, &col2, &col3);
        GPU_link(mat, "tex_box_blend", wnor, col1, col2, col3, blend, &out[0].link, &out[1].link);
        break;
      }
      case SHD_PROJ_SPHERE: {
        /* This projection is known to have a derivative discontinuity.
         * Hide it by turning off mipmapping. */
        sampler_state.disable_filtering_flag(GPU_SAMPLER_FILTERING_MIPMAP);
        GPUNodeLink *gpu_image = GPU_image(mat, ima, iuser, sampler_state);
        GPU_link(mat, "point_texco_remap_square", *texco, texco);
        GPU_link(mat, "point_map_to_sphere", *texco, texco);
        GPU_stack_link(mat, node, gpu_node_name, in, out, gpu_image);
        break;
      }
      case SHD_PROJ_TUBE: {
        /* This projection is known to have a derivative discontinuity.
         * Hide it by turning off mipmapping. */
        sampler_state.disable_filtering_flag(GPU_SAMPLER_FILTERING_MIPMAP);
        GPUNodeLink *gpu_image = GPU_image(mat, ima, iuser, sampler_state);
        GPU_link(mat, "point_texco_remap_square", *texco, texco);
        GPU_link(mat, "point_map_to_tube", *texco, texco);
        GPU_stack_link(mat, node, gpu_node_name, in, out, gpu_image);
        break;
      }
    }
  }

  if (out[0].hasoutput && out[0].link) {
    if (ELEM(ima->alpha_mode, IMA_ALPHA_IGNORE, IMA_ALPHA_CHANNEL_PACKED) ||
        IMB_colormanagement_space_name_is_data(ima->colorspace_settings.name))
    {
      /* Don't let alpha affect color output in these cases. */
      GPU_link(mat, "color_alpha_clear", out[0].link, &out[0].link);
    }
    else {
      /* Output premultiplied alpha depending on alpha socket usage. This makes
       * it so that if we blend the color with a transparent shader using alpha as
       * a factor, we don't multiply alpha into the color twice. And if we do
       * not, then there will be no artifacts from zero alpha areas. */
      if (ima->alpha_mode == IMA_ALPHA_PREMUL) {
        if (out[1].hasoutput) {
          GPU_link(mat, "color_alpha_unpremultiply", out[0].link, &out[0].link);
        }
        else {
          GPU_link(mat, "color_alpha_clear", out[0].link, &out[0].link);
        }
      }
      else {
        if (out[1].hasoutput) {
          GPU_link(mat, "color_alpha_clear", out[0].link, &out[0].link);
        }
        else {
          GPU_link(mat, "color_alpha_premultiply", out[0].link, &out[0].link);
        }
      }
    }
  }

  return true;
}

NODE_SHADER_MATERIALX_BEGIN
#ifdef WITH_MATERIALX
{
  /* Getting node name for Color output. This name will be used for <image> node. */
  std::string image_node_name = node_name("Color");

  NodeItem res = graph_.get_node(image_node_name);
  if (!res.node) {
    res = val(MaterialX::Color4(1.0f, 0.0f, 1.0f, 1.0f));

    Image *image = (Image *)node_->id;
    if (image) {
      NodeTexImage *tex_image = static_cast<NodeTexImage *>(node_->storage);

      std::string image_path = image->id.name;
      if (graph_.export_params.image_fn) {
        Scene *scene = DEG_get_input_scene(graph_.depsgraph);
        Main *bmain = DEG_get_bmain(graph_.depsgraph);
        image_path = graph_.export_params.image_fn(bmain, scene, image, &tex_image->iuser);
      }

      NodeItem vector = get_input_link("Vector", NodeItem::Type::Vector2);
      if (!vector) {
        vector = texcoord_node();
      }
      /* TODO: add math to vector depending of tex_image->projection */

      std::string filtertype;
      switch (tex_image->interpolation) {
        case SHD_INTERP_LINEAR:
          filtertype = "linear";
          break;
        case SHD_INTERP_CLOSEST:
          filtertype = "closest";
          break;
        case SHD_INTERP_CUBIC:
        case SHD_INTERP_SMART:
          filtertype = "cubic";
          break;
        default:
          BLI_assert_unreachable();
      }
      std::string addressmode;
      switch (tex_image->extension) {
        case SHD_IMAGE_EXTENSION_REPEAT:
          addressmode = "periodic";
          break;
        case SHD_IMAGE_EXTENSION_EXTEND:
          addressmode = "clamp";
          break;
        case SHD_IMAGE_EXTENSION_CLIP:
          addressmode = "constant";
          break;
        case SHD_IMAGE_EXTENSION_MIRROR:
          addressmode = "mirror";
          break;
        default:
          BLI_assert_unreachable();
      }

      NodeItem::Type node_type = NodeItem::Type::Color4;
      const char *node_colorspace = nullptr;

      const char *image_colorspace = image->colorspace_settings.name;
      if (IMB_colormanagement_space_name_is_data(image_colorspace)) {
        node_type = NodeItem::Type::Vector4;
      }
      else if (IMB_colormanagement_space_name_is_scene_linear(image_colorspace)) {
        node_colorspace = "lin_rec709";
      }
      else if (IMB_colormanagement_space_name_is_srgb(image_colorspace)) {
        node_colorspace = "srgb_texture";
      }

      res = create_node("image",
                        node_type,
                        {{"texcoord", vector},
                         {"filtertype", val(filtertype)},
                         {"uaddressmode", val(addressmode)},
                         {"vaddressmode", val(addressmode)}});
      res.set_input("file", image_path, NodeItem::Type::Filename);
      res.node->setName(image_node_name);
      if (node_colorspace) {
        res.node->setAttribute("colorspace", node_colorspace);
      }
    }
  }

  if (STREQ(socket_out_->identifier, "Alpha")) {
    res = res[3];
  }
  return res;
}
#endif
NODE_SHADER_MATERIALX_END

}  // namespace nodes::node_shader_tex_image_cc

void register_node_type_sh_tex_image()
{
  namespace file_ns = nodes::node_shader_tex_image_cc;

  static bke::bNodeType ntype;

  sh_node_type_base(&ntype, "ShaderNodeTexImage"_ustr, SH_NODE_TEX_IMAGE);
  ntype.ui_name = "Image Texture";
  ntype.ui_description = "Sample an image file as a texture";
  ntype.enum_name_legacy = "TEX_IMAGE";
  ntype.nclass = NODE_CLASS_TEXTURE;
  ntype.texture_layer_usage = SHADER_NODE_TREE_USAGE_TEXTURE_GENERATOR |
                              SHADER_NODE_TREE_USAGE_MASK_GENERATOR;
  ntype.declare = file_ns::sh_node_tex_image_declare;
  ntype.initfunc = file_ns::node_shader_init_tex_image;
  bke::node_type_storage(
      ntype, "NodeTexImage", node_free_standard_storage, node_copy_standard_storage);
  ntype.gpu_fn = file_ns::node_shader_gpu_tex_image;
  ntype.labelfunc = node_image_label;
  ntype.default_width = bke::NodeWidth::_240;
  ntype.materialx_fn = file_ns::node_shader_materialx;

  bke::node_register_type(ntype);
}

}  // namespace blender
