#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Build assets/nodes/texture_nodes_essentials.blend.

Run from a Blender build:

    blender -b --factory-startup -P tools/utils/make_texture_nodes_essentials.py

The output file ships built-in shader node group assets used to populate the
Material Texture Layers add menu (Generator / Mask submenus). Each group is
tagged via ShaderNodeTree.usage with one or more eShaderNodeTreeUsage flags.
"""

import os
import sys

import bpy


OUTPUT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "assets", "nodes", "texture_nodes_essentials.blend",
)


def reset_blend():
    bpy.ops.wm.read_factory_settings(use_empty=True)


# Catalog UUIDs from assets/blender_assets.cats.txt. Assigning them makes the
# groups show up organized in the Asset Browser and in the shader node
# editor's add menu (which lists node group assets by catalog).
CATALOG_ADJUSTMENTS = "5bf6e94e-e35d-4f6c-beb8-2e356debcaf5"
CATALOG_MASKS = "e259c450-baf3-4b56-ae9e-0bc8d9726e5f"


def mark_node_group_asset(group, description, catalog_id):
    group.asset_mark()
    group.asset_data.description = description
    group.asset_data.catalog_id = catalog_id
    group.asset_generate_preview()


def make_combine_bundle(tree, items, x, y):
    """Create a Combine Bundle node with the given (socket_type, name) items."""
    node = tree.nodes.new('NodeCombineBundle')
    node.location = (x, y)
    for socket_type, name in items:
        node.bundle_items.new(socket_type, name)
    return node


def link_to_output_bundle(tree, src_socket):
    out = next(n for n in tree.nodes if n.bl_idname == 'NodeGroupOutput')
    bundle_input = out.inputs[0]
    tree.links.new(src_socket, bundle_input)


def setup_group(name):
    """Create a fresh shader node tree with one Bundle output, return (group, output_node)."""
    if name in bpy.data.node_groups:
        bpy.data.node_groups.remove(bpy.data.node_groups[name])
    group = bpy.data.node_groups.new(name, 'ShaderNodeTree')
    group.interface.new_socket("Bundle", in_out='OUTPUT', socket_type='NodeSocketBundle')
    out_node = group.nodes.new('NodeGroupOutput')
    out_node.location = (300, 0)
    return group


def setup_mask_group(name):
    """Create a fresh shader node tree with one float Factor output. Pure mask
    generators output a plain float, which the mask add flow feeds directly
    into the layer's Mask Stack."""
    if name in bpy.data.node_groups:
        bpy.data.node_groups.remove(bpy.data.node_groups[name])
    group = bpy.data.node_groups.new(name, 'ShaderNodeTree')
    group.interface.new_socket("Factor", in_out='OUTPUT', socket_type='NodeSocketFloat')
    out_node = group.nodes.new('NodeGroupOutput')
    out_node.location = (300, 0)
    return group


def link_to_output(tree, src_socket, output_name):
    out = next(n for n in tree.nodes if n.bl_idname == 'NodeGroupOutput')
    tree.links.new(src_socket, out.inputs[output_name])


def make_math(tree, operation, x, y, value=None):
    node = tree.nodes.new('ShaderNodeMath')
    node.operation = operation
    node.location = (x, y)
    if value is not None:
        node.inputs[1].default_value = value
    return node


def setup_adjustment_group(name, channel="Base Color"):
    """Create an adjustment group skeleton: Bundle input -> Separate Bundle on
    channel ... Combine Bundle on channel -> Bundle output. Returns
    (group, group_input_node, channel_output_socket, channel_input_socket)."""
    group = setup_group(name)
    # Adjustments take a bundle in and emit a bundle out, so the layer above
    # them in the stack receives a transformed version of the layer below.
    group.interface.new_socket("Bundle", in_out='INPUT', socket_type='NodeSocketBundle')

    in_node = group.nodes.new('NodeGroupInput')
    in_node.location = (-600, 0)

    sep = group.nodes.new('NodeSeparateBundle')
    sep.location = (-400, 0)
    sep.bundle_items.new('RGBA', channel)
    group.links.new(in_node.outputs["Bundle"], sep.inputs["Bundle"])

    combine = make_combine_bundle(group, (("RGBA", channel),), 100, 0)
    link_to_output_bundle(group, combine.outputs["Bundle"])
    return group, in_node, sep.outputs[channel], combine.inputs[channel]


def build_brightness_group():
    # The adjusted channel is "Base Color" so it lines up with the channel
    # vocabulary used by Fill layers and the BSDF inputs.
    group, in_node, color_out, color_in = setup_adjustment_group("Brightness")
    factor = group.interface.new_socket(
        "Factor", in_out='INPUT', socket_type='NodeSocketFloat')
    factor.default_value = 1.0

    mul = group.nodes.new('ShaderNodeMix')
    mul.data_type = 'RGBA'
    mul.blend_type = 'MULTIPLY'
    mul.clamp_factor = True
    mul.location = (-150, 0)
    # Set Factor to 1 so the multiply blend uses the full B input.
    mul.inputs[0].default_value = 1.0
    group.links.new(color_out, mul.inputs[6])  # A (color)
    # Multiply by a uniform float read from the group's Factor input. Wrap it
    # into a Combine Color so the Mix Color B input gets a vector value.
    combine_color = group.nodes.new('ShaderNodeCombineColor')
    combine_color.location = (-380, -200)
    group.links.new(in_node.outputs["Factor"], combine_color.inputs["Red"])
    group.links.new(in_node.outputs["Factor"], combine_color.inputs["Green"])
    group.links.new(in_node.outputs["Factor"], combine_color.inputs["Blue"])
    group.links.new(combine_color.outputs["Color"], mul.inputs[7])  # B (color)

    group.links.new(mul.outputs[2], color_in)  # Result

    group.usage = {'TEXTURE_ADJUSTMENT'}
    mark_node_group_asset(
        group,
        "Multiplies the Base Color channel of the bundle by a uniform Factor. "
        "Useful as a quick brightness/tint over the layer below.",
        CATALOG_ADJUSTMENTS,
    )
    return group


def build_hsv_group():
    group, in_node, color_out, color_in = setup_adjustment_group("HSV")
    for name, default in (("Hue", 0.5), ("Saturation", 1.0), ("Value", 1.0)):
        sock = group.interface.new_socket(name, in_out='INPUT', socket_type='NodeSocketFloat')
        sock.default_value = default

    hsv = group.nodes.new('ShaderNodeHueSaturation')
    hsv.location = (-150, 0)
    group.links.new(in_node.outputs["Hue"], hsv.inputs["Hue"])
    group.links.new(in_node.outputs["Saturation"], hsv.inputs["Saturation"])
    group.links.new(in_node.outputs["Value"], hsv.inputs["Value"])
    group.links.new(color_out, hsv.inputs["Color"])
    group.links.new(hsv.outputs["Color"], color_in)

    group.usage = {'TEXTURE_ADJUSTMENT'}
    mark_node_group_asset(
        group,
        "Shifts hue, saturation and value of the Base Color channel of the "
        "layers below.",
        CATALOG_ADJUSTMENTS,
    )
    return group


def build_color_correction_group():
    group, in_node, color_out, color_in = setup_adjustment_group("Color Correction")
    for name, default in (("Brightness", 0.0), ("Contrast", 0.0), ("Gamma", 1.0)):
        sock = group.interface.new_socket(name, in_out='INPUT', socket_type='NodeSocketFloat')
        sock.default_value = default

    bc = group.nodes.new('ShaderNodeBrightContrast')
    bc.location = (-250, 0)
    group.links.new(color_out, bc.inputs["Color"])
    group.links.new(in_node.outputs["Brightness"], bc.inputs["Bright"])
    group.links.new(in_node.outputs["Contrast"], bc.inputs["Contrast"])

    gamma = group.nodes.new('ShaderNodeGamma')
    gamma.location = (-80, 0)
    group.links.new(bc.outputs["Color"], gamma.inputs["Color"])
    group.links.new(in_node.outputs["Gamma"], gamma.inputs["Gamma"])
    group.links.new(gamma.outputs["Color"], color_in)

    group.usage = {'TEXTURE_ADJUSTMENT'}
    mark_node_group_asset(
        group,
        "Brightness, contrast and gamma correction of the Base Color channel "
        "of the layers below.",
        CATALOG_ADJUSTMENTS,
    )
    return group


def build_udim_tile_group():
    group = setup_mask_group("UDIM Tile")
    tile = group.interface.new_socket("Tile", in_out='INPUT', socket_type='NodeSocketFloat')
    tile.default_value = 1001

    in_node = group.nodes.new('NodeGroupInput')
    in_node.location = (-800, -200)

    coord = group.nodes.new('ShaderNodeTexCoord')
    coord.location = (-800, 100)
    xyz = group.nodes.new('ShaderNodeSeparateXYZ')
    xyz.location = (-620, 100)
    group.links.new(coord.outputs["UV"], xyz.inputs["Vector"])

    floor_x = make_math(group, 'FLOOR', -450, 180)
    group.links.new(xyz.outputs["X"], floor_x.inputs[0])
    floor_y = make_math(group, 'FLOOR', -450, 20)
    group.links.new(xyz.outputs["Y"], floor_y.inputs[0])

    # UDIM number = 1001 + floor(u) + 10 * floor(v).
    row = make_math(group, 'MULTIPLY', -300, 20, value=10.0)
    group.links.new(floor_y.outputs["Value"], row.inputs[0])
    tile_offset = make_math(group, 'ADD', -160, 100)
    group.links.new(floor_x.outputs["Value"], tile_offset.inputs[0])
    group.links.new(row.outputs["Value"], tile_offset.inputs[1])
    number = make_math(group, 'ADD', -20, 100, value=1001.0)
    group.links.new(tile_offset.outputs["Value"], number.inputs[0])

    compare = make_math(group, 'COMPARE', 130, 0)
    compare.inputs[2].default_value = 0.5  # Epsilon: match one tile exactly.
    group.links.new(number.outputs["Value"], compare.inputs[0])
    group.links.new(in_node.outputs["Tile"], compare.inputs[1])
    link_to_output(group, compare.outputs["Value"], "Factor")

    group.usage = {'MASK_GENERATOR'}
    group.default_mask_blend_type = 'ADD'
    mark_node_group_asset(
        group,
        "Mask selecting a single UDIM tile of the UV map by its tile number "
        "(1001, 1002, ...).",
        CATALOG_MASKS,
    )
    return group


def main():
    reset_blend()
    build_brightness_group()
    build_hsv_group()
    build_color_correction_group()
    build_udim_tile_group()

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=OUTPUT_PATH, compress=True)
    print("Wrote: " + OUTPUT_PATH)


if __name__ == "__main__":
    main()
    # Force a clean exit when running via `blender -b -P`.
    sys.exit(0)
