# SPDX-FileCopyrightText: 2026 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

# Tests for multi-layer / multi-pass image support in the shader Image Texture
# node (design phase 15): a multi-layer EXR assigned to the node must expose one
# output socket per pass instead of the fixed Color/Alpha pair.
#
# ./bin/blender --factory-startup --python tests/python/bl_node_multilayer_image.py

import os
import sys
import tempfile
import unittest

import bpy

import numpy as np
import OpenImageIO as oiio


# Distinct constant values per pass, so a render that sampled the wrong pass is
# detectable.
COMBINED_RGBA = (0.1, 0.2, 0.3, 1.0)
DEPTH_VALUE = 0.5

# Value painted into the Base Color pass of a generated multi-layer image, and
# the fill of its Roughness pass (the Principled BSDF input default).
PAINTED_RGBA = (0.25, 0.5, 0.75, 1.0)
ROUGHNESS_VALUE = 0.5


def write_multilayer_exr(filepath, fill=False):
    """Write a multi-layer EXR with a Combined RGBA pass and a Depth pass.

    When ``fill`` is set, each pass is filled with its distinct constant value.
    """
    channels = [
        "ViewLayer.Combined.R", "ViewLayer.Combined.G",
        "ViewLayer.Combined.B", "ViewLayer.Combined.A",
        "ViewLayer.Depth.Z",
    ]
    spec = oiio.ImageSpec(16, 16, len(channels), "float")
    spec.channelnames = channels
    pixels = np.zeros((16, 16, len(channels)), dtype=np.float32)
    if fill:
        pixels[:, :, 0:4] = COMBINED_RGBA
        pixels[:, :, 4] = DEPTH_VALUE
    out = oiio.ImageOutput.create(filepath)
    assert out is not None, "could not create EXR output"
    out.open(filepath, spec)
    out.write_image(pixels)
    out.close()


def read_center_pixel(filepath):
    """Return the center pixel of an image file as a tuple of floats."""
    src = oiio.ImageInput.open(filepath)
    assert src is not None, f"could not open {filepath}"
    spec = src.spec()
    pixels = src.read_image(0, 0, 0, spec.nchannels, "float")
    src.close()
    pixels = np.array(pixels).reshape(spec.height, spec.width, spec.nchannels)
    return tuple(float(c) for c in pixels[spec.height // 2, spec.width // 2])


class MultiLayerImageTextureNodeTest(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.exr_path = os.path.join(self.tmpdir, "multilayer.exr")
        write_multilayer_exr(self.exr_path)

    def _make_node(self, image):
        mat = bpy.data.materials.new("MultiLayerTest")
        mat.use_nodes = True
        node = mat.node_tree.nodes.new("ShaderNodeTexImage")
        node.image = image
        bpy.context.view_layer.update()
        mat.node_tree.update_tag()
        bpy.context.view_layer.update()
        return node

    def test_multi_layer_declares_per_pass_sockets(self):
        image = bpy.data.images.load(self.exr_path)
        # Acquiring a buffer initializes the multi-layer structure.
        _ = image.pixels[0]
        self.assertEqual(image.type, 'MULTILAYER')

        node = self._make_node(image)
        outputs = [s.name for s in node.outputs]

        self.assertIn("Combined", outputs)
        self.assertIn("Depth", outputs)
        # A synthetic Alpha is generated from the RGBA Combined pass.
        self.assertIn("Alpha", outputs)
        self.assertTrue(node.has_layers)

    def test_single_layer_keeps_color_alpha(self):
        # A plain generated image must keep the fixed Color/Alpha outputs.
        image = bpy.data.images.new("Plain", 16, 16)
        node = self._make_node(image)
        outputs = [s.name for s in node.outputs]
        self.assertEqual(outputs, ["Color", "Alpha"])
        self.assertFalse(node.has_layers)


class MultiLayerImageTextureRenderTest(unittest.TestCase):
    """Renders a pass of a multi-layer image. The shader node inliner (design
    phase 16) must expand the multi-pass Image Texture node into single-pass
    copies so the engine samples the requested pass. Both Cycles and the GPU
    engines inline shader node trees, so the path is shared."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.exr_path = os.path.join(self.tmpdir, "multilayer.exr")
        write_multilayer_exr(self.exr_path, fill=True)

    # Socket type for a Separate Bundle item, per pass.
    _BUNDLE_ITEM_TYPE = {"Depth": 'FLOAT', "Combined": 'RGBA'}

    def _render_pass(self, pass_socket_name, engine, via_bundle=False):
        bpy.ops.wm.read_factory_settings(use_empty=False)

        image = bpy.data.images.load(self.exr_path)
        _ = image.pixels[0]
        # Pack the image so Cycles uses the builtin loader, which resolves the
        # pinned pass through the Image.layers catalog.
        image.pack()

        mat = bpy.data.materials.new("RenderTest")
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        nodes.clear()

        tex = nodes.new("ShaderNodeTexImage")
        tex.image = image
        bpy.context.view_layer.update()
        mat.node_tree.update_tag()
        bpy.context.view_layer.update()

        emission = nodes.new("ShaderNodeEmission")
        output = nodes.new("ShaderNodeOutputMaterial")
        links.new(emission.outputs["Emission"], output.inputs["Surface"])

        if via_bundle:
            # Route through the Bundle output and a Separate Bundle node, so the
            # inliner has to expand the bundle into single-pass nodes (phase 15b).
            separate = nodes.new("NodeSeparateBundle")
            separate.bundle_items.new(self._BUNDLE_ITEM_TYPE[pass_socket_name],
                                      pass_socket_name)
            links.new(tex.outputs["Passes"], separate.inputs["Bundle"])
            pass_socket = separate.outputs[pass_socket_name]
        else:
            pass_socket = tex.outputs[pass_socket_name]
        links.new(pass_socket, emission.inputs["Color"])

        cube = bpy.data.objects["Cube"]
        cube.data.materials.clear()
        cube.data.materials.append(mat)

        scene = bpy.context.scene
        scene.render.engine = engine
        if engine == 'CYCLES':
            scene.cycles.samples = 1
        else:
            scene.eevee.taa_render_samples = 1
        scene.render.resolution_x = 16
        scene.render.resolution_y = 16
        scene.render.image_settings.media_type = 'IMAGE'
        scene.render.image_settings.file_format = 'OPEN_EXR'
        scene.view_settings.view_transform = 'Raw'
        out_path = os.path.join(self.tmpdir, f"render_{engine}_{pass_socket_name}.exr")
        scene.render.filepath = out_path
        bpy.ops.render.render(write_still=True)
        return read_center_pixel(out_path)

    def _check_pass(self, pass_socket_name, expected, engine, via_bundle=False):
        rgba = self._render_pass(pass_socket_name, engine, via_bundle)
        # A GPU engine on a headless host without a usable GPU renders black;
        # treat that as "GPU unavailable" and skip rather than fail.
        if engine != 'CYCLES' and max(rgba[:3]) < 1e-4 and max(expected) > 1e-4:
            self.skipTest("GPU engine produced an empty render (no usable GPU)")
        for got, want in zip(rgba[:3], expected[:3]):
            self.assertAlmostEqual(got, want, delta=0.05)

    def test_cycles_render_depth_pass(self):
        # The Depth pass is a constant 0.5; routed through Emission Color it
        # renders as mid-grey.
        self._check_pass("Depth", (DEPTH_VALUE,) * 3, 'CYCLES')

    def test_cycles_render_combined_pass(self):
        self._check_pass("Combined", COMBINED_RGBA, 'CYCLES')

    def test_eevee_render_depth_pass(self):
        self._check_pass("Depth", (DEPTH_VALUE,) * 3, 'BLENDER_EEVEE')

    def test_eevee_render_combined_pass(self):
        self._check_pass("Combined", COMBINED_RGBA, 'BLENDER_EEVEE')

    def test_cycles_render_depth_via_bundle(self):
        # The same pass reached through the Bundle output + a Separate Bundle
        # node (design phase 15b).
        self._check_pass("Depth", (DEPTH_VALUE,) * 3, 'CYCLES', via_bundle=True)

    def test_eevee_render_depth_via_bundle(self):
        self._check_pass("Depth", (DEPTH_VALUE,) * 3, 'BLENDER_EEVEE', via_bundle=True)


class GeneratedMultiLayerImageRenderTest(unittest.TestCase):
    """Renders a pass of a generated multi-layer image, the image behind a Paint
    texture layer. Its layer catalog is authored rather than read from a file,
    so nothing ever loads a multi-layer EXR for it; the pass must still resolve
    for both engines, otherwise Cycles renders it as a missing texture."""

    RES = 16

    def _make_paint_layer(self):
        """Build a material with a Paint texture layer, its image painted."""
        bpy.ops.wm.read_factory_settings(use_empty=False)
        bpy.context.preferences.experimental.use_texture_layers = True

        mat = bpy.data.materials.new("PaintLayer")
        mat.use_nodes = True
        cube = bpy.data.objects["Cube"]
        cube.data.materials.clear()
        cube.data.materials.append(mat)

        with bpy.context.temp_override(material=mat, object=cube, active_object=cube):
            # The default layer bootstraps the stack wired into the BSDF, so the
            # Paint layer's passes are filled with the BSDF input defaults.
            bpy.ops.material.texture_layer_add_default()
            bpy.ops.material.texture_layer_add_paint(
                name="Paint", width=self.RES, height=self.RES)

        image = bpy.data.images["Paint"]
        self.assertEqual(image.source, 'GENERATED')
        self.assertEqual(image.type, 'MULTILAYER')
        # Paint a distinct color into the first pass (Base Color).
        image.pixels = list(PAINTED_RGBA) * (self.RES * self.RES)
        return mat

    def _render_pass(self, pass_socket_name, engine):
        mat = self._make_paint_layer()
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        tex = next(n for n in nodes if n.bl_idname == 'ShaderNodeTexImage')
        emission = nodes.new("ShaderNodeEmission")
        output = next(n for n in nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
        links.new(tex.outputs[pass_socket_name], emission.inputs["Color"])
        links.new(emission.outputs["Emission"], output.inputs["Surface"])

        scene = bpy.context.scene
        scene.render.engine = engine
        if engine == 'CYCLES':
            scene.cycles.samples = 1
        else:
            scene.eevee.taa_render_samples = 1
        scene.render.resolution_x = self.RES
        scene.render.resolution_y = self.RES
        scene.render.image_settings.media_type = 'IMAGE'
        scene.render.image_settings.file_format = 'OPEN_EXR'
        scene.view_settings.view_transform = 'Raw'
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, f"render_{engine}_{pass_socket_name}.exr")
            scene.render.filepath = out_path
            bpy.ops.render.render(write_still=True)
            return read_center_pixel(out_path)

    def _check_pass(self, pass_socket_name, expected, engine):
        rgba = self._render_pass(pass_socket_name, engine)
        if engine != 'CYCLES' and max(rgba[:3]) < 1e-4 and max(expected) > 1e-4:
            self.skipTest("GPU engine produced an empty render (no usable GPU)")
        for got, want in zip(rgba[:3], expected[:3]):
            self.assertAlmostEqual(got, want, delta=0.05)

    def test_cycles_render_painted_pass(self):
        self._check_pass("Base Color", PAINTED_RGBA, 'CYCLES')

    def test_cycles_render_roughness_pass(self):
        # A different pass of the same image, to catch resolving to the wrong one.
        self._check_pass("Roughness", (ROUGHNESS_VALUE,) * 3, 'CYCLES')

    def test_eevee_render_painted_pass(self):
        self._check_pass("Base Color", PAINTED_RGBA, 'BLENDER_EEVEE')

    def test_eevee_render_roughness_pass(self):
        self._check_pass("Roughness", (ROUGHNESS_VALUE,) * 3, 'BLENDER_EEVEE')


def main():
    # Drop arguments meant for Blender so unittest does not choke on them.
    argv = [sys.argv[0]]
    if "--" in sys.argv:
        argv += sys.argv[sys.argv.index("--") + 1:]
    unittest.main(argv=argv)


if __name__ == "__main__":
    main()
