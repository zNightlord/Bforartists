# SPDX-FileCopyrightText: 2026 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
This file does not run anything, its methods are accessed for tests by ``run_blender_setup.py``.
"""
import modules.ui_test_utils as ui


def _author_multilayer_exr(path, layer_names):
    """Write a multi-layer EXR with one "Combined" pass per named layer."""
    from OpenImageIO import ImageSpec, ImageOutput
    import numpy as np

    width = height = 4
    channels = []
    for layer_name in layer_names:
        channels += [layer_name + ".Combined." + c for c in ("R", "G", "B")]

    spec = ImageSpec(width, height, len(channels), "half")
    spec.channelnames = tuple(channels)

    out = ImageOutput.create(path)
    out.open(path, spec)
    out.write_image(np.full((height, width, len(channels)), 0.5, dtype="float16"))
    out.close()


def test_image_layer_name_tracking():
    """A multi-layer image's layer selection follows the layer by name when the
    file is reloaded with a layer inserted ahead of the selected one."""
    import bpy
    import os
    import tempfile

    e, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    exr_path = os.path.join(tmp_dir, "multilayer.exr")
    # Blender sorts EXR layers by name, so "Banana" is the first layer here.
    _author_multilayer_exr(exr_path, ["Banana", "Zebra"])

    # Show the image in an image editor so its ImageUser resolves on redraw.
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    space = area.spaces.active
    image = bpy.data.images.load(exr_path)
    space.image = image
    yield
    yield

    iuser = space.image_user
    t.assertEqual([layer.name for layer in image.layers], ["Banana", "Zebra"])
    # The default selection is the first layer; its name is now recorded.
    t.assertEqual(iuser.multilayer_layer, 0)
    t.assertEqual(iuser.multilayer_layer_name, "Banana")

    # Reload with a layer inserted alphabetically before the selected one.
    _author_multilayer_exr(exr_path, ["Apple", "Banana", "Zebra"])
    image.reload()
    yield
    yield

    t.assertEqual([layer.name for layer in image.layers], ["Apple", "Banana", "Zebra"])
    # The selection tracked "Banana" to its new index instead of staying at 0.
    t.assertEqual(iuser.multilayer_layer_name, "Banana")
    t.assertEqual(iuser.multilayer_layer, 1)


def test_render_result_catalog_copy():
    """A render-result viewer's layer/pass catalog is copied onto Image.layers,
    just like a loaded multi-layer EXR. The copy is refreshed
    centrally when the render finishes, with no image editor involved."""
    import bpy

    e, t, window = ui.test_window()

    # Enable extra passes so the catalog has more than just "Combined".
    view_layer = window.view_layer
    view_layer.use_pass_z = True
    view_layer.use_pass_mist = True

    scene = window.scene
    scene.render.resolution_x = 16
    scene.render.resolution_y = 16
    bpy.ops.render.render()
    yield

    # The catalog is populated by the central sync at render finish; no image
    # editor has been opened.
    render_result = bpy.data.images['Render Result']
    t.assertEqual([layer.name for layer in render_result.layers], ["ViewLayer"])
    pass_names = {p.name for p in render_result.layers[0].passes}
    t.assertIn("Combined", pass_names)
    t.assertIn("Depth", pass_names)
    t.assertIn("Mist", pass_names)


def test_render_result_combined_layer():
    """With a compositor, the render-result catalog gains a synthetic "Composite"
    layer ahead of the real render layers."""
    import bpy

    e, t, window = ui.test_window()

    scene = window.scene
    scene.render.resolution_x = 16
    scene.render.resolution_y = 16
    scene.render.compositor_device = 'CPU'

    # A minimal compositor node group: Render Layers -> Group Output.
    node_group = bpy.data.node_groups.new("Compositor", 'CompositorNodeTree')
    scene.compositing_node_group = node_group
    render_layers = node_group.nodes.new('CompositorNodeRLayers')
    group_output = node_group.nodes.new('NodeGroupOutput')
    node_group.interface.new_socket('Image', in_out='OUTPUT', socket_type='NodeSocketColor')
    node_group.links.new(render_layers.outputs['Image'], group_output.inputs[0])

    bpy.ops.render.render()
    yield

    render_result = bpy.data.images['Render Result']

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield
    area.spaces.active.image = render_result
    yield
    yield

    # The compositor result is exposed as a synthetic first layer.
    t.assertEqual([layer.name for layer in render_result.layers], ["Composite", "ViewLayer"])


def test_udim_gpu_textures():
    """A UDIM tiled image builds its aggregate GPU textures (atlas + tile mapping)
    when displayed, and rebuilds them without errors as the tile set changes."""
    import bpy

    e, t, window = ui.test_window()

    # A generated tiled image with three tiles.
    image = bpy.data.images.new("UDIM", 32, 32, tiled=True)
    for number in (1002, 1011):
        image.tiles.new(tile_number=number)
    t.assertEqual([tile.number for tile in image.tiles], [1001, 1002, 1011])

    # Show it in an image editor; drawing builds the UDIM array + tile-mapping textures.
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    area.spaces.active.image = image
    yield
    yield
    t.assertTrue(image.has_data)

    # Removing a tile invalidates the aggregate textures; the next redraws rebuild them.
    image.tiles.remove(image.tiles[2])
    t.assertEqual([tile.number for tile in image.tiles], [1001, 1002])
    yield
    yield
    t.assertTrue(image.has_data)

    # Adding a tile back likewise triggers a rebuild.
    image.tiles.new(tile_number=1005)
    yield
    yield
    t.assertEqual([tile.number for tile in image.tiles], [1001, 1002, 1005])
    t.assertTrue(image.has_data)


def test_image_cache_gc():
    """An idle image cache buffer is evicted by the time-based garbage collector."""
    import bpy
    import time

    e, t, window = ui.test_window()
    # The image cache garbage collector runs from the 3D viewport redraw.
    area = window.screen.areas[0]
    area.type = 'VIEW_3D'
    yield

    # Make the collector aggressive: collect every second, evict after one second idle.
    system_prefs = bpy.context.preferences.system
    system_prefs.texture_time_out = 1
    system_prefs.texture_collection_rate = 1

    # A generated image whose pixel buffer is cached on first access.
    image = bpy.data.images.new("GCTest", 64, 64)
    _ = image.pixels[0]
    t.assertTrue(image.has_data)

    # Redraw the viewport until the now-idle buffer is collected (a few seconds).
    deadline = time.time() + 20.0
    while image.has_data and time.time() < deadline:
        area.tag_redraw()
        yield

    t.assertFalse(image.has_data, "idle image cache buffer should be evicted by the GC")


def _author_image_sequence(directory, count):
    """Write `count` numbered 8x8 PNG frames; return the path of the first frame."""
    from OpenImageIO import ImageSpec, ImageOutput
    import numpy as np
    import os

    first_path = ""
    for i in range(1, count + 1):
        path = os.path.join(directory, "seq_{:04d}.png".format(i))
        spec = ImageSpec(8, 8, 3, "uint8")
        out = ImageOutput.create(path)
        out.open(path, spec)
        out.write_image(np.full((8, 8, 3), i * 40, dtype="uint8"))
        out.close()
        if i == 1:
            first_path = path
    return first_path


def test_image_sequence_gpu_textures():
    """Stepping through an image sequence in the image editor uploads each frame's
    GPU texture and drops the inactive frames' textures without errors."""
    import bpy
    import os
    import tempfile

    e, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    first_path = _author_image_sequence(tmp_dir, 5)

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    space = area.spaces.active
    image = bpy.data.images.load(first_path)
    image.source = 'SEQUENCE'
    space.image = image

    image_user = space.image_user
    image_user.frame_duration = 5
    image_user.frame_start = 1
    image_user.use_auto_refresh = True

    # Play through the sequence; each redraw uploads the current frame's texture.
    for frame in range(1, 6):
        window.scene.frame_set(frame)
        yield
        yield
    # ...and scrub back, which re-uploads earlier frames.
    for frame in (3, 1):
        window.scene.frame_set(frame)
        yield
        yield

    t.assertTrue(image.has_data)


def _read_exr_channels(path):
    """Return {qualified_channel_name: first_pixel_value} across all subimages of
    an EXR. Channel names are already layer/pass-qualified; multi-view EXRs repeat
    the same channel names across one subimage per view, so the per-subimage
    ``view`` attribute is prepended to keep views distinct. This normalizes over
    single-part vs multi-part layout, so an authored source and a Blender-saved
    copy compare equal."""
    from OpenImageIO import ImageInput
    inp = ImageInput.open(path)
    assert inp is not None, "cannot open " + path
    values = {}
    subimage = 0
    while True:
        if not inp.seek_subimage(subimage, 0):
            break
        spec = inp.spec()
        view = spec.get_string_attribute("view")
        pixels = inp.read_image()
        for ci, name in enumerate(spec.channelnames):
            key = f"{view}/{name}" if view else name
            values[key] = float(pixels[..., ci].reshape(-1)[0])
        subimage += 1
    inp.close()
    return values


def test_multilayer_exr_save_roundtrip():
    """A loaded multi-layer EXR saved back to a multi-layer EXR through the image
    save path preserves every layer, pass, channel and pixel value."""
    import bpy
    import os
    import tempfile
    from OpenImageIO import ImageSpec, ImageOutput
    import numpy as np

    e, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    src = os.path.join(tmp_dir, "src.exr")
    dst = os.path.join(tmp_dir, "dst.exr")

    # Author a multi-layer EXR: two layers, each Combined(RGBA) + Depth(Z), with a
    # distinct constant per channel.
    chan_specs = []
    value = 0.0
    for layer_name in ("Extra", "ViewLayer"):
        for channel in ("R", "G", "B", "A"):
            value += 0.03
            chan_specs.append((f"{layer_name}.Combined.{channel}", value))
        value += 0.03
        chan_specs.append((f"{layer_name}.Depth.Z", value))

    names = [c[0] for c in chan_specs]
    width = height = 4
    spec = ImageSpec(width, height, len(names), "float")
    spec.channelnames = tuple(names)
    pixels = np.zeros((height, width, len(names)), dtype="float32")
    for i, (_, v) in enumerate(chan_specs):
        pixels[:, :, i] = v
    out = ImageOutput.create(src)
    out.open(src, spec)
    out.write_image(pixels)
    out.close()

    # Show the image in an image editor so its ImageUser resolves and data loads.
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield
    image = bpy.data.images.load(src)
    area.spaces.active.image = image
    yield
    yield

    t.assertEqual([layer.name for layer in image.layers], ["Extra", "ViewLayer"])

    # Save it back as a multi-layer EXR through the image save path.
    scene = window.scene
    scene.render.image_settings.media_type = 'MULTI_LAYER_IMAGE'
    scene.render.image_settings.file_format = 'OPEN_EXR_MULTILAYER'
    scene.render.image_settings.color_depth = '32'
    image.save_render(dst, scene=scene)
    yield

    before = _read_exr_channels(src)
    after = _read_exr_channels(dst)
    t.assertEqual(set(before.keys()), set(after.keys()))
    for name in before:
        t.assertAlmostEqual(before[name], after[name], places=3, msg=name)


def test_multiview_multilayer_exr_save_roundtrip():
    """A multi-view multi-layer EXR saved through the image save path preserves
    both views' subimages, channels and pixel values."""
    import bpy
    import os
    import tempfile

    e, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    rendered = os.path.join(tmp_dir, "render.exr")
    dst = os.path.join(tmp_dir, "dst.exr")

    scene = window.scene
    scene.render.resolution_x = 16
    scene.render.resolution_y = 16
    scene.render.use_multiview = True
    scene.render.views_format = 'STEREO_3D'
    scene.view_layers[0].use_pass_z = True

    scene.render.image_settings.media_type = 'MULTI_LAYER_IMAGE'
    scene.render.image_settings.file_format = 'OPEN_EXR_MULTILAYER'
    scene.render.image_settings.color_depth = '32'
    scene.render.image_settings.views_format = 'MULTIVIEW'
    scene.render.filepath = rendered

    bpy.ops.render.render(write_still=True)
    yield

    # Load the rendered multi-view EXR and display it so it resolves.
    image = bpy.data.images.load(rendered)
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    area.spaces.active.image = image
    yield
    yield

    image.save_render(dst, scene=scene)
    yield

    before = _read_exr_channels(rendered)
    after = _read_exr_channels(dst)
    # The saved file must contain both left and right view subimages.
    t.assertTrue([n for n in after if "left" in n.lower()], sorted(after))
    t.assertTrue([n for n in after if "right" in n.lower()], sorted(after))
    # The same per-view channels round-tripped, with matching pixel values.
    t.assertEqual(set(before.keys()), set(after.keys()))
    for name in before:
        t.assertAlmostEqual(before[name], after[name], places=3, msg=name)


def _paint_layer_image(window, t):
    """Build a material with a Paint texture layer, yielding for the UI to catch
    up. Returns its generated multi-layer image and the Image Texture node
    reading it, for use with ``yield from``."""
    import bpy

    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("CatalogMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    bpy.context.active_object.data.materials.append(mat)
    yield

    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_paint(name="Paint", width=8, height=8)
    yield

    image_node = next(n for n in mat.node_tree.nodes if n.bl_idname == 'ShaderNodeTexImage')
    image = image_node.image
    t.assertEqual(image.type, 'MULTILAYER')
    return image, image_node


def test_image_catalog_operators():
    """The layer and pass operators edit the catalog of a generated multi-layer
    image, and the Image Texture node's per-pass outputs follow."""
    import bpy

    _, t, window = ui.test_window()

    image, image_node = yield from _paint_layer_image(window, t)

    # Edit through an image editor, which provides the image and its user.
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    area.spaces.active.image = image
    yield
    iuser = area.spaces.active.image_user

    t.assertEqual(len(image.layers), 1)
    initial_passes = [p.name for p in image.layers[0].passes]

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_add(name="Bump", type='VALUE')
    yield

    t.assertEqual([p.name for p in image.layers[0].passes], initial_passes + ["Bump"])
    added_pass = image.layers[0].passes["Bump"]
    t.assertEqual(added_pass.channels, 1)
    t.assertEqual(added_pass.channel_ids, "X")
    # The new pass is selected, and the node gained a matching output socket.
    t.assertEqual(iuser.multilayer_pass_name, "Bump")
    t.assertIn("Bump", [s.name for s in image_node.outputs])

    # A second pass with the same name is made unique instead of colliding.
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_add(name="Bump", type='COLOR')
    yield
    t.assertEqual([p.name for p in image.layers[0].passes], initial_passes + ["Bump", "Bump.001"])

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_remove()
        bpy.ops.image.pass_remove()
    yield
    t.assertEqual([p.name for p in image.layers[0].passes], initial_passes)
    t.assertNotIn("Bump", [s.name for s in image_node.outputs])

    # A new layer starts with a single Combined pass and becomes the selection
    # of the image editor that added it. The node reads its own layer, so its
    # per-pass outputs are unchanged.
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.layer_add(name="Second")
    yield
    t.assertEqual([layer.name for layer in image.layers], ["", "Second"])
    t.assertEqual([p.name for p in image.layers[1].passes], ["Combined"])
    t.assertEqual(iuser.multilayer_layer_name, "Second")
    t.assertEqual([s.name for s in image_node.outputs], initial_passes + ["Passes"])

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.layer_remove()
    yield
    t.assertEqual([layer.name for layer in image.layers], [""])
    # The last remaining layer and pass cannot be removed.
    with bpy.context.temp_override(window=window, area=area):
        t.assertFalse(bpy.ops.image.layer_remove.poll())


def test_image_catalog_rename():
    """Renaming a pass retargets the selection of every image user reading it,
    and renames the Image Texture node's output socket for that pass."""
    import bpy

    _, t, window = ui.test_window()

    image, image_node = yield from _paint_layer_image(window, t)

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    area.spaces.active.image = image
    yield
    iuser = area.spaces.active.image_user

    image_pass = image.layers[0].passes["Roughness"]
    # Select it in both the image editor and the node, by name.
    iuser.multilayer_pass = [p.name for p in image.layers[0].passes].index("Roughness")
    yield
    t.assertEqual(iuser.multilayer_pass_name, "Roughness")

    image_pass.name = "Gloss"
    yield

    t.assertEqual(image.layers[0].passes[image_pass.name].name, "Gloss")
    t.assertIn("Gloss", [s.name for s in image_node.outputs])
    t.assertNotIn("Roughness", [s.name for s in image_node.outputs])
    # The image editor followed the rename instead of falling back to the first
    # pass.
    t.assertEqual(iuser.multilayer_pass_name, "Gloss")

    # Layer names are made unique, and the layer selection follows too.
    image.layers[0].name = "Body"
    yield
    t.assertEqual(iuser.multilayer_layer_name, "Body")
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.layer_add(name="Body")
    yield
    t.assertEqual([layer.name for layer in image.layers], ["Body", "Body.001"])


def test_image_catalog_read_only():
    """A render result's catalog mirrors live render data, so it is read-only
    and the layer and pass operators do not poll on it."""
    import bpy

    _, t, window = ui.test_window()

    scene = window.scene
    scene.render.resolution_x = 16
    scene.render.resolution_y = 16
    scene.view_layers[0].use_pass_z = True
    bpy.ops.render.render()
    yield

    render_result = bpy.data.images["Render Result"]
    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    area.spaces.active.image = render_result
    yield
    yield

    t.assertEqual([layer.name for layer in render_result.layers], ["ViewLayer"])
    t.assertTrue(render_result.layers[0].is_property_readonly("name"))
    t.assertTrue(render_result.layers[0].passes[0].is_property_readonly("name"))

    with bpy.context.temp_override(window=window, area=area):
        t.assertFalse(bpy.ops.image.layer_add.poll())
        t.assertFalse(bpy.ops.image.pass_add.poll())
        # The panel still lists the catalog, read-only.
        t.assertTrue(bpy.types.IMAGE_PT_image_layers.poll(bpy.context))


def test_image_new_multilayer():
    """image.new with multilayer creates a generated image with an editable
    layer/pass catalog, its single pass filled with the dialog color."""
    import bpy
    import os
    import tempfile

    _, t, window = ui.test_window()

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.new(name="MultiLayer", width=4, height=4, multilayer=True,
                          color=(0.5, 0.25, 0.75, 1.0))
    yield
    yield

    image = bpy.data.images["MultiLayer"]
    t.assertEqual(image.source, 'GENERATED')
    t.assertEqual(image.type, 'MULTILAYER')
    # A single unnamed layer with one color pass, ready to be built up.
    t.assertEqual([layer.name for layer in image.layers], [""])
    t.assertEqual([p.name for p in image.layers[0].passes], ["Combined"])
    t.assertEqual(image.layers[0].passes[0].channels, 4)
    # The catalog is authored, so it is editable and the panel is shown.
    t.assertFalse(image.layers[0].is_property_readonly("name"))
    with bpy.context.temp_override(window=window, area=area):
        t.assertTrue(bpy.ops.image.pass_add.poll())
        t.assertTrue(bpy.types.IMAGE_PT_image_layers.poll(bpy.context))

    # The pass buffer carries the fill color, converted from the gamma-corrected
    # dialog color to the linear space of the float buffers.
    tmp_dir = tempfile.mkdtemp()
    dst = os.path.join(tmp_dir, "multilayer.exr")
    scene = window.scene
    scene.render.image_settings.media_type = 'MULTI_LAYER_IMAGE'
    scene.render.image_settings.file_format = 'OPEN_EXR_MULTILAYER'
    scene.render.image_settings.color_depth = '32'
    image.save_render(dst, scene=scene)
    yield

    from mathutils import Color
    expected = Color((0.5, 0.25, 0.75)).from_srgb_to_scene_linear()
    channels = _read_exr_channels(dst)
    for channel, value in zip(("R", "G", "B"), expected):
        name = next(n for n in channels if n.endswith("." + channel))
        t.assertAlmostEqual(channels[name], value, places=3, msg=name)


def test_image_file_catalog_editing():
    """The catalog of a loaded multi-layer EXR is editable in memory: edits make
    the image modified, keep the pixels of renamed passes, and are written back
    by saving over its own file."""
    import bpy
    import os
    import tempfile
    from OpenImageIO import ImageSpec, ImageOutput
    import numpy as np

    _, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    exr_path = os.path.join(tmp_dir, "editable.exr")

    # Two layers, each Combined(RGB), with a distinct constant per channel so a
    # renamed pass can be told apart from a blank one.
    names, values = [], []
    value = 0.0
    for layer_name in ("Extra", "ViewLayer"):
        for channel in ("R", "G", "B"):
            value += 0.05
            names.append("{}.Combined.{}".format(layer_name, channel))
            values.append(value)
    spec = ImageSpec(4, 4, len(names), "float")
    spec.channelnames = tuple(names)
    pixels = np.zeros((4, 4, len(names)), dtype="float32")
    for i, v in enumerate(values):
        pixels[:, :, i] = v
    out = ImageOutput.create(exr_path)
    out.open(exr_path, spec)
    out.write_image(pixels)
    out.close()

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield
    image = bpy.data.images.load(exr_path)
    area.spaces.active.image = image
    yield
    yield

    t.assertEqual([layer.name for layer in image.layers], ["Extra", "ViewLayer"])
    # A loaded EXR is editable, and starts out matching its file.
    t.assertFalse(image.layers[0].is_property_readonly("name"))
    t.assertFalse(image.is_dirty)
    with bpy.context.temp_override(window=window, area=area):
        t.assertTrue(bpy.ops.image.pass_add.poll())

    # Renaming a pass makes the image modified, since the edit is only in memory.
    image.layers[0].passes["Combined"].name = "Beauty"
    yield
    t.assertEqual([p.name for p in image.layers[0].passes], ["Beauty"])
    t.assertTrue(image.is_dirty)

    # ...and a pass added to the catalog has no data in the file, so it gets a
    # blank buffer of the image's size rather than failing to load.
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_add(name="Mask", type='VALUE', color=(0.0, 0.0, 0.0, 1.0))
    yield
    yield
    t.assertEqual([p.name for p in image.layers[0].passes], ["Beauty", "Mask"])
    t.assertEqual(tuple(image.size), (4, 4))

    # Saving over its own file writes the edited catalog and clears the modified
    # state; the catalog then reloads from the file unchanged.
    scene = window.scene
    scene.render.image_settings.media_type = 'MULTI_LAYER_IMAGE'
    scene.render.image_settings.file_format = 'OPEN_EXR_MULTILAYER'
    scene.render.image_settings.color_depth = '32'
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.save()
    yield
    yield

    t.assertFalse(image.is_dirty)
    t.assertEqual([layer.name for layer in image.layers], ["Extra", "ViewLayer"])
    t.assertEqual([p.name for p in image.layers[0].passes], ["Beauty", "Mask"])

    channels = _read_exr_channels(exr_path)
    # The renamed pass kept the pixels it was loaded with.
    for channel, value in zip(("R", "G", "B"), values[:3]):
        t.assertAlmostEqual(channels["Extra.Beauty." + channel], value, places=3, msg=channel)
    # The added pass was written too, as the blank it was filled with.
    added = [n for n in channels if n.startswith("Extra.Mask.")]
    t.assertTrue(added, sorted(channels))
    for name in added:
        t.assertAlmostEqual(channels[name], 0.0, places=3, msg=name)
    # The untouched layer is unchanged.
    for channel, value in zip(("R", "G", "B"), values[3:]):
        t.assertAlmostEqual(channels["ViewLayer.Combined." + channel], value, places=3, msg=channel)


def test_multifile_image_save_and_load():
    """A multi-layer image saved to a <PASS> path writes one file per pass, named
    by the pass token, and reads them back as its passes."""
    import bpy
    import os
    import tempfile

    _, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.new(name="TextureSet", width=4, height=4, multilayer=True,
                          color=(0.0, 0.0, 0.0, 1.0))
    yield
    image = bpy.data.images["TextureSet"]
    area.spaces.active.image = image
    yield

    image.layers[0].passes[0].name = "BaseColor"
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_add(name="Roughness", type='VALUE', color=(0.5, 0.5, 0.5, 1.0))
        bpy.ops.image.pass_add(name="Normal", type='VECTOR', color=(0.25, 0.5, 1.0, 1.0))
    yield

    # A pass whose file is named differently from the pass itself.
    image.layers[0].passes["Roughness"].token = "rough"
    yield

    # Saving to a <PASS> path writes one file per pass and rebinds the image to
    # that path, keeping its authored catalog.
    image.filepath_raw = os.path.join(tmp_dir, "tex_<PASS>.png")
    image.save()
    yield
    yield

    t.assertEqual(sorted(os.listdir(tmp_dir)),
                  ["tex_BaseColor.png", "tex_Normal.png", "tex_rough.png"])
    t.assertEqual(image.source, 'FILE')
    t.assertEqual(image.type, 'MULTILAYER')
    t.assertTrue(image.filepath.endswith("tex_<PASS>.png"), image.filepath)
    t.assertTrue(image.has_layer_tokens)
    t.assertEqual([p.name for p in image.layers[0].passes],
                  ["BaseColor", "Roughness", "Normal"])

    # Reloading reads every pass back from its own file. Writing the result to a
    # multi-layer EXR is the way to inspect all of them at once.
    image.reload()
    yield
    yield

    dst = os.path.join(tempfile.mkdtemp(), "roundtrip.exr")
    scene = window.scene
    scene.render.image_settings.media_type = 'MULTI_LAYER_IMAGE'
    scene.render.image_settings.file_format = 'OPEN_EXR_MULTILAYER'
    scene.render.image_settings.color_depth = '32'
    image.save_render(dst, scene=scene)
    yield

    channels = _read_exr_channels(dst)
    # A color pass round-trips through the image colorspace...
    t.assertAlmostEqual(channels["BaseColor.R"], 0.0, places=2)
    # ...while value and vector passes hold data, and keep their exact values.
    t.assertAlmostEqual(channels["Roughness.X"], 0.5, places=2)
    for channel, value in zip(("X", "Y", "Z"), (0.25, 0.5, 1.0)):
        t.assertAlmostEqual(channels["Normal." + channel], value, places=2, msg=channel)


def test_multifile_image_detect_layers():
    """The detect operator builds the layer/pass catalog of a multi-file image
    from the files matching its tokenized path, and keeps what still matches."""
    import bpy
    import os
    import tempfile

    _, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()

    # Author the files a texturing pipeline would have produced.
    filenames = ["mat_body_basecolor.png", "mat_body_rough.png", "mat_head_basecolor.png"]
    for filename in filenames:
        source = bpy.data.images.new("src", width=4, height=4)
        source.filepath_raw = os.path.join(tmp_dir, filename)
        source.file_format = 'PNG'
        source.save()
        bpy.data.images.remove(source)

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    image = bpy.data.images.load(os.path.join(tmp_dir, filenames[0]))
    image.filepath = os.path.join(tmp_dir, "mat_<LAYER>_<PASS>.png")
    area.spaces.active.image = image
    yield

    # The panel offers the operator before there is any catalog to show.
    with bpy.context.temp_override(window=window, area=area):
        t.assertTrue(bpy.types.IMAGE_PT_image_layers.poll(bpy.context))
        t.assertTrue(bpy.ops.image.layers_detect.poll())
        bpy.ops.image.layers_detect()
    yield
    yield

    t.assertEqual(image.type, 'MULTILAYER')
    t.assertEqual([layer.name for layer in image.layers], ["body", "head"])
    t.assertEqual([p.name for p in image.layers[0].passes], ["basecolor", "rough"])
    t.assertEqual([p.name for p in image.layers[1].passes], ["basecolor"])
    t.assertEqual(tuple(image.size), (4, 4))

    # Give a pass a name of its own, keeping the file name as its token. A second
    # detect must recognize it by that token and leave it alone, while dropping
    # the passes whose files are gone.
    rough = image.layers[0].passes["rough"]
    rough.token = "rough"
    rough.name = "Roughness"
    rough.type = 'VALUE'
    yield
    t.assertEqual(rough.channels, 1)
    t.assertEqual(rough.channel_ids, "X")

    os.remove(os.path.join(tmp_dir, "mat_head_basecolor.png"))
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.layers_detect()
    yield

    t.assertEqual([layer.name for layer in image.layers], ["body"])
    t.assertEqual([p.name for p in image.layers[0].passes], ["basecolor", "Roughness"])
    t.assertEqual(image.layers[0].passes["Roughness"].channels, 1)


def test_multifile_image_catalog_persists():
    """The catalog of a multi-file image is user data: it is written to the blend
    file and read back as-is, with no re-scan of the directory."""
    import bpy
    import os
    import tempfile

    _, t, window = ui.test_window()

    tmp_dir = tempfile.mkdtemp()
    blend_path = os.path.join(tmp_dir, "multifile.blend")

    area = window.screen.areas[0]
    area.type = 'IMAGE_EDITOR'
    yield

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.new(name="Persisted", width=4, height=4, multilayer=True,
                          color=(0.0, 0.0, 0.0, 1.0))
    yield
    image = bpy.data.images["Persisted"]
    area.spaces.active.image = image
    yield

    image.layers[0].passes[0].name = "BaseColor"
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.image.pass_add(name="Metallic", type='VALUE', color=(0.0, 0.0, 0.0, 1.0))
    yield
    image.layers[0].passes["Metallic"].token = "metal"
    image.filepath_raw = os.path.join(tmp_dir, "persisted_<PASS>.png")
    image.save()
    yield
    yield

    t.assertEqual(sorted(os.listdir(tmp_dir)),
                  ["persisted_BaseColor.png", "persisted_metal.png"])

    bpy.ops.wm.save_as_mainfile(filepath=blend_path)
    yield
    bpy.ops.wm.open_mainfile(filepath=blend_path)
    yield
    _, t, window = ui.test_window()
    yield

    image = bpy.data.images["Persisted"]
    t.assertEqual(image.type, 'MULTILAYER')
    t.assertEqual([p.name for p in image.layers[0].passes], ["BaseColor", "Metallic"])
    metallic = image.layers[0].passes["Metallic"]
    t.assertEqual(metallic.token, "metal")
    t.assertEqual(metallic.type, 'VALUE')
