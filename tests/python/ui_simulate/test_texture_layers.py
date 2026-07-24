# SPDX-FileCopyrightText: 2026 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Verify the Texture Layers experimental Surface panel and shader tree view."""

import modules.ui_test_utils as ui


def _build_shader_tree():
    import bpy
    mat = bpy.data.materials.new("TestMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    mix = nt.nodes.new("ShaderNodeMixShader")
    bsdf1 = nt.nodes.new("ShaderNodeBsdfPrincipled")
    emit = nt.nodes.new("ShaderNodeEmission")
    nt.links.new(bsdf1.outputs["BSDF"], mix.inputs[1])
    nt.links.new(emit.outputs["Emission"], mix.inputs[2])
    nt.links.new(mix.outputs["Shader"], out.inputs["Surface"])

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    return mat, mix, bsdf1, emit


def test_texture_layers_panel():
    """Toggle the experimental flag and verify the panel swap with legacy panels."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    _build_shader_tree()

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    surface_cls = bpy.types.MATERIAL_PT_texture_layers_surface
    legacy_surface_cls = bpy.types.EEVEE_MATERIAL_PT_surface
    legacy_volume_cls = bpy.types.EEVEE_MATERIAL_PT_volume
    legacy_displacement_cls = bpy.types.EEVEE_MATERIAL_PT_displacement

    with bpy.context.temp_override(window=window, area=properties_area):
        t.assertTrue(surface_cls.poll(bpy.context))
        t.assertFalse(legacy_surface_cls.poll(bpy.context))
        t.assertFalse(legacy_volume_cls.poll(bpy.context))
        t.assertFalse(legacy_displacement_cls.poll(bpy.context))

    bpy.context.preferences.experimental.use_texture_layers = False
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        t.assertFalse(surface_cls.poll(bpy.context))
        t.assertTrue(legacy_surface_cls.poll(bpy.context))
        t.assertTrue(legacy_volume_cls.poll(bpy.context))
        t.assertTrue(legacy_displacement_cls.poll(bpy.context))


def _build_shader_tree_with_layers():
    """Build a shader tree with a TextureLayerStack feeding a Principled BSDF
    via a Separate Bundle. The stack contains three layers (top, middle, base);
    only the base layer has a non-blended bundle input."""
    import bpy
    mat = bpy.data.materials.new("LayeredMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    stack.layer_items.new(name="Top")
    stack.layer_items.new(name="Middle")
    stack.layer_items.new(name="Base")

    sep = nt.nodes.new("NodeSeparateBundle")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    return mat, bsdf, stack, sep


def test_texture_layers_stack_node():
    """The Texture Layer Stack node grows its sockets per item, and the last
    item (base) is the only one without Opacity / Mask sockets."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    yield

    t.assertEqual(len(stack.layer_items), 3)
    t.assertEqual([item.name for item in stack.layer_items], ["Top", "Middle", "Base"])

    # The two non-base layers each get Layer + Opacity + Mask sockets; the base
    # layer gets a Layer + Mask socket (no Opacity, but a Mask so it can keep a
    # mask when it becomes the base). The trailing hollow extend socket is what
    # users drag onto in the node editor to add another layer.
    input_socket_names = [s.name for s in stack.inputs]
    t.assertEqual(input_socket_names, ["Top", "Opacity", "Mask",
                                       "Middle", "Opacity", "Mask",
                                       "Base", "Mask",
                                       ""])  # __extend__
    t.assertEqual(stack.inputs[-1].identifier, "__extend__")

    # Selecting the stack node and changing active_index drives which layer
    # the Material panel surfaces.
    mat.node_tree.nodes.active = stack
    stack.active_index = 1
    yield
    t.assertEqual(stack.active_index, 1)

    # Removing an item via RNA shrinks the socket list accordingly.
    stack.layer_items.remove(stack.layer_items[2])
    yield
    t.assertEqual(len(stack.layer_items), 2)
    input_socket_names = [s.name for s in stack.inputs]
    t.assertEqual(input_socket_names, ["Top", "Opacity", "Mask", "Middle", "Mask", ""])

    # Mute flag round-trips through RNA and preserves the socket layout.
    stack.layer_items[0].mute = True
    yield
    t.assertTrue(stack.layer_items[0].mute)


def test_texture_layers_panel_operators():
    """Material panel operators add/remove layers on the active stack."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()

    # Activate the Stack node directly so the panel operators resolve to it
    # (the function also walks BSDF -> Separate Bundle -> Stack, but there is
    # no link from the Separate Bundle to the BSDF in this minimal setup).
    mat.node_tree.nodes.active = stack
    yield

    # Object must be the material's owner so context.material resolves it.
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    initial_count = len(stack.layer_items)

    # Select the middle layer, so we can verify the new layer is inserted
    # directly above it (not appended at the bottom / base position).
    stack.active_index = 1
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add(name="Inserted")
    yield
    t.assertEqual(len(stack.layer_items), initial_count + 1)
    # The new layer takes the active slot, pushing the previously active layer
    # (and everything below it) down; active_index points to the new layer.
    t.assertEqual([it.name for it in stack.layer_items], ["Top", "Inserted", "Middle", "Base"])
    t.assertEqual(stack.active_index, 1)
    t.assertEqual(stack.layer_items[stack.active_index].name, "Inserted")

    # Reorder via the reparent operator (drag-and-drop replacement for up/down).
    src_idx = stack.active_index
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack.name, from_index=src_idx,
            to_stack=stack.name, to_index=src_idx + 1,
        )
    yield
    t.assertEqual(stack.layer_items[stack.active_index].name, "Inserted")

    # Remove the active layer.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_remove()
    yield
    t.assertEqual(len(stack.layer_items), initial_count)


def test_texture_layers_add_fill():
    """material.texture_layer_add_fill inserts a layer wired to a Combine Bundle."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    yield

    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    initial_count = len(stack.layer_items)
    initial_combines = sum(1 for n in mat.node_tree.nodes if n.bl_idname == 'NodeCombineBundle')

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_fill()
    yield

    t.assertEqual(len(stack.layer_items), initial_count + 1)
    combines = [n for n in mat.node_tree.nodes if n.bl_idname == 'NodeCombineBundle']
    t.assertEqual(len(combines), initial_combines + 1)
    fill = combines[-1]
    # Default Fill channels are exposed as Combine Bundle items.
    item_names = [it.name for it in fill.bundle_items]
    t.assertIn("Base Color", item_names)
    t.assertIn("Roughness", item_names)
    # The new layer's Bundle input is connected to the Combine Bundle's output.
    new_item = stack.layer_items[stack.active_index]
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier))
    t.assertTrue(layer_socket.is_linked)
    t.assertEqual(layer_socket.links[0].from_node, fill)


def test_texture_layers_add_paint():
    """material.texture_layer_add_paint inserts a layer wired to an Image
    Texture node with a new generated multi-layer image, whose passes are the
    layer's channels."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    yield

    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    initial_count = len(stack.layer_items)

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_paint()
    yield

    t.assertEqual(len(stack.layer_items), initial_count + 1)
    image_nodes = [n for n in mat.node_tree.nodes if n.bl_idname == 'ShaderNodeTexImage']
    t.assertEqual(len(image_nodes), 1)
    image_node = image_nodes[0]
    image = image_node.image
    t.assertIsNotNone(image)
    t.assertEqual(image.source, 'GENERATED')
    t.assertEqual(image.type, 'MULTILAYER')
    # A single unnamed image layer whose passes are the default channels.
    t.assertEqual(len(image.layers), 1)
    pass_names = [p.name for p in image.layers[0].passes]
    t.assertIn("Base Color", pass_names)
    t.assertIn("Roughness", pass_names)
    # The node exposes per-pass outputs plus the Passes bundle output, which
    # feeds the new layer's Bundle input.
    output_names = [s.name for s in image_node.outputs]
    t.assertIn("Passes", output_names)
    for name in pass_names:
        t.assertIn(name, output_names)
    new_item = stack.layer_items[stack.active_index]
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier))
    t.assertTrue(layer_socket.is_linked)
    t.assertEqual(layer_socket.links[0].from_node, image_node)
    t.assertEqual(layer_socket.links[0].from_socket.name, "Passes")


def test_texture_layers_select_activates_paint_image():
    """Selecting an image-backed layer row in the tree view makes its Image
    Texture node the active paint canvas, so the texture paint slot follows
    the selected layer."""
    import bpy

    e, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    # A material whose stack is fully wired into the BSDF (bootstrapped by the
    # operator), so the tree view lists the layer rows.
    mat = bpy.data.materials.new("PaintSelMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    bpy.context.active_object.data.materials.append(mat)

    # Use the largest area for the tree view, so the click targets are big.
    properties_area = ui.largest_area(window.screen)
    properties_area.type = 'PROPERTIES'
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_paint(name="PaintA", width=8, height=8)
        bpy.ops.material.texture_layer_add_paint(name="PaintB", width=8, height=8)
    yield
    stack = next(n for n in mat.node_tree.nodes
                 if n.bl_idname == 'ShaderNodeTextureLayerStack')

    image_a = bpy.data.images["PaintA"]
    image_b = bpy.data.images["PaintB"]

    # Texture paint mode builds the paint slots from the active paint canvas.
    bpy.ops.object.mode_set(mode='TEXTURE_PAINT')
    yield
    # Dismiss the "Failed to load using Vulkan" fallback dialog that can pop up
    # on headless configurations and would swallow the clicks below.
    e.ret.tap()
    yield

    def active_paint_image():
        m = bpy.data.materials["PaintSelMat"]
        if not m.texture_paint_images:
            return None
        return m.texture_paint_images[m.paint_active_slot]

    # The most recently added Paint layer is the paint target.
    t.assertEqual(active_paint_image(), image_b)

    # Click down the tree view rows (BSDF row, then PaintB / PaintA / Fill);
    # whenever a Paint layer row becomes the active layer, the paint slot must
    # have followed its image. The scan band covers the rows below the panels
    # above the tree (material slots, ID template, Preview and Surface panel
    # headers), without reaching any other widgets.
    layer_images = {"PaintA": image_a, "PaintB": image_b}
    seen = set()
    x = properties_area.x + 200
    area_top = properties_area.y + properties_area.height
    for i, y in enumerate(range(area_top - 245, area_top - 350, -8)):
        # Alternate the click x position so consecutive clicks are never within
        # the double-click distance (a double-click would start renaming).
        e.cursor_position_set(x + (60 if i % 2 else 0), y, move=True)
        yield
        e.leftmouse.tap()
        yield
        if stack.active_index < 0 or stack.active_index >= len(stack.layer_items):
            continue
        name = stack.layer_items[stack.active_index].name
        if name in layer_images:
            seen.add(name)
            t.assertEqual(active_paint_image(), layer_images[name])
    # The scan must have hit both paint layers for the test to mean anything.
    t.assertEqual(seen, {"PaintA", "PaintB"})


def test_texture_layers_paint_undo_keeps_pixels():
    """Painted pass buffers survive a memfile undo that re-reads the Image ID:
    the preserved cache entries are re-keyed onto the re-read pass structs."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    yield

    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_paint(width=8, height=8)
    yield

    image = next(n.image for n in mat.node_tree.nodes if n.bl_idname == 'ShaderNodeTexImage')

    # Paint the first pass (Base Color); this only touches the cached buffer.
    image.pixels = [0.25, 0.5, 0.75, 1.0] * 64
    yield

    # A DNA change to the image plus undo forces the Image ID to be re-read
    # from the memfile, while the buffer cache is preserved across the step.
    with bpy.context.temp_override(window=window):
        bpy.ops.ed.undo_push(message="Before rename")
    image.name = "Renamed"
    with bpy.context.temp_override(window=window):
        bpy.ops.ed.undo_push(message="Rename image")
        bpy.ops.ed.undo()
    yield

    image = bpy.data.images.get("Paint")
    t.assertIsNotNone(image)
    for value, expected in zip(image.pixels[0:4], (0.25, 0.5, 0.75, 1.0)):
        t.assertAlmostEqual(value, expected, places=3)


def test_texture_layers_groups():
    """Adding a Group creates a nested Texture Layer Stack, and the reparent
    operator can move layers in and out of the group (= drag/drop semantics)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    stack_name = stack.name
    initial_count = len(stack.layer_items)
    initial_stacks = sum(
        1 for n in mat.node_tree.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
    )

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_group()
    yield

    # Refresh handles; bpy wrappers can be invalidated when the node array
    # is reallocated by nodes.new().
    stack = mat.node_tree.nodes[stack_name]
    nt = mat.node_tree
    stacks = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(len(stacks), initial_stacks + 1)
    t.assertEqual(len(stack.layer_items), initial_count + 1)
    nested = next(n for n in stacks if n.name != stack_name)
    nested_name = nested.name

    # The new layer's Bundle input is fed by the nested stack.
    group_item = stack.layer_items[stack.active_index]
    group_id = group_item.identifier
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(group_id))
    t.assertTrue(layer_socket.is_linked)
    t.assertEqual(layer_socket.links[0].from_node.name, nested_name)

    # Drag a sibling layer INTO the group: source is some other layer in the
    # outer stack; drop into the nested stack at index 0 (top).
    other_idx = next(
        i for i, it in enumerate(stack.layer_items) if it.identifier != group_id
    )
    other_name = stack.layer_items[other_idx].name
    nested_count_before = len(nested.layer_items)
    outer_count_before = len(stack.layer_items)
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_name, from_index=other_idx,
            to_stack=nested_name, to_index=0,
        )
    yield
    stack = nt.nodes[stack_name]
    nested = nt.nodes[nested_name]
    t.assertEqual(len(nested.layer_items), nested_count_before + 1)
    t.assertEqual(len(stack.layer_items), outer_count_before - 1)
    t.assertEqual(nested.layer_items[0].name, other_name)

    # Drag the moved layer back OUT of the group, placing it at the top.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=nested_name, from_index=0,
            to_stack=stack_name, to_index=0,
        )
    yield
    stack = nt.nodes[stack_name]
    nested = nt.nodes[nested_name]
    t.assertEqual(len(nested.layer_items), nested_count_before)
    t.assertEqual(len(stack.layer_items), outer_count_before)
    t.assertEqual(stack.layer_items[0].name, other_name)

    # Cycle prevention: a group can't be dropped into itself.
    group_idx = next(
        i for i, it in enumerate(stack.layer_items) if it.identifier == group_id
    )
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_name, from_index=group_idx,
            to_stack=nested_name, to_index=0,
        )
    yield
    stack = nt.nodes[stack_name]
    # The cycle move was refused — outer count stays the same.
    t.assertEqual(len(stack.layer_items), outer_count_before)


def _make_tint_adjustment_group():
    """Create a minimal local TEXTURE_ADJUSTMENT-tagged shader node group with
    Bundle in/out: it just passes the input bundle straight to the output. The
    wiring details don't matter — the operator only cares about the first
    Bundle input/output sockets."""
    import bpy
    g = bpy.data.node_groups.new("TestAdjust", "ShaderNodeTree")
    g.usage = {'TEXTURE_ADJUSTMENT'}
    g.interface.new_socket("Bundle", in_out='INPUT', socket_type='NodeSocketBundle')
    g.interface.new_socket("Bundle", in_out='OUTPUT', socket_type='NodeSocketBundle')
    in_node = g.nodes.new('NodeGroupInput')
    out_node = g.nodes.new('NodeGroupOutput')
    g.links.new(in_node.outputs["Bundle"], out_node.inputs["Bundle"])
    g.asset_mark()
    return g


def test_texture_layers_adjustment_wiring():
    """Legacy adjustment wiring, as found in files saved before adjustments
    became closure regions: the layer's Bundle input is fed by the group's
    output and the group's input by the source of the layer below. This
    wiring must remain constructible and structurally intact."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'

    # Wire a Combine Bundle into the bottom (base) layer so the layer below the
    # adjustment has a real bundle source for the adjustment to consume.
    nt = mat.node_tree
    base_item = stack.layer_items[-1]
    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    base_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(base_item.identifier)
    )
    nt.links.new(combine.outputs["Bundle"], base_socket)
    yield

    # Simulate the operator's wiring: insert a layer above the base, link the
    # group's Bundle output to it, and link the group's Bundle input to the
    # bundle source feeding the layer immediately below.
    group = _make_tint_adjustment_group()

    # Place the new layer just above the base layer.
    target_index = len(stack.layer_items) - 1
    new_item = stack.layer_items.new(name="Tint")
    stack.layer_items.move(len(stack.layer_items) - 1, target_index)
    stack.active_index = target_index

    group_node = nt.nodes.new('ShaderNodeGroup')
    group_node.node_tree = group

    new_layer_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier)
    )
    group_out = next(s for s in group_node.outputs if s.type == 'BUNDLE')
    nt.links.new(group_out, new_layer_socket)

    # The adjustment's Bundle input should be wired to the source feeding the
    # layer below — for our setup that's the Combine Bundle output.
    combine_name = combine.name
    group_in = next(s for s in group_node.inputs if s.type == 'BUNDLE')
    combine_out = nt.nodes[combine_name].outputs["Bundle"]
    nt.links.new(combine_out, group_in)
    yield

    # New layer is fed by the adjustment group output.
    new_layer_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier)
    )
    t.assertTrue(new_layer_socket.is_linked)
    t.assertEqual(new_layer_socket.links[0].from_node.name, group_node.name)

    # Adjustment's Bundle input is fed by the Combine Bundle of the base layer.
    group_in = next(s for s in group_node.inputs if s.type == 'BUNDLE')
    t.assertTrue(group_in.is_linked)
    t.assertEqual(group_in.links[0].from_node.name, combine_name)


def test_texture_layers_adjustment_operator():
    """Invoke material.texture_layer_add_node_group with an adjustment asset
    and verify the closure-region wiring: the new layer is CLOSURE-typed and
    its socket is fed by a Closure Output node; the asset's group node sits
    inside the zone, between the paired Closure Input/Output nodes."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'

    # Provide a real bundle source for the base (bottom) layer so the
    # adjustment has something to consume.
    nt = mat.node_tree
    base_item = stack.layer_items[-1]
    base_id = base_item.identifier
    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    combine_name = combine.name
    base_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(base_id)
    )
    nt.links.new(combine.outputs["Bundle"], base_socket)

    # Make the base layer the active one so the new layer is inserted just
    # above it (the operator places new layers above the active position).
    stack.active_index = len(stack.layer_items) - 1

    # Local asset adjustment group — looked up via the LOCAL asset library.
    group = _make_tint_adjustment_group()
    group_name = group.name
    yield

    initial_layers = len(stack.layer_items)
    # Local asset weak ref is "<idtype_name>/<id_name>" — e.g. "NodeTree/TestAdjust".
    import os
    relative_id = "NodeTree" + os.sep + group_name
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node_group(
            'EXEC_DEFAULT',
            asset_library_type='LOCAL',
            asset_library_identifier='',
            relative_asset_identifier=relative_id,
        )
    yield

    del combine_name  # The base layer's source is untouched by the operator.

    # A new CLOSURE-typed layer was inserted, fed by a Closure Output node.
    t.assertEqual(len(stack.layer_items), initial_layers + 1)
    new_item = stack.layer_items[stack.active_index]
    t.assertEqual(new_item.item_type, 'CLOSURE')
    new_layer_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier)
    )
    t.assertEqual(new_layer_socket.type, 'CLOSURE')
    t.assertTrue(new_layer_socket.is_linked)
    closure_out = new_layer_socket.links[0].from_node
    t.assertEqual(closure_out.bl_idname, 'NodeClosureOutput')

    # The zone declares one Bundle input and one Bundle output item.
    t.assertEqual([i.socket_type for i in closure_out.input_items], ['BUNDLE'])
    t.assertEqual([i.socket_type for i in closure_out.output_items], ['BUNDLE'])

    # The asset's group node sits inside the zone: its Bundle output feeds the
    # Closure Output, its Bundle input is fed by the paired Closure Input.
    zone_bundle_in = next(s for s in closure_out.inputs if s.type == 'BUNDLE')
    t.assertTrue(zone_bundle_in.is_linked)
    group_node = zone_bundle_in.links[0].from_node
    t.assertEqual(group_node.bl_idname, 'ShaderNodeGroup')
    t.assertEqual(group_node.node_tree.name, group_name)
    group_in = next(s for s in group_node.inputs if s.type == 'BUNDLE')
    t.assertTrue(group_in.is_linked)
    closure_in = group_in.links[0].from_node
    t.assertEqual(closure_in.bl_idname, 'NodeClosureInput')
    t.assertEqual(closure_in.paired_output, closure_out)


def test_texture_layers_adjustment_inline_roundtrip():
    """Legacy bundle-wired adjustments (files saved before adjustments became
    closure regions) degrade gracefully: the group acts on its wired input
    and the stack still inlines into a well-formed tree."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("InlineMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    stack.layer_items.new(name="Adjusted")
    stack.layer_items.new(name="Base")

    sep = nt.nodes.new("NodeSeparateBundle")
    sep.bundle_items.new('RGBA', "Color")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])
    nt.links.new(sep.outputs["Color"], bsdf.inputs["Base Color"])

    # Base layer: Combine Bundle with red Color.
    base_id = stack.layer_items[-1].identifier
    base_combine = nt.nodes.new('NodeCombineBundle')
    base_combine.bundle_items.new('RGBA', "Color")
    base_combine.inputs["Color"].default_value = (1.0, 0.0, 0.0, 1.0)
    base_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(base_id))
    nt.links.new(base_combine.outputs["Bundle"], base_socket)

    # Adjustment layer: a passthrough TEXTURE_ADJUSTMENT group on top.
    adj_id = stack.layer_items[0].identifier
    adj_group = _make_tint_adjustment_group()
    adj_node = nt.nodes.new('ShaderNodeGroup')
    adj_node.node_tree = adj_group
    adj_layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(adj_id))
    adj_out = next(s for s in adj_node.outputs if s.type == 'BUNDLE')
    nt.links.new(adj_out, adj_layer_socket)
    adj_in = next(s for s in adj_node.inputs if s.type == 'BUNDLE')
    nt.links.new(base_combine.outputs["Bundle"], adj_in)
    yield

    # The wrapper owns the inlined tree's lifetime — keep a reference.
    inline_wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
    inlined = inline_wrapper.node_tree

    # The inlined tree must not contain a Texture Layer Stack node — it should
    # have been fully expanded, with the adjustment's contents inlined too.
    stack_nodes = [n for n in inlined.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(stack_nodes, [])

    # The Material Output's Surface input is connected (basic well-formedness).
    out_inlined = next(n for n in inlined.nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
    surface = out_inlined.inputs["Surface"]
    t.assertTrue(surface.is_linked)
    # Keep the wrapper alive until the assertions are done.
    del inline_wrapper


def test_texture_layer_stack_extend_socket():
    """The Texture Layer Stack node has a hollow trailing extend socket: dragging
    a Bundle link onto it adds a new layer item at the end. The new item
    becomes the base, and the previously-last item gains Opacity/Mask
    sockets — matching the natural "drop into the bottom socket" mental
    model. Names default to "Layer", "Layer.001", ... rather than the source
    socket's name."""
    import bpy

    _, t, window = ui.test_window()

    mat = bpy.data.materials.new("ExtendMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    stack.layer_items.new(name="Layer")
    first_id = stack.layer_items[0].identifier

    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    yield

    # The trailing extend socket is the last input on the stack node.
    extend = next(s for s in stack.inputs if s.identifier == "__extend__")
    t.assertEqual(stack.inputs[-1], extend)
    # handle_dynamic_sockets=True is what triggers the node's insert_link
    # callback (matches what the node editor does on user drag-drop).
    nt.links.new(combine.outputs["Bundle"], extend, handle_dynamic_sockets=True)
    yield

    # A new layer was appended at the end; the previously-only item is now at
    # index 0 and the new item is the base.
    t.assertEqual(len(stack.layer_items), 2)
    t.assertEqual(stack.layer_items[0].identifier, first_id)
    t.assertEqual(stack.active_index, 1)
    new_item = stack.layer_items[1]
    t.assertEqual(new_item.name, "Layer.001")

    new_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier))
    t.assertTrue(new_socket.is_linked)
    t.assertEqual(new_socket.links[0].from_node, combine)

    # The previously-last item (now at items[0]) gains an Opacity socket; the
    # new bottom item (the new base) gets a Layer + Mask socket, but no Opacity.
    ids = [s.identifier for s in stack.inputs]
    t.assertIn("Opacity_" + str(first_id), ids)
    t.assertIn("Mask_" + str(first_id), ids)
    t.assertNotIn("Opacity_" + str(new_item.identifier), ids)
    t.assertIn("Mask_" + str(new_item.identifier), ids)

    # The extend socket is still there at the end, ready for more links.
    t.assertEqual(stack.inputs[-1].identifier, "__extend__")


def test_texture_layer_stack_extend_socket_closure():
    """Dropping a Closure link onto the trailing extend socket creates a
    CLOSURE-typed layer item with a Closure-typed socket. Items created any
    other way default to BUNDLE."""
    import bpy

    _, t, window = ui.test_window()

    mat = bpy.data.materials.new("ExtendClosureMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    stack.layer_items.new(name="Base")
    # Unlinked items are untyped until a link determines their type.
    t.assertEqual(stack.layer_items[0].item_type, 'UNSET')

    closure_in = nt.nodes.new("NodeClosureInput")
    closure_out = nt.nodes.new("NodeClosureOutput")
    t.assertTrue(closure_in.pair_with_output(closure_out))
    closure_out.input_items.new('BUNDLE', "Bundle")
    closure_out.output_items.new('BUNDLE', "Bundle")
    yield

    extend = next(s for s in stack.inputs if s.identifier == "__extend__")
    closure_socket = next(s for s in closure_out.outputs if s.type == 'CLOSURE')
    nt.links.new(closure_socket, extend, handle_dynamic_sockets=True)
    yield

    t.assertEqual(len(stack.layer_items), 2)
    new_item = stack.layer_items[1]
    t.assertEqual(new_item.item_type, 'CLOSURE')
    new_socket = next(
        s for s in stack.inputs if s.identifier == "Layer_" + str(new_item.identifier))
    t.assertEqual(new_socket.type, 'CLOSURE')
    t.assertTrue(new_socket.is_linked)
    t.assertEqual(new_socket.links[0].from_node, closure_out)

    # The extend socket survives for further links.
    t.assertEqual(stack.inputs[-1].identifier, "__extend__")


def test_texture_layer_retype_through_extend_socket():
    """Unlinked layer items settle to UNSET and expose a hollow extend
    socket; dropping a link types the item in place (Bundle or Closure), and
    unlinking resets it so it can be re-typed."""
    import bpy

    _, t, window = ui.test_window()

    mat = bpy.data.materials.new("RetypeMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    item = stack.layer_items.new(name="Layer")
    item_id = item.identifier
    yield

    # An unlinked item settles to UNSET with a hollow socket.
    t.assertEqual(stack.layer_items[0].item_type, 'UNSET')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    t.assertEqual(layer_socket.bl_idname, 'NodeSocketVirtual')

    # Dropping a Bundle link onto the hollow socket types the item in place.
    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    nt.links.new(combine.outputs["Bundle"], layer_socket, handle_dynamic_sockets=True)
    yield
    t.assertEqual(len(stack.layer_items), 1)
    t.assertEqual(stack.layer_items[0].item_type, 'BUNDLE')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    t.assertEqual(layer_socket.type, 'BUNDLE')
    t.assertTrue(layer_socket.is_linked)

    # Unlinking resets the item to UNSET.
    nt.links.remove(layer_socket.links[0])
    yield
    t.assertEqual(stack.layer_items[0].item_type, 'UNSET')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    t.assertEqual(layer_socket.bl_idname, 'NodeSocketVirtual')

    # Re-typing with a Closure link works on the same item.
    closure_in = nt.nodes.new("NodeClosureInput")
    closure_out = nt.nodes.new("NodeClosureOutput")
    t.assertTrue(closure_in.pair_with_output(closure_out))
    closure_out.input_items.new('BUNDLE', "Bundle")
    closure_out.output_items.new('BUNDLE', "Bundle")
    closure_socket = next(s for s in closure_out.outputs if s.type == 'CLOSURE')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    nt.links.new(closure_socket, layer_socket, handle_dynamic_sockets=True)
    yield
    t.assertEqual(len(stack.layer_items), 1)
    t.assertEqual(stack.layer_items[0].item_type, 'CLOSURE')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    t.assertEqual(layer_socket.type, 'CLOSURE')
    t.assertTrue(layer_socket.is_linked)

    # Plain links.new without dynamic-socket handling self-heals through the
    # node's update: the link lands on the hollow socket and the item infers
    # its type from it.
    nt.links.remove(layer_socket.links[0])
    yield
    t.assertEqual(stack.layer_items[0].item_type, 'UNSET')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    nt.links.new(combine.outputs["Bundle"], layer_socket)
    yield
    t.assertEqual(stack.layer_items[0].item_type, 'BUNDLE')
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
    t.assertEqual(layer_socket.type, 'BUNDLE')
    t.assertTrue(layer_socket.is_linked)


def test_texture_layer_convert_to_group():
    """Wrap a single Fill layer in a Texture Layer Group: a new nested
    Texture Layer Stack appears in the tree, the Fill source is moved into
    the nested stack, and the outer layer's Bundle is now fed by the nested
    stack's Result."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    # Wire the active layer to a real Combine Bundle so it has a source to nest.
    nt = mat.node_tree
    stack.active_index = 0
    item = stack.layer_items[0]
    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    combine_name = combine.name
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item.identifier))
    nt.links.new(combine.outputs["Bundle"], layer_socket)

    initial_stacks = sum(
        1 for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
    )
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_convert_to_group()
    yield

    # A new nested stack appeared.
    stacks = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(len(stacks), initial_stacks + 1)
    nested = next(n for n in stacks if n != stack)

    # The outer layer's Bundle is now fed by the nested stack's Result.
    layer_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item.identifier))
    t.assertTrue(layer_socket.is_linked)
    t.assertEqual(layer_socket.links[0].from_node, nested)
    t.assertEqual(layer_socket.links[0].from_socket.name, "Result")

    # The Combine Bundle that used to feed the outer layer now feeds the nested
    # stack's first item — preserving the user's original source.
    t.assertEqual(len(nested.layer_items), 1)
    nested_item = nested.layer_items[0]
    nested_socket = next(
        s for s in nested.inputs if s.identifier == "Layer_" + str(nested_item.identifier)
    )
    t.assertTrue(nested_socket.is_linked)
    t.assertEqual(nested_socket.links[0].from_node.name, combine_name)


def _fully_wired_material():
    """A material whose Principled BSDF is driven by a Texture Layer Stack
    (stack -> Separate Bundle -> BSDF), built via the add-default operator so
    the stack actually reaches the shader output and inlines."""
    import bpy
    mat = bpy.data.materials.new("AdjGrpMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_add_default()
    return mat, ob


def _add_tint_adjustment(mat, ob):
    """Insert a pass-through TEXTURE_ADJUSTMENT asset as a CLOSURE layer via the
    add-node-group operator; returns the new closure item."""
    import bpy
    import os
    stack = next(n for n in mat.node_tree.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    mat.node_tree.nodes.active = stack
    group = _make_tint_adjustment_group()
    relative_id = "NodeTree" + os.sep + group.name
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_add_node_group(
            'EXEC_DEFAULT',
            asset_library_type='LOCAL',
            asset_library_identifier='',
            relative_asset_identifier=relative_id,
        )
    return stack.layer_items[stack.active_index]


def _inline_leftovers(mat):
    import bpy
    bpy.context.view_layer.update()
    w = bpy.types.InlineShaderNodes.from_material(material=mat)
    leftovers = [
        n.bl_idname for n in w.node_tree.nodes
        if n.bl_idname in ('ShaderNodeTextureLayerStack', 'ShaderNodeMaskStack',
                           'NodeClosureInput', 'NodeClosureOutput')
    ]
    surface_linked = next(
        n for n in w.node_tree.nodes if n.bl_idname == 'ShaderNodeOutputMaterial'
    ).inputs["Surface"].is_linked
    del w
    return leftovers, surface_linked


def _closure_input_source(stack):
    """The nested stack's base item socket and whether it is fed by a Closure
    Input node (the adjustment group's hidden input)."""
    base_item = stack.layer_items[-1]
    sock = next(s for s in stack.inputs if s.identifier == "Layer_" + str(base_item.identifier))
    return sock.is_linked and sock.links[0].from_node.bl_idname == 'NodeClosureInput'


def test_texture_layer_add_adjustment_group():
    """Add an empty Adjustment Group: a new closure layer whose zone body is a
    nested stack with a single hidden base fed by the group's Closure Input.
    The material inlines to a well-formed tree (identity pass-through)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, ob = _fully_wired_material()
    nt = mat.node_tree
    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    stack.active_index = 0
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_add_adjustment_group()
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
                 and not _closure_input_source(n))
    nested = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
                  and _closure_input_source(n))
    # Nested has a single hidden base fed by Closure Input.
    t.assertEqual(len(nested.layer_items), 1)
    # The outer group layer is closure-typed, fed by a Closure Output.
    grp_item = next(it for it in stack.layer_items if it.item_type == 'CLOSURE')
    grp_socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(grp_item.identifier))
    t.assertTrue(grp_socket.is_linked)
    t.assertEqual(grp_socket.links[0].from_node.bl_idname, 'NodeClosureOutput')

    leftovers, surface_linked = _inline_leftovers(mat)
    t.assertEqual(leftovers, [])
    t.assertTrue(surface_linked)


def test_texture_layer_convert_adjustment_to_group():
    """Convert an adjustment (closure) layer to an adjustment group: the
    adjustment moves onto a nested stack whose hidden base is fed by the group's
    Closure Input, and the outer layer is re-fed from a Closure Output. The
    material still inlines cleanly."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, ob = _fully_wired_material()
    nt = mat.node_tree
    adj_item = _add_tint_adjustment(mat, ob)
    yield
    t.assertEqual(adj_item.item_type, 'CLOSURE')

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    stack.active_index = next(
        i for i, it in enumerate(stack.layer_items) if it.item_type == 'CLOSURE')
    nt.nodes.active = stack
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_convert_to_group()
    yield

    # A nested stack appeared, fed at its base by a Closure Input, with the
    # adjustment (closure) on top of the base.
    nested = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
                  and _closure_input_source(n))
    t.assertEqual(len(nested.layer_items), 2)
    t.assertEqual(nested.layer_items[0].item_type, 'CLOSURE')

    leftovers, surface_linked = _inline_leftovers(mat)
    t.assertEqual(leftovers, [])
    t.assertTrue(surface_linked)


def test_texture_layer_ungroup_adjustment_group():
    """Ungroup an adjustment group: its adjustment layers flatten back into the
    parent stack, the hidden input base is dropped, and the group's own closure
    zone and nested stack are removed. A single adjustment round-trips back to a
    flat closure layer."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, ob = _fully_wired_material()
    nt = mat.node_tree
    _add_tint_adjustment(mat, ob)
    yield

    def activate_closure_layer():
        stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
                     and not _closure_input_source(n))
        nt.nodes.active = stack
        stack.active_index = next(
            i for i, it in enumerate(stack.layer_items) if it.item_type == 'CLOSURE')

    # Convert the adjustment to a group, then ungroup it again.
    activate_closure_layer()
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_convert_to_group()
    yield
    activate_closure_layer()
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_ungroup()
    yield

    # Back to a single stack with a flat closure adjustment layer; the group's
    # nested stack is gone (no Closure-Input-fed stack remains).
    stacks = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(len(stacks), 1)
    t.assertFalse(_closure_input_source(stacks[0]))
    t.assertTrue(any(it.item_type == 'CLOSURE' for it in stacks[0].layer_items))

    leftovers, surface_linked = _inline_leftovers(mat)
    t.assertEqual(leftovers, [])
    t.assertTrue(surface_linked)


def test_texture_layer_collapse_state_persists():
    """A group layer's tree-view collapse state is stored on the item
    (show_expanded), so it survives undo (which reallocates node pointers) and
    save/load."""
    import bpy
    import os
    import tempfile

    e, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("CollapseMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob_name = ob.name
    ob.data.materials.append(mat)
    with bpy.context.temp_override(material=mat, object=ob):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_group()
    yield

    def group_item():
        m = bpy.data.materials["CollapseMat"]
        st = next(n for n in m.node_tree.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
                  and any((it.name or "").startswith("Group") for it in n.layer_items))
        return next(it for it in st.layer_items if (it.name or "").startswith("Group"))

    # Groups default to expanded; collapse it (as the tree-view arrow would).
    t.assertTrue(group_item().show_expanded)
    group_item().show_expanded = False
    t.assertFalse(group_item().show_expanded)

    # Undo an unrelated layer operation; the collapse must persist.
    m = bpy.data.materials["CollapseMat"]
    ob2 = bpy.data.objects[ob_name]
    st = next(n for n in m.node_tree.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack'
              and any((it.name or "").startswith("Group") for it in n.layer_items))
    st.active_index = len(st.layer_items) - 1
    area = window.screen.areas[0]
    with bpy.context.temp_override(window=window, area=area, material=m, object=ob2):
        bpy.ops.material.texture_layer_add_fill()
    yield
    pos = ui.get_area_center_from_spacetype(window, 'VIEW_3D')
    e.cursor_position_set(*pos, move=True)
    yield
    yield e.ctrl.z()
    t.assertFalse(group_item().show_expanded)

    # Save and reload; the collapse must persist.
    path = os.path.join(tempfile.gettempdir(), "tl_collapse.blend")
    bpy.ops.wm.save_as_mainfile(filepath=path)
    with bpy.context.temp_override(window=window):
        bpy.ops.wm.open_mainfile(filepath=path)
    t.assertFalse(group_item().show_expanded)


def test_texture_layer_ungroup():
    """Ungroup a Texture Layer Group: the nested stack's items are flattened
    into the outer stack at the group's position, and the now-orphaned nested
    stack node is deleted."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    stack_name = stack.name
    nt = mat.node_tree

    # Add a group at the top, then add a Fill INTO the nested stack so we have
    # something to flatten besides the nested stack's empty self.
    stack.active_index = 0
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_group()
    yield
    stack = nt.nodes[stack_name]
    nested = next(
        n for n in nt.nodes
        if n.bl_idname == 'ShaderNodeTextureLayerStack' and n != stack
    )
    nested_name = nested.name

    # Move a sibling layer into the nested stack so it has content.
    # The group sits at index 0 (we activated 0 before adding); pick a sibling.
    other_idx = 1
    sibling_name = stack.layer_items[other_idx].name
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_name, from_index=other_idx,
            to_stack=nested_name, to_index=0,
        )
    yield
    stack = nt.nodes[stack_name]
    nested = nt.nodes[nested_name]
    t.assertEqual(len(nested.layer_items), 1)
    t.assertEqual(nested.layer_items[0].name, sibling_name)

    # The reparent operator switches the active node to the destination stack,
    # so we have to switch back to the outer stack before invoking ungroup —
    # otherwise the ungroup poll resolves "active stack" to the nested one.
    nt.nodes.active = stack
    group_idx = next(
        i for i, it in enumerate(stack.layer_items)
        if "Group" in it.name
    )
    stack.active_index = group_idx
    outer_count_before = len(stack.layer_items)
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_ungroup()
    yield
    stack = nt.nodes[stack_name]

    # The nested stack node is gone, the group layer item is gone, and the
    # nested stack's single item now lives in the outer stack at the group's
    # original position.
    t.assertNotIn(nested_name, [n.name for n in nt.nodes])
    # Outer count: removed 1 (group), added 1 (nested item) → unchanged.
    t.assertEqual(len(stack.layer_items), outer_count_before)
    t.assertEqual(stack.layer_items[group_idx].name, sibling_name)


def test_texture_layer_remove_deletes_closure_zone():
    """Removing a closure-typed layer deletes its closure zone — the paired
    Closure Input/Output nodes and the body node inside — and Convert to
    Group refuses closure layers (a group takes a bundle, not a closure)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    nt = mat.node_tree
    stack_name = stack.name

    # Build a closure zone with a group body and drop its Closure output onto
    # the stack's extend socket, creating a closure-typed layer.
    closure_in = nt.nodes.new("NodeClosureInput")
    closure_out = nt.nodes.new("NodeClosureOutput")
    t.assertTrue(closure_in.pair_with_output(closure_out))
    closure_out.input_items.new('BUNDLE', "Bundle")
    closure_out.output_items.new('BUNDLE', "Bundle")
    group = _make_tint_adjustment_group()
    body = nt.nodes.new('ShaderNodeGroup')
    body.node_tree = group
    nt.links.new(closure_in.outputs["Bundle"],
                 next(s for s in body.inputs if s.type == 'BUNDLE'))
    nt.links.new(next(s for s in body.outputs if s.type == 'BUNDLE'),
                 closure_out.inputs["Bundle"])
    closure_socket = next(s for s in closure_out.outputs if s.type == 'CLOSURE')
    extend = next(s for s in stack.inputs if s.identifier == "__extend__")
    nt.links.new(closure_socket, extend, handle_dynamic_sockets=True)
    yield

    stack = nt.nodes[stack_name]
    closure_idx = len(stack.layer_items) - 1
    t.assertEqual(stack.layer_items[closure_idx].item_type, 'CLOSURE')
    closure_in_name = closure_in.name
    closure_out_name = closure_out.name
    body_name = body.name

    stack.active_index = closure_idx

    # Convert to Group now wraps a closure (adjustment) layer in an adjustment
    # group, so it is available here.
    with bpy.context.temp_override(window=window, area=properties_area):
        t.assertTrue(bpy.ops.material.texture_layer_convert_to_group.poll())

    layers_before = len(stack.layer_items)
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_remove()
    yield

    stack = nt.nodes[stack_name]
    t.assertEqual(len(stack.layer_items), layers_before - 1)
    node_names = [n.name for n in nt.nodes]
    t.assertNotIn(closure_out_name, node_names)
    t.assertNotIn(closure_in_name, node_names)
    t.assertNotIn(body_name, node_names)


def test_texture_layer_remove_deletes_nested():
    """Removing a group layer must also delete its nested stack and any
    nodes that the nested stack owned (e.g. a Combine Bundle feeding one of
    its layers) — they would otherwise be left orphaned in the node tree."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    stack_name = stack.name
    nt = mat.node_tree

    # Add a group, then push a sibling INTO the nested stack so it has content
    # we can verify gets cleaned up.
    stack.active_index = 0
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_group()
    yield
    stack = nt.nodes[stack_name]
    nested = next(
        n for n in nt.nodes
        if n.bl_idname == 'ShaderNodeTextureLayerStack' and n != stack
    )
    nested_name = nested.name

    # Move a sibling into the nested stack so it gets a real layer item.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_name, from_index=1,
            to_stack=nested_name, to_index=0,
        )
    yield
    stack = nt.nodes[stack_name]
    nested = nt.nodes[nested_name]

    # Wire a Combine Bundle into the nested stack's first (moved-in) item.
    nested_item = nested.layer_items[0]
    combine = nt.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Color")
    combine_name = combine.name
    nested_socket = next(
        s for s in nested.inputs if s.identifier == "Layer_" + str(nested_item.identifier)
    )
    nt.links.new(combine.outputs["Bundle"], nested_socket)
    yield

    # Reparent set the active node to the nested stack — flip it back.
    nt.nodes.active = stack
    group_idx = next(
        i for i, it in enumerate(stack.layer_items)
        if "Group" in it.name
    )
    stack.active_index = group_idx
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_remove()
    yield
    stack = nt.nodes[stack_name]

    # The nested stack node was deleted, AND the Combine Bundle that only
    # fed the nested stack was deleted too.
    node_names = {n.name for n in nt.nodes}
    t.assertNotIn(nested_name, node_names)
    t.assertNotIn(combine_name, node_names)


def test_mask_stack_node():
    """The Mask Stack node grows float Mask + Opacity sockets per item, mirroring
    Texture Layer Stack but for floats."""
    import bpy

    _, t, window = ui.test_window()

    mat = bpy.data.materials.new("MaskStackMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    stack = nt.nodes.new("ShaderNodeMaskStack")
    stack.mask_items.new(name="Top")
    stack.mask_items.new(name="Mid")
    stack.mask_items.new(name="Base")
    yield

    t.assertEqual(len(stack.mask_items), 3)
    t.assertEqual([item.name for item in stack.mask_items], ["Top", "Mid", "Base"])

    # Every item gets Mask + Opacity (Float), including the base: its Opacity
    # blends between no mask and the base mask.
    input_socket_names = [s.name for s in stack.inputs]
    t.assertEqual(input_socket_names, ["Top", "Opacity",
                                       "Mid", "Opacity",
                                       "Base", "Opacity",
                                       ""])  # __extend__
    t.assertTrue(all(
        s.type == 'VALUE'
        for s in stack.inputs
        if s.identifier != "__extend__"
    ))
    t.assertEqual(stack.outputs["Result"].type, 'VALUE')

    # Identifier-based socket lookup matches the convention from texture stack.
    base = stack.mask_items[2]
    base_sock = next(s for s in stack.inputs if s.identifier == "Mask_" + str(base.identifier))
    t.assertEqual(base_sock.name, "Base")


def test_mask_stack_extend_socket():
    """Dragging a float link onto the Mask Stack's hollow extend socket adds a
    new mask layer at the end (matches Texture Layer Stack behavior)."""
    import bpy

    _, t, window = ui.test_window()

    mat = bpy.data.materials.new("MaskExtendMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    stack = nt.nodes.new("ShaderNodeMaskStack")
    stack.mask_items.new(name="Mask")
    first_id = stack.mask_items[0].identifier

    value_node = nt.nodes.new('ShaderNodeValue')
    value_node.outputs[0].default_value = 0.5
    yield

    extend = next(s for s in stack.inputs if s.identifier == "__extend__")
    nt.links.new(value_node.outputs[0], extend, handle_dynamic_sockets=True)
    yield

    t.assertEqual(len(stack.mask_items), 2)
    new_item = stack.mask_items[1]
    t.assertEqual(new_item.name, "Mask.001")
    t.assertEqual(stack.active_index, 1)
    new_sock = next(s for s in stack.inputs if s.identifier == "Mask_" + str(new_item.identifier))
    t.assertTrue(new_sock.is_linked)
    # The previous-only item is now at index 0 and should have an Opacity socket.
    # The new item is the base, which also has an Opacity socket (base masks
    # blend between no mask and the mask).
    ids = [s.identifier for s in stack.inputs]
    t.assertIn("Opacity_" + str(first_id), ids)
    t.assertIn("Opacity_" + str(new_item.identifier), ids)


def test_mask_stack_inline_roundtrip():
    """Wiring a Mask Stack into a Texture Layer's Mask socket and inlining the
    material expands both stacks. The inlined tree must not contain either
    stack node anymore — both have been fully expanded."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("MaskInlineMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    stack.layer_items.new(name="Top")
    stack.layer_items.new(name="Base")

    sep = nt.nodes.new("NodeSeparateBundle")
    sep.bundle_items.new('RGBA', "Color")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])
    nt.links.new(sep.outputs["Color"], bsdf.inputs["Base Color"])

    # Top + Base get a Combine Bundle each so the stack composites a real result.
    top_id = stack.layer_items[0].identifier
    base_id = stack.layer_items[-1].identifier
    for item_id, color in ((top_id, (0.0, 1.0, 0.0, 1.0)), (base_id, (1.0, 0.0, 0.0, 1.0))):
        cb = nt.nodes.new('NodeCombineBundle')
        cb.bundle_items.new('RGBA', "Color")
        cb.inputs["Color"].default_value = color
        layer_sock = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item_id))
        nt.links.new(cb.outputs["Bundle"], layer_sock)

    # Wire a Mask Stack into the top layer's Mask socket. Two mask items: a
    # constant 0.5 base and a constant 1.0 top (multiply blend on top).
    mask_stack = nt.nodes.new("ShaderNodeMaskStack")
    mask_stack.mask_items.new(name="Top")
    mask_stack.mask_items.new(name="Base")
    mask_stack.mask_items[1].blend_type = 'MIX'
    mask_stack.mask_items[0].blend_type = 'MULTIPLY'
    base_mask = mask_stack.mask_items[-1]
    top_mask = mask_stack.mask_items[0]
    base_mask_sock = next(
        s for s in mask_stack.inputs if s.identifier == "Mask_" + str(base_mask.identifier)
    )
    base_mask_sock.default_value = 0.5
    top_mask_sock = next(
        s for s in mask_stack.inputs if s.identifier == "Mask_" + str(top_mask.identifier)
    )
    top_mask_sock.default_value = 1.0

    layer_mask_sock = next(
        s for s in stack.inputs if s.identifier == "Mask_" + str(top_id)
    )
    nt.links.new(mask_stack.outputs["Result"], layer_mask_sock)
    yield

    inline_wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
    inlined = inline_wrapper.node_tree

    # Both stack node types must be fully expanded — neither should appear in
    # the inlined tree.
    stack_nodes = [
        n for n in inlined.nodes
        if n.bl_idname in {'ShaderNodeTextureLayerStack', 'ShaderNodeMaskStack'}
    ]
    t.assertEqual(stack_nodes, [])

    # The Material Output's Surface input is connected (basic well-formedness).
    out_inlined = next(n for n in inlined.nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
    t.assertTrue(out_inlined.inputs["Surface"].is_linked)
    del inline_wrapper


def test_texture_layer_add_white_mask():
    """Adding a White Mask creates a Mask Stack feeding the active layer's
    Mask socket, with one mask item set to a constant 1.0 float."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    nt = mat.node_tree
    # Active a non-base layer: White Mask poll rejects the base.
    stack.active_index = 0
    item = stack.layer_items[0]

    initial_mask_stacks = sum(
        1 for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack'
    )
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_white_mask()
    yield

    # A Mask Stack appeared, fed into the layer's Mask socket.
    mask_stacks = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack']
    t.assertEqual(len(mask_stacks), initial_mask_stacks + 1)
    new_mask_stack = mask_stacks[-1]
    layer_mask_sock = next(
        s for s in stack.inputs if s.identifier == "Mask_" + str(item.identifier)
    )
    t.assertTrue(layer_mask_sock.is_linked)
    t.assertEqual(layer_mask_sock.links[0].from_node, new_mask_stack)

    # The new mask item exists with default value 1.0.
    t.assertEqual(len(new_mask_stack.mask_items), 1)
    new_mask_item = new_mask_stack.mask_items[0]
    new_sock = next(
        s for s in new_mask_stack.inputs
        if s.identifier == "Mask_" + str(new_mask_item.identifier)
    )
    t.assertEqual(new_sock.default_value, 1.0)

    # Adding a second mask reuses the same Mask Stack.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_black_mask()
    yield
    mask_stacks2 = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack']
    t.assertEqual(len(mask_stacks2), initial_mask_stacks + 1)
    t.assertEqual(len(new_mask_stack.mask_items), 2)


def _make_mask_generator_group():
    """A minimal MASK_GENERATOR-tagged shader node group exposing a single
    Float output (a constant 0.5). This is what mask assets look like in the
    library — they output a single mask value, not a bundle."""
    import bpy
    g = bpy.data.node_groups.new("TestMaskGen", "ShaderNodeTree")
    g.usage = {'MASK_GENERATOR'}
    g.interface.new_socket("Mask", in_out='OUTPUT', socket_type='NodeSocketFloat')
    in_node = g.nodes.new('NodeGroupInput')
    out_node = g.nodes.new('NodeGroupOutput')
    val = g.nodes.new('ShaderNodeValue')
    val.outputs[0].default_value = 0.5
    g.links.new(val.outputs[0], out_node.inputs["Mask"])
    g.asset_mark()
    _ = in_node  # group input is a placeholder — masks generate, they don't transform
    return g


def test_texture_layer_add_mask_asset():
    """Adding a mask asset via the Mask menu must route through the active
    layer's mask stack — not insert a new layer fed at the Bundle socket."""
    import bpy
    import os

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    mat.node_tree.nodes.active = stack
    bpy.context.view_layer.objects.active = bpy.context.active_object
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'

    nt = mat.node_tree
    # Pick a non-base layer.
    stack.active_index = 0
    item = stack.layer_items[0]

    group = _make_mask_generator_group()
    group_name = group.name
    yield

    initial_layers = len(stack.layer_items)
    initial_mask_stacks = sum(
        1 for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack'
    )
    relative_id = "NodeTree" + os.sep + group_name
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node_group(
            'EXEC_DEFAULT',
            asset_library_type='LOCAL',
            asset_library_identifier='',
            relative_asset_identifier=relative_id,
            as_mask=True,
        )
    yield

    # No new layer was inserted into the texture layer stack.
    t.assertEqual(len(stack.layer_items), initial_layers)
    # A Mask Stack was created and fed into the active layer's Mask socket.
    mask_stacks = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack']
    t.assertEqual(len(mask_stacks), initial_mask_stacks + 1)
    new_mask_stack = mask_stacks[-1]
    layer_mask_sock = next(
        s for s in stack.inputs if s.identifier == "Mask_" + str(item.identifier)
    )
    t.assertTrue(layer_mask_sock.is_linked)
    t.assertEqual(layer_mask_sock.links[0].from_node, new_mask_stack)
    # The Mask Stack has a new mask item, and the asset's group node is wired
    # into that item's Mask socket.
    t.assertEqual(len(new_mask_stack.mask_items), 1)
    new_mask_item = new_mask_stack.mask_items[0]
    new_mask_sock = next(
        s for s in new_mask_stack.inputs
        if s.identifier == "Mask_" + str(new_mask_item.identifier)
    )
    t.assertTrue(new_mask_sock.is_linked)
    src = new_mask_sock.links[0].from_node
    t.assertEqual(src.bl_idname, 'ShaderNodeGroup')
    t.assertEqual(src.node_tree.name, group_name)


def test_adjustment_closure_inline():
    """Inline a stack with a closure adjustment in the middle of multiple
    layers. The closure evaluates with the accumulator (Base composited with
    Mid) as its bundle input, feeding the accumulator through directly — no
    duplicate partial-accumulator subtree is emitted, so every Mix node in
    the inlined tree is reachable from the output (no orphans)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("ClosureMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")

    sep = nt.nodes.new("NodeSeparateBundle")
    sep.bundle_items.new('RGBA', "Color")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])
    nt.links.new(sep.outputs["Color"], bsdf.inputs["Base Color"])

    # A pass-through closure zone for the adjustment layer.
    closure_in = nt.nodes.new("NodeClosureInput")
    closure_out = nt.nodes.new("NodeClosureOutput")
    t.assertTrue(closure_in.pair_with_output(closure_out))
    closure_out.input_items.new('BUNDLE', "Bundle")
    closure_out.output_items.new('BUNDLE', "Bundle")
    nt.links.new(closure_in.outputs["Bundle"], closure_out.inputs["Bundle"])
    closure_socket = next(s for s in closure_out.outputs if s.type == 'CLOSURE')
    yield

    # Build [Adjustment (closure), Mid, Base] via the extend socket.
    def extend_drop(from_socket):
        extend = next(s for s in stack.inputs if s.identifier == "__extend__")
        nt.links.new(from_socket, extend, handle_dynamic_sockets=True)

    extend_drop(closure_socket)
    for color in ((0.0, 1.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0)):
        cb = nt.nodes.new('NodeCombineBundle')
        cb.bundle_items.new('RGBA', "Color")
        cb.inputs["Color"].default_value = color
        extend_drop(next(s for s in cb.outputs if s.type == 'BUNDLE'))
    yield

    t.assertEqual([i.item_type for i in stack.layer_items],
                  ['CLOSURE', 'BUNDLE', 'BUNDLE'])

    inline_wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
    inlined = inline_wrapper.node_tree

    # Stack must be fully expanded and the tree well-formed.
    stack_nodes = [n for n in inlined.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(stack_nodes, [])
    out_inlined = next(n for n in inlined.nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
    t.assertTrue(out_inlined.inputs["Surface"].is_linked)

    # One Mix per composed channel: Mid over Base, then the closure's
    # (pass-through) result over the accumulator.
    visited = set()
    reachable_mixes = [0]

    def walk(node):
        if node in visited:
            return
        visited.add(node)
        if node.bl_idname == 'ShaderNodeMix':
            reachable_mixes[0] += 1
        for inp in node.inputs:
            for link in inp.links:
                walk(link.from_node)

    walk(out_inlined)
    t.assertEqual(reachable_mixes[0], 2)

    # No orphaned subtrees: every Mix node emitted during inlining is
    # reachable from the material output.
    total_mixes = sum(1 for n in inlined.nodes if n.bl_idname == 'ShaderNodeMix')
    t.assertEqual(total_mixes, reachable_mixes[0])
    del inline_wrapper


def test_texture_layer_closure_region_inline():
    """Closure-typed layers: the inliner evaluates each closure layer's body
    with the running accumulator as its bundle input. One closure zone is
    shared by two layers of the same stack — each layer must get its own
    evaluation context (a shared context would collide in the value cache and
    feed the second layer the first layer's accumulator)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("ClosureRegionMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")

    sep = nt.nodes.new("NodeSeparateBundle")
    sep.bundle_items.new('RGBA', "Color")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])
    nt.links.new(sep.outputs["Color"], bsdf.inputs["Base Color"])

    # A single pass-through closure zone, used by two layers.
    closure_in = nt.nodes.new("NodeClosureInput")
    closure_out = nt.nodes.new("NodeClosureOutput")
    t.assertTrue(closure_in.pair_with_output(closure_out))
    closure_out.input_items.new('BUNDLE', "Bundle")
    closure_out.output_items.new('BUNDLE', "Bundle")
    nt.links.new(closure_in.outputs["Bundle"], closure_out.inputs["Bundle"])
    closure_socket = next(s for s in closure_out.outputs if s.type == 'CLOSURE')
    yield

    # Build the stack top-to-bottom via the extend socket:
    # [Adjust2 (closure), Adjust1 (closure), Mid, Base].
    def extend_drop(from_socket):
        extend = next(s for s in stack.inputs if s.identifier == "__extend__")
        nt.links.new(from_socket, extend, handle_dynamic_sockets=True)

    extend_drop(closure_socket)
    extend_drop(closure_socket)
    for color in ((0.0, 1.0, 0.0, 1.0), (1.0, 0.0, 0.0, 1.0)):
        cb = nt.nodes.new('NodeCombineBundle')
        cb.bundle_items.new('RGBA', "Color")
        cb.inputs["Color"].default_value = color
        extend_drop(next(s for s in cb.outputs if s.type == 'BUNDLE'))
    yield

    t.assertEqual(len(stack.layer_items), 4)
    t.assertEqual(stack.layer_items[0].item_type, 'CLOSURE')
    t.assertEqual(stack.layer_items[1].item_type, 'CLOSURE')
    t.assertEqual(stack.layer_items[2].item_type, 'BUNDLE')
    t.assertEqual(stack.layer_items[3].item_type, 'BUNDLE')

    inline_wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
    inlined = inline_wrapper.node_tree

    # Stack must be fully expanded and the tree well-formed.
    stack_nodes = [n for n in inlined.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack']
    t.assertEqual(stack_nodes, [])
    out_inlined = next(n for n in inlined.nodes if n.bl_idname == 'ShaderNodeOutputMaterial')
    t.assertTrue(out_inlined.inputs["Surface"].is_linked)

    # Composition emits one Mix per channel per composed layer: Mid over Base,
    # then each closure layer's (pass-through) result over the accumulator.
    visited = set()
    mix_count = [0]

    def walk(node):
        if node in visited:
            return
        visited.add(node)
        if node.bl_idname == 'ShaderNodeMix':
            mix_count[0] += 1
        for inp in node.inputs:
            for link in inp.links:
                walk(link.from_node)

    walk(out_inlined)
    t.assertGreaterEqual(
        mix_count[0],
        3,
        f"Inlined tree has only {mix_count[0]} Mix nodes feeding the BSDF — "
        "expected at least 3 (Mid over Base plus one per closure layer).",
    )
    del inline_wrapper


def test_shader_node_tree_usage():
    """ShaderNodeTree exposes a usage bitflag (Material / Group / Generator / etc)."""
    import bpy
    nt = bpy.data.node_groups.new("usage_test", "ShaderNodeTree")
    # Default usage is empty (no flags). Asset enumeration code treats this as
    # "not exposed via the texture-layer add menus".
    assert nt.usage == set()
    nt.usage = {'TEXTURE_GENERATOR'}
    assert nt.usage == {'TEXTURE_GENERATOR'}
    nt.usage = {'TEXTURE_GENERATOR', 'MASK_GENERATOR'}
    assert nt.usage == {'TEXTURE_GENERATOR', 'MASK_GENERATOR'}
    bpy.data.node_groups.remove(nt)
    yield  # generator test stub


def test_node_editor_stack_item_operators():
    """The generic NODE_OT item operators manage stack items from the node
    editor, and the sidebar draws the items list for the active stack node."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    nt = mat.node_tree

    # Turn the properties editor into a shader node editor showing the
    # material's tree (it follows the active object's material).
    area = ui.get_window_area_by_type(window, 'PROPERTIES')
    area.type = 'NODE_EDITOR'
    area.spaces.active.tree_type = 'ShaderNodeTree'
    area.spaces.active.show_region_ui = True
    yield

    t.assertEqual(area.spaces.active.edit_tree, nt)

    nt.nodes.active = stack
    stack.select = True
    yield  # Let the sidebar draw the items list once (uiList + color prop).

    stack.active_index = 2
    with bpy.context.temp_override(window=window, area=area):
        # Add inserts after the active item, initialized from it.
        bpy.ops.node.texture_layer_stack_item_add(show_dialog=False)
        t.assertEqual(len(stack.layer_items), 4)
        t.assertEqual(stack.active_index, 3)
        t.assertEqual(stack.layer_items[3].name, "Base.001")

        # Move the new item up, then remove it.
        bpy.ops.node.texture_layer_stack_item_move(direction='UP')
        t.assertEqual(stack.active_index, 2)
        t.assertEqual(stack.layer_items[2].name, "Base.001")
        bpy.ops.node.texture_layer_stack_item_remove()
        t.assertEqual(len(stack.layer_items), 3)
        t.assertEqual([i.name for i in stack.layer_items], ["Top", "Middle", "Base"])

    # The socket-type color shown in the list is exposed on layer items.
    t.assertEqual(len(stack.layer_items[0].color), 4)

    # Mask stack: same generic operators under their own idnames.
    mask_stack = nt.nodes.new("ShaderNodeMaskStack")
    mask_stack.mask_items.new(name="MaskA")
    mask_stack.mask_items.new(name="MaskB")
    nt.nodes.active = mask_stack
    yield

    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.node.mask_stack_item_add(show_dialog=False)
        t.assertEqual(len(mask_stack.mask_items), 3)
        bpy.ops.node.mask_stack_item_remove()
        t.assertEqual(len(mask_stack.mask_items), 2)


def test_texture_layer_selection_sync():
    """Selection is synced both ways between the layers and the node editor: a
    layer points the active node and the node selection at the nodes making it
    up, and selecting one of those nodes makes its layer active."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("SyncMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    bpy.context.active_object.data.materials.append(mat)
    nt = mat.node_tree

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_node(type='ShaderNodeTexNoise', channel="Roughness")
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    noise = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexNoise')
    noise_combine = next(l.to_node for l in nt.links if l.from_node == noise)
    fill = next(n for n in nt.nodes
                if n.bl_idname == 'NodeCombineBundle' and n != noise_combine)

    # The added layer is active, and points at its own nodes rather than at the
    # stack node: both the generator and the Combine Bundle it feeds.
    t.assertEqual(stack.active_index, 0)
    t.assertEqual(nt.nodes.active, noise_combine)
    t.assertEqual({n.name for n in nt.nodes if n.select}, {noise.name, noise_combine.name})

    # Selecting a node of the other layer makes that layer active, so the layer
    # operators act on it: removing takes the Fill layer, not the Noise one.
    # (Names, since removing nodes invalidates the references above.)
    fill_name, noise_name, combine_name = fill.name, noise.name, noise_combine.name
    nt.nodes.active = fill
    yield
    with bpy.context.temp_override(window=window, area=properties_area):
        t.assertTrue(bpy.ops.material.texture_layer_remove.poll())
        bpy.ops.material.texture_layer_remove()
    yield
    t.assertEqual(len(stack.layer_items), 1)
    node_names = [n.name for n in nt.nodes]
    t.assertIn(noise_name, node_names)
    t.assertNotIn(fill_name, node_names)

    # The layer left behind is active and selected in its place.
    t.assertEqual(nt.nodes.active.name, combine_name)
    t.assertEqual({n.name for n in nt.nodes if n.select}, {noise_name, combine_name})


def test_texture_layer_node_placement_order():
    """Layer nodes are placed in the node editor in the same order as the
    layers and the stack's item sockets, so their links do not cross: each
    layer's nodes sit below those of the layers above it, and inserting in the
    middle pushes the layers below down instead of overlapping them."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("PlacementMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    bpy.context.active_object.data.materials.append(mat)
    nt = mat.node_tree
    bsdf = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeBsdfPrincipled')

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_fill()
        bpy.ops.material.texture_layer_add_node(type='ShaderNodeTexNoise', channel="Roughness")
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')

    # The bootstrapped Fill covers the BSDF's channels in the BSDF's own input
    # order, so the links from the Separate Bundle do not cross either.
    fill = next(n for n in nt.nodes if n.bl_idname == 'NodeCombineBundle' and n.label == "Fill")
    channels = [i.name for i in fill.bundle_items]
    t.assertEqual(channels, [s.name for s in bsdf.inputs if s.name in set(channels)])

    # The bootstrapped chain reads left to right without nodes overlapping:
    # Combine Bundle, stack, Separate Bundle, BSDF.
    separate = next(n for n in nt.nodes if n.bl_idname == 'NodeSeparateBundle')
    chain = [fill, stack, separate, bsdf]
    for left, right in zip(chain, chain[1:]):
        t.assertLess(left.location[0] + max(left.width, 140.0), right.location[0],
                     "%s overlaps %s" % (left.name, right.name))

    def source_y(index):
        item = stack.layer_items[index]
        socket = next(s for s in stack.inputs if s.identifier == "Layer_" + str(item.identifier))
        return socket.links[0].from_node.location[1] if socket.is_linked else None

    def check_order():
        ys = [source_y(i) for i in range(len(stack.layer_items))]
        ys = [y for y in ys if y is not None]
        for upper, lower in zip(ys, ys[1:]):
            # Ordered top to bottom, and spaced enough not to overlap.
            t.assertGreater(upper - lower, 100.0, str(ys))

    check_order()

    # Insert in the middle: the layer below has to make room.
    stack.active_index = 1
    with bpy.context.temp_override(window=window, area=properties_area):
        nt.nodes.active = stack
        bpy.ops.material.texture_layer_add_fill()
    yield

    t.assertEqual(len(stack.layer_items), 4)
    check_order()


def test_texture_layer_selection_sync_from_node_editor():
    """Clicking a layer's node in the node editor makes its layer active."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("SyncClickMat")
    mat.use_nodes = True
    bpy.ops.mesh.primitive_cube_add()
    bpy.context.active_object.data.materials.append(mat)
    nt = mat.node_tree

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
        bpy.ops.material.texture_layer_add_fill()
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    base_fill = next(l.from_node for l in nt.links
                     if l.to_node == stack and l.to_socket.identifier ==
                     "Layer_" + str(stack.layer_items[1].identifier))

    # The added layer (index 0) is active; the base layer's Combine Bundle is
    # the click target below.
    t.assertEqual(stack.active_index, 0)

    area = ui.largest_area(window.screen)
    area.type = 'NODE_EDITOR'
    area.spaces.active.tree_type = 'ShaderNodeTree'
    yield
    t.assertEqual(area.spaces.active.edit_tree, nt)

    region = next(r for r in area.regions if r.type == 'WINDOW')
    with bpy.context.temp_override(window=window, area=area, region=region):
        bpy.ops.node.view_all()
    yield

    x, y = region.view2d.view_to_region(
        base_fill.location[0] + 10.0, base_fill.location[1] - 10.0, clip=False)
    with bpy.context.temp_override(window=window, area=area, region=region):
        bpy.ops.node.select(location=(x, y))
    yield

    t.assertEqual(nt.nodes.active, base_fill)
    t.assertEqual(stack.active_index, 1)


def test_texture_layer_add_builtin_node():
    """Built-in shader nodes with a texture-layer usage trait can be added as
    generator layers (wrapped in a Combine Bundle keyed to a BSDF channel) and
    as masks (their float output feeds the layer's mask stack)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf, stack, sep = _build_shader_tree_with_layers()
    nt = mat.node_tree
    sep.bundle_items.new('RGBA', "Color")
    nt.links.new(sep.outputs["Color"], bsdf.inputs["Base Color"])
    nt.nodes.active = stack
    stack.active_index = 0

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node(type='ShaderNodeTexNoise', channel='Roughness')

    t.assertEqual(len(stack.layer_items), 4)
    t.assertEqual(stack.layer_items[stack.active_index].name, "Noise Texture")

    # The layer's source is a Combine Bundle keyed to the chosen channel; the
    # item type comes from the matching BSDF input (Roughness is a float), and
    # the node's float output feeds it.
    combine = next(
        n for n in nt.nodes
        if n.bl_idname == 'NodeCombineBundle' and len(n.bundle_items) == 1)
    t.assertEqual(combine.bundle_items[0].name, "Roughness")
    t.assertEqual(combine.bundle_items[0].socket_type, 'FLOAT')
    noise = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexNoise')
    in_link = next(l for l in nt.links if l.to_node == combine)
    t.assertEqual(in_link.from_node, noise)
    t.assertEqual(in_link.from_socket.name, "Factor")
    t.assertTrue(any(l.from_node == combine and l.to_node == stack for l in nt.links))

    # As mask: no new layer; the node's float output feeds a mask stack on the
    # active (non-base) layer.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node(type='ShaderNodeTexWave', as_mask=True)

    t.assertEqual(len(stack.layer_items), 4)
    mask_stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack')
    t.assertEqual([i.name for i in mask_stack.mask_items][0], "Wave Texture")
    wave = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexWave')
    t.assertTrue(any(
        l.from_node == wave and l.to_node == mask_stack and l.from_socket.name == "Factor"
        for l in nt.links))

    # A node with no float value of its own drives the mask from its Color, not
    # from the Alpha that accompanies it.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node(type='ShaderNodeTexImage', as_mask=True)
    image = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexImage')
    link = next(l for l in nt.links if l.from_node == image)
    t.assertEqual(link.to_node, mask_stack)
    t.assertEqual(link.from_socket.name, "Color")

    # The Paint mask, which adds an Image Texture of its own, matches it.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_paint_mask()
    paint_image = [n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexImage'][-1]
    link = next(l for l in nt.links if l.from_node == paint_image)
    t.assertEqual(link.to_node, mask_stack)
    t.assertEqual(link.from_socket.name, "Color")


def test_texture_layer_mask_remove():
    """The mask remove operator deletes the active mask item with its source
    nodes, and deletes the Mask Stack itself when the last item goes."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    nt = mat.node_tree
    nt.nodes.active = stack
    stack.active_index = 0

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_white_mask()
        bpy.ops.material.texture_layer_add_paint_mask()

    mask_stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack')
    t.assertEqual([i.name for i in mask_stack.mask_items], ["Paint", "White"])
    t.assertTrue(any(n.bl_idname == 'ShaderNodeTexImage' for n in nt.nodes))

    # Let the tree view draw once with mask child rows present.
    yield

    # Remove the active (Paint) mask: its image source node goes with it.
    nt.nodes.active = mask_stack
    mask_stack.active_index = 0
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_mask_remove()
    t.assertEqual([i.name for i in mask_stack.mask_items], ["White"])
    t.assertFalse(any(n.bl_idname == 'ShaderNodeTexImage' for n in nt.nodes))

    # Removing the last mask deletes the now-empty Mask Stack node, so the
    # layer's Mask input falls back to its default.
    mask_stack.active_index = 0
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_mask_remove()
    t.assertFalse(any(n.bl_idname == 'ShaderNodeMaskStack' for n in nt.nodes))
    mask_socket = next(s for s in stack.inputs if s.identifier.startswith("Mask_"))
    t.assertFalse(mask_socket.is_linked)


def test_texture_layer_mask_reparent():
    """The reparent operator also reorders masks within a Mask Stack and moves
    them between mask stacks, but never across stack kinds."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    nt = mat.node_tree
    nt.nodes.active = stack

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    def layer_mask_source(layer_index):
        socket = [s for s in stack.inputs if s.identifier.startswith("Mask_")][layer_index]
        return socket.links[0].from_node if socket.is_linked else None

    # Masks on the Top layer: [Paint, White]; a White mask on Middle. Adding a
    # mask activates it, so selecting the next layer means pointing the active
    # node back at the stack, as activating a layer row does.
    with bpy.context.temp_override(window=window, area=properties_area):
        stack.active_index = 0
        bpy.ops.material.texture_layer_add_white_mask()
        bpy.ops.material.texture_layer_add_paint_mask()
        nt.nodes.active = stack
        stack.active_index = 1
        bpy.ops.material.texture_layer_add_white_mask()

    stack_a = layer_mask_source(0)
    stack_b = layer_mask_source(1)
    t.assertEqual([i.name for i in stack_a.mask_items], ["Paint", "White"])
    t.assertEqual([i.name for i in stack_b.mask_items], ["White"])

    with bpy.context.temp_override(window=window, area=properties_area):
        # Reorder within one mask stack.
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_a.name, from_index=0, to_stack=stack_a.name, to_index=1)
        t.assertEqual([i.name for i in stack_a.mask_items], ["White", "Paint"])

        # Move the Paint mask to the other layer's mask stack; its image
        # source node moves with it.
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack_a.name, from_index=1, to_stack=stack_b.name, to_index=0)
        t.assertEqual([i.name for i in stack_a.mask_items], ["White"])
        t.assertEqual([i.name for i in stack_b.mask_items], ["Paint", "White"])
        image = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTexImage')
        t.assertTrue(any(l.from_node == image and l.to_node == stack_b for l in nt.links))

        # Mask <-> texture-layer moves are refused.
        result = bpy.ops.material.texture_layer_reparent(
            from_stack=stack_a.name, from_index=0, to_stack=stack.name, to_index=0)
        t.assertEqual(result, {'CANCELLED'})
        t.assertEqual([i.name for i in stack_a.mask_items], ["White"])
        t.assertEqual(len(stack.layer_items), 3)


def test_texture_layer_default_mask_blend():
    """Mask sources can declare a default blend mode: node group assets via
    ShaderNodeTree.default_mask_blend_type (baked into asset metadata), and
    built-in nodes via their type (AO defaults to Multiply)."""
    import bpy
    import os

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    nt = mat.node_tree
    nt.nodes.active = stack
    stack.active_index = 0

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'

    # RNA round-trip; Mix is the default.
    g = bpy.data.node_groups.new("BlendMaskGen", "ShaderNodeTree")
    t.assertEqual(g.default_mask_blend_type, 'MIX')
    g.usage = {'MASK_GENERATOR'}
    g.default_mask_blend_type = 'MULTIPLY'
    g.interface.new_socket("Mask", in_out='OUTPUT', socket_type='NodeSocketFloat')
    out_node = g.nodes.new('NodeGroupOutput')
    val = g.nodes.new('ShaderNodeValue')
    val.outputs[0].default_value = 0.5
    g.links.new(val.outputs[0], out_node.inputs["Mask"])
    g.asset_mark()  # Bakes usage + default blend into the asset metadata.
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        # Base mask first: the default blend only applies to non-base items.
        bpy.ops.material.texture_layer_add_white_mask()
        bpy.ops.material.texture_layer_add_node_group(
            'EXEC_DEFAULT',
            asset_library_type='LOCAL',
            asset_library_identifier='',
            relative_asset_identifier="NodeTree" + os.sep + g.name,
            as_mask=True,
        )

    mask_stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeMaskStack')
    t.assertEqual(mask_stack.mask_items[0].name, g.name)
    t.assertEqual(mask_stack.mask_items[0].blend_type, 'MULTIPLY')
    t.assertEqual(mask_stack.mask_items[1].blend_type, 'MIX')

    # Built-in node default: AO masks multiply over the masks below.
    nt.nodes.active = stack
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node(
            type='ShaderNodeAmbientOcclusion', as_mask=True)
    t.assertEqual(mask_stack.mask_items[0].name, "Ambient Occlusion")
    t.assertEqual(mask_stack.mask_items[0].blend_type, 'MULTIPLY')


def test_texture_layer_add_node_group_at_index():
    """The asset add operator accepts an explicit target stack and index, the
    backend of dropping an asset between rows of the layers tree view."""
    import bpy
    import os

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, _bsdf, stack, _sep = _build_shader_tree_with_layers()
    nt = mat.node_tree
    nt.nodes.active = stack
    stack.active_index = 0

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'

    g = bpy.data.node_groups.new("DropGen", "ShaderNodeTree")
    g.usage = {'TEXTURE_GENERATOR'}
    g.interface.new_socket("Bundle", in_out='OUTPUT', socket_type='NodeSocketBundle')
    out_node = g.nodes.new('NodeGroupOutput')
    combine = g.nodes.new('NodeCombineBundle')
    combine.bundle_items.new('RGBA', "Base Color")
    g.links.new(combine.outputs["Bundle"], out_node.inputs[0])
    g.asset_mark()
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node_group(
            'EXEC_DEFAULT',
            asset_library_type='LOCAL',
            asset_library_identifier='',
            relative_asset_identifier="NodeTree" + os.sep + g.name,
            to_stack=stack.name,
            to_index=2,
        )

    t.assertEqual(len(stack.layer_items), 4)
    t.assertEqual(stack.layer_items[2].name, g.name)
    t.assertEqual(stack.active_index, 2)
    # The group node feeds the inserted layer's Bundle socket.
    layer_sock = next(
        s for s in stack.inputs
        if s.identifier == "Layer_" + str(stack.layer_items[2].identifier))
    t.assertTrue(layer_sock.is_linked)
    t.assertEqual(layer_sock.links[0].from_node.node_tree, g)


def test_texture_layer_channels_inline():
    """Per-layer channel toggles: disabled channels are dropped during stack
    composition, so no Mix node is emitted for them and the layers below show
    through unchanged."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("ChannelsMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    stack = nt.nodes.new("ShaderNodeTextureLayerStack")
    sep = nt.nodes.new("NodeSeparateBundle")
    sep.bundle_items.new('RGBA', "Base Color")
    nt.links.new(stack.outputs["Result"], sep.inputs["Bundle"])
    nt.links.new(sep.outputs["Base Color"], bsdf.inputs["Base Color"])

    def extend_drop(from_socket):
        extend = next(s for s in stack.inputs if s.identifier == "__extend__")
        nt.links.new(from_socket, extend, handle_dynamic_sockets=True)

    # Two layers, each contributing Base Color + Roughness.
    for _i in range(2):
        cb = nt.nodes.new('NodeCombineBundle')
        cb.bundle_items.new('RGBA', "Base Color")
        cb.bundle_items.new('FLOAT', "Roughness")
        extend_drop(next(s for s in cb.outputs if s.type == 'BUNDLE'))
    yield

    top = stack.layer_items[0]
    t.assertTrue(top.channel_enabled("Roughness"))
    t.assertTrue(top.channel_enabled("Base Color"))

    def count_mixes():
        wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
        count = sum(1 for n in wrapper.node_tree.nodes if n.bl_idname == 'ShaderNodeMix')
        del wrapper
        return count

    # Both channels enabled: one Mix per channel of the top layer.
    t.assertEqual(count_mixes(), 2)

    # Disabling Roughness on the top layer drops its Mix; the base layer's
    # Roughness passes through.
    top.channel_set_enabled("Roughness", False)
    t.assertFalse(top.channel_enabled("Roughness"))
    t.assertEqual(count_mixes(), 1)

    # Re-enabling restores the contribution.
    top.channel_set_enabled("Roughness", True)
    t.assertEqual(count_mixes(), 2)

    # Toggles survive copying the material (item deep copy).
    mat_copy = mat.copy()
    top.channel_set_enabled("Roughness", False)
    copy_stack = next(
        n for n in mat_copy.node_tree.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    t.assertTrue(copy_stack.layer_items[0].channel_enabled("Roughness"))
    yield


def _build_principled_material():
    import bpy
    mat = bpy.data.materials.new("ChannelToggleMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    return mat, bsdf


def test_texture_channel_toggle():
    """MATERIAL_OT_texture_channel_toggle routes a BSDF input through the
    stack and back, preserving the value in both directions and reusing the
    stack across a full disable."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree
    bsdf.inputs["Roughness"].default_value = 0.77

    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    nt.nodes.active = bsdf
    yield  # Draw the Surface panel: BSDF details with the channel toggles.

    def stacks():
        return [n for n in nt.nodes if n.bl_idname == "ShaderNodeTextureLayerStack"]

    # Enable Roughness: bootstraps stack + Combine/Separate Bundle wiring.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
    yield

    t.assertEqual(len(stacks()), 1)
    sep = next(n for n in nt.nodes if n.bl_idname == "NodeSeparateBundle")
    comb = next(n for n in nt.nodes if n.bl_idname == "NodeCombineBundle")
    t.assertTrue(sep.define_signature)
    t.assertEqual([i.name for i in sep.bundle_items], ["Roughness"])
    t.assertTrue(bsdf.inputs["Roughness"].is_linked)

    # The base layer was seeded with the BSDF's value, so the evaluated
    # shader is unchanged by the toggle.
    rough_in = next(s for s in comb.inputs if s.name == "Roughness")
    t.assertAlmostEqual(rough_in.default_value, 0.77, places=5)
    wrapper = bpy.types.InlineShaderNodes.from_material(material=mat)
    inline_bsdf = next(
        n for n in wrapper.node_tree.nodes if n.bl_idname == 'ShaderNodeBsdfPrincipled')
    t.assertAlmostEqual(inline_bsdf.inputs["Roughness"].default_value, 0.77, places=5)
    del wrapper

    # Enable a second channel: same stack, new Separate Bundle item.
    bsdf.inputs["Metallic"].default_value = 0.33
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Metallic"].identifier)
    yield
    t.assertEqual(len(stacks()), 1)
    t.assertEqual([i.name for i in sep.bundle_items], ["Roughness", "Metallic"])

    # Draw the layer's Channels panel with the material-level state.
    stack = stacks()[0]
    nt.nodes.active = stack
    stack.active_index = 0
    yield

    # Disable Metallic after editing the base value: the edited value is
    # copied back to the BSDF socket, and the channel leaves every layer's
    # Combine Bundle since it is no longer part of the stack. The value is
    # preserved on the BSDF socket, so re-enabling would re-seed it.
    metal_in = next(s for s in comb.inputs if s.name == "Metallic")
    metal_in.default_value = 0.9
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Metallic"].identifier)
    yield
    t.assertFalse(bsdf.inputs["Metallic"].is_linked)
    t.assertAlmostEqual(bsdf.inputs["Metallic"].default_value, 0.9, places=5)
    t.assertEqual([i.name for i in sep.bundle_items], ["Roughness"])
    t.assertFalse(any(s.name == "Metallic" for s in comb.inputs))

    # Disable the last channel, re-enable: the existing stack is reused.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
        t.assertFalse(bsdf.inputs["Roughness"].is_linked)
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
    yield
    t.assertEqual(len(stacks()), 1)
    t.assertTrue(bsdf.inputs["Roughness"].is_linked)


def test_texture_channel_manual_link():
    """Manually linking the Separate Bundle to a BSDF input (as a node editor
    drag does) completes the channel: an extend-socket link creates the item
    named after the target and seeds the base layer; name-mismatched links
    are cross-wiring and leave other nodes alone."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree
    with bpy.context.temp_override(material=mat):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
    yield

    sep = next(n for n in nt.nodes if n.bl_idname == "NodeSeparateBundle")
    comb = next(n for n in nt.nodes if n.bl_idname == "NodeCombineBundle")

    # Extend-socket drag onto Metallic: item is created and base layer seeded.
    bsdf.inputs["Metallic"].default_value = 0.4
    nt.links.new(bsdf.inputs["Metallic"], sep.outputs["__extend__"],
                 handle_dynamic_sockets=True)
    yield
    t.assertEqual([i.name for i in sep.bundle_items], ["Roughness", "Metallic"])
    t.assertTrue(bsdf.inputs["Metallic"].is_linked)
    metal_in = next(s for s in comb.inputs if s.name == "Metallic")
    t.assertAlmostEqual(metal_in.default_value, 0.4, places=5)

    # Cross-wiring an existing item into a differently named input does not
    # touch the base layer.
    rough_out = next(s for s in sep.outputs if s.name == "Roughness")
    nt.links.new(bsdf.inputs["IOR"], rough_out, handle_dynamic_sockets=True)
    yield
    t.assertTrue(bsdf.inputs["IOR"].is_linked)
    t.assertFalse(any(s.name == "IOR" for s in comb.inputs))
    t.assertEqual([i.name for i in sep.bundle_items], ["Roughness", "Metallic"])


def test_texture_channel_combine_sync():
    """Bundle signatures propagate through the Texture Layer Stack: linking an
    empty Combine Bundle to a layer syncs its items from the BSDF's Separate
    Bundle. Every layer carries the full channel set (per-layer contribution is
    the layer's own toggle), so Sync restores a removed channel."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree
    with bpy.context.temp_override(material=mat):
        for name in ("Roughness", "Metallic"):
            bpy.ops.material.texture_channel_toggle(
                node=bsdf.name, socket=bsdf.inputs[name].identifier)
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == "ShaderNodeTextureLayerStack")

    area = ui.get_window_area_by_type(window, 'PROPERTIES')
    area.type = 'NODE_EDITOR'
    area.spaces.active.tree_type = 'ShaderNodeTree'
    yield
    t.assertEqual(area.spaces.active.edit_tree, nt)

    # A new layer starts untyped with a hollow socket; linking a fresh
    # Combine Bundle types it, and Sync fills the combine's items from the
    # BSDF's Separate Bundle, through the stack.
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.material.texture_layer_add()
    yield
    layer_socket = next(
        s for s in stack.inputs if s.type == 'CUSTOM' and s.identifier != "__extend__")
    comb2 = nt.nodes.new("NodeCombineBundle")
    nt.links.new(layer_socket, comb2.outputs["Bundle"], handle_dynamic_sockets=True)
    nt.nodes.active = comb2
    yield
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.node.sockets_sync(node_name=comb2.name)
    t.assertEqual([i.name for i in comb2.bundle_items], ["Roughness", "Metallic"])
    t.assertEqual([i.socket_type for i in comb2.bundle_items], ['FLOAT', 'FLOAT'])

    # A removed channel is out of sync now: every layer carries the full
    # channel set, so Sync restores it.
    comb2.bundle_items.remove(comb2.bundle_items["Metallic"])
    yield
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.node.sockets_sync(node_name=comb2.name)
    t.assertEqual([i.name for i in comb2.bundle_items], ["Roughness", "Metallic"])

    # A type conflict is out of sync: Sync rewrites to the channel signature.
    comb2.bundle_items["Roughness"].socket_type = 'RGBA'
    yield
    with bpy.context.temp_override(window=window, area=area):
        bpy.ops.node.sockets_sync(node_name=comb2.name)
    t.assertEqual([i.name for i in comb2.bundle_items], ["Roughness", "Metallic"])
    t.assertEqual([i.socket_type for i in comb2.bundle_items], ['FLOAT', 'FLOAT'])


def test_texture_layer_fill_textured_channels():
    """A new Fill layer carries the channels the material textures, falling
    back to the default channel set when none are enabled."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree
    with bpy.context.temp_override(material=mat):
        for name in ("Base Color", "Roughness"):
            bpy.ops.material.texture_channel_toggle(
                node=bsdf.name, socket=bsdf.inputs[name].identifier)
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == "ShaderNodeTextureLayerStack")
    nt.nodes.active = stack
    with bpy.context.temp_override(material=mat):
        bpy.ops.material.texture_layer_add_fill()
    yield

    fill = next(n for n in nt.nodes
                if n.bl_idname == "NodeCombineBundle" and n.label == "Fill")
    t.assertEqual([i.name for i in fill.bundle_items], ["Base Color", "Roughness"])
    t.assertEqual([i.socket_type for i in fill.bundle_items], ['RGBA', 'FLOAT'])


def test_texture_layer_mask_selected_ops():
    """With a mask row selected (active Mask Stack node + active mask index,
    as tree view activation sets), the layer operators act on the mask: the
    generic remove button removes that mask, and mask adds insert into the
    same stack above it."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree

    with bpy.context.temp_override(material=mat):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
        bpy.ops.material.texture_layer_add()  # Non-base layer with a Mask socket.
        bpy.ops.material.texture_layer_add_white_mask()
        bpy.ops.material.texture_layer_add_black_mask()
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == "ShaderNodeTextureLayerStack")
    mask_stack = next(n for n in nt.nodes if n.bl_idname == "ShaderNodeMaskStack")
    t.assertEqual([i.name for i in mask_stack.mask_items], ["Black", "White"])
    layer_count = len(stack.layer_items)

    # Select the second mask, as clicking its tree row does.
    nt.nodes.active = mask_stack
    mask_stack.active_index = 1
    yield

    # Adding a mask with a mask selected inserts into the same stack, above
    # the selected one.
    with bpy.context.temp_override(material=mat):
        t.assertTrue(bpy.ops.material.texture_layer_add_white_mask.poll())
        bpy.ops.material.texture_layer_add_white_mask()
    t.assertEqual([i.name for i in mask_stack.mask_items], ["Black", "White.001", "White"])

    # The generic remove operator removes the selected mask, not the layer.
    with bpy.context.temp_override(material=mat):
        t.assertTrue(bpy.ops.material.texture_layer_remove.poll())
        bpy.ops.material.texture_layer_remove()
    yield
    t.assertEqual([i.name for i in mask_stack.mask_items], ["Black", "White"])
    t.assertEqual(len(stack.layer_items), layer_count)

    # Asset/menu polls also accept the selected mask.
    with bpy.context.temp_override(material=mat):
        t.assertTrue(bpy.ops.material.texture_layer_add_paint_mask.poll())


def test_texture_layer_reorder_rebuilds_sockets():
    """Reordering layers within a stack rebuilds its per-item sockets: the
    layer moved off the base position gains Opacity/Mask sockets and the new
    base loses them (regression: the same-stack reparent skipped the update
    tag, leaving the previous order's sockets)."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat, bsdf = _build_principled_material()
    nt = mat.node_tree
    with bpy.context.temp_override(material=mat):
        bpy.ops.material.texture_channel_toggle(
            node=bsdf.name, socket=bsdf.inputs["Roughness"].identifier)
        bpy.ops.material.texture_layer_add_node(
            type="ShaderNodeTexNoise", channel="Roughness")
    yield

    stack = next(n for n in nt.nodes if n.bl_idname == "ShaderNodeTextureLayerStack")
    t.assertEqual([i.name for i in stack.layer_items], ["Noise Texture", "Fill"])
    t.assertEqual([s.name for s in stack.inputs],
                  ["Noise Texture", "Opacity", "Mask", "Fill", "Mask", ""])

    # Same-stack drag & drop reorder: Fill on top, Noise becomes the base.
    with bpy.context.temp_override(material=mat):
        bpy.ops.material.texture_layer_reparent(
            from_stack=stack.name, from_index=0, to_stack=stack.name, to_index=1)
    yield
    t.assertEqual([i.name for i in stack.layer_items], ["Fill", "Noise Texture"])
    t.assertEqual([s.name for s in stack.inputs],
                  ["Fill", "Opacity", "Mask", "Noise Texture", "Mask", ""])


def test_texture_layers_add_generator_routes_channel():
    """Adding a generator that targets a channel not yet routed through the
    stack must also wire that channel out through the Separate Bundle to the
    BSDF, otherwise the new layer feeds the stack but never reaches the shader."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("GenMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    bsdf_name = bsdf.name

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    bpy.context.view_layer.objects.active = ob
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    # Route only Metallic through a stack (bootstraps the stack + Separate Bundle
    # feeding the BSDF); Base Color stays a plain, unrouted value.
    metallic = next(s for s in bsdf.inputs if s.name == "Metallic")
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_channel_toggle(node=bsdf_name, socket=metallic.identifier)
    yield

    bsdf = nt.nodes[bsdf_name]
    t.assertFalse(next(s for s in bsdf.inputs if s.name == "Base Color").is_linked)

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    nt.nodes.active = stack

    # Add a built-in Checker Texture generator driving the Base Color channel.
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_node(
            'EXEC_DEFAULT', type="ShaderNodeTexChecker", channel="Base Color")
    yield

    # Base Color is now routed through the stack: fed by the Separate Bundle.
    bsdf = nt.nodes[bsdf_name]
    base_color = next(s for s in bsdf.inputs if s.name == "Base Color")
    t.assertTrue(base_color.is_linked)
    t.assertEqual(base_color.links[0].from_node.bl_idname, 'NodeSeparateBundle')


def test_texture_layers_add_default_and_teardown():
    """The (+) default add on a BSDF with no stack creates a stack + Fill layer
    covering Base Color, Roughness, Metallic and Normal, all routed to the BSDF.
    Removing that single layer tears the whole stack (and Separate Bundle) down
    again."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("DefMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    bsdf_name = bsdf.name

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    bpy.context.view_layer.objects.active = ob
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    def has(idname):
        return any(n.bl_idname == idname for n in nt.nodes)

    t.assertFalse(has('ShaderNodeTextureLayerStack'))
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
    yield

    t.assertTrue(has('ShaderNodeTextureLayerStack'))
    bsdf = nt.nodes[bsdf_name]
    for chan in ("Base Color", "Roughness", "Metallic", "Normal"):
        sock = next(s for s in bsdf.inputs if s.name == chan)
        t.assertTrue(sock.is_linked, chan)
        t.assertEqual(sock.links[0].from_node.bl_idname, 'NodeSeparateBundle')

    stack = next(n for n in nt.nodes if n.bl_idname == 'ShaderNodeTextureLayerStack')
    nt.nodes.active = stack
    stack.active_index = 0
    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_remove()
    yield

    t.assertFalse(has('ShaderNodeTextureLayerStack'))
    t.assertFalse(has('NodeSeparateBundle'))
    bsdf = nt.nodes[bsdf_name]
    t.assertFalse(next(s for s in bsdf.inputs if s.name == "Base Color").is_linked)


def test_texture_layers_fill_inherits_hide_value():
    """A Fill layer's Normal channel input inherits the BSDF socket's
    hide-value flag, so it does not show an editable vector; Roughness (with a
    visible value) does not."""
    import bpy

    _, t, window = ui.test_window()
    bpy.context.preferences.experimental.use_texture_layers = True

    mat = bpy.data.materials.new("NrmMat")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    bpy.ops.mesh.primitive_cube_add()
    ob = bpy.context.active_object
    ob.data.materials.append(mat)
    bpy.context.view_layer.objects.active = ob
    bpy.context.object.active_material = mat
    properties_area = ui.get_window_area_by_type(window, 'PROPERTIES')
    properties_area.spaces.active.context = 'MATERIAL'
    yield

    with bpy.context.temp_override(window=window, area=properties_area):
        bpy.ops.material.texture_layer_add_default()
    yield

    combine = next(n for n in nt.nodes if n.bl_idname == 'NodeCombineBundle')
    normal_in = next(s for s in combine.inputs if s.name == "Normal")
    t.assertTrue(normal_in.hide_value)
    rough_in = next(s for s in combine.inputs if s.name == "Roughness")
    t.assertFalse(rough_in.hide_value)
