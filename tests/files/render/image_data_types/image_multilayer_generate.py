# SPDX-FileCopyrightText: 2026 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

# Regenerates the multi-layer image texture render test assets:
#   textures/image_multilayer.exr    - a multi-part EXR with BaseColor / Roughness
#                                       / Normal passes, each spatially varying so
#                                       its effect on the shaded result is visible.
#   image_multilayer_passes.blend     - Image Texture node with the per-pass output
#                                       sockets wired into a Principled BSDF.
#   image_multilayer_bundle.blend     - same, routed through the node's Bundle
#                                       output and a Separate Bundle node.
#
# Run with:
#   blender --background --factory-startup --python image_multilayer_generate.py
#
# Afterwards regenerate the reference renders with BLENDER_TEST_UPDATE=1.

import math
import os

import bpy
import numpy as np
import OpenImageIO as oiio

DIR = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.join(DIR, "textures", "image_multilayer.exr")
RES = 128


def hsv_to_rgb(h, s, v):
    i = np.floor(h * 6).astype(int) % 6
    f = h * 6 - np.floor(h * 6)
    p, q, t = v * (1 - s), v * (1 - f * s), v * (1 - (1 - f) * s)
    r = np.choose(i, [v, q, p, p, t, v])
    g = np.choose(i, [t, v, v, q, p, p])
    b = np.choose(i, [p, p, t, v, v, q])
    return r, g, b


def write_exr():
    ys, xs = np.mgrid[0:RES, 0:RES].astype(np.float32)
    u, v = xs / (RES - 1), ys / (RES - 1)
    # BaseColor: horizontal hue sweep, clearly colorful.
    r, g, b = hsv_to_rgb(u, 1.0, 0.85)
    base = np.stack([r, g, b, np.ones_like(r)], -1)
    # Roughness: vertical gradient, smooth (top) to rough (bottom); the specular
    # highlights sharpen towards the top.
    rough = (0.05 + 0.9 * v)[..., None]
    # Normal: egg-carton bumps as a tangent-space normal map (flat = 0.5,0.5,1.0).
    freq, amp = 4.0, 0.3
    nx = amp * np.sin(2 * np.pi * freq * u) * np.cos(2 * np.pi * freq * v)
    ny = amp * np.cos(2 * np.pi * freq * u) * np.sin(2 * np.pi * freq * v)
    nz = np.sqrt(np.clip(1 - nx**2 - ny**2, 0, 1))
    normal = np.stack([nx * 0.5 + 0.5, ny * 0.5 + 0.5, nz * 0.5 + 0.5], -1)

    # Multi-part EXR: one part (subimage) per pass, matching how Blender writes
    # multi-layer render results. A single-part multi-channel EXR is not used
    # because the passes are then mapped positionally and get scrambled.
    parts = [(["ViewLayer.BaseColor.R", "ViewLayer.BaseColor.G", "ViewLayer.BaseColor.B",
               "ViewLayer.BaseColor.A"], "ViewLayer.BaseColor", base),
             (["ViewLayer.Roughness.V"], "ViewLayer.Roughness", rough),
             (["ViewLayer.Normal.X", "ViewLayer.Normal.Y", "ViewLayer.Normal.Z"],
              "ViewLayer.Normal", normal)]
    specs = []
    for chans, name, _ in parts:
        spec = oiio.ImageSpec(RES, RES, len(chans), "float")
        spec.channelnames = chans
        spec.attribute("name", name)
        specs.append(spec)
    os.makedirs(os.path.dirname(TEX), exist_ok=True)
    out = oiio.ImageOutput.create(TEX)
    out.open(TEX, tuple(specs))  # declare the multi-part file
    out.write_image(parts[0][2].astype(np.float32))
    for i in range(1, len(parts)):
        out.open(TEX, specs[i], "AppendSubimage")
        out.write_image(parts[i][2].astype(np.float32))
    out.close()


def build_scene():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.samples = 32
    scene.cycles.use_denoising = False
    scene.render.resolution_x = 128
    scene.render.resolution_y = 128
    scene.render.resolution_percentage = 100
    scene.view_settings.view_transform = 'Standard'

    # Dim world fill so shadowed bump sides are not fully black.
    world = bpy.data.worlds.new("World")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = (0.15, 0.15, 0.15, 1.0)

    # Plane facing the camera, filling the frame.
    bpy.ops.mesh.primitive_plane_add(size=2.6)
    plane = bpy.context.active_object
    plane.rotation_euler = (math.radians(90), 0, math.radians(90))

    cam_data = bpy.data.cameras.new("Camera")
    cam_data.lens = 35
    cam = bpy.data.objects.new("Camera", cam_data)
    cam.location = (2.958, 0, 0)
    cam.rotation_euler = (math.radians(90), 0, math.radians(90))
    bpy.context.collection.objects.link(cam)
    scene.camera = cam

    # Point light in front of the plane: lights the whole face and its specular
    # reflection reveals the roughness gradient and the normal bumps.
    light_data = bpy.data.lights.new("Lamp", 'POINT')
    light_data.energy = 100.0
    light = bpy.data.objects.new("Lamp", light_data)
    light.location = (4.076, 1.005, 2.904)
    bpy.context.collection.objects.link(light)
    return plane


def make_material(plane, mode):
    img = bpy.data.images.load(TEX)
    # Acquiring a buffer builds the multi-layer pass catalog.
    _ = img.pixels[0]
    mat = bpy.data.materials.new("MultiLayer")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled")
    nt.links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])
    nmap = nt.nodes.new("ShaderNodeNormalMap")
    nmap.uv_map = "UVMap"
    nt.links.new(nmap.outputs["Normal"], bsdf.inputs["Normal"])
    tex = nt.nodes.new("ShaderNodeTexImage")
    tex.image = img
    # Rebuild the node declaration so the per-pass output sockets appear.
    bpy.context.view_layer.update()
    mat.node_tree.update_tag()
    bpy.context.view_layer.update()

    if mode == 'passes':
        nt.links.new(tex.outputs["BaseColor"], bsdf.inputs["Base Color"])
        nt.links.new(tex.outputs["Roughness"], bsdf.inputs["Roughness"])
        nt.links.new(tex.outputs["Normal"], nmap.inputs["Color"])
    else:  # bundle: extract each pass from the node's Bundle output.
        sep = nt.nodes.new("NodeSeparateBundle")
        for socket_type, name in (('RGBA', "BaseColor"), ('FLOAT', "Roughness"),
                                  ('VECTOR', "Normal")):
            sep.bundle_items.new(socket_type, name)
        nt.links.new(tex.outputs["Passes"], sep.inputs["Bundle"])
        nt.links.new(sep.outputs["BaseColor"], bsdf.inputs["Base Color"])
        nt.links.new(sep.outputs["Roughness"], bsdf.inputs["Roughness"])
        nt.links.new(sep.outputs["Normal"], nmap.inputs["Color"])

    plane.data.materials.clear()
    plane.data.materials.append(mat)


def main():
    write_exr()
    for mode in ('passes', 'bundle'):
        plane = build_scene()
        make_material(plane, mode)
        bpy.ops.wm.save_as_mainfile(
            filepath=os.path.join(DIR, f"image_multilayer_{mode}.blend"), relative_remap=True)


main()
