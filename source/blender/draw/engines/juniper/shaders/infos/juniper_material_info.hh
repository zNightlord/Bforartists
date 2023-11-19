/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "gpu_shader_create_info.hh"
#include "../../juniper_defines.hh"

GPU_SHADER_INTERFACE_INFO(juniper_surf_iface, "interp")
        .smooth(Type::VEC3, "P")
        .smooth(Type::VEC3, "N");

GPU_SHADER_CREATE_INFO(juniper_surface_mesh)
        .vertex_source("juniper_surf_test_vert.glsl")
        .vertex_in(0, Type::VEC3, "pos")
        .vertex_in(1, Type::VEC3, "nor")
        .vertex_out(juniper_surf_iface)
        .typedef_source("juniper_defines.hh")
        .typedef_source("juniper_shader_shared.hh")
        .additional_info("juniper_material_lib")
        .additional_info("draw_modelmat_new", "draw_resource_id_varying", "draw_view");

GPU_SHADER_CREATE_INFO(juniper_light_data)
        .push_constant(Type::IVEC4, "light_groups")
        .storage_buf(SH_BUF_LIGHTDATA_SLOT, Qualifier::READ, "LightData", "light_buf[]")
        .storage_buf(SH_BUF_LIGHTMETA_SLOT, Qualifier::READ, "LightMeta", "light_meta_buf");

GPU_SHADER_CREATE_INFO(juniper_material_lib)
        .sampler(SH_TEX_DEPTH_SLOT, ImageType::FLOAT_2D, "in_depth_tx")
        .sampler(SH_TEX_NORMAL_SLOT, ImageType::FLOAT_2D, "in_normal_tx");

GPU_SHADER_CREATE_INFO(juniper_surface_prepass)
        .define("JUNIPER_PREPASS")
        .fragment_source("juniper_surf_prepass_frag.glsl")
        .additional_info("juniper_surface_mesh", "juniper_light_data")
        .fragment_out(0, Type::VEC3, "out_normal");

GPU_SHADER_CREATE_INFO(juniper_surface_test)
        .define("JUNIPER_SHADING_PASS")
        .additional_info("juniper_surface_mesh")
        .fragment_source("juniper_surf_test_frag.glsl")
        .fragment_out(0, Type::VEC4, "out_color")
        .additional_info("juniper_light_data");

GPU_SHADER_CREATE_INFO(juniper_film_composite)
        .fragment_out(0, Type::VEC4, "out_color")
        .fragment_source("juniper_film_composite.glsl")
        .sampler(SH_TEX_NORMAL_SLOT, ImageType::FLOAT_2D, "normal_buffer")
        .sampler(SH_TEX_DEPTH_SLOT, ImageType::FLOAT_2D, "depth_buffer")
        .sampler(SH_TEX_COLOR_SLOT, ImageType::FLOAT_2D, "color_buffer")
        .typedef_source("juniper_shader_shared.hh")
        .additional_info("juniper_light_data", "draw_fullscreen", "draw_view")
        .do_static_compilation(true)
        .depth_write(DepthWrite::ANY);