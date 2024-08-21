/* SPDX-License-Identifier: GPL-2.0-or-later */
 
#include "gpu_shader_create_info.hh"



GPU_SHADER_INTERFACE_INFO(npr_prepass_iface, "")
    .flat(Type::UINT, "id")
    // .smooth(Type::VEC3, "normal")
    // .smooth(Type::VEC4, "tangent")
    ;
 
GPU_SHADER_CREATE_INFO(npr_prepass_mesh)
    .vertex_source("npr_prepass_vert.glsl")
    .vertex_in(0, Type::VEC3, "pos")
    // .vertex_in(1, Type::VEC3, "nor")
    // .vertex_in(2, Type::VEC4, "tan")
    .vertex_out(npr_prepass_iface)
    .fragment_source("npr_prepass_frag.glsl")
    // .fragment_out(0, Type::UINT, "out_id")
    // .fragment_out(1, Type::VEC3, "out_normal")
    // .fragment_out(2, Type::VEC4, "out_tangent")
    .additional_info("draw_modelmat_new", "draw_view", "draw_resource_handle_new")
    .do_static_compilation(true);

GPU_SHADER_CREATE_INFO(npr_deferred_pass)
    .fragment_source("npr_deferred_frag.glsl")
    // .sampler(0, ImageType::UINT_2D, "id_tx")
    // .sampler(1, ImageType::FLOAT_2D, "normal_tx")
    .sampler(0, ImageType::DEPTH_2D, "depth_tx")
    .sampler(1, ImageType::FLOAT_2D, "contour_tx")
    .fragment_out(0, Type::VEC4, "out_color")
    .fragment_out(1, Type::VEC4, "out_line_color")
    .fragment_out(2, Type::FLOAT, "out_line_data")
    .additional_info("draw_fullscreen", "draw_view")
    .do_static_compilation(true);
