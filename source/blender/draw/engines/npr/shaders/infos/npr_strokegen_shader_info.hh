/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "bnpr_defines.hh"
#include "gpu_shader_create_info.hh"

/* -------------------------------------------------------------------- */
/** \name Tutorial
 * \{ */

/* For details, see "gpu_shader_create_info.hh" */
// GPU_SHADER_CREATE_INFO(bnpr_strokegen_test)
// .do_static_compilation(true)

/* -------------------------------------------------------------------- */
/* Name of other infos to recursively merge with this one.
 * No data slot must overlap otherwise we throw an error. */
// .additional_info("eevee_shared", "draw_view", "draw_view_culling")

/* -------------------------------------------------------------------- */
/* Macros */
// .define("DOF_BOKEH_TEXTURE", "false")

/* -------------------------------------------------------------------- */
/** Resources bindings points
// .uniform_buf(6, "DepthOfFieldData", "dof_buf")
// .storage_buf(0, Qualifier::READ_WRITE, "LightCullingData", "light_cull_buf")
// .storage_buf(1, Qualifier::READ, "LightData", "in_light_buf[]")
// .storage_buf(2, Qualifier::WRITE, "LightData", "out_light_buf[]")
// .sampler(0, ImageType::FLOAT_2D, "downsample_tx")
// .image(0, GPU_RGBA16F, Qualifier::READ_WRITE, ImageType::FLOAT_2D, "inout_color_lod0_img")
/*          eGPUTextureFormat                     ImageType  */

/* -------------------------------------------------------------------- */
/** Comptue shader
// .local_group_size(CULLING_SELECT_GROUP_SIZE) /* <== from "bnpr_defines.hh" */
// .compute_source("eevee_light_culling_select_comp.glsl");

/* -------------------------------------------------------------------- */
// Vertex & Fragment shader
// .vertex_in(0, Type::VEC3, "pos")
// .vertex_in(1, Type::VEC3, "nor")
// .vertex_source("eevee_geom_mesh_vert.glsl")
//
// .fragment_out(0, Type::VEC4, "out_radiance", DualBlend::SRC_0)
// .fragment_out(1, Type::VEC4, "out_transmittance", DualBlend::SRC_1)
// .fragment_source("eevee_surf_forward_frag.glsl")
//
/* In order to use .vertex_out for vs output,
 * we firstly need to define an interface:
// GPU_SHADER_INTERFACE_INFO(interface_info_name, "interp")
// .smooth(Type::VEC3, "a") /* smooth: conventional interpolation for fragments */
// .flat(Type::VEC3, "b"); /* flat: no interpolation, instead, use attribute from a "provoking
// vertex" */ .noperspective(Type::VEC3, "c") /* interpolation, without perspective correction */
/* Then use the interface to declare a .vertex_out: */
// .vertex_out(interface_info_name)

/** \} */

/* -------------------------------------------------------------------- */
/** \shared shader infos
 * \{ */

/** GPU Scan / Segscan --------------------
 */
GPU_SHADER_CREATE_INFO(bnpr_scan_uint_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uint")
    .define("SCAN_OP", "u32_add")
    .define("SCAN_ZERO_VAL", "0u")
    .define("SCAN_FUNCTION_TAG", "_u32_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_uvec2_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uvec2")
    .define("SCAN_OP", "uvec2_add")
    .define("SCAN_ZERO_VAL", "uvec2(0u, 0u)")
    .define("SCAN_FUNCTION_TAG", "uvec2_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_uvec3_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uvec3")
    .define("SCAN_OP", "uvec3_add")
    .define("SCAN_ZERO_VAL", "uvec3(0u, 0u, 0u)")
    .define("SCAN_FUNCTION_TAG", "uvec3_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_uvec4_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uvec4")
    .define("SCAN_OP", "uvec4_add")
    .define("SCAN_ZERO_VAL", "uvec4(0u, 0u, 0u, 0u)")
    .define("SCAN_FUNCTION_TAG", "uvec4_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_float_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "float")
    .define("SCAN_OP", "f32_add")
    .define("SCAN_ZERO_VAL", ".0f")
    .define("SCAN_FUNCTION_TAG", "_f32_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_test_inputs).additional_info("bnpr_scan_uvec3_add");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_inputs)
    .additional_info("bnpr_scan_uvec3_add")
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanTest");

/** GPU 1D Segmented & Looped Convolution --------------------
 */
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test)
    .typedef_source("bnpr_shader_shared.hh")
    .define("LOOPCONV1D_TAG", "build_patch")
    .define("LOOPCONV1D_MAX_RADIUS", NPR_SEGLOOPCONV1D_CONV_RADIUS_STR)
    .define("DATA_TYPE_LOOPCONV1D", "uint");

GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test_build_patch)
    /* input functions related to buffer load/store can only be defined in a .glsl file
     * see npr_strokegen_segloopconv1d_inputs_lib.glsl */
    .define("FUNC_DEVICE_STORE_LOOPCONV1D_PATCH_ID", "func_device_store_loopconv1d_patch_id");

GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test_convolution)
    /* input functions related to buffer load/store can only be defined in a .glsl file
     * see npr_strokegen_segloopconv1d_inputs_lib.glsl */
    .define("FUNC_DEVICE_LOAD_LOOPCONV1D_PATCH_ID", "func_device_load_loopconv1d_patch_id")
    .define("FUNC_DEVICE_LOAD_LOOPCONV1D_DATA", "func_device_load_loopconv1d_data");

/** GPU Compaction --------------------
 */
GPU_SHADER_CREATE_INFO(npr_compaction_off)
    .define("COMPACTION_LIB_EXCLUDE", "1");

/** GPU List Ranking --------------------
 */

/* -------------------------------------------------------------------- */
/** \Geometry extraction from GPUBatch(es)
 * \{ */
/* Collect & Transform Contour Edges */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_boostrap)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE_BOOSTRAP_GEOM_EXTRACT", "1")
    .define("DECODE_IBO_EXCLUDE",       "1") /* exclude ibo decode shader */
    .additional_info("npr_compaction_off")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .local_group_size(32)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_extract)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE__GEOM_EXTRACT", "1")
    .define("GLOBAL_COUNTER", "ssbo_bnpr_mesh_pool_counters_.num_contour_edges")
    .define("DECODE_IBO_EXCLUDE",       "1") /* exclude ibo decode shader */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("CP_TAG", "contour_edge")
    .storage_buf(0, Qualifier::READ, "uint", "buf_ibo[]")
    .storage_buf(1, Qualifier::READ, "float", "buf_vbo[]") 
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "buf_strokegen_mesh_pool[]")
    .storage_buf(3, Qualifier::READ, "ObjectMatrices", "drw_matrix_buf[]")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_edge_to_contour_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_wedge_flooding_pointers_out_[]")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_ib_fmt_u16")
    .push_constant(Type::INT, "pcs_num_edges")
    .push_constant(Type::INT, "pcs_num_ib_offset")
    .push_constant(Type::INT, "pcs_rsc_handle")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")
    .push_constant(Type::INT, "pcs_dbg_wedge_flooding_")
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_extract_ibo_16bits)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_extract")
    .define("_KERNEL_MULTICOMPILE__INDEX_BUFFER_16BIT", "1");

GPU_SHADER_CREATE_INFO(strokegen_fill_dispatch_args_per_contour_edge)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_EDGE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(2, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_contour_edge_")
    .push_constant(Type::INT, "pc_per_contour_edge_dispatch_group_size_") /* group size of dispatched kernel */
    .push_constant(Type::INT, "pc_dispatch_for_all_edges_") 
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_extract_calc_contour_edge_raster_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("_KERNEL_MULTICOMPILE_CALC_CONTOUR_EDGE_RASTER_DATA", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "buf_strokegen_mesh_pool[]")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_edge_to_contour_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_") 
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

/* Collect Mesh Verts */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_collect_verts)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE__COMPACT_VBO", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_meshbatch_ibo_[]")
    .storage_buf(1, Qualifier::READ, "SSBO_Data_PosNor", "ssbo_meshbatch_vbo_[]") /* encoded posnor vbo */
    .storage_buf(2, Qualifier::WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(3, Qualifier::READ, "ObjectMatrices", "drw_matrix_buf[]")
    .storage_buf(4, Qualifier::WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(5, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_rsc_handle_")
    .push_constant(Type::INT, "pcs_meshbatch_num_verts_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

/* Collect Mesh Edge Adjacency */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_collect_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("_KERNEL_MULTICOMPILE__COMPACT_EDGE_ADJ_IBO", "1")
    .define("DECODE_IBO_INCLUDE", "1")

    .storage_buf(0, Qualifier::READ, "uint", "IBO_BUF[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_edge_to_vert_[]") /* encoded posnor vbo */
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_extract_collect_edges_16bits)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_extract_collect_edges")
    .define("_KERNEL_MULTICOMPILE__INDEX_BUFFER_16BIT", "1");

/* Draw Mesh Contour Edges */
GPU_SHADER_CREATE_INFO(bnpr_geom_fill_draw_args_contour_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h") /*DrawCommand*/
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE",                   "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS",  "1")
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_mesh_pool_draw_args_")
    .local_group_size(GROUP_SIZE_FILL_ARGS)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_INTERFACE_INFO(bnpr_v2f_geom_draw_contour_edges, "")
    /* .flat(Type::UINT, "id") */
    .smooth(Type::VEC4, "color") 
    .smooth(Type::VEC3, "normal")
    .smooth(Type::VEC4, "tangent");

GPU_SHADER_CREATE_INFO(bnpr_geom_draw_contour_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .additional_info("draw_modelmat_new", "draw_view", "draw_resource_handle_new")
    
    .storage_buf(0, Qualifier::READ, "uint", "buf_strokegen_mesh_pool[]")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_contour_edge_rank_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_contour_edge_list_len_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_contour_edge_list_head_[]")
    .push_constant(Type::VEC2, "pcs_screen_size_inv_")

    .vertex_source("npr_strokegen_mesh_contour_vert.glsl")
    .vertex_in(0, Type::VEC3, "pos")
    .vertex_in(1, Type::VEC3, "nor")
    .vertex_in(2, Type::VEC4, "tan")
    .vertex_out(bnpr_v2f_geom_draw_contour_edges)
    .fragment_source("npr_strokegen_mesh_contour_frag.glsl")
    .fragment_out(0, Type::VEC4, "out_col")
    .fragment_out(1, Type::VEC3, "out_normal")
    .fragment_out(2, Type::VEC4, "out_tangent");

/* Collect Contour Pixels */
GPU_SHADER_CREATE_INFO(bnpr_compress_contour_pixels)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__COMPRESS_CONTOUR_IMG", "1")

    .image(0, GPU_RGBA32F, Qualifier::READ, ImageType::FLOAT_2D, "tex2d_contour_image_")
    .image(1, GPU_R32UI, Qualifier::READ_WRITE, ImageType::UINT_2D, "tex2d_contour_pix_marks_")
    .push_constant(Type::VEC2, "pcs_screen_size_")

    .local_group_size(GROUP_SIZE_X_CONTOUR_PIXEL_COMPRESS, GROUP_SIZE_Y_CONTOUR_PIXEL_COMPRESS)
    .compute_source("npr_strokegen_pixel_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_compress_contour_pixels_dbg)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .additional_info("bnpr_compress_contour_pixels")
    .image(2, GPU_R32UI, Qualifier::WRITE, ImageType::UINT_2D, "tex2d_contour_pix_marks_dbg_")
    .define("_KERNEL_MULTICOMPILE__COMPRESS_CONTOUR_IMG_VALIDATE", "1"); 
/** \} */



/* -------------------------------------------------------------------- */
/** \GPU Meshing
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_verts)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__VERT_MERGE", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_vert_spatial_map_headers_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_vert_spatial_map_payloads_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_vert_merged_id_[]")
    .storage_buf(3, Qualifier::READ, "float", "ssbo_vbo_full_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .push_constant(Type::INT, "pcs_hash_map_size_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_gpu_meshing_compute.glsl");
;
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_verts_bootstrap)     
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_merge_verts")
    .define("_KERNEL_MULTICOMPILE__VERT_MERGE_BOOTSTRAP", "1");
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_verts_spatial_hashing)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_merge_verts")
    .define("_KERNEL_MULTICOMPILE__VERT_MERGE_BUILD_HASHMAP", "1");
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_verts_deduplicate)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_merge_verts")
    .define("_KERNEL_MULTICOMPILE__VERT_MERGE_DEDUPLICATE", "1");


/*
 * uint pcs_hash_map_size_
 * int pcs_edge_count_              note: all edges, NOT just contours 
 *
 * uint ssbo_vert_merged_id_[]
 * uint ssbo_edge_index_map_headers_[]      note: reuse vertex hash buffer
 * uint ssbo_edge_spatial_map_payloads_[]   note: reuse vertex payload buffer
 * uint ssbo_edge_to_vert_[]
 * uint ssbo_edge_to_edges_[]
*/
GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_adjecency)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_vert_merged_id_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_edge_index_map_headers_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_edge_spatial_map_payloads_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(6, Qualifier::WRITE, "uint", "ssbo_wedge_flooding_pointers_out_[]")
    .storage_buf(7, Qualifier::READ, "float", "ssbo_vbo_full_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_hash_map_size_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_gpu_meshing_compute.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_edges_bootstrap)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_adjecency")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_BOOTSTRAP", "1");
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_edges_hashing)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_adjecency")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_HASHING", "1");
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_edges_fill)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_adjecency")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_FIND_ADJ", "1");


GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("GLOBAL_COUNTER", "ssbo_bnpr_mesh_pool_counters_.num_filtered_edges")
    .define("CP_TAG", "filtered_edge")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_wedge_flooding_pointers_in_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_wedge_flooding_pointers_out_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(3, Qualifier::WRITE, "uint", "ssbo_filtered_edge_to_edge_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_iter)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_last_iter)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER", "1")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER", "1");


/*
 * uint pcs_vert_id_offset_;
 * uint pcs_vert_count_;
 * uint ssbo_wedge_flooding_pointers_in_[];
 * uint ssbo_edge_to_edges_[];
 * uint ssbo_edge_to_vert_[];
 * uint ssbo_filtered_vert_to_vert_[];
 */
GPU_SHADER_CREATE_INFO(bnpr_meshing_compact_filtered_verts)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE__COMPACT_FILTERED_VERTS", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("GLOBAL_COUNTER", "ssbo_bnpr_mesh_pool_counters_.num_filtered_verts")
    .define("CP_TAG", "filtered_vert")

    .storage_buf(0, Qualifier::WRITE, "uint", "ssbo_filtered_vert_to_vert_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_wedge_flooding_pointers_in_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_vert_merged_id_[]")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");






/* In: 
 * int pc_meshing_dispatch_group_size_
 * ssbo_bnpr_mesh_pool_counters_
 * ssbo_indirect_dispatch_args_per_filtered_edge_
 * ssbo_indirect_dispatch_args_per_filtered_vert_
*/
GPU_SHADER_CREATE_INFO(strokegen_meshing_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_filtered_edge_")
    .storage_buf(2, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_filtered_vert_")
    .push_constant(Type::INT, "pc_meshing_dispatch_group_size_") 

    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_fill_dispatch_args_per_filtered_edge) 
    .do_static_compilation(true)
    .additional_info("strokegen_meshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_EDGE", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_fill_dispatch_args_per_filtered_vert)
    .do_static_compilation(true)
    .additional_info("strokegen_meshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__MESHING__PER_FILTERED_VERT", "1");


/*
 * uint ssbo_bnpr_mesh_pool_counters_.num_filtered_edges/verts
 * uint ssbo_edge_to_edges_[];
 * uint ssbo_vert_to_edge_list_header_[]; 
 * uint ssbo_filtered_edge_to_edge_[];
 * uint ssbo_filtered_vert_to_vert_[];
 * uint ssbo_vert_quadric_data_in_[];  <- (reuse) ssbo_bnpr_mesh_pool_
 * uint ssbo_vert_quadric_data_out_[]; <- (reuse) ssbo_mesh_buffer_reuse_0_
 * uint ssbo_edge_quadric_data_[];     <- (reuse) ssbo_vert_merged_id_
 * float ssbo_vbo_full_[]; rw
 * float ssbo_filtered_normal_vert_[];  <- (reuse) ssbo_mesh_buffer_reuse_1_
 * float ssbo_filtered_normal_edge_[];  <- (reuse) ssbo_vert_quadric_data_in_ at iter#0
 * ubo_view_matrices_
 */
GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING_", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("COMPACTION_LIB_EXCLUDE") /* Remove compaction code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("QUADRICS_FILTERING_INCLUDE", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_filtered_edge_to_edge_[]")
    .storage_buf(5, Qualifier::READ, "uint", "ssbo_filtered_vert_to_vert_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "float", "ssbo_filtered_normal_vert_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "float", "ssbo_filtered_normal_edge_in_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "float", "ssbo_filtered_normal_edge_out_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_edge_quadric_data_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "uint", "ssbo_vert_quadric_data_in_[]")
    .storage_buf(12, Qualifier::READ_WRITE, "uint", "ssbo_vert_quadric_data_out_[]")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_") 
    .push_constant(Type::INT, "pcs_filtered_quadric_type_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_edge_id_offset_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_vert_id_offset_")
    .push_constant(Type::FLOAT, "pcs_quadric_deviation_")
    .push_constant(Type::FLOAT, "pcs_geodist_deviation_")
    .push_constant(Type::FLOAT, "pcs_positiion_regularization_scale_")
    .push_constant(Type::INT, "pcs_test_bilateral_filtering_")
    
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_edge_normal_)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_mesh_filtering_")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__EDGE_NORMAL", "1")
    .push_constant(Type::INT, "pcs_first_iter_")
    .push_constant(Type::INT, "pcs_use_normal_filtering_");

GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_edge_quadric_)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_mesh_filtering_")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__EDGE_QUADRIC", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_init_vert_quadric_)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_mesh_filtering_")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC", "1")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC_INIT", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_diffuse_vert_quadric_)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_mesh_filtering_")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__VERT_QUADRIC", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_mesh_filtering_move_verts_)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_mesh_filtering_")
    .define("_KERNEL_MULTICOMPILE__MESH_FILTERING__MOVE_VERTS", "1");

/** \} */



/* -------------------------------------------------------------------- */
/** \Bare minumum compute shader test
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_strokegen_test_xxx)
    .do_static_compilation(true)
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "buf_test[]")
    .storage_buf(1, Qualifier::READ, "uint", "buf_ibo[]")
    .local_group_size(GROUP_SIZE_STROKEGEN_TEST) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_test_comp.glsl");
/** \} */

/* -------------------------------------------------------------------- */
/** \GPU (Segmented)Scan Test
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_scan_test_upsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_scan_test_inputs")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_UPSWEEP", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_in_scan_data_buf_[]")
    .storage_buf(1, Qualifier::READ_WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_out_scan_data_buf_[]")
    .storage_buf(2, Qualifier::WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_scan_test_aggregate)
    .do_static_compilation(true)
    .additional_info("bnpr_scan_test_inputs")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_AGGREGATE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_AGGRG) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_scan_test_dwsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_scan_test_inputs")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_DWSWEEP", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_out_scan_data_buf_[]")
    .storage_buf(1, Qualifier::READ, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_upsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_UPSWEEP", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_in_scan_data_buf_[]")
    .storage_buf(1, Qualifier::READ_WRITE, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_out_scan_data_buf_[]")
    .storage_buf(2, Qualifier::WRITE, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_aggregate)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_AGGREGATE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_AGGRG) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_dwsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_DWSWEEP", "1")
    .storage_buf(0, Qualifier::READ_WRITE, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_out_scan_data_buf_[]")
    .storage_buf(1, Qualifier::READ, BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");
/** \} */

/* -------------------------------------------------------------------- */
/** \GPU 1D Segmented Looped Convolution
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_segloopconv1D_test_build_patch)
    .do_static_compilation(true)
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_segloopconv1d_patch_table_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_debug_segloopconv1d_data_[]")
    .uniform_buf(0, "UBData_SegLoopConv1D", "ubo_segloopconv1d_")
    .additional_info("npr_segloopconv1D_test")
    .additional_info("npr_segloopconv1D_test_build_patch")
    .local_group_size(GROUP_SIZE_SEGLOOPCONV1D_TEST)
    .compute_source("npr_strokegen_segloopconv1d_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_segloopconv1D_test_convolution)
    .do_static_compilation(true)
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION", "1")
    .storage_buf(0, Qualifier::READ, "uint", "ssbo_segloopconv1d_patch_table_[]")
    .storage_buf(1, Qualifier::READ_WRITE, NPR_SEGLOOPCONV1D_TEST_DATA_TYPE_STR, "ssbo_in_segloopconv1d_data_[]")
    .storage_buf(2, Qualifier::READ_WRITE, NPR_SEGLOOPCONV1D_TEST_DATA_TYPE_STR, "ssbo_out_segloopconv1d_data_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_debug_segloopconv1d_data_[]")
    .uniform_buf(0, "UBData_SegLoopConv1D", "ubo_segloopconv1d_")
    .additional_info("npr_segloopconv1D_test")
    .additional_info("npr_segloopconv1D_test_convolution")
    .local_group_size(GROUP_SIZE_SEGLOOPCONV1D_TEST)
    .compute_source("npr_strokegen_segloopconv1d_test_comp.glsl");
/** \} */


/* -------------------------------------------------------------------- */
/** \GPU List Ranking Test
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_list_ranking_setup_input_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SETUP", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_links_in_[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_list_ranking_links_out_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_custom_") /* 0:=testing pass, 1:=custom pass */
    
    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.h") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__LIST_RANKING_ANCHORS", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_list_ranking_splice_counters_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_list_ranking_indirect_dispatch_args_per_anchor")
    .storage_buf(3, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_list_ranking_indirect_dispatch_args_per_spliced")
    .storage_buf(4, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .push_constant(Type::INT, "pc_listranking_custom_") /* 0:=testing pass, 1:=custom pass */
    .push_constant(Type::INT, "pc_listranking_indirect_arg_granularity_") /* 0 := anchors, 1:= spliced nodes */
    .push_constant(Type::INT, "pc_listranking_counter_buffer_slot_id_") /* which counter to use */
    .push_constant(Type::INT, "pc_listranking_dispatch_group_size_") /* group size of dispatched kernel */
    
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_tagging)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_TAGGING", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_tags_in_[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_list_ranking_tags_out_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_to_node_in_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_anchor_to_node_out_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(6, Qualifier::WRITE, "uint", "ssbo_list_ranking_splice_counters_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(8, Qualifier::WRITE, "uint", "ssbo_list_ranking_addressing_counters_[]") 
    .storage_buf(9, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")
    .push_constant(Type::INT, "pc_listranking_tagging_iter_")
    .push_constant(Type::INT, "pc_num_splice_iters_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_compact_anchors)
    .do_static_compilation(true)
    .typedef_source("draw_shader_shared.h")
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_COMPACT_ANCHORS", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_tags_in_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_to_node_in_[]")
    .storage_buf(4, Qualifier::WRITE, "uint", "ssbo_list_ranking_anchor_to_node_out_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_list_ranking_spliced_node_id_out_[]")
    .storage_buf(6, Qualifier::WRITE, "uint", "ssbo_list_ranking_anchor_to_next_anchor_[]")
    .storage_buf(7, Qualifier::WRITE, "uint", "ssbo_list_ranking_node_to_anchor_out_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_splice_counters_[]")
    .storage_buf(10, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")
    .push_constant(Type::INT, "pc_num_splice_iters_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_splice_out_nodes)
    .do_static_compilation(true)
    .typedef_source("draw_shader_shared.h")
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SPLICE_NODES", "1")
    
    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_tags_in_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_to_node_in_[]")
    .storage_buf(4, Qualifier::WRITE, "uint", "ssbo_list_ranking_anchor_to_node_out_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_spliced_node_id_in_[]")
    .storage_buf(6, Qualifier::WRITE, "uint", "ssbo_list_ranking_anchor_to_next_anchor_[]")
    .storage_buf(7, Qualifier::WRITE, "uint", "ssbo_list_ranking_node_to_anchor_out_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_splice_counters_[]")
    .storage_buf(10, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_sublist_pointer_jumping)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_SUBLIST_POINTER_JUMPING", "1")
    .define("MAX_NUM_JUMPS_BNPR_LIST_RANK_TEST", "26")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_list_ranking_node_to_anchor_in_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_to_node_in_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_addressing_counters_[]") 
    .storage_buf(8, Qualifier::WRITE, "uint", "ssbo_list_ranking_serialized_topo_[]")
    .storage_buf(9, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")
    .push_constant(Type::INT, "pc_listranking_jumping_iter_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_looped_pointer_jumping)
    .do_static_compilation(true)
    .additional_info("strokegen_list_ranking_test_sublist_pointer_jumping")
    .define("_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD", "1")
    .push_constant(Type::INT, "pc_listranking_ranking_pass_with_broken_loops_");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_mark_loop_head_tail)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_BREAK_CIRCLES")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_to_node_in_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[]")
    .storage_buf(2, Qualifier::WRITE, "uint", "ssbo_list_ranking_per_anchor_sublist_jumping_info_out_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_list_ranking_node_to_anchor_in_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_list_ranking_anchor_counters_[]")
    .storage_buf(7, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");
    
GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_relinking)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_RELINKING", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_spliced_node_id_in_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_links_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_list_ranking_splice_counters_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_list_ranking_per_anchor_sublist_jumping_info_in_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_list_ranking_serialized_topo_[]")
    .storage_buf(6, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_relink_iter_")
    .push_constant(Type::INT, "pc_listranking_num_relink_iters_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_looped_relinking)
    .do_static_compilation(true)
    .additional_info("strokegen_list_ranking_test_relinking")
    .define("_KERNEL_MULTICOMPILE__LIST_RANKING_SUBLIST_POINTER_JUMPING__FIND_LOOP_HEAD", "1");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_output)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_OUTPUT", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_ranks_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_list_ranking_serialized_topo_[]")
    .storage_buf(2, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
    .storage_buf(3, Qualifier::WRITE, "uint", "ssbo_list_ranking_output_ranks_[]")
    .storage_buf(4, Qualifier::WRITE, "uint", "ssbo_list_ranking_output_list_len_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_list_ranking_output_list_addr_[]")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

/** \} */
