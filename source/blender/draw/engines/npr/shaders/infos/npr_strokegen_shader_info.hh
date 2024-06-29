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

/* -------------------------------------------------------------------- */
/** \Geometry extraction from GPUBatch(es)
 * \{ */
/* Collect & Transform Contour Edges */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_boostrap)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
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
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__GEOM_EXTRACT", "1")
    .define("DECODE_IBO_EXCLUDE",       "1") /* exclude ibo decode shader */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("TOPO_DIAGONOSIS_INCLUDE", "1")
    .storage_buf(0, Qualifier::READ, "uint", "buf_ibo[]")
    .storage_buf(1, Qualifier::READ, "float", "buf_vbo[]") 
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_temp_data_[]")
    .storage_buf(3, Qualifier::READ, "ObjectMatrices", "drw_matrix_buf[]")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(5, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_contour_[]")
    .storage_buf(7, Qualifier::READ, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]") 
    .storage_buf(9, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(10, Qualifier::WRITE, "uint", "ssbo_face_to_vert_draw_depth_[]")
    /* debugging */
    .define("INCLUDE_DEBUG_LINE_CONFIG", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .storage_buf(11, Qualifier::READ, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(12, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(13, Qualifier::WRITE, "uint", "ssbo_dbg_lines_[]")
    .storage_buf(14, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    /* ---------------- */
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_") 
    .push_constant(Type::INT, "pcs_ib_fmt_u16")
    .push_constant(Type::INT, "pcs_num_edges")
    .push_constant(Type::INT, "pcs_num_ib_offset")
    .push_constant(Type::INT, "pcs_rsc_handle")
    .push_constant(Type::INT, "pcs_edge_visualize_mode_")
    .push_constant(Type::INT, "pcs_chain_interpo_contour_")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    .push_constant(Type::FLOAT, "pcs_dbg_geom_scale_")
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
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_EDGE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(2, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_contour_edge_")
    .push_constant(Type::INT, "pc_per_contour_edge_dispatch_group_size_") /* group size of dispatched kernel */
    .push_constant(Type::INT, "pc_dispatch_for_all_edges_") 
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_fill_dispatch_args_per_contour_vert)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_VERT", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_contour_vert_")
    .push_constant(Type::INT, "pc_per_contour_vert_dispatch_group_size_") /* group size of dispatched kernel */
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_fill_dispatch_args_per_contour_fragment)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_FRAG", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(2, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_bnpr_mesh_contour_frag_dispatch_args_")
    .push_constant(Type::INT, "pc_dispatch_for_all_frags_") 
    .push_constant(Type::INT, "pc_per_contour_frag_dispatch_group_size_") /* group size of dispatched kernel */
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_fill_dispatch_args_per_contour_2d_sample)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__PER_CONTOUR_2D_SAMPLE", "1")
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_bnpr_mesh_contour_2d_sample_dispatch_args_")
    .push_constant(Type::INT, "pc_per_contour_2d_sample_dispatch_group_size_") /* group size of dispatched kernel */
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_extract_mesh_contour_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("_KERNEL_MULTICOMPILE__EXTRACT_MESH_CONTOUR_DATA", "1")
    .define("VE_CIRCULATOR_INCLUDE", "1")
    .define("INCLUDE_VERTEX_CURV_MAX", "1")
    .define("USE_CONTOUR_TRANSFER_DATA_BUFFER", "1")
    .define("INCLUDE_LOAD_STORE_CONTOUR_RASTER_DATA", "1")
    
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_contour_temp_data_[]")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_contour_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]") 
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]") 
    .storage_buf(8, Qualifier::READ, "float", "ssbo_vbo_full_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_") 
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_contour_edge_transfer_data_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "uint", "ssbo_vcurv_max_[]")
    .storage_buf(12, Qualifier::READ_WRITE, "uint", "ssbo_contour_raster_data_[]")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_geom_extract_comp.glsl");


GPU_SHADER_CREATE_INFO(strokegen_build_contour_fragments)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("INCLUDE_LOAD_STORE_CONTOUR_RASTER_DATA", "1")
    
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_frag_to_contour_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_raster_data_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_frag_raster_data_[]")
    .storage_buf(4, Qualifier::WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_contour_segmentation_")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_tree_scan_input_contour_fragment_idmapping_[]")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .sampler(0, ImageType::FLOAT_2D, "tex_remeshed_surf_depth_")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_")
    .push_constant(Type::VEC2, "pcs_screen_size_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_geom_extract_comp.glsl"); 

GPU_SHADER_CREATE_INFO(strokegen_clear_frag_to_contour_idmapping)
    .do_static_compilation(true)
    .additional_info("strokegen_build_contour_fragments")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__CLEAR_BUFFER", "1");

GPU_SHADER_CREATE_INFO(strokegen_prep_segscan_frag_to_contour_idmapping)
    .do_static_compilation(true)
    .additional_info("strokegen_build_contour_fragments")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__SETUP_SEGSCAN", "1");

GPU_SHADER_CREATE_INFO(strokegen_finish_segscan_frag_to_contour_idmapping)
    .do_static_compilation(true)
    .additional_info("strokegen_build_contour_fragments")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__IDMAPPING__FINISH_SEGSCAN", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_frag_visibility_test)
    .do_static_compilation(true)
    .additional_info("strokegen_build_contour_fragments")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__VISIBILITY_TEST", "1")
    .push_constant(Type::FLOAT, "pcs_visibility_thresh_");

GPU_SHADER_CREATE_INFO(strokegen_visibility_split_contour_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("INCLUDE_LOAD_STORE_CONTOUR_RASTER_DATA", "1")
    
    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_frag_to_contour_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_raster_data_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_frag_raster_data_[]")
    .storage_buf(4, Qualifier::WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_contour_segmentation_")
#define NUM_SSBO_strokegen_visibility_split_contour_edges 5
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_geom_extract_comp.glsl"); 

GPU_SHADER_CREATE_INFO(strokegen_contour_frag_setup_visibility_segmentation)
    .do_static_compilation(true)
    .additional_info("strokegen_visibility_split_contour_edges")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__SETUP_SEGSCAN", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_visibility_split_contour_edges
    .storage_buf(SSBO_OFFSET + 0, Qualifier::WRITE, "uint", "ssbo_tree_scan_input_contour_visibility_split_0_[]")
    .storage_buf(SSBO_OFFSET + 1, Qualifier::WRITE, "uint", "ssbo_tree_scan_input_contour_visibility_split_1_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_contour_visibility_split_step_0)
    .do_static_compilation(true)
    .additional_info("strokegen_visibility_split_contour_edges")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_0", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_visibility_split_contour_edges
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "uint", "ssbo_contour_visibility_split_info_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_contour_visibility_split_step_1)
    .do_static_compilation(true)
    .additional_info("strokegen_visibility_split_contour_edges")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_1", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_visibility_split_contour_edges
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "uint", "ssbo_contour_visibility_split_info_[]")
    .storage_buf(SSBO_OFFSET + 1, Qualifier::READ_WRITE, "uint", "ssbo_frag_seg_head_to_visibility_split_contour_[]")
    .storage_buf(SSBO_OFFSET + 2, Qualifier::READ, "uint", "ssbo_tree_scan_output_contour_visibility_split_0_[]")
    .storage_buf(SSBO_OFFSET + 3, Qualifier::READ, "uint", "ssbo_tree_scan_output_contour_visibility_split_1_[]")
    .storage_buf(SSBO_OFFSET + 4, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_contour_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_contour_visibility_split_step_2)
    .do_static_compilation(true)
    .additional_info("strokegen_visibility_split_contour_edges")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2_3", "1")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2", "1")
    .define("USE_CONTOUR_TRANSFER_DATA_BUFFER", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_visibility_split_contour_edges
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "uint", "ssbo_contour_visibility_split_info_[]")
    .storage_buf(SSBO_OFFSET + 1, Qualifier::READ_WRITE, "uint", "ssbo_frag_seg_head_to_visibility_split_contour_[]")
    .storage_buf(SSBO_OFFSET + 2, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(SSBO_OFFSET + 3, Qualifier::READ_WRITE, "uint", "ssbo_contour_edge_transfer_data_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_contour_visibility_split_step_3)
    .do_static_compilation(true)
    .additional_info("strokegen_visibility_split_contour_edges")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_2_3", "1")
    .define("_KERNEL_MULTICOMPILE__PROCESS_CONTOUR_FRAGMENTS__SPLIT_VISIBILITY__GENERATE_NEW_CONTOURS_STEP_3", "1")
    .define("USE_CONTOUR_TRANSFER_DATA_BUFFER", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_visibility_split_contour_edges
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "uint", "ssbo_contour_visibility_split_info_[]")
    .storage_buf(SSBO_OFFSET + 1, Qualifier::READ_WRITE, "uint", "ssbo_frag_seg_head_to_visibility_split_contour_[]")
    .storage_buf(SSBO_OFFSET + 2, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(SSBO_OFFSET + 3, Qualifier::READ_WRITE, "uint", "ssbo_contour_edge_transfer_data_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_fill_cotour_edge_ranking_inputs)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__FILL_LIST_RANKING_INPUTS_FOR_CONTOUR_EDGES", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_serialize_contour_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    .define("USE_CONTOUR_TRANSFER_DATA_BUFFER", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_contour_edge_rank_in_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_contour_edge_list_len_in_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_list_ranking_list_head_info_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_edge_transfer_data_[]")
    .storage_buf(4, Qualifier::WRITE, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(7, Qualifier::WRITE, "uint", "ssbo_contour_snake_vpos_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(11, Qualifier::WRITE, "UBData_SegLoopConv1D", "ssbo_segloopconv1d_info_")
    .storage_buf(12, Qualifier::WRITE, "uint", "ssbo_in_segloopconv1d_data_[]")
    .storage_buf(13, Qualifier::READ, "uint", "ssbo_list_ranking_addressing_counters_[]")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_serialize_contour_edges_pass_0)
    .do_static_compilation(true)
    .additional_info("strokegen_serialize_contour_edges")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_SERIALIZATION__PASS_0", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_segmentation)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    .define("INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "scan_data_buf_0_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "scan_data_buf_1_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "scan_output_buf_0_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "scan_output_buf_1_[]")
    .storage_buf(8, Qualifier::WRITE, "uint", "ssbo_contour_snake_seg_rank_[]")
    .storage_buf(9, Qualifier::WRITE, "uint", "ssbo_contour_snake_seg_len_[]")
    .storage_buf(10, Qualifier::WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_contour_segmentation_")
    .storage_buf(11, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_setup_contour_segmentation)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_segmentation")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__SETUP", "1");

GPU_SHADER_CREATE_INFO(strokegen_finish_contour_segmentation)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_segmentation")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_SEGMENTATION__FINISH", "1");


GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    .define("INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_resample_raster_data_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_tree_scan_input_2d_resampler_accumulate_curvlen_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_arc_len_param_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_tree_scan_input_2d_resampler_alloc_samples_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_start_sample_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_tree_scan_input_2d_resample_contour_idmapping_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_2d_sample_to_contour_[]")
    .storage_buf(8, Qualifier::READ, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(9, Qualifier::READ, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(10, Qualifier::READ, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
    .storage_buf(12, Qualifier::READ, "uint", "ssbo_contour_snake_vpos_[]")
    .storage_buf(13, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(14, Qualifier::READ_WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_2d_resampler_")
#define NUM_SSBO_strokegen_contour_2d_resample 15

    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")

    .push_constant(Type::VEC2, "pcs_screen_size_")
    .push_constant(Type::FLOAT, "pcs_sample_rate_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_arclen_parameterization)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__PREP_ARCLEN_PARAM", "1")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_"); 
    
GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_alloc_samples)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__ALLOC_SAMPLES", "1")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_"); 

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_alloc_samples_finish)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__ALLOC_SAMPLES_FINISH", "1")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_"); 

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_idmapping_clear_buffer)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__CLEAR_BUFFER", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_idmapping_setup_segscan)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__SETUP_SEGSCAN", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_idmapping_finish)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_resample")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__IDMAPPING__FINISH", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_sample_eval)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    .define("INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD", "1")
    .define("USE_CONTOUR_2D_SAMPLE_GEOMETRY_BUFFER", "1")
    .define("USE_CONTOUR_2D_SAMPLE_TOPOLOGY_BUFFER", "1")
    
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_arc_len_param_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_to_start_sample_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_2d_sample_to_contour_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(7, Qualifier::READ, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(8, Qualifier::READ, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(9, Qualifier::READ, "uint", "ssbo_contour_snake_seg_rank_[]")
    .storage_buf(10, Qualifier::READ, "uint", "ssbo_contour_snake_seg_len_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
#define NUM_SSBO_strokegen_contour_2d_sample_eval 12

    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_")

    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")

    .push_constant(Type::VEC2, "pcs_screen_size_")
    .push_constant(Type::FLOAT, "pcs_sample_rate_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_eval_position)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_sample_eval")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_POSITION", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_contour_2d_sample_eval
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_resample_raster_data_[]")
#undef SSBO_OFFSET
    ;

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_eval_topo_step_0)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_sample_eval")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__STEP_0", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_contour_2d_sample_eval
    .storage_buf(SSBO_OFFSET + 0, Qualifier::READ_WRITE, "UBData_SegLoopConv1D", "ssbo_segloopconv1d_info_")
#undef SSBO_OFFSET
    .push_constant(Type::INT, "pcs_segment_by_seg_");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_segmentation_prep_seg_tails)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_sample_eval")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SETMENTATION__PREP_SEGTAILS", "1")
    .push_constant(Type::INT, "pcs_segment_by_seg_");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_segmentation_setup_segscan)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_sample_eval")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__SETUP_SEGSCAN", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_contour_2d_sample_eval
    .storage_buf(SSBO_OFFSET + 0, Qualifier::WRITE, "uint", "ssbo_tree_scan_input_2d_sample_segmentation_0_[]")
    .storage_buf(SSBO_OFFSET + 1, Qualifier::WRITE, "uint", "ssbo_tree_scan_input_2d_sample_segmentation_1_[]")
    .storage_buf(SSBO_OFFSET + 2, Qualifier::READ_WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_2d_resampler_")
#undef SSBO_OFFSET
    .push_constant(Type::INT, "pcs_segment_by_seg_");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_eval_topo_remove_fake_corners)
    .do_static_compilation(true)
    .additional_info("strokegen_contour_2d_sample_eval")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__REMOVE_FAKE_CORNERS", "1");

GPU_SHADER_CREATE_INFO(strokegen_contour_2d_resample_eval_topo_finish)
    .do_static_compilation(true)

    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY", "1")
    .define("_KERNEL_MULTICOMPILE__CONTOUR_EDGES_2D_RESAMPLE__EVALUATE_TOPOLOGY__SEGMENTATION__FINISH_SEGSCAN", "1")

    .define("USE_CONTOUR_2D_SAMPLE_GEOMETRY_BUFFER", "1")
    .define("USE_CONTOUR_2D_SAMPLE_TOPOLOGY_BUFFER", "1")
    
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_tree_scan_input_2d_sample_segmentation_0_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_tree_scan_output_2d_sample_segmentation_0_[]")
    .storage_buf(5, Qualifier::READ, "uint", "ssbo_tree_scan_input_2d_sample_segmentation_1_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_tree_scan_output_2d_sample_segmentation_1_[]")
    
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_")

    .push_constant(Type::INT, "pcs_segment_by_seg_")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    .push_constant(Type::FLOAT, "pcs_sample_rate_")
    
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_calc_contour_2d_curve_render_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CALC_CONTOUR_2D_STROKE_RENDER_DATA", "1")
    .define("USE_CONTOUR_2D_SAMPLE_GEOMETRY_BUFFER", "1")
    .define("USE_CONTOUR_2D_SAMPLE_TOPOLOGY_BUFFER", "1")
    .define("USE_CONTOUR_2D_STROKE_MESH_BUFFER", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_stroke_mesh_pool_[]")

    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_")

    .push_constant(Type::VEC2, "pcs_screen_size_")
    .push_constant(Type::FLOAT, "pcs_stroke_width_")
    
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");


GPU_SHADER_CREATE_INFO(strokegen_calc_contour_edges_render_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__CALC_CONTOUR_EDGES_RENDER_DATA", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    .define("INCLUDE_CONTOUR_CURVE_TOPOLOGY_LOAD", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_contour_snake_vpos_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
    .storage_buf(5, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
#define NUM_SSBO_strokegen_calc_contour_edges_render_data 6
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) 
    .compute_source("npr_strokegen_contour_processing.glsl");

GPU_SHADER_CREATE_INFO(strokegen_calc_contour_edges_draw_data)
    .do_static_compilation(true)
    .additional_info("strokegen_calc_contour_edges_render_data")
    .define("_KERNEL_MULTICOMPILE__CALC_CONTOUR_EDGES_DRAW_DATA", "1")
#define SSBO_OFFSET NUM_SSBO_strokegen_calc_contour_edges_render_data
    .storage_buf(SSBO_OFFSET + 0, Qualifier::WRITE, "uint", "buf_strokegen_mesh_pool[]"); 
#undef SSBO_OFFSET

/* Collect Mesh Verts */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_collect_verts)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE__COPY_VBO", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_meshbatch_ibo_[]")
    .storage_buf(1, Qualifier::READ, "float", "ssbo_meshbatch_vbo_[]") /* encoded posnor vbo */
    .storage_buf(2, Qualifier::WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(3, Qualifier::READ, "ObjectMatrices", "drw_matrix_buf[]")
    .storage_buf(4, Qualifier::WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_prev_")
    .storage_buf(5, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_rsc_handle_")
    .push_constant(Type::INT, "pcs_meshbatch_num_verts_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

/* Collect Mesh Edge Adjacency */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract_collect_edges)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("_KERNEL_MULTICOMPILE__COPY_EDGE_ADJ_IBO", "1")
    .define("DECODE_IBO_INCLUDE", "1")

    .storage_buf(0, Qualifier::READ, "uint", "IBO_BUF[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_edge_to_vert_[]") /* encoded posnor vbo */
    .push_constant(Type::INT, "pcs_edge_count_")

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
    .typedef_source("draw_shader_shared.hh") /*DrawCommand*/
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE",                   "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS",  "1")
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS__CONTOUR_EDGES",  "1")
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_contour_mesh_draw_args_")
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
    .define("_KERNEL_MULTICOMPILE_DRAW_CONTOUR__EDGES", "1")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    
    .storage_buf(0, Qualifier::READ, "uint", "buf_strokegen_mesh_pool[]")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(4, Qualifier::READ, "uint", "ssbo_contour_snake_list_head_[]")
    .storage_buf(5, Qualifier::READ, "uint", "ssbo_contour_to_contour_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_contour_snake_seg_rank_[]")
    .storage_buf(7, Qualifier::READ, "uint", "ssbo_contour_snake_seg_len_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]")
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

/* Draw Mesh 2D-Resampled Contour  Curves */
GPU_SHADER_CREATE_INFO(bnpr_geom_fill_draw_args_contour_2d_samples)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /*DrawCommand*/
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE",                   "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS",  "1")
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS__CONTOUR_SAMPLES",  "1")
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_2d_sample_draw_args_")
    .local_group_size(GROUP_SIZE_FILL_ARGS)
    .compute_source("npr_strokegen_geom_extract_comp.glsl"); 

GPU_SHADER_INTERFACE_INFO(bnpr_v2f_geom_draw_contour_2d_samples, "")
    .smooth(Type::VEC4, "color") 
    .smooth(Type::VEC3, "normal")
    .smooth(Type::VEC4, "tangent");

GPU_SHADER_CREATE_INFO(bnpr_geom_draw_contour_2d_samples)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .additional_info("draw_modelmat_new", "draw_view", "draw_resource_handle_new")
    .define("_KERNEL_MULTICOMPILE_DRAW_CONTOUR__2D_SAMPLES", "1")
    .define("USE_CONTOUR_2D_STROKE_MESH_BUFFER", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_stroke_mesh_pool_[]")
    .push_constant(Type::VEC2, "pcs_screen_size_") 

    .vertex_source("npr_strokegen_mesh_contour_vert.glsl")
    .vertex_in(0, Type::VEC3, "pos")
    .vertex_in(1, Type::VEC3, "nor")
    .vertex_in(2, Type::VEC4, "tan")
    .vertex_out(bnpr_v2f_geom_draw_contour_edges)
    .fragment_source("npr_strokegen_mesh_contour_frag.glsl")
    .fragment_out(0, Type::VEC4, "out_col")
    .fragment_out(1, Type::VEC3, "out_normal")
    .fragment_out(2, Type::VEC4, "out_tangent");

/* Draw remeshed surface ----------------------------------------- */
GPU_SHADER_CREATE_INFO(bnpr_geom_fill_draw_args_remeshed_surface)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /*DrawCommand*/
    .additional_info("npr_compaction_off") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS",  "1")
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_contour_mesh_draw_args_")
    .local_group_size(GROUP_SIZE_FILL_ARGS)
    .compute_source("npr_strokegen_geom_extract_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_geom_fill_draw_args_remeshed_surface_depth)
    .additional_info("bnpr_geom_fill_draw_args_remeshed_surface")
    .do_static_compilation(true)
    .define("_KERNEL_MULTICOMPILE_FILL_DRAW_ARGS_DEPTH", "1"); 

GPU_SHADER_INTERFACE_INFO(bnpr_v2f_geom_draw_remeshed_surface_depth, "")
    .flat(Type::UINT, "id")
    .smooth(Type::VEC4, "color") 
    .smooth(Type::VEC3, "normal")
    .smooth(Type::VEC4, "tangent");

GPU_SHADER_CREATE_INFO(bnpr_geom_draw_remeshed_surface_depth)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .additional_info("draw_modelmat_new", "draw_view", "draw_resource_handle_new")
    .define("INCLUDE_CONTOUR_FLAGS_LOAD_STORE", "1")
    
    .storage_buf(0, Qualifier::READ, "uint", "ssbo_face_to_vert_draw_depth_[]")
    .storage_buf(1, Qualifier::READ, "float", "ssbo_vbo_full_[]") 
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")

    .vertex_source("npr_strokegen_remeshed_surface_vert.glsl")
    .vertex_in(0, Type::VEC3, "pos")
    .vertex_in(1, Type::VEC3, "nor")
    .vertex_in(2, Type::VEC4, "tan")
    .vertex_out(bnpr_v2f_geom_draw_remeshed_surface_depth)
    .fragment_source("npr_strokegen_remeshed_surface_frag.glsl")
    .fragment_out(0, Type::VEC4, "out_col");



/* Collect Contour Pixels --------------------------------- */
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
    .define("VERT_FLAGS_INCLUDED", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_vert_spatial_map_headers_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_vert_spatial_map_payloads_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_vert_merged_id_[]")
    .storage_buf(3, Qualifier::READ, "float", "ssbo_vbo_full_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .push_constant(Type::INT, "pcs_hash_map_size_")
    .push_constant(Type::INT, "pcs_vert_count_")

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
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_vert_merged_id_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_edge_index_map_headers_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_edge_spatial_map_payloads_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(6, Qualifier::READ, "float", "ssbo_vbo_full_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    /* Clear args */
    .storage_buf(9, Qualifier::WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_in_")
    .storage_buf(10, Qualifier::WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(11, Qualifier::WRITE, "SSBOData_StrokeGenEdgeSplitCounters", "ssbo_edge_split_counters_[]")
    .storage_buf(12, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_vert_debug_draw_args_")
    
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_hash_map_size_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

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
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_edges_hashing_finish)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_adjecency")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_HASHING_FINISH", "1");
GPU_SHADER_CREATE_INFO(bnpr_meshing_merge_edges_fill)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_adjecency")
    .define("_KERNEL_MULTICOMPILE__EDGE_ADJACENCY_FIND_ADJ", "1");



/* -------------------------------------------------------------------- */
/** \ Edge Selection
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("GLOBAL_COUNTER", "ssbo_bnpr_mesh_pool_counters_.num_filtered_edges")
    .define("CP_TAG", "filtered_edge")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_wedge_flooding_pointers_in_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_wedge_flooding_pointers_out_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(5, Qualifier::WRITE, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(7, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_output_edge_to_selected_edge_")
    .push_constant(Type::INT, "pcs_output_selected_edge_to_edge_")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_") 

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_iter)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_iter_init)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__INIT", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_last_iter_output_flags)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER", "1")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER", "1")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER__OUTPUT_FLAGS", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_wedge_flooding_last_iter_compaction)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_wedge_flooding")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER", "1")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER", "1")
    .define("_KERNEL_MULTICOMPILE__WEDGE_FLOODING__ITER__LAST_ITER__OUTPUT_CAMPACTION", "1");


GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_selection)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__SELECT_EDGES", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1") /* Remove compaction code */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_selection_mark_selection_border)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_selection")
    .define("_KERNEL_MULTICOMPILE__SELECT_EDGES__MARK_BORDERS", "1");




GPU_SHADER_CREATE_INFO(strokegen_meshing_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
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
/** \} */




GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("VE_CIRCULATOR_INCLUDE", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")
    .define("USE_DYNAMESH_VERT_SELECTION_INDEXING", "1")
    .define("INCLUDE_VERTEX_POSITION", "1")
    .define("INCLUDE_VERTEX_NORMAL", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("VERT_FLAGS_INCLUDED", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(5, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(6, Qualifier::READ, "uint", "ssbo_selected_vert_to_vert_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_vnor_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
#define NUM_SSBO_bnpr_meshing_surf_filtering_ 11

    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_") 

    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_geom_filtering_comp.glsl");


GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vnor_filtering)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_surf_filtering_")
    .push_constant(Type::INT, "pcs_vnor_filtering_iter_")
#define SSBO_OFFSET NUM_SSBO_bnpr_meshing_surf_filtering_
    .storage_buf(SSBO_OFFSET+0, Qualifier::READ_WRITE, "uint", "ssbo_vnor_temp_in_[]")
    .storage_buf(SSBO_OFFSET+1, Qualifier::READ_WRITE, "uint", "ssbo_vnor_temp_out_[]")
#undef SSBO_OFFSET
    .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__VNOR_FILTERING", "1");


GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vpos_filtering)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_surf_filtering_")
    .define("INCLUDE_VERTEX_REMESH_LEN", "1")
#define SSBO_OFFSET NUM_SSBO_bnpr_meshing_surf_filtering_
    .storage_buf(SSBO_OFFSET, Qualifier::READ_WRITE, "uint", "ssbo_vpos_temp_[]")
    .storage_buf(SSBO_OFFSET+1, Qualifier::READ_WRITE, "uint", "ssbo_vtx_remesh_len_[]")
#undef SSBO_OFFSET
    .push_constant(Type::INT, "pcs_vpos_filtering_iter_")
    .push_constant(Type::INT, "pcs_num_vpos_filtering_iters_")
    .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__VPOS_FILTERING", "1");


GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vquadric_common)
    .do_static_compilation(true)
    .define("QUADRICS_FILTERING_INCLUDE", "1")
    .additional_info("bnpr_meshing_surf_filtering_")
#define SSBO_OFFSET NUM_SSBO_bnpr_meshing_surf_filtering_
    .storage_buf(SSBO_OFFSET+0, Qualifier::READ_WRITE, "uint", "ssbo_vert_quadric_data_in_[]")
    .storage_buf(SSBO_OFFSET+1, Qualifier::READ_WRITE, "uint", "ssbo_vert_quadric_data_out_[]")
#undef SSBO_OFFSET
    .push_constant(Type::INT, "pcs_vq_filtering_iter_")
    .push_constant(Type::INT, "pcs_filtered_quadric_type_")
    .push_constant(Type::FLOAT, "pcs_quadric_deviation_")
    .push_constant(Type::FLOAT, "pcs_geodist_deviation_")
    .push_constant(Type::FLOAT, "pcs_position_regularization_scale_"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vquadric_diffusion)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_surf_filtering_vquadric_common")
    .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__VQUADRICS_DIFFUSION", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_quadric_vpos_filtering)
    .do_static_compilation(true)
    .define("QUADRICS_FILTERING_INCLUDE", "1")
    .additional_info("bnpr_meshing_surf_filtering_vquadric_common")
    .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING", "1")
    // .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__QUADRIC_VPOS_FILTERING__CONSTRAINED_SOLVE", "1")
    ;


GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_subd_common)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_surf_filtering_")
#define SSBO_OFFSET NUM_SSBO_bnpr_meshing_surf_filtering_
    .storage_buf(SSBO_OFFSET, Qualifier::READ_WRITE, "uint", "ssbo_vpos_subd_[]")
    .storage_buf(SSBO_OFFSET+1, Qualifier::READ_WRITE, "uint", "ssbo_epos_subd_[]")
    .push_constant(Type::INT, "pcs_subdiv_type_")
    .push_constant(Type::INT, "pcs_loop_subd_enable_crease_"); 
#undef SSBO_OFFSET 

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_subdiv_vpos_smoothing)
     .do_static_compilation(true)
     .additional_info("bnpr_meshing_surf_filtering_subd_common")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_subdiv_vpos_smoothing_finish)
     .do_static_compilation(true)
     .additional_info("bnpr_meshing_surf_filtering_subd_common")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__SUBDIV_VPOS_FILTERING_FINISH", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_subdiv_edge_points)
     .do_static_compilation(true)
     .additional_info("bnpr_meshing_surf_filtering_subd_common")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__EDGES", "1")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__LOOP_SUBD_EDGE_POINTS", "1"); 


GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vcurv_smoothing)
     .do_static_compilation(true)
     .additional_info("bnpr_meshing_surf_filtering_")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING", "1")
     .define("INCLUDE_VERTEX_CURV_MAX", "1")
     .define("INCLUDE_VERTEX_REMESH_LEN", "1")
#define SSBO_OFFSET NUM_SSBO_bnpr_meshing_surf_filtering_
    .storage_buf(SSBO_OFFSET, Qualifier::READ_WRITE, "uint", "ssbo_vcurv_max_[]")
    .storage_buf(SSBO_OFFSET+1, Qualifier::READ_WRITE, "uint", "ssbo_vcurv_max_temp_[]")
    .storage_buf(SSBO_OFFSET+2, Qualifier::READ_WRITE, "uint", "ssbo_vtx_remesh_len_[]")
#undef SSBO_OFFSET
    .push_constant(Type::INT, "pcs_vcurv_smooth_iter_");

GPU_SHADER_CREATE_INFO(bnpr_meshing_surf_filtering_vcurv_smoothing_output_remesh_len)
     .do_static_compilation(true)
     .additional_info("bnpr_meshing_surf_filtering_vcurv_smoothing")
     .define("_KERNEL_MULTICOMPILE__SURF_FILTERING__VCURVE_SMOOTHING__OUTPUT_REMESH_LEN", "1");
/** \} */


/* -------------------------------------------------------------------- */
/** \ Vertex Selection
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_select_verts)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("_KERNEL_MULTICOMPILE__SELECT_VERTS", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_selected_vert_to_vert_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(9, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")

    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::IVEC4, "pcs_vertex_selection_slots_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_wedge_flooding_compute.glsl");

GPU_SHADER_CREATE_INFO(strokegen_select_verts_from_selected_edges)
    .do_static_compilation(true)
    .additional_info("strokegen_select_verts")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("_KERNEL_MULTICOMPILE__SELECT_VERTS__FROM_SELECTED_EDGES", "1")
    .define("_KERNEL_MULTICOMPILE__SELECT_VERTS__FROM_SELECTED_EDGES_MAIN", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1");

GPU_SHADER_CREATE_INFO(strokegen_expand_verts_from_selected_edges)
    .do_static_compilation(true)
    .additional_info("strokegen_select_verts")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("_KERNEL_MULTICOMPILE__SELECT_VERTS__FROM_SELECTED_EDGES", "1")
    .define("_KERNEL_MULTICOMPILE__EXPAND_VERTS__FROM_SELECTED_EDGES", "1")
    .push_constant(Type::IVEC4, "pcs_vertex_selection_slots_out_")
    ;

GPU_SHADER_CREATE_INFO(strokegen_compact_selected_verts)
    .do_static_compilation(true)
    .additional_info("strokegen_select_verts")

    .define("_KERNEL_MULTICOMPILE__COMPACT_SELECTED_VERTS", "1")
    .define("GLOBAL_COUNTER", "ssbo_bnpr_mesh_pool_counters_.num_filtered_verts")
    .define("CP_TAG", "selected_vert")

    .push_constant(Type::INT, "pcs_vertex_select_all_slots_"); 


/** \} */

/* -------------------------------------------------------------------- */
/** \ Fill Dispatch Args for Remeshed Topology
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING", "1")
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_split_edge)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_EDGE", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenEdgeSplitCounters", "ssbo_edge_split_counters_[]")
    .storage_buf(1, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_split_edge_")
    .push_constant(Type::INT, "pcs_edge_split_dispatch_group_size_")
    .push_constant(Type::INT, "pcs_split_iter_"); 

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_collapsed_edge)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_COLLAPSED_EDGE", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenEdgeCollapseCounters", "ssbo_edge_collapse_counters_[]")
    .storage_buf(1, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_collapsed_edge_")
    .push_constant(Type::INT, "pcs_edge_collapse_dispatch_group_size_")
    .push_constant(Type::INT, "pcs_collapse_iter_"); 

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_flip_edge)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_FLIP_EDGE", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenEdgeFlipCounters", "ssbo_edge_flip_counters_[]")
    .storage_buf(1, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_flip_edge_")
    .push_constant(Type::INT, "pcs_edge_flip_dispatch_group_size_")
    .push_constant(Type::INT, "pcs_flip_iter_"); 

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_split_face)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__REMESHING__PER_SPLIT_FACE", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenFaceSplitCounters", "ssbo_face_split_counters_[]")
    .storage_buf(1, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_split_face_")
    .push_constant(Type::INT, "pcs_face_split_dispatch_group_size_")
    .push_constant(Type::INT, "pcs_split_iter_"); 

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_remeshed_edge)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_EDGE", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_remeshed_edges_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_only_selected_elems_")
    .push_constant(Type::INT, "pcs_remeshed_edges_dispatch_group_size_"); 

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_dispatch_args_per_remeshed_vert)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_dispatch_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS_PER_REMESHED_VERT", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_")
    .storage_buf(1, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::WRITE, "DispatchCommand", "ssbo_indirect_dispatch_args_per_remeshed_verts_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_only_selected_elems_")
    .push_constant(Type::INT, "pcs_remeshed_verts_dispatch_group_size_"); 

/** \} */


/* -------------------------------------------------------------------- */
/** \ Fill Draw Args for Remeshed Topology
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_draw_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DRAW_ARGS__REMESHING", "1")
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_remeshing_fill_draw_args_dbg_lines)
    .do_static_compilation(true)
    .additional_info("strokegen_remeshing_fill_draw_args")
    .define("_KERNEL_MULTICOMPILE__FILL_DRAW_ARGS__REMESHING__DBG_LINES", "1")
    .define("INCLUDE_DEBUG_LINE_CONFIG", "1")
    
    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::WRITE, "DrawCommand", "ssbo_bnpr_vert_debug_draw_args_")
    .push_constant(Type::INT, "pcs_line_type_");
/** \} */


/* -------------------------------------------------------------------- */
/** \ Edge Split
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    /* topo lib multicompile macros */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")
    .define("INCLUDE_VERTEX_REMESH_LEN", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_in_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenEdgeSplitCounters", "ssbo_edge_split_counters_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_per_edge_split_info_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_per_split_edge_info_[]")
    .storage_buf(11, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(12, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(13, Qualifier::READ_WRITE, "uint", "ssbo_vtx_remesh_len_[]")
    .storage_buf(14, Qualifier::READ_WRITE, "uint", "ssbo_epos_subd_[]")
    .storage_buf(15, Qualifier::READ_WRITE, "uint", "ssbo_vnor_[]")

    .push_constant(Type::INT, "pcs_split_iter_")
    .push_constant(Type::FLOAT, "pcs_remesh_edge_len_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    /* test contour vert insertion */
    .define("INCLUDE_VERTEX_NORMAL", "1")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_split_mode_")
    /**/

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_edge_split_comp.glsl");
;

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split_init)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_split")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT_INIT", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split_compact)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_split")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT_COMPACT", "1")
    .define("GLOBAL_COUNTER", "ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges_pass_1")
    .define("CP_TAG", "split_select_long_edges");

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split_exclude_border)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_split")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXCLUDE_BORDER", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split_resolve_conflict)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_split")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT_RESOLVE_CONFLICT", "1")
    .define("GLOBAL_COUNTER", "ssbo_edge_split_counters_[pcs_split_iter_].num_split_edges")
    .define("CP_TAG", "split_resolve_conflict"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_split_execute)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_split")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_SPLIT_EXECUTE", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 
/** \} */



/* -------------------------------------------------------------------- */
/** \ Edge Collapse
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    /* topo lib multicompile macros */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("VE_CIRCULATOR_INCLUDE", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_in_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenEdgeCollapseCounters", "ssbo_edge_collapse_counters_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_per_edge_collapse_info_in_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_per_edge_collapse_info_out_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "uint", "ssbo_per_collapse_edge_info_[]")
    .storage_buf(12, Qualifier::READ_WRITE, "uint", "ssbo_per_vert_collapse_wedge_id_[]") 
    .storage_buf(13, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(14, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .push_constant(Type::INT, "pcs_collapse_iter_")
    .push_constant(Type::FLOAT, "pcs_remesh_edge_len_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_edge_collapse_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_init)
    .do_static_compilation(true) 
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_INIT", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_compact) 
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_COMPACT", "1")
    .define("INCLUDE_VERTEX_REMESH_LEN", "1")
    .define("GLOBAL_COUNTER", "ssbo_edge_collapse_counters_[pcs_collapse_iter_].num_collapsed_edges_pass_1")
    .define("CP_TAG", "collapse_select_short_edges"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_resolve_conflict_0)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_0", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_resolve_conflict_1)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_1", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_resolve_conflict_2)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_RESOLVE_CONFLICT__PASS_2", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_collapse_execute)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_collapse")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_COLLAPSE_EXECUTE", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("INCLUDE_VERTEX_REMESH_LEN", "1"); 
/** \} */

/* -------------------------------------------------------------------- */
/** \ Edge Flip
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    /* topo lib multicompile macros */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("VE_CIRCULATOR_INCLUDE", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_in_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenEdgeFlipCounters", "ssbo_edge_flip_counters_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_per_edge_flip_info_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_per_flip_edge_info_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_vertex_edge_flip_info_[]")
    .storage_buf(11, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(12, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(13, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .push_constant(Type::INT, "pcs_flip_opti_goal_type_") 
    .push_constant(Type::INT, "pcs_flip_iter_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_edge_flip_comp.glsl")
;

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip_init)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_flip")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP_INIT", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip_compact)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_flip")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP_COMPACT", "1")
    .define("GLOBAL_COUNTER", "ssbo_edge_flip_counters_[pcs_flip_iter_].num_flip_edges_pass_1")
    .define("CP_TAG", "flip_inital_selection");

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip_validate)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_flip")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP_VALIDATE", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip_resolve_conflict)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_flip")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP_RESOLVE_CONFLICT", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_edge_flip_execute)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_edge_flip")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP", "1")
    .define("_KERNEL_MULTICOMPILE__EDGE_FLIP_EXECUTE", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1"); 
/** \} */



/* -------------------------------------------------------------------- */
/** \ Face Split
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_meshing_face_split)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    /* topo lib multicompile macros */
    .define("DECODE_IBO_EXCLUDE", "1") /* Remove ibo code */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")

    .storage_buf(0, Qualifier::READ, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(1, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(2, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_in_")
    .storage_buf(3, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(4, Qualifier::READ_WRITE, "SSBOData_StrokeGenFaceSplitCounters", "ssbo_face_split_counters_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(9, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .storage_buf(11, Qualifier::READ_WRITE, "uint", "ssbo_per_face_split_info_[]")
    .push_constant(Type::INT, "pcs_split_iter_")
    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_face_split_comp.glsl"); 

GPU_SHADER_CREATE_INFO(bnpr_meshing_face_split_init)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_face_split")
    .define("_KERNEL_MULTICOMPILE__FACE_SPLIT_INIT", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_face_split_work_generation)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_face_split")
    .define("_KERNEL_MULTICOMPILE__FACE_SPLIT_WORK_GEN", "1");

GPU_SHADER_CREATE_INFO(bnpr_meshing_face_split_execute)
    .do_static_compilation(true)
    .additional_info("bnpr_meshing_face_split")
    .define("_KERNEL_MULTICOMPILE__FACE_SPLIT_EXECUTE", "1");

/** \} */



/* -------------------------------------------------------------------- */
/** \Surface Geometry Analysis
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh")
    /* topo lib multicompile macros */
    .define("WINGED_EDGE_TOPO_INCLUDE", "1")
    .define("VERT_WEDGE_LIST_TOPO_INCLUDE", "1")
    .define("VERT_FLAGS_INCLUDED", "1")
    .define("VE_CIRCULATOR_INCLUDE", "1")
    .define("USE_DYNAMESH_EDGE_SELECTION_INDEXING", "1")
    .define("USE_DYNAMESH_VERT_SELECTION_INDEXING", "1")
    .define("INCLUDE_VERTEX_POSITION", "1")
    .define("INCLUDE_DEBUG_LINE_CONFIG", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "SSBOData_StrokeGenDynamicMeshCounters", "ssbo_dyn_mesh_counters_out_")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .storage_buf(2, Qualifier::READ, "uint", "ssbo_selected_edge_to_edge_[]")
    .storage_buf(3, Qualifier::READ, "uint", "ssbo_selected_vert_to_vert_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_vert_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_edge_to_edges_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_vert_to_edge_list_header_[]")
    .storage_buf(7, Qualifier::READ_WRITE, "uint", "ssbo_vert_flags_[]")
    .storage_buf(8, Qualifier::READ_WRITE, "float", "ssbo_vbo_full_[]")
    .storage_buf(9, Qualifier::WRITE, "uint", "ssbo_dbg_lines_[]")
#define NUM_SSBO_BASE 10

    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")

    .push_constant(Type::INT, "pcs_edge_count_")
    .push_constant(Type::INT, "pcs_vert_count_")
    .push_constant(Type::INT, "pcs_rsc_handle")
    .push_constant(Type::INT, "pcs_output_dbg_geom_")
    .push_constant(Type::FLOAT, "pcs_dbg_geom_scale_")
    .push_constant(Type::INT, "pcs_only_selected_verts_")

    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT)
    .compute_source("npr_strokegen_geom_analysis_comp.glsl"); 

/* Estimate order 0 vertex attributes */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_0_vert_attrs)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1");
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_0_vert_normal)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_0_vert_attrs")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__NORMAL", "1")
    .define("INCLUDE_VERTEX_NORMAL", "1")
    .storage_buf(NUM_SSBO_BASE, Qualifier::READ_WRITE, "uint", "ssbo_vnor_[]")
    .push_constant(Type::INT, "pcs_output_vertex_facing_flag_");

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_0_vert_normal_voroarea)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_0_vert_normal")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA", "1")
    .define("INCLUDE_VERTEX_VORONOI_AREA", "1")
    .storage_buf(NUM_SSBO_BASE + 1, Qualifier::READ_WRITE, "uint", "ssbo_varea_[]");

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_0_vert_normal_topoflags)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_0_vert_normal")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER", "1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .storage_buf(NUM_SSBO_BASE + 1, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]");

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_0_vert_normal_voroarea_topoflags)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_0_vert_normal")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__AREA", "1")
    .define("INCLUDE_VERTEX_VORONOI_AREA", "1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__BORDER", "1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_0__CREASE", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .storage_buf(NUM_SSBO_BASE + 1, Qualifier::READ_WRITE, "uint", "ssbo_varea_[]")
    .storage_buf(NUM_SSBO_BASE + 2, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]");


/* Estimate order 1 vertex attributes */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1", "1")
    .define("INCLUDE_VERTEX_NORMAL", "1")
    .define("INCLUDE_VERTEX_VORONOI_AREA", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .storage_buf(NUM_SSBO_BASE + 0u, Qualifier::READ_WRITE, "uint", "ssbo_vnor_[]")
    .storage_buf(NUM_SSBO_BASE + 1u, Qualifier::READ_WRITE, "uint", "ssbo_varea_[]")
    .storage_buf(NUM_SSBO_BASE + 2u, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]");
#define NUM_SSBO_1 ((NUM_SSBO_BASE + 3u))

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1_main)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__MAIN", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1_vert_curv_pass_0)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURV_PER_FACE", "1")
    .storage_buf(NUM_SSBO_1 + 0, Qualifier::READ_WRITE, "uint", "ssbo_edge_vtensors_[]");

GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1_main_curvature)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_1_main")

    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__MAIN", "1")
    .define("INCLUDE_VERTEX_CURV_TENSOR", "1")
    .define("INCLUDE_VERTEX_CURV_MAX", "1")
    .storage_buf(NUM_SSBO_1 + 0, Qualifier::READ_WRITE, "uint", "ssbo_edge_vtensors_[]")
    .storage_buf(NUM_SSBO_1 + 1, Qualifier::READ_WRITE, "uint", "ssbo_vcurv_pdirs_k1k2_[]")
    .storage_buf(NUM_SSBO_1 + 2, Qualifier::READ_WRITE, "uint", "ssbo_vcurv_max_[]")
    .push_constant(Type::INT, "pcs_output_curv_tensors_"); 

/* Rusinkiewicz's Curvature Estimator
 * Outputs both Curvature and Tensor
 * Requires 2 passes, one per-edge, one per-vertex 
 * https://github.com/Forceflow/trimesh2/blob/main/libsrc/TriMesh_curvature.cc */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1_main_curvature_rusinkiewicz)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_1_main_curvature")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__CURVTENSOR", "1"); 

/* JacquesOlivierLachaud's Curvature Estimator
 * Outputs both Curvature and Tensor, single pass
 * https://github.com/dcoeurjo/CorrectedNormalCurrent/blob/master/CorrectedNormalCurrentFormula.h
 * https://dgtal-team.github.io/doc-nightly/moduleCurvatureMeasures.html
 * https://github.com/CGAL/cgal/issues/7063 */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_order_1_main_curvature_jacques)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis_order_1_main_curvature")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING", "1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVTENSOR", "1")
    .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1_GRAD_VDOTN", "1")
    // .define("INCLUDE_VERTEX_RADIAL_NORMAL", "1")
    
    /* Only for solving Gaussian & Mean curvatures, but seems not correct, 
     * need fix if we want this cheaper version */
    // .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_2RING", "1")
    // .define("_KERNEL_MULTICOMPILE__CALC_VERT_ATTRS_ORDER_1__INTERPO_CURVATURE", "1")

    .push_constant(Type::INT, "pcs_output_maxcurv_with_cusp_function_")
    ; 

#undef NUM_SSBO_1

/* Estimate order 0 edge attributes */
GPU_SHADER_CREATE_INFO(bnpr_geom_analysis_feature_edges)
    .do_static_compilation(true)
    .additional_info("bnpr_geom_analysis")
    .define("_KERNEL_MULTICOMPILE__CALC_FEATURE_EDGES", "1")
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1")
    .define("EDGE_FLAGS_INCLUDED", "1")
    .storage_buf(NUM_SSBO_BASE, Qualifier::READ_WRITE, "uint", "ssbo_edge_flags_[]")
    .push_constant(Type::INT, "pcs_only_selected_edges_");

#undef NUM_SSBO_BASE


/** \} */



/* -------------------------------------------------------------------- */
/** \Draw Debug Primitives
 * \{ */
GPU_SHADER_INTERFACE_INFO(bnpr_v2f_geom_draw_normal_lines, "")
    /* .flat(Type::UINT, "id") */
    .smooth(Type::VEC4, "color") ;

GPU_SHADER_CREATE_INFO(bnpr_geom_draw_debug_lines)
    .do_static_compilation(true)
    .typedef_source("bnpr_shader_shared.hh")
    .additional_info("draw_modelmat_new", "draw_view", "draw_resource_handle_new")
    .define("INCLUDE_DEBUG_LINE_CONFIG", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_dbg_lines_[]")
    .storage_buf(1, Qualifier::READ_WRITE, "SSBOData_StrokeGenMeshPoolCounters", "ssbo_bnpr_mesh_pool_counters_")
    .uniform_buf(0, "ViewMatrices", "ubo_view_matrices_")
    .push_constant(Type::INT, "pcs_line_type_")

    .vertex_source("npr_strokegen_debug_lines_vert.glsl")
    .vertex_in(0, Type::VEC3, "pos")
    .vertex_in(1, Type::VEC3, "nor")
    .vertex_in(2, Type::VEC4, "tan")
    .vertex_out(bnpr_v2f_geom_draw_normal_lines)
    .fragment_source("npr_strokegen_debug_lines_frag.glsl")
    .fragment_out(0, Type::VEC4, "out_col")
/*     .fragment_out(1, Type::VEC3, "out_normal")
    .fragment_out(2, Type::VEC4, "out_tangent") */;
/** \} */








/* -------------------------------------------------------------------- */
/** \Bare minumum compute shader test
 * \{ */
// GPU_SHADER_CREATE_INFO(bnpr_strokegen_test_xxx)
//     .do_static_compilation(true)
//     .storage_buf(0, Qualifier::READ_WRITE, "uint", "buf_test[]")
//     .storage_buf(1, Qualifier::READ, "uint", "buf_ibo[]")
//     .local_group_size(GROUP_SIZE_STROKEGEN_TEST) /* <== from "bnpr_defines.hh" */
//     .compute_source("npr_strokegen_test_comp.glsl");
/** \} */

/** GPU Scan / Segscan --------------------
 */
GPU_SHADER_CREATE_INFO(bnpr_scan_uint_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uint")
    .define("SCAN_OP", "u32_add")
    .define("SCAN_ZERO_VAL", "0u")
    .define("SCAN_FUNCTION_TAG", "_u32_add");
GPU_SHADER_CREATE_INFO(bnpr_scan_uint_min)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uint")
    .define("SCAN_OP", "u32_min")
    .define("SCAN_ZERO_VAL", "0x3fffffffu") // 2 bit(s) occupied by segscan hf
    .define("SCAN_FUNCTION_TAG", "_u32_min");
GPU_SHADER_CREATE_INFO(bnpr_segscan_uint_storage)
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanType_uint")
    .define("SEGSCAN_STRUCT_TYPE_ENCODED", "uint")
    .define("SEGSCAN_ENCODE", "segscan_uint_hf_encode")
    .define("SEGSCAN_DECODE", "segscan_uint_hf_decode"); 

GPU_SHADER_CREATE_INFO(bnpr_scan_uvec2_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uvec2")
    .define("SCAN_OP", "uvec2_add")
    .define("SCAN_ZERO_VAL", "uvec2(0u, 0u)")
    .define("SCAN_FUNCTION_TAG", "uvec2_add");
GPU_SHADER_CREATE_INFO(bnpr_segscan_uvec2_storage)
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanType_uvec2")
    .define("SEGSCAN_STRUCT_TYPE_ENCODED", "uvec2")
    .define("SEGSCAN_ENCODE", "segscan_uvec2_hf_encode")
    .define("SEGSCAN_DECODE", "segscan_uvec2_hf_decode"); 

// uvec3 is not supported since std430 layout requires 8-byte or 16-byte alignment
// GPU_SHADER_CREATE_INFO(bnpr_scan_uvec3_add) 
//     .typedef_source("bnpr_shader_shared.hh")
//     .define("SCAN_DATA_TYPE", "uvec3")
//     .define("SCAN_OP", "uvec3_add")
//     .define("SCAN_ZERO_VAL", "uvec3(0u, 0u, 0u)")
//     .define("SCAN_FUNCTION_TAG", "uvec3_add");

GPU_SHADER_CREATE_INFO(bnpr_scan_uvec4_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "uvec4")
    .define("SCAN_OP", "uvec4_add")
    .define("SCAN_ZERO_VAL", "uvec4(0u, 0u, 0u, 0u)")
    .define("SCAN_FUNCTION_TAG", "uvec4_add");
GPU_SHADER_CREATE_INFO(bnpr_segscan_uvec4_storage)
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanType_uvec4")
    .define("SEGSCAN_STRUCT_TYPE_ENCODED", "uvec4")
    .define("SEGSCAN_ENCODE", "segscan_uvec4_hf_encode")
    .define("SEGSCAN_DECODE", "segscan_uvec4_hf_decode");    

GPU_SHADER_CREATE_INFO(bnpr_scan_float_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "float")
    .define("SCAN_OP", "f32_add")
    .define("SCAN_ZERO_VAL", ".0f")
    .define("SCAN_FUNCTION_TAG", "_f32_add");
GPU_SHADER_CREATE_INFO(bnpr_segscan_float_storage)
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanType_float")
    .define("SEGSCAN_STRUCT_TYPE_ENCODED", "uvec2")
    .define("SEGSCAN_ENCODE", "segscan_float_hf_encode")
    .define("SEGSCAN_DECODE", "segscan_float_hf_decode"); 

GPU_SHADER_CREATE_INFO(bnpr_scan_vec3_add)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SCAN_DATA_TYPE", "vec3")
    .define("SCAN_OP", "vec3_add")
    .define("SCAN_ZERO_VAL", "vec3(.0f, .0f, .0f)")
    .define("SCAN_FUNCTION_TAG", "_vec3_add");
GPU_SHADER_CREATE_INFO(bnpr_segscan_vec3_storage)
    .define("SEGSCAN_STRUCT_TYPE", "SSBOData_SegScanType_vec3")
    .define("SEGSCAN_STRUCT_TYPE_ENCODED", "uvec4")
    .define("SEGSCAN_ENCODE", "segscan_vec3_hf_encode")
    .define("SEGSCAN_DECODE", "segscan_vec3_hf_decode"); 


GPU_SHADER_CREATE_INFO(bnpr_scan_test_inputs)
    .additional_info("bnpr_scan_float_add")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_TEST", "1");

GPU_SHADER_CREATE_INFO(bnpr_scan_uint_add_inputs)
    .additional_info("bnpr_scan_uint_add")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_INFO_SSBO", "1");


GPU_SHADER_CREATE_INFO(bnpr_segscan_test_inputs)
    .additional_info("bnpr_scan_uvec4_add")
    .additional_info("bnpr_segscan_uvec4_storage")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_INFO_SSBO", "0")

    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_TEST")
    .push_constant(Type::INT, "pcs_segscan_test_random_seed_");

GPU_SHADER_CREATE_INFO(bnpr_segscan_uint_add_inputs)
    .additional_info("bnpr_scan_uint_add")
    .additional_info("bnpr_segscan_uint_storage")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_INFO_SSBO", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_segscan_uint_min_inputs)
    .additional_info("bnpr_scan_uint_min")
    .additional_info("bnpr_segscan_uint_storage")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_INFO_SSBO", "1"); 

GPU_SHADER_CREATE_INFO(bnpr_segscan_float_add_inputs)
    .additional_info("bnpr_scan_float_add")
    .additional_info("bnpr_segscan_float_storage")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_INFO_SSBO", "1"); 


/** GPU Compaction --------------------
 */
GPU_SHADER_CREATE_INFO(npr_compaction_off)
    .define("COMPACTION_LIB_EXCLUDE_DEFAULT_CODEGEN", "1");

/** GPU List Ranking --------------------
 */



/* -------------------------------------------------------------------- */
/** \GPU (Segmented)Scan Test
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_scan_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__SCAN", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_")
    .storage_buf(1, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_scan_dispatch_args_")
    
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

#define GPU_SHADER_CREATE_INFO__TREE_SCAN(name, input_shader_info, scan_data_type)                     \
GPU_SHADER_CREATE_INFO(name##_upsweep)                                                                 \
    .do_static_compilation(true)                                                                       \
    .additional_info(#input_shader_info)                                                               \
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_UPSWEEP", "1")                                           \
                                                                                                       \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type, "bnpr_in_scan_data_buf_[]")                \
    .storage_buf(1, Qualifier::READ_WRITE, #scan_data_type, "bnpr_out_scan_data_buf_[]")               \
    .storage_buf(2, Qualifier::WRITE, #scan_data_type, "bnpr_scan_block_sum_buf_[]")                   \
    .storage_buf(3, Qualifier::READ_WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_")               \
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")                                    \
                                                                                                       \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_SWEEP)                                                      \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                              \
                                                                                                       \
GPU_SHADER_CREATE_INFO(name##_aggregate)                                                               \
    .do_static_compilation(true)                                                                       \
    .additional_info(#input_shader_info)                                                               \
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_AGGREGATE", "1")                                         \
                                                                                                       \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type, "bnpr_scan_block_sum_buf_[]")              \
                                                                                                       \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_AGGRG)                                                      \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                              \
                                                                                                       \
GPU_SHADER_CREATE_INFO(name##_dwsweep)                                                                 \
    .do_static_compilation(true)                                                                       \
    .additional_info(#input_shader_info)                                                               \
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_DWSWEEP", "1")                                           \
                                                                                                       \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type, "bnpr_out_scan_data_buf_[]")               \
    .storage_buf(1, Qualifier::READ, #scan_data_type, "bnpr_scan_block_sum_buf_[]")                    \
                                                                                                       \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_SWEEP)                                                      \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                              \


GPU_SHADER_CREATE_INFO__TREE_SCAN(bnpr_scan_test, bnpr_scan_test_inputs, SCAN_DATA_TYPE)

GPU_SHADER_CREATE_INFO__TREE_SCAN(bnpr_scan_uint_add, bnpr_scan_uint_add_inputs, SCAN_DATA_TYPE)


#define GPU_SHADER_CREATE_INFO__TREE_SEGSCAN(name, input_shader_info, scan_data_type_str)     \
GPU_SHADER_CREATE_INFO(name##_upsweep)                                                        \
    .do_static_compilation(true)                                                              \
    .additional_info(#input_shader_info)                                                      \
    .define("IS_TREE_SEG_SCAN", "1")                                                          \
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_UPSWEEP", "1")                              \
                                                                                              \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type_str, "bnpr_in_scan_data_buf_[]")   \
    .storage_buf(1, Qualifier::READ_WRITE, #scan_data_type_str, "bnpr_out_scan_data_buf_[]")  \
    .storage_buf(2, Qualifier::WRITE, #scan_data_type_str, "bnpr_scan_block_sum_buf_[]")      \
    .storage_buf(3, Qualifier::READ_WRITE, "UBData_TreeScan", "ssbo_tree_scan_infos_")        \
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")                           \
                                                                                              \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_SWEEP)                                             \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                     \
                                                                                              \
GPU_SHADER_CREATE_INFO(name##_aggregate)                                                      \
    .do_static_compilation(true)                                                              \
    .additional_info(#input_shader_info)                                                      \
    .define("IS_TREE_SEG_SCAN", "1")                                                          \
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_AGGREGATE", "1")                            \
                                                                                              \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type_str, "bnpr_scan_block_sum_buf_[]") \
                                                                                              \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_AGGRG)                                             \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                     \
                                                                                              \
GPU_SHADER_CREATE_INFO(name##_dwsweep)                                                        \
    .do_static_compilation(true)                                                              \
    .additional_info(#input_shader_info)                                                      \
    .define("IS_TREE_SEG_SCAN", "1")                                                          \
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_DWSWEEP", "1")                              \
                                                                                              \
    .storage_buf(0, Qualifier::READ_WRITE, #scan_data_type_str, "bnpr_out_scan_data_buf_[]")  \
    .storage_buf(1, Qualifier::READ, #scan_data_type_str, "bnpr_scan_block_sum_buf_[]")       \
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")                           \
                                                                                              \
    .local_group_size(GROUP_SIZE_BNPR_SCAN_SWEEP)                                             \
    .compute_source("npr_strokegen_scan_test_comp.glsl");                                     \


GPU_SHADER_CREATE_INFO__TREE_SEGSCAN(bnpr_segscan_test, bnpr_segscan_test_inputs, SEGSCAN_STRUCT_TYPE_ENCODED)

GPU_SHADER_CREATE_INFO__TREE_SEGSCAN(bnpr_segscan_uint_add, bnpr_segscan_uint_add_inputs, SEGSCAN_STRUCT_TYPE_ENCODED)

GPU_SHADER_CREATE_INFO__TREE_SEGSCAN(bnpr_segscan_uint_min, bnpr_segscan_uint_min_inputs, SEGSCAN_STRUCT_TYPE_ENCODED)

GPU_SHADER_CREATE_INFO__TREE_SEGSCAN(bnpr_segscan_float_add, bnpr_segscan_float_add_inputs, SEGSCAN_STRUCT_TYPE_ENCODED)
/** \} */


/* -------------------------------------------------------------------- */
/** \GPU 1D Segmented Looped Convolution
 * \{ */
GPU_SHADER_CREATE_INFO(strokegen_segloopconv1d_fill_dispatch_args)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS", "1")
    .define("_KERNEL_MULTICOMPILE__FILL_DISPATCH_ARGS__SEGLOOPCONV1D", "1")

    .storage_buf(0, Qualifier::READ_WRITE, "UBData_SegLoopConv1D", "ssbo_segloopconv1d_info_")
    .storage_buf(1, Qualifier::READ_WRITE, "DispatchCommand", "ssbo_segloopconv1d_dispatch_args_")
    .push_constant(Type::INT, "pc_segloopconv1d_dispatch_group_size_")
    
    .local_group_size(32)
    .compute_source("npr_strokegen_fill_indirect_args_comp.glsl");

GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SEGLOOPCONV1D_USE_ADVANCED_INPUT", "1")
    .define("LOOPCONV1D_TAG", "build_patch")
    .define("LOOPCONV1D_MAX_RADIUS", NPR_TEST_SEGLOOPCONV1D_CONV_RADIUS_STR)
    .define("DATA_TYPE_LOOPCONV1D", "float")
    .define("_KERNEL_MULTI_COMPILE__SEGLOOPCONV1D_INFO_SSBO", "0")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__TEST", "1");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test_build_patch)
    .additional_info("npr_segloopconv1D_test")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_debug_segloopconv1d_data_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_test_convolution)
    .additional_info("npr_segloopconv1D_test")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_debug_segloopconv1d_data_[]");

GPU_SHADER_CREATE_INFO(npr_segloopconv1D_seg_denoising)
    .typedef_source("bnpr_shader_shared.hh")
    .define("LOOPCONV1D_TAG", "seg_denoising")
    .define("LOOPCONV1D_MAX_RADIUS", "8")
    .define("DATA_TYPE_LOOPCONV1D", "uint")
    .define("_KERNEL_MULTI_COMPILE__SEGLOOPCONV1D_INFO_SSBO", "1")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__SEG_DENOISING", "1");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_seg_denoising_build_patch)
    .additional_info("npr_segloopconv1D_seg_denoising")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_list_len_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_flags_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_seg_denoising_convolution)
    .additional_info("npr_segloopconv1D_seg_denoising")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_rank_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_snake_list_len_[]");

GPU_SHADER_CREATE_INFO(npr_segloopconv1D_2d_sample_processing)
    .typedef_source("bnpr_shader_shared.hh")
    .define("SEGLOOPCONV1D_USE_ADVANCED_INPUT", "1")
    .define("_KERNEL_MULTI_COMPILE__SEGLOOPCONV1D_INFO_SSBO", "1")
    .define("LOOPCONV1D_TAG", "2d_sample_processing")
    .define("LOOPCONV1D_MAX_RADIUS", "32")
    .define("USE_CONTOUR_2D_SAMPLE_TOPOLOGY_BUFFER", "1")
    .define("USE_CONTOUR_2D_SAMPLE_GEOMETRY_BUFFER", "1")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_PROCESSING", "1")
    .push_constant(Type::VEC2, "pcs_screen_size_")
    .image(0, GPU_RGBA32F, Qualifier::WRITE, ImageType::FLOAT_2D, "tex2d_contour_dbg_");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_2d_samples_build_patch)
    .additional_info("npr_segloopconv1D_2d_sample_processing")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0", "1")
    .define("DATA_TYPE_LOOPCONV1D", "vec2")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_corner_detection_convolution_step_0)
    .additional_info("npr_segloopconv1D_2d_sample_processing")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_0", "1")
    .define("DATA_TYPE_LOOPCONV1D", "vec2")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_corner_detection_convolution_step_1)
    .additional_info("npr_segloopconv1D_2d_sample_processing")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CORNER_DETECTION__STEP_1", "1")
    .define("DATA_TYPE_LOOPCONV1D", "float")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_2d_sample_calc_tangent_curvature_convolution)
    .additional_info("npr_segloopconv1D_2d_sample_processing")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2DSAMPLE_CALC_TANGENT_CURVATURE", "1")
    .define("DATA_TYPE_LOOPCONV1D", "vec2")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]");
GPU_SHADER_CREATE_INFO(npr_segloopconv1D_2d_sample_visibility_denoising)
    .additional_info("npr_segloopconv1D_2d_sample_processing")
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION__2D_SAMPLE_VISIBILITY_DENOISING", "1")
    .define("DATA_TYPE_LOOPCONV1D", "uint")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_topology_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_contour_2d_sample_geometry_[]");

#define GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_BUILD_PATCH(name, input_shader_info_build_patch) \
GPU_SHADER_CREATE_INFO(strokegen_segloopconv1D_##name##_build_patch)                     \
    .do_static_compilation(true)                                                         \
    .additional_info(#input_shader_info_build_patch)                                     \
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_BUILD_PATCH_TABLE", "1")                    \
    .storage_buf(0, Qualifier::READ_WRITE, "uint", "ssbo_segloopconv1d_patch_table_[]")  \
    .storage_buf(1, Qualifier::READ, "UBData_SegLoopConv1D", "ssbo_segloopconv1d_info_") \
    .uniform_buf(0, "UBData_SegLoopConv1D", "ubo_segloopconv1d_")                        \
                                                                                         \
    .local_group_size(GROUP_SIZE_SEGLOOPCONV1D_TEST)                                     \
    .compute_source("npr_strokegen_segloopconv1D_test_comp.glsl");                       \

GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_BUILD_PATCH(test, npr_segloopconv1D_test_build_patch)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_BUILD_PATCH(seg_denoising, npr_segloopconv1D_seg_denoising_build_patch)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_BUILD_PATCH(samples_2d, npr_segloopconv1D_2d_samples_build_patch)


#define GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(name, input_shader_info_conv) \
GPU_SHADER_CREATE_INFO(strokegen_segloopconv1D_##name##_convolution)                     \
    .do_static_compilation(true)                                                         \
    .define("_KERNEL_MULTICOMPILE__1DSEGLOOP_CONVOLUTION", "1")                          \
    .additional_info(#input_shader_info_conv)                                            \
    .storage_buf(0, Qualifier::READ, "uint", "ssbo_segloopconv1d_patch_table_[]")        \
    .storage_buf(1, Qualifier::READ_WRITE, "uint", "ssbo_in_segloopconv1d_data_[]")      \
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "ssbo_out_segloopconv1d_data_[]")     \
    .storage_buf(3, Qualifier::READ, "UBData_SegLoopConv1D", "ssbo_segloopconv1d_info_") \
    .uniform_buf(0, "UBData_SegLoopConv1D", "ubo_segloopconv1d_")                        \
                                                                                         \
    .local_group_size(GROUP_SIZE_SEGLOOPCONV1D_TEST)                                     \
    .compute_source("npr_strokegen_segloopconv1D_test_comp.glsl");                       \

GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(test, npr_segloopconv1D_test_convolution)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(seg_denoising, npr_segloopconv1D_seg_denoising_convolution)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(corner_detection_step_0, npr_segloopconv1D_corner_detection_convolution_step_0)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(corner_detection_step_1, npr_segloopconv1D_corner_detection_convolution_step_1)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(calc_2d_sample_tangent_curv, npr_segloopconv1D_2d_sample_calc_tangent_curvature_convolution)
GPU_SHADER_CREATE_INFO__SEGLOOPCONV1D_CONV(denoise_2d_sample_visibility, npr_segloopconv1D_2d_sample_visibility_denoising)
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
    .typedef_source("draw_shader_shared.hh") /* Always needed for indirect args */
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
    .typedef_source("draw_shader_shared.hh")
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
    .typedef_source("draw_shader_shared.hh")
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
    .storage_buf(7, Qualifier::WRITE, "uint", "ssbo_list_ranking_serialized_topo_[]")
    .storage_buf(8, Qualifier::READ, "SSBOData_ListRankingInputs", "ssbo_list_ranking_inputs_")
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
    .storage_buf(3, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_addressing_counters_[]")
    .storage_buf(4, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_output_ranks_[]")
    .storage_buf(5, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_output_list_len_[]")
    .storage_buf(6, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_list_head_info_[]")
    .push_constant(Type::INT, "pcs_contour_edge_linking_output_pass_type_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_output_pass_0)
    .do_static_compilation(true)
    .additional_info("strokegen_list_ranking_test_output")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_OUTPUT__PASS_0", "1");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_output_pass_1)
    .do_static_compilation(true)
    .additional_info("strokegen_list_ranking_test_output")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_OUTPUT__PASS_1", "1");


/** \} */
