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

/** GPU List Ranking --------------------
 */

/* -------------------------------------------------------------------- */
/** \Geometry extraction from GPUBatch(es)
 * \{ */
GPU_SHADER_CREATE_INFO(bnpr_geom_extract)
    .typedef_source("bnpr_shader_shared.hh")
    .do_static_compilation(true)
    .do_static_compilation(true)
    .storage_buf(0, Qualifier::READ, "uint", "buf_ibo[]")
    .storage_buf(1, Qualifier::READ, "vec3", "buf_vbo[]")
    .storage_buf(2, Qualifier::READ_WRITE, "uint", "buf_strokegen_mesh_pool[]")
    .push_constant(Type::BOOL, "pcs_ib_fmt_u16")
    .push_constant(Type::UINT, "pcs_num_verts")
    .push_constant(Type::UINT, "pcs_num_ib_offset")
    .local_group_size(GROUP_SIZE_STROKEGEN_GEOM_EXTRACT) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_geom_extract_comp.glsl");
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
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SCAN_TEST_DATA_TYPE_STR,
                 "bnpr_in_scan_data_buf_[]")
    .storage_buf(1,
                 Qualifier::READ_WRITE,
                 BNPR_SCAN_TEST_DATA_TYPE_STR,
                 "bnpr_out_scan_data_buf_[]")
    .storage_buf(2, Qualifier::WRITE, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_scan_test_aggregate)
    .do_static_compilation(true)
    .additional_info("bnpr_scan_test_inputs")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_AGGREGATE", "1")
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SCAN_TEST_DATA_TYPE_STR,
                 "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_AGGRG) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_scan_test_dwsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_scan_test_inputs")
    .define("_KERNEL_MULTI_COMPILE__TREE_SCAN_DWSWEEP", "1")
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SCAN_TEST_DATA_TYPE_STR,
                 "bnpr_out_scan_data_buf_[]")
    .storage_buf(1, Qualifier::READ, BNPR_SCAN_TEST_DATA_TYPE_STR, "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_upsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_UPSWEEP", "1")
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_in_scan_data_buf_[]")
    .storage_buf(1,
                 Qualifier::READ_WRITE,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_out_scan_data_buf_[]")
    .storage_buf(2,
                 Qualifier::WRITE,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_scan_block_sum_buf_[]")
    .uniform_buf(0, "UBData_TreeScan", "ubo_bnpr_tree_scan_infos_")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_SWEEP) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_aggregate)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_AGGREGATE", "1")
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_scan_block_sum_buf_[]")
    .local_group_size(GROUP_SIZE_BNPR_SCAN_TEST_AGGRG) /* <== from "bnpr_defines.hh" */
    .compute_source("npr_strokegen_scan_test_comp.glsl");

GPU_SHADER_CREATE_INFO(bnpr_segscan_test_dwsweep)
    .do_static_compilation(true)
    .additional_info("bnpr_segscan_test_inputs")
    .define("IS_TREE_SEG_SCAN", "1")
    .define("_KERNEL_MULTI_COMPILE__TREE_SEG_SCAN_DWSWEEP", "1")
    .storage_buf(0,
                 Qualifier::READ_WRITE,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_out_scan_data_buf_[]")
    .storage_buf(1,
                 Qualifier::READ,
                 BNPR_SEG_SCAN_TEST_STRUCT_TYPE_STR,
                 "bnpr_scan_block_sum_buf_[]")
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
    .storage_buf(1,
                 Qualifier::READ_WRITE,
                 NPR_SEGLOOPCONV1D_TEST_DATA_TYPE_STR,
                 "ssbo_in_segloopconv1d_data_[]")
    .storage_buf(2,
                 Qualifier::READ_WRITE,
                 NPR_SEGLOOPCONV1D_TEST_DATA_TYPE_STR,
                 "ssbo_out_segloopconv1d_data_[]")
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
    .storage_buf(10, Qualifier::READ_WRITE, "uint", "ssbo_list_ranking_debug_[]")
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
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_splice_iter_")
    .push_constant(Type::INT, "pc_listranking_jumping_iter_")
    .push_constant(Type::INT, "pc_listranking_ranking_pass_with_broken_loops_")

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
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    .push_constant(Type::INT, "pc_listranking_relink_iter_")
    .push_constant(Type::INT, "pc_listranking_num_relink_iters_")

    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

GPU_SHADER_CREATE_INFO(strokegen_list_ranking_test_upload_cpu_test_data)
    .do_static_compilation(true)
    .typedef_source("bnpr_defines.hh")
    .typedef_source("bnpr_shader_shared.hh")
    .define("_KERNEL_MULTICOMPILE__TEST_LIST_RANKING_FILL_CPU_DATA", "1")

    .storage_buf(0, Qualifier::READ, "uint", "ssbo_list_ranking_links_in_[]")
    .storage_buf(1, Qualifier::WRITE, "uint", "ssbo_list_ranking_links_out_[]")
    .uniform_buf(0, "UBData_ListRanking", "ubo_list_ranking_splicing_")
    
    .local_group_size(GROUP_SIZE_BNPR_LIST_RANK_TEST)
    .compute_source("npr_strokegen_list_ranking_test_comp.glsl");

/** \} */
