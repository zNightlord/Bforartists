/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup bnpr
 */

#pragma once

#include "draw_pass.hh"
#include "draw_manager.hh"

#include "npr_strokegen_shader.hh"
#include "bnpr_shader_shared.hh"
#include "npr_strokegen_buffer_pool.hh"
#include "npr_strokegen_texture_pool.hh"
#include "strokegen_mesh_raster_pass.hh"

#include <random> 


namespace blender::npr::strokegen {
class Instance;

struct SurfaceDebugContext {
  enum DbgLineType {
    vnor = 0, // to be deprecated
    general = 2,
    edges = 1 // to be deprecated
  }; // match to shader defines
  bool dbg_lines; 
  bool dbg_vert_normal;
  bool dbg_vert_curv;
  bool dbg_vert_contour_grad; 
  bool dbg_edges; 

  float dbg_line_length;
};

class StrokegenMeshComputePass : public draw::PassSimple {
public:
  StrokegenMeshComputePass(const char *name = "unkown strokegen compute pass") : draw::PassSimple(name) {}
};

class StrokeGenPassModule // similar to "LineDrawingRenderPass"
{
public:
  int get_num_passes_extract_geom() { return curr_mesh_id_extract_geom + 1; }
  int get_num_passes_remeshed_surf_depth() { return curr_mesh_id_surf_depth + 1; }

 private:
  /** Compute Passes */
  draw::PassSimple pass_comp_test = {"Strokegen Compute Test"};

  int curr_mesh_id_extract_geom; 
  std::array<StrokegenMeshComputePass, 1024> pass_extract_geom_arr;
  StrokegenMeshComputePass &pass_extract_geom() { return pass_extract_geom_arr[curr_mesh_id_extract_geom]; }
  draw::PassSimple pass_process_contours = {"StrokeGen Process Contours"}; 
  draw::PassSimple pass_compress_contour_pixels = {"Generate Contour Pixel Mask"}; 

  draw::PassSimple pass_scan_test = {"Bnpr GPU Blelloch Scan Test"};
  draw::PassSimple pass_segscan_test = {"Bnpr GPU Blelloch SegScan Test"};
  draw::PassSimple pass_conv1d_test = {"Test GPU 1d conv on circular segments"};
  draw::PassSimple pass_listranking_test = {"List ranking test"};

  /** Draw Passes */
  StrokegenMeshRasterPass pass_draw_contour_edges = {"Draw Contour Edges"}; // Inherited from draw::PassMain
  StrokegenMeshRasterPass pass_draw_contour_2d_samples = {"Draw 2D Contour Curves"}; // Inherited from draw::PassMain
  int curr_mesh_id_surf_depth; 
  std::array<StrokegenMeshRasterPass, 1024> pass_draw_remeshed_surface_depth_;  
  StrokegenMeshRasterPass pass_draw_debug_lines_ = {"Draw Debug Lines"}; // Inherited from draw::PassMain 

  /** Dependent Modules */
  StrokeGenShaderModule &shaders_;
  GPUBufferPoolModule   &buffers_;
  GPUTexturePoolModule  &textures_;

  /** Sync states */
  bool first_frame; 
  int strokegen_frame_id;
  int strokegen_obj_id; 


public:
  StrokeGenPassModule(
      StrokeGenShaderModule &strokegen_shaders,
      GPUBufferPoolModule &strokegen_buffers,
      GPUTexturePoolModule &strokegen_textures
      )
    : shaders_(strokegen_shaders),
      buffers_(strokegen_buffers),
      textures_(strokegen_textures),
      first_frame(true), 
      strokegen_frame_id(0)
  {
    ;
  }

  ~StrokeGenPassModule() {}



  /* -------------------------------------------------------------------- */
  /** \name Passes
     * \{ */
  enum eType
  {
    SCAN_TEST = 0,
    SEGSCAN_TEST,
    SEGLOOPCONV_TEST,
    LIST_RANKING_TEST,

    GEOM_EXTRACTION,
    CONTOUR_PROCESS, 
    INDIRECT_DRAW_CONTOUR_EDGES,
    INDIRECT_DRAW_CONTOUR_2D_SAMPLES,
    INDIRECT_DRAW_REMESHED_DEPTH,
    COMPRESS_CONTOUR_PIXELS,

    INDIRECT_DRAW_DBG_VNOR,

  };

  PassSimple& get_compute_pass(eType passType, int pass_id = 0);
  PassMain &get_render_pass(eType passType, int pass_id = 0);
  void prepare_validation_passes(int frame_counter);
  void init_surface_depth_passes();
  /** \} */



  /* -------------------------------------------------------------------- */
  /** \name Sync with Draw Module
     * \{ */
  void on_begin_sync(int frame_counter);
  void on_end_sync(); 
  /** \} */



  SurfaceDebugContext surf_dbg_ctx; 


  /* -------------------------------------------------------------------- */
  /** \name Rebuild Render Passes
     * \{ */
  bool boostrap_before_extract_first_batch;
  int num_total_mesh_verts;
  int get_vtxid_offset() { return 0; }
  int num_total_mesh_edges;
  int get_edgeid_offset() { return 0; }
  int num_total_mesh_tris;

  bool has_strokegen_enabled_mesh; 
  void init_mesh_extraction_passes();
  void append_per_mesh_pass(
      Object* ob,
      gpu::Batch* gpu_batch_line_adj,
      gpu::Batch* gpu_batch_surf,
      ResourceHandle& rsc_handle,
      const DRWView* drw_view
    );


  // ---------------------------------------------------------------------------
  void append_subpass_cpy_vbo(gpu::Batch* gpu_batch_surf, int batch_resource_index, int num_verts);
  void append_subpass_cpy_line_adj_ibo(gpu::Batch* gpu_batch_line_adj,
                                         gpu::GPUIndexBufType ib_type,
                                         int num_added_edges);

  // ---------------------------------------------------------------------------
  void append_subpass_meshing_merge_verts(int num_verts_in, bool debug = false);
  void append_subpass_meshing_wedge_adjacency(int num_edges_in, int num_verts_in, bool debug = false);


  // ---------------------------------------------------------------------------
  struct EdgeFloodingOptions {
   bool compact_edges;
   bool output_selected_to_edge; // has effect only when .compact_edges==true
   bool output_edge_to_selected; // has effect only when .compact_edges==true
  };
  void append_subpass_select_remeshed_edges(int num_edges, int num_verts, EdgeFloodingOptions options);
  void append_subpass_mark_selection_border_edges(int num_edges, int num_verts);


  enum GPUMeshQuadricFilter {
    GeomNormalPlane = 0,
    ViewNormalPlane = 1,
    ViewTangentPlane = 2,
    ViewBinormalPlane = 3
  };
  enum ContourType {
    Interpolated = 1,
    Raw = 2
  };
  struct GPURemeshingParameters {
    int num_filtering_iters; // deprecate
    int max_num_remesh_dbg_iters;
     
    int num_quadric_diffusion_iters; // deprecate
    float quadric_deviation; // deprecate
    float geodist_deviation; // deprecate
    int seconds_sync_view_mat;
    int num_vtx_smooth_iters; 

    int num_edge_flooding_iters;
    float position_regularization_scale;
    GPUMeshQuadricFilter alternate_filter_0;
    GPUMeshQuadricFilter alternate_filter_1;
    int edge_visualize_mode;
    int contour_mode; // 0: no visuals, 1: interpo contour, 2: raw contour
    bool visualize_contour_edges;
    int iters_test_subdiv; 

    float remeshing_targ_edge_len;
    int remeshing_split_iters; 
    int remeshing_collapse_iters;
    int remeshing_flip_iters;
    int remeshing_iters;
    int remeshing_delaunay_flip_iters;

    int subdiv_type;
    bool subdiv_use_crease;

    bool denoise_cusp_segmentation;
    bool cusp_eval_opti;
    float visibility_thresh;

    int dbg_matching_line_mode; 
    int dbg_history_trace_steps;
    int dbg_history_trace_passes;
    int dbg_ndv_grad_mode;
  } meshing_params;


  // ---------------------------------------------------------------------------
  struct SelectVertsFromEdgesContext {
    // each vertex has 4 selection flag/slots, -1 if non-active
    int4 active_selection_slots; 

    bool expand_selection; // expand selection to any vertex connected to selected edge
    int4 input_selection_slots_for_expansion;
    int4 output_selection_slots_for_expansion;  

    bool compact_verts; // output indexing table for verts selected in active slots
    bool compact_all_slots_selected; // compact when all slots are selected
    int4 selection_slots_for_compation; 
  };
  void append_subpass_select_verts_from_selected_edges(SelectVertsFromEdgesContext ctx, int num_edges, int num_verts);


  enum EdgeSplitMode {
    // Match to shader macros
    // #define EDGE_SPLIT_LONG_EDGES 0u
    LongEdge = 0u,
    // #define EDGE_SPLIT_CONTOUR_EDGES 1u
    InterpContour = 1u, 
    // #define EDGE_SPLIT_LOOP_SUBDIV 2u
    LoopSubdivSplit = 2u, 
  }; 
  enum EdgeFlipOptiGoal {
    // Match to shader macros
    // #define EDGE_FLIP_OPTI_VALENCE 0u
    Valence = 0,
    // #define EDGE_FLIP_OPTI_DELAUNAY 1u
    Delaunay = 1,
    // #define EDGE_FLIP_OPTI_SQRT3_SUBDIV 2u
    SqrtSubdiv = 2,
    // #define EDGE_FLIP_OPTI_SQRT3_SUBDIV 3u
    LoopSubdivFlip = 3, 
  };
  void append_subpass_split_edges(EdgeSplitMode mode, int iter_split, int num_edges, int num_verts);
  void append_subpass_collapse_edges(int iter_remesh, int iter_collapse, int num_edges, int num_verts);
  void append_subpass_flip_edges(EdgeFlipOptiGoal opti_goal, int iter_flip, int num_edges, int num_verts);
  void append_subpass_split_faces(int iter_split, int num_edges, int num_verts); 
  void append_subpass_fill_dispatch_args_remeshed_edges_(int num_static_edges, bool only_selected_edges);
  void append_subpass_fill_dispatch_args_remeshed_verts_(int num_static_verts, bool only_selected_elems_);


  // ---------------------------------------------------------------------------
  // Surface Analysis
  struct SurfaceAnalysisContext {

    // Vertex attributes
    bool order_0_only_selected;  // only calculate order-0 attrs on selected elems
    bool calc_vert_normal;
    bool output_vertex_facing_flag; // mark facing by dot(n, view_dir) 
    GPUStorageBuf *ssbo_vnor_; 
    bool calc_vert_voronoi_area;
    GPUStorageBuf *ssbo_varea_;
    bool calc_vert_topo_flags;
    // note: due to shader variant issues (I'm lazy to generate all the variants),
    // calc_vert_topo_flags conflicts against calc_vert_voronoi_area

    bool order_1_only_selected;
    bool order_1_only_contour; // make sure order-0 attrs for neighboring verts are properly computed
    bool calc_vert_curvature;
    enum CurvatureEstimator { Rusinkiewicz = 0, Jacques } curvature_estimator;
    bool output_curvature_tensors;
    bool output_maxcurv_with_cusp_function; 
    GPUStorageBuf *ssbo_vcurv_tensor_; 
    GPUStorageBuf *ssbo_vcurv_pdirs_k1k2_; 
    GPUStorageBuf *ssbo_edge_vtensors_; // temp buffer holding partial tensors

    bool calc_vert_contour_grad; 
    GPUStorageBuf *ssbo_vgrad_contour_; 

    // Edge attrs
    bool calc_feature_edges; 
    bool only_selected_edges; 


    SurfaceAnalysisContext()
      : order_0_only_selected(false),
        calc_vert_normal(false),
        output_vertex_facing_flag(false), 
        ssbo_vnor_(nullptr),
        calc_vert_voronoi_area(false),
        ssbo_varea_(nullptr),
        calc_vert_topo_flags(false), 
        order_1_only_selected(false),
        order_1_only_contour(false), 
        calc_vert_curvature(false),
        curvature_estimator(Jacques),
        output_curvature_tensors(false),
        output_maxcurv_with_cusp_function(false),
        ssbo_vcurv_tensor_(nullptr),
        ssbo_vcurv_pdirs_k1k2_(nullptr),
        ssbo_edge_vtensors_(nullptr),
        calc_vert_contour_grad(false),
        ssbo_vgrad_contour_(nullptr), 
        calc_feature_edges(false),
        only_selected_edges(false)
    {
    }

    void set_calc_vert_normal(bool val, bool output_facing_flags)
    {
      calc_vert_normal = val;
      output_vertex_facing_flag = output_facing_flags; 
    }
    void set_calc_vert_voronoi_area(bool val) { calc_vert_voronoi_area = val; }
    void set_calc_vert_topo_flags(bool val) { calc_vert_topo_flags = val; }

    void set_calc_vert_contour_grad(bool val)
    {
      calc_vert_contour_grad = val; 
    }

    void set_calc_vert_curvature(bool val, CurvatureEstimator algo, bool output_tensors, bool output_cusp_and_maxcurv) 
    { 
      calc_vert_curvature = val;
      curvature_estimator = algo; 
      if (calc_vert_curvature) {
        set_calc_vert_normal(true, output_vertex_facing_flag); 
        set_calc_vert_voronoi_area(true);
        output_curvature_tensors = output_tensors;
        output_maxcurv_with_cusp_function = output_cusp_and_maxcurv; 
      } 
    }

    void set_calc_feature_edges(bool val, bool only_selected)
    {
      calc_feature_edges = true;
      only_selected_edges = only_selected; 
    }
  };
  void GetSurfaceAnalysisContext_InitPass(SurfaceAnalysisContext &surf_analysis_ctx) const;
  void GetSurfaceAnalysisContext_VertexRelocationPass(SurfaceAnalysisContext &surf_analysis_ctx) const;
  void GetSurfaceAnalysisContext_ContourInsertionPass(SurfaceAnalysisContext &surf_analysis_ctx) const;
  void GetSurfaceAnalysisContext_CuspDetectionPass(SurfaceAnalysisContext &surf_analysis_ctx) const;
  void GetSurfaceAnalysisContext_CurvatureForAdaptiveRemeshing(SurfaceAnalysisContext &surf_analysis_ctx) const;

  void append_subpasses_estimate_curvature_for_adaptive_remeshing(ResourceHandle& rsc_handle,
                                          int num_edges,
                                          int num_verts, bool output_dbg_lines = false);
  void append_subpass_init_temporal_records(int num_edges,
                                            SurfaceAnalysisContext surf_analysis_ctx_contour);
  void append_subpass_reconstruct_last_frame_temporal_records();
  void append_subpasses_sqrt_subdiv(int num_edges, int num_verts);
  void bind_rsc_for_loop_subd_tree_processing(StrokegenMeshComputePass::PassBase<DrawCommandBuf>& sub, int num_edges, int& out_ssbo_offset);
  void append_subpass_build_loop_subd_tree_upwards_for_face_edges(int num_edges);
  void append_subpass_build_loop_subd_tree_downwards_init(int num_edges);
  void append_subpass_build_loop_subd_tree_downwards_(int num_edges);
  void append_subpasses_loop_subdiv(int num_edges, int num_verts);

  void append_subpass_surf_geom_analysis(
      ResourceHandle& rsc_handle, int num_verts,
      int num_edges,
      const SurfaceAnalysisContext& ctx, const SurfaceDebugContext& dbg_options
      ); 

  void rebuild_pass_dbg_geom_drawcall(SurfaceDebugContext dbg_ctx);


  // ---------------------------------------------------------------------------
  // Mesh Filtering
  void bind_rsc_for_bnpr_meshing_surf_filtering_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf>& sub,
      int num_verts,
      int num_edges,
      int& out_num_ssbo);

  enum VertexRelocationMode {
    TangentialSmoothing = 0, 
    QuadricFiltering = 1,
    Sqrt3SubdivSmoothCache = 2, 
    Sqrt3SubdivSmoothApply = 3,
    LoopSubdivSmoothCache = 4,
    LoopSubdivSmoothApply = 5,
  };
  void append_subpass_vertex_relocation(VertexRelocationMode mode,
                                        int num_edges,
                                        int num_verts,
                                        int num_vnor_filter_iters,
                                        bool only_selected_verts);
  void append_subpass_vertex_curv_smoothing(
    int iter_smooth, int num_verts, int num_edges,
    bool only_selected_verts = false, bool output_remesh_len = false
  );


  // ---------------------------------------------------------------------------
  void append_subpass_extract_contour_edges(gpu::Batch* gpu_batch_line_adj,
                                            ResourceHandle& rsc_handle,
                                            gpu::Batch* edge_batch,
                                            int num_edges,
                                            gpu::GPUIndexBufType ib_type,
                                            int edge_visualize_mode,
                                            int contour_visualize_mode);

  void rebuild_pass_process_contours();
  void append_subpass_fill_dispatch_args_contour_edges(PassSimple& pass, bool all_contour_edges);
  void append_subpass_fill_dispatch_args_contour_verts(PassSimple &pass);
  void append_subpass_fill_dispatch_args_contour_frags(PassSimple &pass, bool all_contour_frags);
  void append_subpass_fill_dispatch_args_contour_2d_samples(PassSimple &pass);
  void append_subpass_fill_dispatch_args_temporal_records_(int obj_id, int rec_frame_id);
  void append_subpass_setup_contour_edge_data();

  // ---------------------------------------------------------------------------
  void append_pass_remeshed_surface_depth_drawcall();

  // ---------------------------------------------------------------------------
  void append_subpass_contour_edges_soft_rasterization();
  void append_subpass_visibility_split_contour_edges(); // split edges based on soft rasterization

  // ---------------------------------------------------------------------------
  void append_subpass_fill_contour_list_ranking_inputs(); 
  
  // ---------------------------------------------------------------------------
  void append_subpass_serialize_contour_edges();
  void append_subpass_contour_segmentation();

  void bind_rsc_for_contour_2d_resample_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf>& sub, float2 screen_res, float pcs_sample_rate, int& out_ssbo_offset);
  void bind_rsc_for_contour_2d_sample_evaluation_(
      draw::detail::Pass<DrawCommandBuf>::PassBase<DrawCommandBuf> &sub, float2 screen_res, float pcs_sample_rate, int &out_ssbo_offset); 
  void append_subpass_contour_arclen_parameterization(float2 screen_res, float sample_rate);
  void append_subpass_contour_generate_2d_samples(float2 screen_res, float sample_rate); 
  void append_subpass_2d_sample_segmentation(float2 screen_res, float sample_rate, bool is_segmentation_by_curve_pass);
  void append_subpass_calc_contour_edges_draw_data(); 

  // ---------------------------------------------------------------------------
  void rebuild_pass_contour_edge_drawcall();
  void rebuild_pass_contour_2d_samples_drawcall();
  void rebuild_pass_compress_contour_pixels(bool debug = false);



  // ---------------------------------------------------------------------------
  bool test_scan;
  struct ScanSettings {
    bool is_validation_shader;
    int frame_counter;

    bool use_indirect_dispatch;
    GPUStorageBuf *ssbo_scan_infos_;

    GPUStorageBuf *ssbo_in_scan_data_;
    GPUStorageBuf *ssbo_out_scan_data_;
    GPUStorageBuf *ssbo_scan_block_sum_;

    eShaderType shader_upsweep;
    eShaderType shader_aggregate;
    eShaderType shader_dwsweep;
  };
  void append_subpass_scan(ScanSettings scan_settings, PassSimple& pass);
  void append_subpass_segscan(ScanSettings scan_settings,
                                 PassSimple& pass);
  struct SegLoopConv1DSettings {
    bool is_validation_shader;

    bool use_indirect_dispatch;
    GPUStorageBuf *ssbo_segloopconv1d_info_;

    GPUStorageBuf *ssbo_in_segloopconv1d_data_;
    GPUStorageBuf *ssbo_out_segloopconv1d_data_;
    GPUStorageBuf *ssbo_segloopconv1d_patch_table_;

    eShaderType shader_build_patch_table;
    bool lazy_dispatch; // skip shader_build_patch_table and reuse ssbo_segloopconv1d_patch_table_ 
    eShaderType shader_convolution;
  };
  bool test_segloopconv; 
  void append_subpass_segloopconv1d(SegLoopConv1DSettings settings, PassSimple& pass);

  bool test_list_ranking; 
  bool test_looped_pass_list_ranking;
  enum ListRankingPassUsage { // match to shader defines
    TestListRanking = 0,      // OUTPUT_PASS_TYPE__TEST
    ContourEdgeLinking = 1,   // OUTPUT_PASS_TYPE__CONTOUR_EDGE_LINKING
  };
  void append_subpass_list_ranking(ListRankingPassUsage passType, PassSimple& pass_listranking, bool looped_pass_list_ranking);
  void rebuild_pass_list_ranking_pointer_jumping(PassSimple &pass_listranking,
                                                 int num_splice_iters,
                                                 const int num_jump_iters,
                                                 int jumping_info_offset,
                                                 bool loop_breaking_pass,
                                                 bool loop_ranking_pass);
  void rebuild_pass_list_ranking_fill_args(PassSimple& pass_listranking, bool per_anchor, bool per_spliced, int splicing_or_relinking_iter, int group_size_x, bool custom_pass);
  void print_list_ranking_nodes(int head_node_id, uint* computed_ranks, uint* computed_topo, uint* computed_links) const;
  /** \} */




  /* -------------------------------------------------------------------- */
  /** \name Scan Tester
   * \{ */
  template<typename T>
  void validate_pass_scan_test(bool (*equals)(const T &, const T &));
  template<typename T>
  bool validate_inter_block_exclusive_scan(
      const T *bufferInputVals,
      const T *bufferPrefixSum,
      bool (*equals)(const T &, const T &),
      uint num_scan_items,
      uint blk_size,
      uint numBlocks
      );
  template<typename T>
  bool validate_exclusive_scan(
      const T *bufferInputVals,
      const T *bufferPrefixSum,
      bool (*equals)(const T &, const T &),
      uint num_scan_items
      );
  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name Segment Scan Validation
   * \{ */
  template<typename T, typename TD>
  void validate_segscan(
    T (*func_decode)(TD), 
    bool (*func_equals)(const T&, const T&),
    uint (*func_get_hf)(const T&),
    T (*func_scan_op)(const T&, const T&),
    T zero_val,
    bool inclusive = false
  );
  template<typename T, typename TD>
  bool validate_segscan_internal(
    const TD *input,
    const TD *output,
    const int blk_size, const int num_elems,
    T (*func_decode)(TD), 
    uint (*func_get_hf)(const T&),
    bool (*func_equals)(const T &, const T &),
    T (*func_scan_op)(const T&, const T&),
    T zero_value,
    bool inclusive = false
  );
  /** \} */


  /* -------------------------------------------------------------------- */
  /** \name 1D Segmented Looped Convolution Validation
   * \{ */
  template<typename T>
  void validate_segloopconv1d(
    bool (*func_equals)(const T&, const T&),
    T (*func_conv_op)(const T&, const T&)
  );

  template<typename T>
  bool validate_segloopconv1d_internal(
    const T *input,
    const T *output,
    const int *seg_topo,
    int num_elems,
    int conv_radius,
    bool (*func_equals)(const T &, const T &),
    T (*func_conv_op)(const T&, const T&)
  );
  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name List Ranking Validation
   * \{ */
  bool validate_list_ranking() const;

  /** \} */


private:

};


template<typename T>
void StrokeGenPassModule::validate_pass_scan_test(bool (*equals)(const T &, const T &))
{
  SSBO_BnprScanData &buf_scan_inputs = buffers_.ssbo_in_scan_data_;
  buf_scan_inputs.read();
  T *data_scan_inputs = reinterpret_cast<T *>(buf_scan_inputs.data());

  SSBO_BnprScanData &buf_scan_output = buffers_.ssbo_out_scan_data_;
  buf_scan_output.read();
  T *data_scan_output = reinterpret_cast<T *>(buf_scan_output.data());

  bool valid_inter_block_scan = StrokeGenPassModule::validate_inter_block_exclusive_scan<T>(
      data_scan_inputs,
      data_scan_output,
      equals,

      buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items,
      GROUP_SIZE_BNPR_SCAN_SWEEP * 2u,
      buffers_.ubo_bnpr_tree_scan_infos_.num_thread_groups
      );
  if (!valid_inter_block_scan)
    fprintf(stderr, "bnpr: error: INTER-BLOCK scan test failed");

  bool valid_global_scan = StrokeGenPassModule::validate_exclusive_scan<T>(
      data_scan_inputs,
      data_scan_output,
      equals,
      buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items
      );
  if (!valid_global_scan)
    fprintf(stderr, "bnpr: error: GLOBAL scan test failed");
}


template<typename T>
bool StrokeGenPassModule::validate_inter_block_exclusive_scan(
    const T *const bufferInputVals,
    const T *const bufferPrefixSum,
    bool (*equals)(const T &, const T &),
    uint num_scan_items,
    uint blk_size,
    uint numBlocks
)
{
  Vector<int> failedElems(0);

  for (uint blk_id = 0u; blk_id < numBlocks; blk_id++) {
    for (uint blk_offset = 0u; blk_offset < blk_size - 1u; blk_offset++) {
      // index might go out of bound
      uint index = blk_id * blk_size + blk_offset;
      if (index >= num_scan_items - 2u)
        break;

      if (false == equals(
              bufferPrefixSum[index] + bufferInputVals[index],
              bufferPrefixSum[index + 1]
              )) {
        failedElems.append(index);
      }
    }
  }

  if (!failedElems.is_empty()) {
    return false;
  }

  return true;
}


template<typename T>
bool StrokeGenPassModule::validate_exclusive_scan(const T *bufferInputVals,
                                                  const T *bufferPrefixSum,
                                                  bool (* equals)(const T &, const T &),
                                                  uint num_scan_items)
{
  return validate_inter_block_exclusive_scan(
      bufferInputVals,
      bufferPrefixSum,
      equals,
      num_scan_items,
      num_scan_items,
      1
      );
}

template<typename T, typename TD>
void StrokeGenPassModule::validate_segscan(
  T (*func_decode)(TD), 
  bool (*func_equals)(const T&, const T&),
  uint (*func_get_hf)(const T&),
  T (*func_scan_op)(const T&, const T&),
  T zero_val,
  bool inclusive)
{
  SSBO_BnprScanData &buf_segscan_inputs = buffers_.ssbo_in_scan_data_;
  buf_segscan_inputs.read();
  TD *data_segscan_inputs = reinterpret_cast<TD *>(buf_segscan_inputs.data());

  SSBO_BnprScanData &buf_segscan_output = buffers_.ssbo_out_scan_data_;
  buf_segscan_output.read();
  TD *data_segscan_output = reinterpret_cast<TD *>(buf_segscan_output.data());

  bool succ = validate_segscan_internal<T, TD>(
    data_segscan_inputs,
    data_segscan_output,
    GROUP_SIZE_BNPR_SCAN_SWEEP * 2u,
    buffers_.ubo_bnpr_tree_scan_infos_.num_scan_items,
    func_decode, func_get_hf, func_equals, func_scan_op, zero_val,
    inclusive
  );

  if (!succ)
    fprintf(stderr, "strokegen error: segment scan test failed");
}

template<typename T, typename TD>
bool StrokeGenPassModule::validate_segscan_internal(
  const TD *input,
  const TD *output,
  const int blk_size,
  const int num_elems,
  T (*func_decode)(TD), 
  uint (*func_get_hf)(const T &),
  bool (*func_equals)(const T &, const T &),
  T (*func_scan_op)(const T&, const T&),
  T zero_value,
  const bool inclusive)
{
  T maxErrorFound = zero_value;
  T currScanValue = zero_value;

  for (int i = 0; i < num_elems; i++)
  {
    T input_t = func_decode(input[i]); 
    T output_t = func_decode(output[i]);

    if (func_get_hf(input_t)) {
      currScanValue = zero_value;
    }

    if (inclusive) {
      // Inclusive scan
      currScanValue = func_scan_op(currScanValue, input_t);
    }

    if (false == func_equals(currScanValue, output_t)) {
      return false;
    }

    if (!inclusive) {
      // Exclusive scan
      currScanValue = func_scan_op(currScanValue, input_t);
    }
  }

  return true;
}

template <typename T>
void StrokeGenPassModule::validate_segloopconv1d(
  bool(* func_equals)(const T&, const T&),
  T(* func_conv_op)(const T&, const T&))
{
  SSBO_SegLoopConv1DData &buf_conv_inputs = buffers_.ssbo_in_segloopconv1d_data_;
  buf_conv_inputs.read();
  T *data_conv_inputs = reinterpret_cast<T *>(buf_conv_inputs.data());

  SSBO_SegLoopConv1DData &buf_conv_output = buffers_.ssbo_out_segloopconv1d_data_;
  buf_conv_output.read();
  T *data_conv_output = reinterpret_cast<T *>(buf_conv_output.data());

  SSBO_SegLoopConvDebugData &buf_dbg_data = buffers_.ssbo_debug_segloopconv1d_data_;
  buf_dbg_data.read();
  int* data_conv_debug = reinterpret_cast<int *>(buf_dbg_data.data());

  bool succ = validate_segloopconv1d_internal<T>(
    data_conv_inputs, data_conv_output, data_conv_debug,
    buffers_.ubo_segloopconv1d_.num_conv_items,
    NPR_TEST_SEGLOOPCONV1D_CONV_RADIUS,
    func_equals, func_conv_op
  );

  if (!succ)
    fprintf(stderr, "strokegen error: segloopconv1d test failed");
}

template <typename T>
bool StrokeGenPassModule::validate_segloopconv1d_internal(
  const T* input,
  const T* output,
  const int* seg_topo,
  int num_elems,
  int conv_radius,
  bool(* func_equals)(const T&, const T&),
  T(* func_conv_op)(const T&, const T&)
)
{
  T currOutput;

  for (int i = 0; i < num_elems; ++i)
  {
    const int seg_head = seg_topo[i * 4];
    const int seg_tail = seg_topo[i * 4 + 1];

    currOutput = input[i];
    { // convolution
      int neigh_id = i;
      for (int d = 1; d <= conv_radius; ++d)
      {
        neigh_id--;
        if (neigh_id < seg_head) neigh_id = seg_tail;

        currOutput = func_conv_op(currOutput, input[neigh_id]);
      }

      neigh_id = i;
      for (int d = 1; d <= conv_radius; ++d)
      {
        neigh_id++;
        if (neigh_id > seg_tail) neigh_id = seg_head;

        currOutput = func_conv_op(currOutput, input[neigh_id]);
      }
    }

    if (!func_equals(currOutput, output[i]))
      return false;
  }

  return true;
}
}
