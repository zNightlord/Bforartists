/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <optional>
#include <utility>

#include <fmt/format.h>

#include "BLI_array.hh"
#include "BLI_index_mask.hh"
#include "BLI_math_base.h"
#include "BLI_math_base.hh"
#include "BLI_math_geom.h"
#include "BLI_math_matrix.hh"
#include "BLI_math_numbers.hh"
#include "BLI_math_vector.h"
#include "BLI_math_vector.hh"
#include "BLI_multi_value_map.hh"
#include "BLI_offset_indices.hh"
#include "BLI_span.hh"
#include "BLI_task.hh"
#include "BLI_vector.hh"

#include "BKE_curveprofile.h"
#include "BKE_mesh.hh"
#include "BKE_mesh_mapping.hh"

#include "GEO_mesh_bevel.hh"

#define DEBUG_BEVEL
#ifdef DEBUG_BEVEL
#  include "draw_debug.hh"
#endif

namespace blender::geometry {

/**
 * The un-transformed 2D storage of profile vertex locations. Also, for non-custom profiles
 * this serves as a cache for the results of the expensive calculation of u parameter values to
 * get even spacing on superellipse for current BevelParams seg and pro_super_r.
 */
struct ProfileSpacing {
  /** The profile's seg+1 x values. */
  Array<double> xvals;
  /** The profile's seg+1 y values. */
  Array<double> yvals;
  /** The profile's segments_power_2+1 x values. */
  Array<double> xvals_2;
  /** The profile's segments_power_2+1 y values, (seg_2 = power of 2 >= seg). */
  Array<double> yvals_2;
  /** The power of two greater than or equal to the number of segments. */
  int segments_power_2;
  /** How far "out" the profile is, used at the start of subdivision. */
  float fullness;
};

/** Different kinds of patterns for mesh fragments that will replace bevel-involved vertices. */
enum class MeshKind {
  None,         /* No mesh. */
  Line,         /* A line of connected vertices. */
  TerminalPoly, /* A simple polygon, one side is a terminal beveled edge. */
  Adj,          /* "Adjacent edges" mesh pattern. */
  TriFan,       /* A simple polygon - fan filled. */
  Cutoff,       /* A triangulated face at the end of each profile. */
};

class MeshPattern {
 public:
  MeshKind kind;
  int num_anchors;
  int num_segs;

  MeshPattern() : kind(MeshKind::None), num_anchors(0), num_segs(0) {}
  MeshPattern(MeshKind kind, int num_anchors, int num_segs)
      : kind(kind), num_anchors(num_anchors), num_segs(num_segs)
  {
  }

  /** Return a 4-tuple with the number of vertices, edges, faces, corners.  */
  int4 num_elements() const;

  /** Return the vertex index of the give anchor vertex. */
  int anchor_vert(int anchor_index) const;
};

/** Helper for keeping track of angle kind. */
enum class AngleKind {
  /** Angle less than 180 degrees. */
  Smaller = -1,
  /** 180 degree angle. */
  Straight = 0,
  /** Angle greater than 180 degrees. */
  Larger = 1,
};

/** Return what range[i+1] would when wrapping around at end. */
static inline int index_range_next(const IndexRange range, const int index)
{
  return (index >= range.size() - 1) ? range.first() : range[index + 1];
}

/** Return what range[i-1] would when wraping around at the beginning. */
static inline int index_range_prev(const IndexRange range, const int index)
{
  return (index <= 0) ? range.last() : range[index - 1];
}

/** Mesh information, including some topology information. */
class MeshInfo {
 public:
  const Mesh &mesh;

  MeshInfo(const Mesh &mesh);

  /** For each Mesh vert v, vert_edges[v] is a Span of edge indices incident on v. */
  GroupedSpan<int> vert_edges() const
  {
    return vert_edges_;
  }

  /** For each Mesh edge e, edge_faces[e] is a Span of face indices incident on e. */
  GroupedSpan<int> edge_faces() const
  {
    return edge_faces_;
  };

 private:
  IndexMaskMemory memory_;
  Array<int> vert_to_edge_offsets_;
  Array<int> vert_to_edge_indices_;
  Array<int> edge_to_face_map_offsets_;
  Array<int> edge_to_face_map_indices_;
  GroupedSpan<int> vert_edges_;
  GroupedSpan<int> edge_faces_;
};

/** Bevel parameters and state.
 * A note on some terminology used in the state variables:
 *
 * A "beveled vertex" is, when we are vertex-beveling, a vertex selected by the params
 * and is eligible for beveling (attached to at least two edges).
 * A "beveled edge" is, when we are edge-beveling, an edge selected by the params
 * and is eligibile for beveling (attached to exactly two faces).
 *
 * A bevvert is a "bevel involved vertex". This means:
 *    If we are doing a vertex-only bevel, then it is a beveled vertex.
 *    If we are doing an edge bevel, then it is attached to a beveled edge (not, not a
 *    bevel-involved edge, just a beveled edge).`
 *
 * A bevedge is a "bevel involved edge". This means:
 *    If we are doing a vertex-only bevel, then it is attached to a beveled vertex.
 *    If we are doing an edge bevel, then it is a beveled edge, or another edge
 *    attached to a bevvert (i.e., attached to a vert where a beveled edge is attached).
 *    We can tell which edges are actually beveled because they will have a nonzero
 *    bevel weight (which is a property of bevedges).
 *
 * There is a separate indexing scheme for bevverts and bevedges, each going from zero
 * to the total number of them minus one. There are member functions to go back and
 * forth between the two indexing systems.
 *
 * There are also newverts, newedges, newfaces, and newcorners, all of which will be
 * added to the argument Mesh as part of making the new beveled mesh.
 * During construction, these all have indices starting at zero and going to the number
 * of such elements, minus 1. When we need to mix new and old elements, we offset
 * the new element indices by the size of the corresponding elements in the source Mesh.
 *
 * Both bevverts and bevedges will have their corresponding verts and edges replaced by
 * new mesh fragments, where each fragment is made up of newvertes, newedges, etc.
 * The associated new elements form a particular patterns that, in the case of a bevvert,
 * can have several kinds (see #MeshKind, above).
 */
class BevelState {
 public:
  /** Input parameters. */
  BevelParameters params;
  /** Input Mesh with some topology information.. */
  MeshInfo mesh_info;
  /** Should we attempt to adjust offsets? */
  bool offset_adjust = false;
  /** Number of bevel-involved verts. */
  int bevverts_num = 0;
  /** Number of bevel-involved edges. */
  int bevedges_num = 0;
  /** Number of bevel-involved faces. */
  int bevfaces_num = 0;
  /** Number of new verts. */
  int newverts_num = 0;
  /** Number of new edges. */
  int newedges_num = 0;
  /** Number of new faces. */
  int newfaces_num = 0;
  /** Number of new corners. */
  int newcorners_num = 0;
  /** Parameter values for evenly spaced profile points. */
  ProfileSpacing pro_spacing;
  /** Parameter values for evenly spaced profile points for the miter profiles. */
  ProfileSpacing pro_spacing_miter;
  /** Profile shape parameter converted to a superellipse exponent. */
  float pro_super_r = 2;

  /** Construct initial state from input. */
  BevelState(const Mesh &src_mesh,
             const IndexMask &selection,
             const BevelParameters &bevel_params);

  void initialize_profile_data();

  void order_bevedges();

  void set_bevvert_mesh_topology();

  void set_bevedge_widths();

  void build_vertex_meshes();

  /** Indices of Mesh vertices that have corresponding bevverts.
   */
  const IndexMask &bevverts_mask() const
  {
    return bevverts_mask_;
  }

  /** Indices of Mesh edges that have corresponding bevedges. */
  const IndexMask &bevedges_mask() const
  {
    return bevedges_mask_;
  }

  /** Map from Mesh vertex index to bevverts index. */
  Span<int> vert_bevverts() const
  {
    return vert_bevverts_;
  }

  /** Map from Mesh edge index to bevedges index. */
  Span<int> edge_bevedges() const
  {
    return edge_bevedges_;
  }

  /* The following are indexed by a bevvert index, 0 to bevverts_num - 1. */

  /** The Mesh vert associated with a given bevvert. */
  Span<int> bevvert_mesh_verts() const
  {
    return bevvert_mesh_verts_;
  }

  /** The weight that mulitplies the vertex bevel amount for a given bevvert. */
  Span<float> bevvert_weights() const
  {
    return bevvert_weights_;
  }

  /** The bevedges associated with a given bevvert.
   * The are sorted to be CCW around the vertex, as much as possible. */
  GroupedSpan<int> bevvert_bevedges() const
  {
    return bevvert_bevedges_;
  }

  /** Corresponding to each bevedge in bevvert_bevvedges, this gives
   * the common face between that edge and the next one cyclically,
   * if there is such an edge, else it will be -1. */
  GroupedSpan<int> bevvert_faces() const
  {
    return bevvert_faces_;
  }

  /** What pattern of Mesh will replace the given bevvert. */
  Span<MeshPattern> bevvert_meshpatterns() const
  {
    return bevvert_meshpatterns_;
  }

  /** The range of newvert indices used in the mesh pattern for the given bevvert. */
  OffsetIndices<int> bevvert_newverts() const
  {
    return bevvert_newverts_;
  }

  /** The range of newedge indices used in the mesh pattern for the given bevvert. */
  OffsetIndices<int> bevvert_newedges() const
  {
    return bevvert_newedges_;
  }

  /** The range of newface indices used in the mesh pattern for the given bevvert. */
  OffsetIndices<int> bevvert_newfaces() const
  {
    return bevvert_newfaces_;
  }

  /* The following are indexed by a bevedge index, 0 to bevvedges_num - 1. */

  /** The Mesh edge associated with the given bevedge. */
  Span<int> bevedge_mesh_edges() const
  {
    return bevedge_mesh_edges_;
  }

  /** The weight that mulitplies the edge bevel amount for a given bevedge. */
  Span<float> bevedge_weights() const
  {
    return bevedge_weights_;
  }

  /** Left and right edge specs for each bevedge.
   * In order, the left and right sides (looking down edge) at the end of the origin end
   * of the bevedge, followed by the same for the destination end.
   * For edge beveling with type Offset, Width, or Depth: the width of each side of the centerline.
   * For edge beveling with type Percent or Absolute: the distance along the previous and next
   * edges, either absolutely or as a fraction.
   * For vertex beveling, the 0th spec is the distance to offset the on-edge vertex.
   */
  Span<float4> bevedge_widths() const
  {
    return bevedge_widths_;
  }

  /** Near and far end indices in the vmesh index space where the given bevedge is attached.
   * If the edge is beveled, this is where the left side of the beveled edge is attached. */
  Span<int2> bevedge_attach_verts() const
  {
    return bevedge_attach_verts_;
  }

  /** The newedges that are needed to form the edge polygons for the given bevvedge. */
  OffsetIndices<int> bevedge_newedges() const
  {
    return bevedge_newedges_;
  }

  /** The newfaces that are  the edge polygons for the given bevvedge. */
  OffsetIndices<int> bevedge_newfaces() const
  {
    return bevedge_newedges_;
  }

  /* The following are indexed by a newvert index, 0 to newverts_num - 1. */

  /** The cooredinates of the given newvert. */
  Span<float3> newvert_positions() const
  {
    return newvert_positions_;
  }

  /* The following are indexed by a newedge index, 0 to newedges_num - 1. */

  /** The start and end  indices of vertices for the given newedge. If the indices are in the
   * range of vertices in the source mesh, they come from that; else they are newvert indices
   * offset by the number of vertices in the source mesh. */
  Span<int2> newedge_vertpairs() const
  {
    return newedge_vertpairs_;
  }

  /* The following are indexed by a newface index, 0 to newfaces_num - 1. */

  /** Returns a contiguous chunk of face corners, represented as an #IndexRange, like a #Mesh
   * faces(). The corners are newcorner indices. */
  OffsetIndices<int> newface_faces() const
  {
    return newface_faces_;
  }

  /* The following are indexed by a newcorner index, 0 to newcorners_num - 1. */

  /** Analog of a #Mesh corner_verts() function. The indices are either in the
   * range of vertices of the source mesh, or they are newvert indices offset by
   * the number of vertices in the source mesh. */
  Span<int> newcorner_verts() const
  {
    return newcorner_verts_;
  }

  /** Analog of a #Mesh corner_edges() function. The indices are either in the
   * range of edges of the source mesh, or they are newedge indices offset by
   * the number of edges in the source mesh. */
  Span<int> newcorner_edges() const
  {
    return newcorner_edges_;
  }

  /** Is the edge corresponding to bevedge \a be beveled? */
  bool bevedge_is_beveled(const int be) const
  {
    return bevedge_weights_[be] > 0.0f ? true : false;
  }

  /** Return the next bevedge position around bevvert \a bv after \a edge_pos. */
  int next_edge_pos(const int bv, const int edge_pos) const
  {
    return edge_pos == bevvert_bevedges_[bv].size() - 1 ? 0 : edge_pos + 1;
  }

  /** Return the previous bevedge position around bevvert \a bv before \a edge_pos. */
  int prev_edge_pos(const int bv, const int edge_pos) const
  {
    return edge_pos == 0 ? bevvert_bevedges_[bv].size() - 1 : edge_pos + 1;
  }

  /** Return the next face after the bevedge in position \a edge_pos of bevvert \a bv. */
  int face_next(const int bv, const int edge_pos) const
  {
    return bevvert_faces_[bv][edge_pos];
  }

  /** Return the previous face before the bevedge in position \a edge_pos of bevvert \a bv. */
  int face_prev(const int bv, const int edge_pos) const
  {
    return bevvert_faces_[bv][prev_edge_pos(bv, edge_pos)];
  }

  /** Return the end of the bevedge \a be edge that the bevvert \a bv is
   * on: 0 if it is at the start of the edge, 1 if it is at the end.
   */
  int bevedge_vert_end(const int bv, const int be) const
  {
    return bevvert_mesh_verts_[bv] == mesh_info.mesh.edges()[bevedge_mesh_edges_[be]][0] ? 0 : 1;
  }

  /** Return the bevedge poisition of the last bevedge attached to newvert, or -1 if none. */
  int last_attached_bevedge_pos(const int bv, const int newvert) const;

 private:
  IndexMaskMemory memory_;
  IndexMask bevedges_mask_;
  IndexMask bevverts_mask_;
  Array<int> vert_bevverts_;
  Array<int> edge_bevedges_;
  Array<int> bevvert_mesh_verts_;
  Array<float> bevvert_weights_;
  GroupedSpan<int> bevvert_bevedges_;
  Array<int> bevvert_bevedges_offsets_;
  Array<int> bevvert_bevedges_indices_;
  GroupedSpan<int> bevvert_faces_;
  Array<int> bevvert_faces_indices_;
  Array<MeshPattern> bevvert_meshpatterns_;
  OffsetIndices<int> bevvert_newverts_;
  Array<int> bevvert_newverts_offsets_;
  OffsetIndices<int> bevvert_newedges_;
  Array<int> bevvert_newedges_offsets_;
  OffsetIndices<int> bevvert_newfaces_;
  Array<int> bevvert_newfaces_offsets_;
  Array<int> bevedge_mesh_edges_;
  Array<float> bevedge_weights_;
  Array<float4> bevedge_widths_;
  Array<int2> bevedge_attach_verts_;
  Array<int2> newedge_vertpairs_;
  OffsetIndices<int> bevedge_newedges_;
  OffsetIndices<int> bevedge_newfaces_;
  Array<float3> newvert_positions_;
  OffsetIndices<int> newface_faces_;
  Array<int> newcorner_verts_;
  Array<int> newcorner_edges_;
};

MeshInfo::MeshInfo(const Mesh &mesh) : mesh(mesh)
{
  /* Calculate the vert_to_edge topology map. */
  vert_edges_ = bke::mesh::build_vert_to_edge_map(
      mesh.edges(), mesh.verts_num, vert_to_edge_offsets_, vert_to_edge_indices_);
  /* Calculate the edge_to_face topology map. */
  edge_faces_ = bke::mesh::build_edge_to_face_map(mesh.faces(),
                                                  mesh.corner_edges(),
                                                  mesh.edges_num,
                                                  edge_to_face_map_indices_,
                                                  edge_to_face_map_offsets_);
}

/** Functions for debug printing. */

template<typename T> static void print_span(Span<T> span, const char *label)
{
  fmt::print("{}:", label);
  for (const int i : span.index_range()) {
    if (i % 10 == 0) {
      fmt::print("\n[{}] ", i);
    }
    fmt::print("{} ", span[i]);
  }
  fmt::println("");
}

static void print_float3(const float3 &v)
{
  fmt::print("({},{},{})", v[0], v[1], v[2]);
}

[[maybe_unused]] static void print_float3_span(Span<float3> span, const char *label)
{
  fmt::print("{}:", label);
  for (const int i : span.index_range()) {
    if (i % 10 == 0) {
      fmt::print("\n[{}] ", i);
    }
    print_float3(span[i]);
    fmt::print(" ");
  }
  fmt::println("");
}

static void print_indexmask(const IndexMask &index_mask, const char *label)
{
  fmt::print("{}:", label);
  index_mask.foreach_index([&](const int v, const int mask) {
    if (mask % 10 == 0) {
      fmt::print("\n[{}] ", mask);
    }
    fmt::print("{} ", v);
  });
  fmt::println("");
}

static void print_groupedspan(const GroupedSpan<int> &groupedspan, const char *label)
{
  fmt::println("{}:", label);
  for (const int i : groupedspan.index_range()) {
    fmt::print("[{}] ", i);
    for (int v : groupedspan[i]) {
      fmt::print("{} ", v);
    }
    fmt::println("");
  }
}

static void print_meshpattern(const MeshPattern &pat)
{
  static const char *kind_names[] = {"None", "Line", "TerminalPoly", "Adj", "TriFan", "Cutoff"};
  fmt::println("{} anchors={} segs={}", kind_names[int(pat.kind)], pat.num_anchors, pat.num_segs);
}

static void print_anchor_newvert_positions(const BevelState &bs, const char *label)
{
  fmt::println("{}", label);
  for (const int bv : IndexRange(bs.bevverts_num)) {
    const MeshPattern &pat = bs.bevvert_meshpatterns()[bv];
    fmt::print("bv {}:", bv);
    for (const int a : IndexRange(pat.num_anchors)) {
      const int apos = pat.anchor_vert(a);
      const int nv_index = bs.bevvert_newverts()[bv][apos];
      float3 co = bs.newvert_positions()[nv_index];
      fmt::print(" [{}]({}.{},{})", apos, co[0], co[1], co[2]);
    }
    fmt::println("");
  }
}

static void print_bevedge_attach_verts(const BevelState &bs, const char *label)
{
  fmt::println("{}", label);
  for (const int bv : IndexRange(bs.bevverts_num)) {
    fmt::println("bv {}: ", bv);
    Span<int> edges = bs.bevvert_bevedges()[bv];
    for (const int epos : edges.index_range()) {
      const int be = edges[epos];
      const int end = bs.bevedge_vert_end(bv, be);
      const int nv = bs.bevedge_attach_verts()[be][end];
      fmt::print(" [{}]: be={} nv={}", epos, be, nv);
    }
    fmt::println("");
  }
}

[[maybe_unused]] static void dump_bevel_state(const BevelState &bs, const char *label)
{
  bool vertex_only = bs.params.affect_type == BevelAffect::Vertices;
  fmt::println("\nBevelState {}", label);
  print_indexmask(bs.bevverts_mask(), "bevverts_mask");
  print_indexmask(bs.bevedges_mask(), "bevedges_mask");
  print_span(bs.vert_bevverts(), "vert_bevverts");
  print_span(bs.edge_bevedges(), "edge_bevedges");
  print_span(bs.bevvert_mesh_verts(), "bevvert_mesh_verts");
  if (vertex_only) {
    print_span(bs.bevvert_weights(), "bevvert_weights");
  }
  print_span(bs.bevedge_mesh_edges(), "bevedge_mesh_edges");
  if (!vertex_only) {
    print_span(bs.bevedge_weights(), "bevedge_weights");
  }
  print_groupedspan(bs.bevvert_bevedges(), "bevvert_bevedges");
  print_groupedspan(bs.bevvert_faces(), "bevvert_faces");
  if (bs.bevedge_attach_verts().size() > 0) {
    print_bevedge_attach_verts(bs, "bevedge_attach_verts");
  }
  if (bs.newvert_positions().size() > 0) {
    print_anchor_newvert_positions(bs, "anchor newvert positions");
  }
}

static void draw_bevvert(int bv, const BevelState &bs)
{
  constexpr uint life = draw::drw_debug_persistent_lifetime;
  constexpr float pntsize = 0.04f;
  const float4 orig_vert_col = {1, 0, 0, 1};
  const float4 first_bndv_col = {0, 0, 1, 1};
  const float4 bndv_col = {0, 1, 0, 1};
  const MeshPattern &pat = bs.bevvert_meshpatterns()[bv];
  float3 bevvert_co = bs.mesh_info.mesh.vert_positions()[bs.bevvert_mesh_verts()[bv]];
  draw::drw_debug_point(bevvert_co, pntsize, orig_vert_col, life);
  for (const int a : IndexRange(pat.num_anchors)) {
    const int apos = pat.anchor_vert(a);
    const int nv_index = bs.bevvert_newverts()[bv][apos];
    float3 co = bs.newvert_positions()[nv_index];
    draw::drw_debug_point(co, pntsize, a == 0 ? first_bndv_col : bndv_col, life);
  }
}

[[maybe_unused]] static void draw_all_bevverts(const BevelState &bs)
{
  for (const int bv : IndexRange(bs.bevverts_num)) {
    draw_bevvert(bv, bs);
  }
}

namespace geom {
/** Functions and data for geometric calculations. */

constexpr double bevel_epsilon_d = 1e-6;
constexpr float bevel_epsilon = 1e-6f;
constexpr float bevel_epsilon_sq = 1e-12f;
constexpr float bevel_epsilon_big = 1e-4f;
constexpr float bevel_epsilon_big_sq = 1e-8f;
constexpr float bevel_epsilon_ang = DEG2RADF(2.0f);
constexpr float bevel_small_ang = DEG2RADF(10.0f);
/** Difference in dot products that corresponds to 10 degree difference between vectors. */
const float bevel_small_ang_dot = (1.0f - math::cos(bevel_small_ang));
/** Difference in dot products that corresponds to 2.0 degree difference between vectors. */
const float bevel_epsilon_ang_dot = (1.0f - cosf(bevel_epsilon_ang));
constexpr float bevel_max_adjust_pct = 10.0f;
constexpr float bevel_max_auto_adjust_pct = 300.0f;
constexpr double bevel_match_spec_weight = 0.2;

/** Return the direction of the edge leaving bevvert \a bv and position \a epos. */
static float3 edge_dir(const int bv, const int epos, const BevelState &bs)
{
  const Mesh &mesh = bs.mesh_info.mesh;
  int be = bs.bevvert_bevedges()[bv][epos];
  int2 e = mesh.edges()[bs.bevedge_mesh_edges()[be]];
  float3 dir = mesh.vert_positions()[e[1]] - mesh.vert_positions()[e[0]];
  return bs.bevedge_vert_end(bv, be) == 0 ? dir : -dir;
}

/** Return the length of the edge for bevvert \a bv and position \a epos. */
static float edge_length(const int bv, const int epos, const BevelState &bs)
{
  return math::length(edge_dir(bv, epos, bs));
}

/** Return the angle going counterclockwise when viewed from the vertex normal side,
 * from the edge in position \a e1_pos to the edge in position \a e2_pos of bevedges
 * connected to Bevvert \a bv.
 */
static float angle_between_edges(const int bv,
                                 const int e1_pos,
                                 const int e2_pos,
                                 const BevelState &bs)
{
  const Mesh &mesh = bs.mesh_info.mesh;
  const int v = bs.bevvert_mesh_verts()[bv];
  float3 dir1 = math::normalize(edge_dir(bv, e1_pos, bs));
  float3 dir2 = math::normalize(edge_dir(bv, e2_pos, bs));
  float ang = angle_normalized_v3v3(dir1, dir2);
  /* Angles are in [0,pi]. Need to compare cross product with normal to see if it is reflex. */
  const int face_next = bs.face_next(bv, e1_pos);
  const int face_prev = bs.face_prev(bv, e2_pos);
  float3 normal = face_next != -1 ?
                      mesh.face_normals()[face_next] :
                      (face_prev != -1 ? mesh.face_normals()[face_prev] : mesh.vert_normals()[v]);
  if (math::dot(math::cross(dir1, dir2), normal) < 0.0f) {
    /* Angle is reflex. */
    ang = 2 * math::numbers::pi - ang;
  }
  return ang;
}

/**
 * \return True if d1 and d2 are parallel or nearly parallel.
 */
static bool nearly_parallel_normalized(const float3 d1, const float3 d2)
{
  return compare_ff(math::abs(math::dot(d1, d2)), 1.0f, bevel_epsilon_ang_dot);
}

/** Return whether the angle is less than, equal to, or larger than 180 degrees when viewed from
 * the positive normal side of faces (if available) or vertex. */
[[maybe_unused]] static AngleKind edges_angle_kind(const int bv,
                                                   const int e1_pos,
                                                   const int e2_pos,
                                                   const BevelState &bs)
{
  float ang = angle_between_edges(bv, e1_pos, e2_pos, bs);
  if (math::abs(ang - math::numbers::pi) < bevel_epsilon_ang) {
    return AngleKind::Straight;
  }
  else if (ang < math::numbers::pi) {
    return AngleKind::Smaller;
  }
  return AngleKind::Larger;
}

/** Return the angle between the two faces adjacent to the bevedge that
 * is at position \a edge_pos of bevvert \a bv.
 * If there are not two, return 0. */
static float edge_face_angle(const int bv, const int edge_pos, const BevelState &bs)
{
  const int face_next = bs.face_next(bv, edge_pos);
  const int face_prev = bs.face_prev(bv, edge_pos);
  if (face_prev == -1 || face_next == -1) {
    return 0.0f;
  }
  const float3 &norm_prev = bs.mesh_info.mesh.face_normals()[face_prev];
  const float3 &norm_next = bs.mesh_info.mesh.face_normals()[face_next];
  return math::numbers::pi - angle_normalized_v3v3(norm_prev, norm_next);
}

/** Find the center axis of the given bevel vert, for use in vertex beveling.
 * Note; Don't use the vertex normal -- it can give unwanted results.
 */
static float3 bevvert_axis(const int bv, const BevelState &bs)
{
  float3 vert_axis(0.0f, 0.0f, 0.0f);
  Span<int> edges = bs.bevvert_bevedges()[bv];
  for (const int be : edges) {
    float3 edir = edge_dir(bv, be, bs);
    if (bs.bevedge_vert_end(bv, be) == 0) {
      edir = -edir;
    }
    vert_axis += math::normalize(edir);
  }
  return vert_axis;
}

/** Assuming \a co is on the edge defined by \a be and \a epos,
 * find out if is on or between those two vertices. If it is, return -1.
 * Otherwise return the index of the vertex that is the end that is nearer \a co.
 */
static bool is_outside_edge(
    const int bv, const int epos, const float3 co, const BevelState &bs, float3 *r_closer)
{
  const Mesh &mesh = bs.mesh_info.mesh;
  int be = bs.bevvert_bevedges()[bv][epos];
  int2 edge = mesh.edges()[bs.bevedge_mesh_edges()[be]];
  float3 l1 = mesh.vert_positions()[edge[0]];
  float3 l2 = mesh.vert_positions()[edge[1]];
  float lenu;
  float3 u = math::normalize_and_get_length(l2 - l1, lenu);
  float lambda = math::dot(u, co - l1);
  if (lambda <= -bevel_epsilon_big * lenu) {
    *r_closer = l1;
    return true;
  }
  if (lambda > (1.0 + bevel_epsilon_big) * lenu) {
    *r_closer = l2;
    return true;
  }
  return false;
}

/** \a co should be approximately on the plane between the edges in positions \a e1_pos
 * and \a e2_pos at bevvert \a bv. They should share face \a face (which cannot be -1).
 * Is it between those edges, sweeping CCW?
 */
static bool point_between_edges(const int bv,
                                const int e1_pos,
                                const int e2_pos,
                                const float3 co,
                                const int face,
                                const BevelState &bs)
{
  BLI_assert(face != -1);
  const Mesh &mesh = bs.mesh_info.mesh;
  float3 dir1 = -math::normalize(edge_dir(bv, e1_pos, bs));
  float3 dir2 = -math::normalize(edge_dir(bv, e2_pos, bs));
  const int v = bs.bevvert_mesh_verts()[bv];
  float3 dirco = math::normalize(mesh.vert_positions()[v] - co);

  float ang11 = angle_normalized_v3v3(dir1, dir2);
  float ang1co = angle_normalized_v3v3(dir1, dirco);
  /* Angles are in [0,pi]. Need to compare cross product with normal to see if they are reflex.
   */
  float3 norm = math::cross(dir1, dir2);
  float3 face_norm = mesh.face_normals()[face];
  if (math::dot(norm, face_norm) < 0.0f) {
    ang11 = float(math::numbers::pi * 2.0) - ang11;
  }
  norm = math::cross(dir1, dirco);
  if (math::dot(norm, face_norm) < 0.0f) {
    ang1co = float(math::numbers::pi * 2.0) - ang1co;
  }
  return (ang11 - ang1co > -bevel_epsilon_ang);
}

/** Return coordinates of a point a distance \a dist along the edge specified by \a bv and \a
 * epos. However, clamp it to be no more than the length of the edge.
 */
static float3 slide_dist(const int bv, const int epos, const float dist, const BevelState &bs)
{
  float len;
  float3 normalized_dir = math::normalize_and_get_length(edge_dir(bv, epos, bs), len);
  float d = dist;
  if (d >= len) {
    /* Don't go quite all the way to the other end. */
    d = len - float(50.0 - bevel_epsilon_d);
  }
  int origin_v = bs.bevvert_mesh_verts()[bv];
  return bs.mesh_info.mesh.vert_positions()[origin_v] + d * normalized_dir;
}

/** Return the coordinates of the offset by \a dist of the origin end of the edge specified
 * by \a bv and \a epos,
 * on the left side if \a left is true, else the right side.
 * Do the offset in the plane that is on that side if it exists, else choose an arbitrary
 * plane normal.
 */
static float3 offset_bevedge(
    const int bv, const int epos, const float dist, bool left, const BevelState &bs)
{
  const Mesh &mesh = bs.mesh_info.mesh;
  float3 dir = math::normalize(edge_dir(bv, epos, bs));
  float3 normal;
  int face = left ? bs.face_prev(bv, epos) : bs.face_next(bv, epos);
  if (face != -1) {
    normal = mesh.face_normals()[face];
  }
  else {
    /* An arbitrary direction not the same as dir. */
    if (math::abs(dir[0]) < math::abs(dir[1])) {
      normal = {1.0f, 0.0f, 0.0f};
    }
    else {
      normal = {0.0f, 1.0f, 0.0f};
    }
  }
  float3 fdir = math::normalize(left ? math::cross(dir, normal) : math::cross(normal, dir));
  int origin_v = bs.bevvert_mesh_verts()[bv];
  return mesh.vert_positions()[origin_v] + dist * fdir;
}

/** When the offset type is Percent or Absolute, fill in the coordinates
 * of the lines whose intersection defines the boundary point between e1 and e2 with common
 * vert v, as defined in the parameters of offset_meet.
 */
static void offset_meet_lines_percent_or_absolute(const int /*bv*/,
                                                  const int /*e1_pos*/,
                                                  const int /*e2_pos*/,
                                                  MutableSpan<float3> /*pts*/,
                                                  const BevelState & /*bs*/)
{
  /* Get points the specified distance along each leg.
   * The legs we need are:
   *   e0 : the next edge around e1->fnext (==f1) after e1.
   *   e3 : the prev edge around e2->fprev (==f2) before e2.
   *   e4 : the previous edge around f1 before e1 (may be e2).
   *   e5 : the next edge around f2 after e2 (may be e1).
   */
  /* TODO */
  BLI_assert(false);
}

/** Return the meeting point between the offset edges for edges at positions \a e1_pos and \a
 * e2_pos. Do this in \a face, if that is not -1. Offset edge is on right of both edges, where e1
 * enters bv and e2 leave it. When offsets are equal, the new point is on the edge bisector, with
 * length offset/sin(angle/2), but if the offsets are not equal (we allow for because the bevel
 * modifier has edge weights that may lead to different offsets) then the meeting point can be
 * found by intersecting offset lines.
 *
 * \a edges_between: If this is true, there are edges between e1 and e2 in CCW order so they
 * do not share a common face. We want the meeting point to be on an existing face so it
 * should be dropped onto one of the intermediate faces, if possible.
 * \a e_in_plane: If we need to drop from the calculated offset lines to one of the faces,
 * we do not want to drop onto the 'in plane' face, so if this is not -1 skip this edge's faces.
 */
static float3 offset_meet(const int bv,
                          const int e1_pos,
                          const int e2_pos,
                          const float e1_spec_r,
                          const float e2_spec_l,
                          const int face,
                          const bool edges_between,
                          const int e_in_plane_pos,
                          const BevelState &bs)
{
  const Mesh &mesh = bs.mesh_info.mesh;
  const int v = bs.bevvert_mesh_verts()[bv];
  const float3 bevvert_pos = mesh.vert_positions()[v];
  /* Get direction along e1 into bv and along e2 from bv. */
  const float3 dir1 = -edge_dir(bv, e1_pos, bs);
  const float3 dir2 = edge_dir(bv, e2_pos, bs);
  Span<int> bevedges = bs.bevvert_bevedges()[bv];
  const int e1_next_pos = index_range_next(bevedges.index_range(), e1_pos);
  const int e2_prev_pos = index_range_prev(bevedges.index_range(), e2_pos);
  const float3 dir1n = -edge_dir(bv, e1_next_pos, bs);
  const float3 dir2p = edge_dir(bv, e2_prev_pos, bs);
  float3 vert_norm = mesh.vert_normals()[v];
  float ang = angle_v3v3(dir1, dir2);
  if (ang < bevel_epsilon_ang) {
    /* Special case: e1 and e2 are parallel, or nearly so;
     * put offset point perp to both, from bevvert.v on a suitable plane.
     * This code used to just use offset and dir1, but that makes for visible errors
     * on a circle with > 200 sides, which trips this "nearly perp" code (see #61214).
     * so use the average of the two, and the offset formula for angle bisector.
     * If offsets are different, we're out of luck:
     * Use the max of the two (so get consistent looking results if the same situation
     * arises elsewhere in the object but with opposite roles for e1 and e2.
     */
    float3 norm_v;
    if (face != -1) {
      norm_v = mesh.face_normals()[face];
    }
    else {
      /* Get average of face normals between e and e2. */
      norm_v = {0.0f, 0.0f, 0.0f};
      int fcount = 0;
      for (int i = e1_pos; i != e2_pos; i = (i + 1) % bevedges.size()) {
        const int face = bs.face_next(bv, i);
        if (face != -1) {
          norm_v += mesh.face_normals()[face];
          fcount++;
        }
      }
      if (fcount == 0) {
        norm_v = vert_norm;
      }
      else {
        norm_v = norm_v / float(fcount);
      }
    }
    float3 norm_perp1 = math::normalize(math::cross(dir1 + dir2, norm_v));
    float d = math::max(e1_spec_r, e2_spec_l) / math::cos(ang / 2.0);
    return bevvert_pos + d * norm_perp1;
  }
  if (math::abs(ang - math::numbers::pi) < bevel_epsilon_ang) {
    /* Special case: e1 and e2 are anti-parallel, so bevel is into a zero-area face.
     * Just make the offset point on the common line, at offset distance from v.
     */
    return slide_dist(bv, e2_pos, math::max(e1_spec_r, e2_spec_l), bs);
  }
  /* Get normal to plane where meet point should be, using cross product instead of
   * face's normal, in case f is non-planar.
   * Except: sometimes locally there can be a small angle between dir1 and dir2 that leads
   * to a normal that is actually almost perpendicular to the face normal;
   * in this case it looks wrong to use the local (cross-product) normal, so use the face normal
   * if the angle between dir1 and dir2 is smallish.
   * If e1-v-e2 is a reflex angle (viewed from vertex normal side), need to flip.
   * Use face's normal to figure out which side to look at angle from, as even if face is
   * non-planar, this will be more accurate than the vertex normal.
   */
  float3 norm_v1;
  float3 norm_v2;
  if (face != -1 && ang < bevel_small_ang) {
    norm_v1 = mesh.face_normals()[face];
    norm_v2 = norm_v1;
  }
  else if (!edges_between) {
    norm_v1 = math::normalize(math::cross(dir2, dir1));
    if (math::dot(norm_v1, face != -1 ? mesh.face_normals()[face] : vert_norm) < 0.0f) {
      norm_v1 = -norm_v1;
    }
    norm_v2 = norm_v1;
  }
  else {
    /* Separate faces; get face norms at corners for each separately. */
    norm_v1 = math::normalize(math::cross(dir1n, dir1));
    int f = bs.face_next(bv, e1_pos);
    if (math::dot(norm_v1, f != -1 ? mesh.face_normals()[f] : vert_norm) < 0.0f) {
      norm_v1 = -norm_v1;
    }
    norm_v2 = math::normalize(math::cross(dir2, dir2p));
    f = bs.face_prev(bv, e2_pos);
    if (math::dot(norm_v2, f != -1 ? mesh.face_normals()[f] : vert_norm) < 0.0f) {
      norm_v2 = -norm_v2;
    }
  }
  /* Get vectors perp to each edge, perp to norm_v, and pointing into face. */
  float3 norm_perp1 = math::normalize(math::cross(dir1, norm_v1));
  float3 norm_perp2 = math::normalize(math::cross(dir2, norm_v2));
  /* Get points on two lines to intersect in order to find the meet point. */
  Array<float3, 4> off_pts(4);
  if (ELEM(bs.params.offset_type, BevelOffsetType::Percent, BevelOffsetType::Absolute)) {
    offset_meet_lines_percent_or_absolute(bv, e1_pos, e2_pos, off_pts.as_mutable_span(), bs);
  }
  else {
    off_pts[0] = bevvert_pos + e1_spec_r * norm_perp1;
    off_pts[1] = off_pts[0] + dir1;
    off_pts[2] = bevvert_pos + e2_spec_l * norm_perp2;
    off_pts[3] = off_pts[2] + dir2;
  }
  /* Intersect the offset lines. */
  float3 meetco;
  float3 isect2;
  int isect_kind = isect_line_line_v3(
      off_pts[0], off_pts[1], off_pts[2], off_pts[3], meetco, isect2);
  if (isect_kind == 0) {
    /* Lines are collinear: we already tested for this, but this used a different epsilon. */
    return off_pts[0]; /* Just do something. */
  }
  /* The lines intersect, but is it at a reasonable place?
   * One problem to check: if one of the offsets is 0, then we don't want an intersection
   * that is outside that edge itself. This can happen if angle between them is > 180 degrees,
   * or if the offset amount is > the edge length.
   */
  float3 vcloser;
  if (e1_spec_r == 0.0f && is_outside_edge(bv, e1_pos, meetco, bs, &vcloser)) {
    meetco = vcloser;
  }
  if (e2_spec_l == 0.0f && is_outside_edge(bv, e2_pos, meetco, bs, &vcloser)) {
    meetco = vcloser;
  }
  if (edges_between && e1_spec_r > 0.0f && e2_spec_l > 0.0f) {
    /* Try to drop meetco to a face between e1 and e2. */
    if (isect_kind == 2) {
      /* Lines didn't meet in 3d: get average of meetco and isect2. */
      mid_v3_v3v3(meetco, meetco, isect2);
    }

    for (int i = e1_pos; i != e2_pos; i = (i + 1) % bevedges.size()) {
      const int fnext = bs.face_next(bv, i);
      if (fnext == -1) {
        continue;
      }
      const float3 fnext_norm = mesh.face_normals()[fnext];
      float4 plane;
      plane_from_point_normal_v3(plane, bevvert_pos, fnext_norm);
      float3 dropco;
      closest_to_plane3_normalized_v3(dropco, plane, meetco);
      /* Don't drop to the faces next to the in plane edge. */
      if (e_in_plane_pos != -1) {
        const int face = bs.face_next(bv, e_in_plane_pos);
        if (face != -1) {
          ang = angle_v3v3(fnext_norm, mesh.face_normals()[face]);
          if (math::abs(ang) < bevel_small_ang ||
              math::abs(ang - math::numbers::pi) < bevel_small_ang)
          {
            continue;
          }
        }
      }
      if (point_between_edges(bv, i, (i + 1) % bevedges.size(), dropco, fnext, bs)) {
        meetco = dropco;
        break;
      }
    }
  }

  return meetco;
}

/* This was changed from 0.25f to fix bug #86768.
 * Original bug #44961 remains fixed with this value.
 * Update: changed again from 0.0001f to fix bug #95335.
 * Original two bugs remained fixed.
 */
const float good_angle = 0.1f;

/** Calculate the meeting point between the right offset edge  of the edges a position \a
 * e1_pos and the left offset edge of the edge at position \a e2_pos, where \a e1 precedes \a
 * e2 in CCW order. Assume that at most of the specified widths is non-zero. It is possible
 * that no such meeting point exists (if the angle between is reflex), in which case false is
 * returned.  It is also possible that the angle is so close to a straight angle that the
 * position of the meet point would be a big spike, so return false in that case too.
 * Otherwise, return the meeting point in *r_co and the angle between in *r_angle and return
 * true.
 */
static bool try_offset_meet_edge(const int bv,
                                 const int e1_pos,
                                 const int e2_pos,
                                 const float e1_spec_r,
                                 const float e2_spec_l,
                                 float3 *r_co,
                                 float *r_angle,
                                 const BevelState &bs)
{
  BLI_assert(e1_spec_r == 0.0f || e2_spec_l == 0.0f);
  float ang = angle_between_edges(bv, e1_pos, e2_pos, bs);
  if (math::abs(ang) < good_angle || math::numbers::pi - ang < good_angle) {
    return false;
  }
  if (r_angle) {
    *r_angle = ang;
  }
  if (r_co == nullptr) {
    return true;
  }
  const Mesh &mesh = bs.mesh_info.mesh;
  float3 co = mesh.vert_positions()[bs.bevvert_mesh_verts()[bv]];
  float sinang = math::sin(ang);
  BLI_assert(sinang != 0.0f);
  if (e1_spec_r == 0.0f) {
    float3 dir1 = math::normalize(edge_dir(bv, e1_pos, bs));
    *r_co = co + (e2_spec_l / sinang) * dir1;
  }
  else {
    float3 dir2 = math::normalize(edge_dir(bv, e2_pos, bs));
    *r_co = co + (e1_spec_r / sinang) * dir2;
  }
  return true;
}

/** Return true if it will look good to put the meeting point where try_offset_on_edge_between
 * would put it. This means that neither side sees a reflex angle or too close to a straight
 * angle.
 */
static bool good_slide(const int bv,
                       const int e1_pos,
                       const int e2_pos,
                       const int eslide_pos,
                       const float e1_spec_r,
                       const float e2_spec_l,
                       const BevelState &bs)
{
  return try_offset_meet_edge(bv, e1_pos, eslide_pos, e1_spec_r, 0.0f, nullptr, nullptr, bs) &&
         try_offset_meet_edge(bv, eslide_pos, e2_pos, 0.0f, e2_spec_l, nullptr, nullptr, bs);
}

/** Calculate the best place for a meeting point for the offsets from edges at \a e1_pos and \a
 * e2_pos on the in-between edge at \a eon_pos. Viewed from the vertex normal side, the CCW
 * order of these edges is e1, eon, e2. Return true if we placed in *r_co as compromise between
 * where two edges met. If we did, put the ratio of sines of angles in *r_sin_ratio too.
 * However, if the offset_type is Percent or Absolute, we just slide along eon by the specified
 * amount.
 */
static bool try_offset_on_edge_between(const int bv,
                                       const int e1_pos,
                                       const int e2_pos,
                                       const int eon_pos,
                                       const float e1_spec_r,
                                       const float e2_spec_l,
                                       float3 *r_co,
                                       float *r_sin_ratio,
                                       const BevelState &bs)
{
  BLI_assert(bs.bevedge_is_beveled(bs.bevvert_bevedges()[bv][e1_pos]) &&
             bs.bevedge_is_beveled(bs.bevvert_bevedges()[bv][e2_pos]) &&
             !bs.bevedge_is_beveled(bs.bevvert_bevedges()[bv][eon_pos]));
  BLI_assert(r_co != nullptr);
  float ang1, ang2;
  float3 meet1, meet2;
  bool ok1 = try_offset_meet_edge(bv, e1_pos, eon_pos, e1_spec_r, 0.0f, &meet1, &ang1, bs);
  bool ok2 = try_offset_meet_edge(bv, eon_pos, e2_pos, 0.0f, e2_spec_l, &meet2, &ang2, bs);
  if (r_sin_ratio != nullptr) {
    *r_sin_ratio = ang1 == 0.0f ? 1.0f : math::sin(ang2) / math::sin(ang1);
  }
  if (ELEM(bs.params.offset_type, BevelOffsetType::Percent, BevelOffsetType::Absolute)) {
    float3 v_co = bs.mesh_info.mesh.vert_positions()[bs.bevvert_mesh_verts()[bv]];
    float eon_len;
    float3 eon_unit_dir = math::normalize_and_get_length(edge_dir(bv, eon_pos, bs), eon_len);
    if (bs.params.offset_type == BevelOffsetType::Percent) {
      /* TODO: use weights, if required. */
      *r_co = v_co + (bs.params.offset / 100.0f) * eon_len * eon_unit_dir;
    }
    else {
      *r_co = v_co + bs.params.offset * eon_unit_dir;
    }
    return true;
  }
  if (ok1 && ok2) {
    /* The two sides will likely lead to different meet points.
     * Compromise on the midpint between them.
     */
    mid_v3_v3v3(*r_co, meet1, meet2);
    return true;
  }
  else if (ok1 && !ok2) {
    *r_co = meet1;
  }
  else if (!ok1 && ok2) {
    *r_co = meet2;
  }
  else {
    /* Neither offset line met eon.
     * This should only happen if all three lines are on top of each other.
     */
    *r_co = slide_dist(bv, eon_pos, e1_spec_r, bs);
  }
  return false;
}

/** Return true if the specified edge is between two faces with 180 degres between their normals.
 */
static bool bevedge_on_plane(const int bv, const int epos, const BevelState &bs)
{
  const int face_prev = bs.face_prev(bv, epos);
  const int face_next = bs.face_next(bv, epos);
  const Mesh &mesh = bs.mesh_info.mesh;
  if (face_prev != -1 && face_next != -1) {
    float dot = math::dot(mesh.face_normals()[face_prev], mesh.face_normals()[face_next]);
    return math::abs(dot + 1.0f) <= bevel_epsilon_big ||
           math::abs(dot - 1.0f) <= bevel_epsilon_big;
  }
  return false;
}

/** Return the closest point on the line (a, b) to the given bevedge specified by \a bv and \a
 * epos. */
static float3 project_to_edge(
    const int bv, const int epos, const float3 &a, const float3 &b, const BevelState &bs)
{
  float3 co1, co2;
  const int be = bs.bevvert_bevedges()[bv][epos];
  const Mesh &mesh = bs.mesh_info.mesh;
  const int2 mesh_e = mesh.edges()[bs.bevedge_mesh_edges()[be]];
  const float3 e1 = mesh.vert_positions()[mesh_e[0]];
  const float3 e2 = mesh.vert_positions()[mesh_e[1]];
  if (!isect_line_line_v3(e1, e2, a, b, co1, co2)) {
    return e1;
  }
  return co1;
}

}  // end namespace geom

namespace profile {
/** Structures, functions and constants related to superellipse profiles. */

/**
 * Profile specification:
 * The profile is a path defined with start, middle, and end control points projected onto a
 * plane (plane_no is normal, plane_co is a point on it) via lines in a given direction (proj_dir).
 *
 * Many interesting profiles are in family of superellipses:
 *     (abs(x/a))^r + abs(y/b))^r = 1
 * r==2 => ellipse; r==1 => line; r < 1 => concave; r > 1 => bulging out.
 * Special cases: let r==0 mean straight-inward, and r==4 mean straight outward.
 *
 * After the parameters are all set, the actual profile points are calculated and pointed to
 * by prof_co. We also may need profile points for a higher resolution number of segments
 * for the subdivision while making the ADJ vertex mesh pattern, and that goes in prof_co_2.
 */
struct Profile {
  /** Superellipse r parameter. */
  float super_r;
  /** Height for profile cutoff face sides. */
  float height;
  /** Start control point for profile. */
  float3 start;
  /** Mid control point for profile. */
  float3 middle;
  /** End control point for profile. */
  float3 end;
  /** Normal of plane to project to. */
  float3 plane_no;
  /** Coordinate on plane to project to. */
  float3 plane_co;
  /** Direction of projection line. */
  float3 proj_dir;
  /** segments+1 profile coordinates. */
  Array<float3, 10> prof_co;
  /** Like prof_co, but for seg power of 2 >= seg. */
  Array<float3, 10> prof_co_2;
  /** Mark a special case so the these parameters aren't reset with others. */
  bool special_params;
};

/** Holds the profiles for each anchor vertex in a vertex mesh pattern.
 * Given an inline capacity to make the need for allocates rare.
 */
typedef Array<Profile, 20> AnchorProfiles;

/* Values for super_r to give particular special shapes. */
constexpr float pro_square_r = 1e4f;
constexpr float pro_circle_r = 2.0f;
constexpr float pro_line_r = 1.0f;
constexpr float pro_square_in_r = 0.0f;

/**
 * Get the coordinate on the superellipse (x^r + y^r = 1), at parameter value x
 * (or, if !rbig, mirrored (y=x)-line).
 * rbig should be true if r > 1.0 and false if <= 1.0.
 * Assume r > 0.0.
 */
static double superellipse_co(double x, float r, bool rbig)
{
  BLI_assert(r > 0.0f);

  /* If r<1, mirror the superellipse function by (y=x)-line to get a numerically stable range
   * Possible because of symmetry, later mirror back. */
  double dr = r;
  if (rbig) {
    return math::pow((1.0 - math::pow(x, dr)), (1.0 / dr));
  }
  return 1.0 - math::pow((1.0 - math::pow(1.0 - x, dr)), (1.0 / dr));
}

/** Find xnew > x0 so that distance((x0,y0), (xnew, ynew)) = dtarget.
 * False position Illinois method used because the function is somewhat linear
 * -> linear interpolation converges fast.
 * Assumes that the gradient is always between 1 and -1 for x in [x0, x0+dtarget].
 */
static double find_superellipse_chord_endpoint(double x0, double dtarget, float r, bool rbig)
{
  double y0 = superellipse_co(x0, r, rbig);
  const double tol = 1e-13; /* accumulates for many segments so use low value. */
  const int maxiter = 10;

  /* For gradient between -1 and 1, xnew can only be in [x0 + sqrt(2)/2*dtarget, x0 + dtarget].
   */
  double xmin = x0 + math::numbers::sqrt2 / 2.0 * dtarget;
  xmin = std::min(xmin, 1.0);
  double xmax = x0 + dtarget;
  xmax = std::min(xmax, 1.0);
  double ymin = superellipse_co(xmin, r, rbig);
  double ymax = superellipse_co(xmax, r, rbig);

  /* NOTE: using distance**2 (no sqrt needed) does not converge that well. */
  double dmaxerr = math::sqrt(math::pow((xmax - x0), 2.0) + math::pow((ymax - y0), 2.0)) - dtarget;
  double dminerr = math::sqrt(math::pow((xmin - x0), 2.0) + math::pow((ymin - y0), 2.0)) - dtarget;

  double xnew = xmax - dmaxerr * (xmax - xmin) / (dmaxerr - dminerr);
  bool lastupdated_upper = true;

  for (int iter = 0; iter < maxiter; iter++) {
    double ynew = superellipse_co(xnew, r, rbig);
    double dnewerr = math::sqrt(math::pow((xnew - x0), 2.0) + math::pow((ynew - y0), 2.0)) -
                     dtarget;
    if (fabs(dnewerr) < tol) {
      break;
    }
    if (dnewerr < 0) {
      xmin = xnew;
      ymin = ynew;
      dminerr = dnewerr;
      if (!lastupdated_upper) {
        xnew = (dmaxerr / 2 * xmin - dminerr * xmax) / (dmaxerr / 2 - dminerr);
      }
      else {
        xnew = xmax - dmaxerr * (xmax - xmin) / (dmaxerr - dminerr);
      }
      lastupdated_upper = false;
    }
    else {
      xmax = xnew;
      ymax = ynew;
      dmaxerr = dnewerr;
      if (lastupdated_upper) {
        xnew = (dmaxerr * xmin - dminerr / 2 * xmax) / (dmaxerr - dminerr / 2);
      }
      else {
        xnew = xmax - dmaxerr * (xmax - xmin) / (dmaxerr - dminerr);
      }
      lastupdated_upper = true;
    }
  }
  return xnew;
}

/**
 * This search procedure to find equidistant points (x,y) in the first
 * superellipse quadrant works for every superellipse exponent but is more
 * expensive than known solutions for special cases.
 * Call the point on superellipse that intersects x=y line mx.
 * For r>=1 use only the range x in [0,mx] and mirror the rest along x=y line,
 * for r<1 use only x in [mx,1]. Points are initially spaced and iteratively
 * repositioned to have the same distance.
 */
static void find_even_superellipse_chords_general(int seg,
                                                  float r,
                                                  MutableSpan<double> xvals,
                                                  MutableSpan<double> yvals)
{
  const int smoothitermax = 10;
  const double error_tol = 1e-7;
  int imax = (seg + 1) / 2 - 1; /* Ceiling division - 1. */

  bool seg_odd = seg % 2;

  bool rbig;
  double mx;
  if (r > 1.0f) {
    rbig = true;
    mx = math::pow(0.5, 1.0 / r);
  }
  else {
    rbig = false;
    mx = 1 - math::pow(0.5, 1.0 / r);
  }

  /* Initial positions, linear spacing along x axis. */
  for (int i = 0; i <= imax; i++) {
    xvals[i] = i * mx / seg * 2;
    yvals[i] = superellipse_co(xvals[i], r, rbig);
  }
  yvals[0] = 1;

  /* Smooth distance loop. */
  for (int iter = 0; iter < smoothitermax; iter++) {
    double sum = 0.0;
    double dmin = 2.0;
    double dmax = 0.0;
    /* Update distances between neighbor points. Store the highest and
     * lowest to see if the maximum error to average distance (which isn't
     * known yet) is below required precision. */
    for (int i = 0; i < imax; i++) {
      double d = math::sqrt(math::pow((xvals[i + 1] - xvals[i]), 2.0) +
                            math::pow((yvals[i + 1] - yvals[i]), 2.0));
      sum += d;
      dmax = std::max(d, dmax);
      dmin = std::min(d, dmin);
    }
    /* For last distance, weight with 1/2 if seg_odd. */
    double davg;
    if (seg_odd) {
      sum += math::numbers::sqrt2 / 2 * (yvals[imax] - xvals[imax]);
      davg = sum / (imax + 0.5);
    }
    else {
      sum += math::sqrt(math::pow((xvals[imax] - mx), 2.0) + math::pow((yvals[imax] - mx), 2.0));
      davg = sum / (imax + 1.0);
    }
    /* Max error in tolerance? -> Quit. */
    bool precision_reached = true;
    if (dmax - davg > error_tol) {
      precision_reached = false;
    }
    if (dmin - davg < error_tol) {
      precision_reached = false;
    }
    if (precision_reached) {
      break;
    }

    /* Update new coordinates. */
    for (int i = 1; i <= imax; i++) {
      xvals[i] = find_superellipse_chord_endpoint(xvals[i - 1], davg, r, rbig);
      yvals[i] = superellipse_co(xvals[i], r, rbig);
    }
  }

  /* Fill remaining. */
  if (!seg_odd) {
    xvals[imax + 1] = mx;
    yvals[imax + 1] = mx;
  }
  for (int i = imax + 1; i <= seg; i++) {
    yvals[i] = xvals[seg - i];
    xvals[i] = yvals[seg - i];
  }

  if (!rbig) {
    for (int i = 0; i <= seg; i++) {
      double temp = xvals[i];
      xvals[i] = 1.0 - yvals[i];
      yvals[i] = 1.0 - temp;
    }
  }
}

/**
 * Find equidistant points `(x0,y0), (x1,y1)... (xn,yn)` on the superellipse
 * function in the first quadrant. For special profiles (linear, arc,
 * rectangle) the point can be calculated easily, for any other profile a more
 * expensive search procedure must be used because there is no known closed
 * form for equidistant parametrization.
 * `xvals` and `yvals` should be size `n+1`.
 */
static void find_even_superellipse_chords(int n,
                                          float r,
                                          MutableSpan<double> xvals,
                                          MutableSpan<double> yvals)
{
  bool seg_odd = n % 2;
  int n2 = n / 2;

  /* Special cases. */
  if (r == pro_line_r) {
    /* Linear spacing. */
    for (int i = 0; i <= n; i++) {
      xvals[i] = double(i) / n;
      yvals[i] = 1.0 - double(i) / n;
    }
    return;
  }
  if (r == pro_circle_r) {
    double temp = M_PI_2 / n;
    /* Angle spacing. */
    for (int i = 0; i <= n; i++) {
      xvals[i] = math::sin(i * temp);
      yvals[i] = math::cos(i * temp);
    }
    return;
  }
  if (r == pro_square_in_r) {
    /* n is even, distribute first and second half linear. */
    if (!seg_odd) {
      for (int i = 0; i <= n2; i++) {
        xvals[i] = 0.0;
        yvals[i] = 1.0 - double(i) / n2;
        xvals[n - i] = yvals[i];
        yvals[n - i] = xvals[i];
      }
    }
    /* n is odd, so get one corner-cut chord. */
    else {
      double temp = 1.0 / (n2 + math::numbers::sqrt2 / 2.0);
      for (int i = 0; i <= n2; i++) {
        xvals[i] = 0.0;
        yvals[i] = 1.0 - double(i) * temp;
        xvals[n - i] = yvals[i];
        yvals[n - i] = xvals[i];
      }
    }
    return;
  }
  if (r == pro_square_r) {
    /* n is even, distribute first and second half linear. */
    if (!seg_odd) {
      for (int i = 0; i <= n2; i++) {
        xvals[i] = double(i) / n2;
        yvals[i] = 1.0;
        xvals[n - i] = yvals[i];
        yvals[n - i] = xvals[i];
      }
    }
    /* n is odd, so get one corner-cut chord. */
    else {
      double temp = 1.0 / (n2 + math::numbers::sqrt2 / 2);
      for (int i = 0; i <= n2; i++) {
        xvals[i] = double(i) * temp;
        yvals[i] = 1.0;
        xvals[n - i] = yvals[i];
        yvals[n - i] = xvals[i];
      }
    }
    return;
  }
  /* For general case use the more expensive search algorithm. */
  find_even_superellipse_chords_general(n, r, xvals, yvals);
}

/**
 * Find the profile's "fullness," which is the fraction of the space it takes up way from the
 * boundvert's centroid to the original vertex for a non-custom profile, or in the case of a
 * custom profile, the average "height" of the profile points along its centerline.
 */
static float find_profile_fullness(BevelState *bs)
{
  int nseg = bs->params.segments;

  /* Precalculated fullness for circle profile radius and more common low seg values. */
  constexpr int circle_fullness_segs = 11;
  static const float circle_fullness[circle_fullness_segs] = {
      0.0f,   /* nsegs == 1 */
      0.559f, /* 2 */
      0.642f, /* 3 */
      0.551f, /* 4 */
      0.646f, /* 5 */
      0.624f, /* 6 */
      0.646f, /* 7 */
      0.619f, /* 8 */
      0.647f, /* 9 */
      0.639f, /* 10 */
      0.647f, /* 11 */
  };

  float fullness;
  if (bs->params.profile_type == BevelProfileType::Custom) {
    /* Set fullness to the average "height" of the profile's sampled points. */
    fullness = 0.0f;
    for (int i = 0; i < nseg; i++) { /* Don't use the end points. */
      fullness += float(bs->pro_spacing.xvals[i] + bs->pro_spacing.yvals[i]) / (2.0f * nseg);
    }
  }
  else {
    /* An offline optimization process found fullness that led to closest fit to sphere as
     * a function of r and ns (for case of cube corner). */
    if (bs->pro_super_r == pro_line_r) {
      fullness = 0.0f;
    }
    else if (bs->pro_super_r == pro_circle_r && nseg > 0 && nseg <= circle_fullness_segs) {
      fullness = circle_fullness[nseg - 1];
    }
    else {
      /* Linear regression fit found best linear function, separately for even/odd segs. */
      if (nseg % 2 == 0) {
        fullness = 2.4506f * bs->params.profile - 0.00000300f * nseg - 0.6266f;
      }
      else {
        fullness = 2.3635f * bs->params.profile + 0.000152f * nseg - 0.6060f;
      }
    }
  }
  return fullness;
}

/**
 * Fill matrix r_mat so that a point in the sheared parallelogram with corners
 * va, vmid, vb (and the 4th that is implied by it being a parallelogram)
 * is the result of transforming the unit square by multiplication with r_mat.
 * If it can't be done because the parallelogram is degenerate, return false,
 * else return true.
 * Method:
 * Find vo, the origin of the parallelogram with other three points va, vmid, vb.
 * Also find vd, which is in direction normal to parallelogram and 1 unit away
 * from the origin.
 * The quarter circle in first quadrant of unit square will be mapped to the
 * quadrant of a sheared ellipse in the parallelogram, using a matrix.
 * The matrix mat is calculated to map:
 *    (0,1,0) -> va
 *    (1,1,0) -> vmid
 *    (1,0,0) -> vb
 *    (0,1,1) -> vd
 * We want M to make M*A=B where A has the left side above, as columns
 * and B has the right side as columns - both extended into homogeneous coords.
 * So M = B*(Ainverse).  Doing Ainverse by hand gives the code below.
 */
static bool make_unit_square_map(const float3 va,
                                 const float3 vmid,
                                 const float3 vb,
                                 float4x4 &r_mat)
{
  const float3 va_vmid = vmid - va;
  const float3 vb_vmid = vmid - vb;

  if (math::is_zero(va_vmid) || math::is_zero(vb_vmid)) {
    return false;
  }

  if (math::abs(angle_v3v3(va_vmid, vb_vmid) - math::numbers::pi) <= geom::bevel_epsilon_ang) {
    return false;
  }

  const float3 vo = va - vb_vmid;
  const float3 vddir = math::normalize(math::cross(vb_vmid, va_vmid));
  const float3 vd = vo + vddir;

  /* The cols of m are: `vmid - va, vmid - vb, vmid + vd - va -vb, va + vb - vmid`;
   * Blender transform matrices are stored such that `m[i][*]` is `i-th` column;
   * the last elements of each col remain as they are in unity matrix. */
  const float3 col0 = vmid - va;
  const float3 col1 = vmid - vb;
  const float3 col2 = vmid + vd - va - vb;
  const float3 col3 = va + vb - vmid;
  r_mat[0] = float4(col0[0], col0[1], col0[2], 0.0f);
  r_mat[1] = float4(col1[0], col1[1], col1[2], 0.0f);
  r_mat[2] = float4(col2[0], col2[1], col2[2], 0.0f);
  r_mat[3] = float4(col3[0], col3[1], col3[2], 1.0f);

  return true;
}

/**
 * Helper for #calculate_profiles that builds the 3D locations for the segments
 * and the higher power of 2 segments.
 */
static void calculate_profile_segments(const Profile &profile,
                                       const float4x4 map,
                                       const bool use_map,
                                       const bool reversed,
                                       const int ns,
                                       Span<double> xvals,
                                       Span<double> yvals,
                                       MutableSpan<float3> r_prof_co)
{
  /* Iterate over the vertices along the boundary arc. */
  for (int k = 0; k <= ns; k++) {
    float3 co;
    if (k == 0) {
      co = profile.start;
    }
    else if (k == ns) {
      co = profile.end;
    }
    else {
      if (use_map) {
        const float3 p = reversed ? float3(yvals[ns - k], xvals[ns - k], 0.0f) :
                                    float3(xvals[k], yvals[k], 0.0f);
        /* Do the 2D->3D transformation of the profile coordinates. */
        co = math::transform_point(map, p);
      }
      else {
        co = math::interpolate(profile.start, profile.end, float(k) / float(ns));
      }
    }
    /* Finish the 2D->3D transformation by projecting onto the final profile plane. */
    float3 &prof_co_k = r_prof_co[k];
    if (!math::is_zero(profile.proj_dir)) {
      float3 co2 = co + profile.proj_dir;
      if (!isect_line_plane_v3(prof_co_k, co, co2, profile.plane_co, profile.plane_no)) {
        /* Shouldn't happen. */
        prof_co_k = co;
      }
    }
    else {
      prof_co_k = co;
    }
  }
}

/**
 * Fills the ProfileSpacing struct with the 2D coordinates for the profile's vertices.
 * The superellipse used for multi-segment profiles does not have a closed-form way
 * to generate evenly spaced points along an arc. We use an expensive search procedure
 * to find the parameter values that lead to bp->seg even chords.
 * We also want spacing for a number of segments that is a power of 2 >= bp->seg (but at least
 * 4). Use doubles because otherwise we cannot come close to float precision for final results.
 *
 * \param pro_spacing: The struct to fill. Changes depending on whether there needs
 * to be a separate miter profile.
 */
static void set_profile_spacing(BevelState *bs, ProfileSpacing *pro_spacing, bool custom)
{
  int segments = bs->params.segments;

  if (segments <= 1) {
    /* Only 1 segment, we don't need any profile information. */
    pro_spacing->segments_power_2 = 0;
    return;
  }

  int segments_power_2 = std::max(power_of_2_max_i(bs->params.segments), 4);

  /* Sample the seg_2 segments used during vertex mesh subdivision. */
  bs->pro_spacing.segments_power_2 = segments_power_2;
  if (segments_power_2 == segments) {
    pro_spacing->xvals_2 = pro_spacing->xvals;
    pro_spacing->yvals_2 = pro_spacing->yvals;
  }
  else {
    pro_spacing->xvals_2 = Array<double>(segments_power_2 + 1);
    pro_spacing->yvals_2 = Array<double>(segments_power_2 + 1);
    if (custom) {
      /* Make sure the curve profile widget's sample table is full of the segments_power_2
       * samples.
       */
      BKE_curveprofile_init((CurveProfile *)bs->params.custom_profile, short(segments_power_2));

      /* Copy segment locations into the profile spacing struct. */
      for (const int i : IndexRange(segments_power_2 + 1)) {
        pro_spacing->xvals_2[i] = double(bs->params.custom_profile->segments[i].y);
        pro_spacing->yvals_2[i] = double(bs->params.custom_profile->segments[i].x);
      }
    }
    else {
      find_even_superellipse_chords(
          segments_power_2, bs->pro_super_r, pro_spacing->xvals_2, pro_spacing->yvals_2);
    }
  }

  /* Sample the input number of segments. */
  pro_spacing->xvals = Array<double>(segments + 1);
  pro_spacing->yvals = Array<double>(segments + 1);
  if (custom) {
    /* Make sure the curve profile's sample table is full. */
    if (bs->params.custom_profile->segments_len != segments ||
        !bs->params.custom_profile->segments)
    {
      BKE_curveprofile_init((CurveProfile *)bs->params.custom_profile, short(segments));
    }

    /* Copy segment locations into the profile spacing struct. */
    for (const int i : IndexRange(segments + 1)) {
      pro_spacing->xvals[i] = double(bs->params.custom_profile->segments[i].y);
      pro_spacing->yvals[i] = double(bs->params.custom_profile->segments[i].x);
    }
  }
  else {
    find_even_superellipse_chords(
        segments, bs->pro_super_r, pro_spacing->xvals, pro_spacing->yvals);
  }
}

static void calculate_profiles(const int bv, AnchorProfiles &profiles, const BevelState &bs)
{
  if (bs.params.segments == 1) {
    /* Profiles are unnecessary for 1-segment bevels. */
    return;
  }
  const MeshPattern &pat = bs.bevvert_meshpatterns()[bv];
  Span<int> edges = bs.bevvert_bevedges()[bv];
  IndexRange newverts = bs.bevvert_newverts()[bv];
  const float3 bv_v_pos = bs.mesh_info.mesh.vert_positions()[bs.bevvert_mesh_verts()[bv]];
  for (const int a : IndexRange(pat.num_anchors)) {
    Profile &profile = profiles[a];
    /* First fill in the profile parameters. */
    const int apos = pat.anchor_vert(a);
    const int anext = (a + 1) % pat.num_anchors;
    const int aposnext = pat.anchor_vert(anext);
    profile.start = bs.newvert_positions()[newverts[apos]];
    profile.end = bs.newvert_positions()[newverts[aposnext]];
    /* Fallback will be linear interopolation bewtween start and end. */
    profile.middle = math::interpolate(profile.start, profile.end, 0.5f);
    profile.super_r = pro_line_r;
    profile.plane_co = float3(0.0, 0.0, 0.0);
    profile.plane_no = float3(0.0, 0.0, 0.0);
    profile.proj_dir = float3(0.0, 0.0, 0.0);
    if (bs.params.affect_type == BevelAffect::Vertices) {
      profile.middle = bv_v_pos;
      continue;
    }
    const int be_pos = bs.last_attached_bevedge_pos(bv, apos);
    if (be_pos == -1) {
      /* TODO: this shouldn't happen. */
      BLI_assert(false);
      continue;
    }
    const int be = edges[be_pos];
    if (bs.bevedge_is_beveled(be)) {
      /* Projection direction is along be towards bv. */
      profile.proj_dir = -math::normalize(geom::edge_dir(bv, be_pos, bs));
      profile.middle = geom::project_to_edge(bv, be_pos, profile.start, profile.end, bs);
      /* Usual plane to project to is the one containing start, middle, and end. */
      const float3 d1 = math::normalize(profile.middle - profile.start);
      const float3 d2 = math::normalize(profile.middle - profile.end);
      if (!geom::nearly_parallel_normalized(d1, d2)) {
        /* Usual plane is fine. */
        profile.plane_no = math::normalize(math::cross(d1, d2));
        profile.plane_co = profile.start;
        profile.super_r = bs.pro_super_r;
      }
      else {
        /* It seems that the beveled edge is coplanar with the two anchor verts.
         * We want to make that plane the profile plane, if possible.
         */
        const int be_prev_pos = be_pos == 0 ? edges.size() - 1 : be_pos - 1;
        const int be_next_pos = be_pos == edges.size() - 1 ? 0 : be_pos + 1;
        const int be_prev = edges[be_prev_pos];
        const int be_next = edges[be_next_pos];
        if (bs.bevedge_is_beveled(be_prev) && bs.bevedge_is_beveled(be_next) && be_prev != be_next)
        {
          const float3 d3 = math::normalize(geom::edge_dir(bv, be_prev_pos, bs));
          const float3 d4 = math::normalize(geom::edge_dir(bv, be_next_pos, bs));
          if (!geom::nearly_parallel_normalized(d3, d4)) {
            const float3 co3 = profile.start + d3;
            const float3 co4 = profile.end + d4;
            float3 meetco, isect2;
            if (isect_line_line_v3(profile.start, co3, profile.end, co4, meetco, isect2)) {
              profile.middle = meetco;
              profile.super_r = bs.pro_super_r;
            }
          }
        }
        else {
          profile.middle = bv_v_pos;
          profile.super_r = bs.pro_super_r;
        }
        if (profile.super_r != pro_line_r) {
          const float3 d5 = math::normalize(profile.middle - profile.start);
          const float3 d6 = math::normalize(profile.middle - profile.end);
          if (!geom::nearly_parallel_normalized(d5, d6)) {
            profile.plane_no = math::normalize(math::cross(d5, d6));
            profile.plane_co = profile.start;
          }
        }
      }
    }
    /* TODO: miters */

    /* Now that the profile parameters are set, we can calculate the positions
     */
    profile.prof_co.reinitialize(bs.params.segments + 1);
    profile.prof_co_2.reinitialize(bs.pro_spacing.segments_power_2 + 1);
    float4x4 map;
    bool use_map = bs.params.profile_type == BevelProfileType::Superellipse &&
                   profile.super_r != pro_line_r;
    if (use_map) {
      use_map = make_unit_square_map(profile.start, profile.middle, profile.end, map);
    }
    /* TODO: cutoff method. */
    /* Calculate the 3D locations for the profile points. */
    calculate_profile_segments(profile,
                               map,
                               use_map,
                               false,
                               bs.params.segments,
                               bs.pro_spacing.xvals,
                               bs.pro_spacing.yvals,
                               profile.prof_co.as_mutable_span());
    if (bs.params.segments != bs.pro_spacing.segments_power_2) {
      calculate_profile_segments(profile,
                                 map,
                                 use_map,
                                 false,
                                 bs.pro_spacing.segments_power_2,
                                 bs.pro_spacing.xvals_2,
                                 bs.pro_spacing.yvals_2,
                                 profile.prof_co_2.as_mutable_span());
    }
    else {
      std::copy(profile.prof_co.begin(), profile.prof_co.end(), profile.prof_co_2.begin());
    }
  }
}

/**
 * Find the point on given profile at parameter \a at_index which goes from 0 to \a nseg as
 * the profile moves from `profile.start` to `profile.end`.
 * We assume that nseg is either the profile's seg number or a power of 2 less than
 * or equal to the power of 2 >= seg.
 * In the latter case, we subsample the profile.prof_co_2, which will not necessarily
 * give equal spaced chords, but is in fact more what is desired by the cubic subdivision
 * method used to make the adj pattern.
 */
static float3 get_profile_point(const Profile &profile, const int at_index, const int nseg)
{
  const int profile_nseg = profile.prof_co.size() - 1;
  if (profile_nseg == 1) {
    return at_index == 0 ? profile.start : profile.end;
  }
  if (nseg == profile_nseg) {
    return profile.prof_co[at_index];
  }
  BLI_assert(is_power_of_2_i(nseg) && nseg < profile.prof_co_2.size());
  const int subsample_spacing = (profile.prof_co_2.size() - 1) / nseg;
  return profile.prof_co_2[at_index * subsample_spacing];
}

}  // end namespace profile

/** Hole Filling Mesh Pattern. */

namespace adj {
/** The "Adj" pattern is the only really complicated pattern.
 * The pattern for nv anchor verts and ns segments is a set of
 * concentric rings of vertices, which form rings of faces.
 * The outer most ring has a skeleton of nv anchor verts, and
 * between successive anchors there are ns edge segments and
 * hence (ns - 2) non-anchor verts (sometimes we'll call them
 * "span" verts) between the successive anchor verts.
 *
 * Each successive inner ring has 2 less verts than the next
 * outward ring.  If ns is odd, we end up with an nv-gon as
 * the innermost ring; if ns is even, we end up with a single
 * vertex in the innermost ring.
 *
 * The edges that go between the verts of the same ring are
 * called "ring edges". There are also edges that go between
 * two successive rings, called "cross-ring" edges. The span
 * vertices are connected one-to-one between the rings, while
 * an anchor vertex has edges to the just-preceding and
 * just-succeeding span vert in the next outer ring.
 * The ring edges and cross ring edges together form rings
 * of faces, which are all quads except possibly the center
 * face, which is an nv-gon ns is odd.
 *
 * The verts, edges, and faces are indexed by starting at the
 * innermost ring and going counterclockwise from the 0th
 * anchor vertex, and continuing with successive rings.
 * The edges index the ring edges first, then the cross ring
 * edges.
 *
 * TODO: link to external doc that has a picture of all this.
 */

static inline int odd(int i)
{
  return i % 2 == 1;
}

static inline int v_ringstart_odd(int r, int nv)
{
  return nv * r * r;
}

static inline int v_ringlen_odd(int r, int nv)
{
  return nv * (2 * r + 1);
}

static inline int v_num_rings_odd(int ns)
{
  return (ns - 1) / 2 + 1;
}

static inline int v_ringstart_even(int r, int nv)
{
  return r == 0 ? 0 : 1 + nv * (r - 1) * r;
}

static inline int v_ringlen_even(int r, int nv)
{
  return r == 0 ? 1 : 2 * nv * r;
}

static inline int v_num_rings_even(int ns)
{
  return ns / 2 + 1;
}

static int v_ringstart(int r, int nv, int ns)
{
  return odd(ns) ? v_ringstart_odd(r, nv) : v_ringstart_even(r, nv);
}

static int v_ringlen(int r, int nv, int ns)
{
  return odd(ns) ? v_ringlen_odd(r, nv) : v_ringlen_even(r, nv);
}

static int v_num_rings(int ns)
{
  return odd(ns) ? v_num_rings_odd(ns) : v_num_rings_even(ns);
}

static int v_total_verts(int nv, int ns)
{
  return v_ringstart(v_num_rings(ns), nv, ns);
}

static inline int f_ringstart_odd(int r, int nv)
{
  return v_ringstart_even(r, nv);
}

static inline int f_ringlen_odd(int r, int nv)
{
  return v_ringlen_even(r, nv);
}

static inline int f_num_rings_odd(int ns)
{
  return (ns - 1) / 2 + 1;
}

static inline int f_ringstart_even(int r, int nv)
{
  return v_ringstart_odd(r, nv);
}

static inline int f_ringlen_even(int r, int nv)
{
  return v_ringlen_odd(r, nv);
}

static inline int f_num_rings_even(int ns)
{
  return ns / 2;
}

static int f_ringstart(int r, int nv, int ns)
{
  return odd(ns) ? f_ringstart_odd(r, nv) : f_ringstart_even(r, nv);
}

static int f_ringlen(int r, int nv, int ns)
{
  return odd(ns) ? f_ringlen_odd(r, nv) : f_ringlen_even(r, nv);
}

static int f_num_rings(int ns)
{
  return odd(ns) ? f_num_rings_odd(ns) : f_num_rings_even(ns);
}

static int f_total_faces(int nv, int ns)
{
  return f_ringstart(f_num_rings(ns), nv, ns);
}

static inline int e_ringstart_odd(int r, int nv)
{
  return v_ringstart_odd(r, nv);
}

static inline int e_ringlen_odd(int r, int nv)
{
  return v_ringlen_odd(r, nv);
}

static inline int e_num_rings_odd(int ns)
{
  return (ns - 1) / 2 + 1;
}

static inline int e_crossring_start_odd(int r, int nv, int ns)
{
  return e_ringstart_odd(e_num_rings_odd(ns), nv) + nv * r * (r + 1);
}

static inline int e_crossring_len_odd(int r, int nv)
{
  return 2 * nv * (r + 1);
}

static inline int e_num_crossrings_odd(int ns)
{
  return ns > 1 ? (ns - 1) / 2 : 0;
}

static inline int e_ringstart_even(int r, int nv)
{
  return v_ringstart_even(r + 1, nv) - 1;
}

static inline int e_ringlen_even(int r, int nv)
{
  return v_ringlen_even(r + 1, nv);
}

static inline int e_num_rings_even(int ns)
{
  return ns / 2;
}

static inline int e_crossring_start_even(int r, int nv, int ns)
{
  return e_ringstart_even(e_num_rings_even(ns), nv) + nv * r * r;
}

static inline int e_crossring_len_even(int r, int nv)
{
  return nv * (2 * r + 1);
}

static inline int e_num_crossrings_even(int ns)
{
  return ns / 2;
}

static int e_ringstart(int r, int nv, int ns)
{
  return odd(ns) ? e_ringstart_odd(r, nv) : e_ringstart_even(r, nv);
}

static int e_ringlen(int r, int nv, int ns)
{
  return odd(ns) ? e_ringlen_odd(r, nv) : e_ringlen_even(r, nv);
}

static int e_num_rings(int ns)
{
  return odd(ns) ? e_num_rings_odd(ns) : e_num_rings_even(ns);
}

static int e_crossring_start(int r, int nv, int ns)
{
  return odd(ns) ? e_crossring_start_odd(r, nv, ns) : e_crossring_start_even(r, nv, ns);
}

static int e_crossring_len(int r, int nv, int ns)
{
  return odd(ns) ? e_crossring_len_odd(r, nv) : e_crossring_len_even(r, nv);
}

static int e_num_crossrings(int ns)
{
  return odd(ns) ? e_num_crossrings_odd(ns) : e_num_crossrings_even(ns);
}

static int e_num_edges(int nv, int ns)
{
  return e_crossring_start(e_num_crossrings(ns), nv, ns);
}

static int face_ring(int f, int nv, int ns)
{
  if (odd(ns)) {
    return f < 1 ? 0 : int(((nv + math::sqrt(nv * nv + 4 * nv * (f - 1))) / (2 * nv)));
  }
  else {
    return f < 0 ? 0 : int(math::sqrt(f / nv));
  }
}

static int vertex_ring(int v, int nv, int ns)
{
  if (odd(ns)) {
    return v < 0 ? 0 : int(math::sqrt(v / nv));
  }
  else {
    return v < 1 ? 0 : int(((nv + math::sqrt(nv * nv + 4 * nv * (v - 1))) / (2 * nv)));
  }
}

static inline int v_anchor_div(int r, int nv, int ns)
{
  return v_ringlen(r, nv, ns) / nv;
}

static inline int f_anchor_div(int r, int nv, int ns)
{
  return f_ringlen(r, nv, ns) / nv;
}

[[maybe_unused]] static int3 v_ring_anchor_offset(int v, int nv, int ns)
{
  int ring = vertex_ring(v, nv, ns);
  int ring_offset = v - v_ringstart(ring, nv, ns);
  int anchor_index = ring_offset / v_anchor_div(ring, nv, ns);
  int offset = ring_offset % v_anchor_div(ring, nv, ns);
  return int3(ring, anchor_index, offset);
}

static int3 f_ring_anchor_offset(int f, int nv, int ns)
{
  int ring = face_ring(f, nv, ns);
  int ring_offset = f - f_ringstart(ring, nv, ns);
  int anchor_index = ring_offset / f_anchor_div(ring, nv, ns);
  int offset = ring_offset % f_anchor_div(ring, nv, ns);
  return int3(ring, anchor_index, offset);
}

static int rao_to_vert(int r, int a, int o, int nv, int ns)
{
  int delta = (v_anchor_div(r, nv, ns) * a + o) % v_ringlen(r, nv, ns);
  return v_ringstart(r, nv, ns) + delta;
}

[[maybe_unused]] static int rao_to_face(int r, int a, int o, int nv, int ns)
{
  int delta = (f_anchor_div(r, nv, ns) * a + o) % f_ringlen(r, nv, ns);
  return f_ringstart(r, nv, ns) + delta;
}

static int face_outer_vertex_ring(int f, int nv, int ns)
{
  return odd(ns) ? face_ring(f, nv, ns) : face_ring(f, nv, ns) + 1;
}

static std::pair<int4, int4> face_vertices_and_edges(int f, int nv, int ns)
{
  if (f == 0 and odd(ns)) {
    /* The center ngon has vertex and edge indices 0, ..., nv-1,
     * and we expect the caller to know this as this API doesn't let
     * us return an ngon. */
    return std::pair<int4, int4>(int4(0, 1, 2, 3), int4(0, 1, 2, 3));
  }
  auto [fr, fa, fo] = f_ring_anchor_offset(f, nv, ns);
  int outer_r = face_outer_vertex_ring(f, nv, ns);
  int inner_r = outer_r - 1;
  int vr = outer_r;
  int va = fa;
  int vo = fo;
  int vroot = rao_to_vert(vr, va, vo, nv, ns);
  int vafter = vroot + 1;
  int vopposite, vbefore;
  if (fo == 0) {
    vopposite = rao_to_vert(inner_r, va, vo, nv, ns);
    vbefore = va > 0 ? vroot - 1 : vroot + v_ringlen(outer_r, nv, ns) - 1;
  }
  else {
    vopposite = rao_to_vert(inner_r, va, vo, nv, ns);
    vbefore = vopposite > v_ringstart(inner_r, nv, ns) ?
                  vopposite - 1 :
                  v_ringstart(inner_r, nv, ns) + v_ringlen(inner_r, nv, ns) - 1;
  }
  int4 vertices = int4(vroot, vafter, vopposite, vbefore);
  int e0 = odd(ns) ? vroot : vroot - 1;
  int f_ring_pos = f - f_ringstart(fr, nv, ns);
  int e_cross_r = odd(ns) ? fr - 1 : fr;
  int e_outer_r = fr;
  int e_inner_r = e_outer_r - 1;
  int e1 = e_crossring_start(e_cross_r, nv, ns) + f_ring_pos;
  int e2, e3;
  if (fo == 0) {
    e2 = fa > 0 ? e1 - 1 : e1 + e_crossring_len(e_cross_r, nv, ns) - 1;
    e3 = e0 > e_ringstart(e_outer_r, nv, ns) ?
             e0 - 1 :
             e_ringstart(e_outer_r, nv, ns) + e_ringlen(e_outer_r, nv, ns) - 1;
  }
  else {
    if (odd(ns)) {
      e2 = vopposite > e_ringstart(e_inner_r, nv, ns) ?
               vopposite - 1 :
               e_ringstart(e_inner_r, nv, ns) + e_ringlen(e_inner_r, nv, ns) - 1;
    }
    else {
      e2 = vopposite - 1 > e_ringstart(e_inner_r, nv, ns) ?
               vopposite - 2 :
               e_ringstart(e_inner_r, nv, ns) + e_ringlen(e_inner_r, nv, ns) - 1;
    }
    e3 = e1 - 1;
  }
  int4 edges = int4(e0, e1, e2, e3);
  return std::pair<int4, int4>(vertices, edges);
}

[[maybe_unused]] static void print_adj_pattern(int nv, int ns)
{
  fmt::println("\nnv = {}, ns = {}", nv, ns);
  int nvrings = v_num_rings(ns);
  fmt::println("verts");
  fmt::println("v_num_rings = {}", nvrings);
  for (const int r : IndexRange(nvrings)) {
    fmt::println("r = {}: start = {}, len = {}", r, v_ringstart(r, nv, ns), v_ringlen(r, nv, ns));
  }
  fmt::println("total verts = {}", v_total_verts(nv, ns));
  fmt::println("faces");
  int nfrings = f_num_rings(ns);
  fmt::println("f_num_rings = {}", nfrings);
  for (const int r : IndexRange(nfrings)) {
    fmt::println("r = {}: start = {}, len = {}", r, f_ringstart(r, nv, ns), f_ringlen(r, nv, ns));
  }
  fmt::println("total faces = {}", f_total_faces(nv, ns));
  fmt::println("edges");
  int nerings = e_num_rings(ns);
  fmt::println("e_num_rings = {}", nerings);
  for (const int r : IndexRange(nerings)) {
    fmt::println("r = {}: start = {}, len = {}", r, e_ringstart(r, nv, ns), e_ringlen(r, nv, ns));
  }
  int necrossrings = e_num_crossrings(ns);
  fmt::println("e_num_crossrings = {}", necrossrings);
  for (const int r : IndexRange(necrossrings)) {
    fmt::println("r = {}: start = {}, len = {}",
                 r,
                 e_crossring_start(r, nv, ns),
                 e_crossring_len(r, nv, ns));
  }
  fmt::println("total edges = {}", e_num_edges(nv, ns));
  fmt::println("");
  for (const int f : IndexRange(f_total_faces(nv, ns))) {
    auto [fv, fe] = face_vertices_and_edges(f, nv, ns);
    fmt::println("f = {}, fverts = ({},{},{},{}), fedges = ({},{},{},{})",
                 f,
                 fv[0],
                 fv[1],
                 fv[2],
                 fv[3],
                 fe[0],
                 fe[1],
                 fe[2],
                 fe[3]);
  }
}

/** A structure to hold the vertex positions of an Adj pattern. */
struct AdjVerts {
  /** The coordinates of vertices arranged according to indexing given above.
   * The inline size 125 will accommodate anchors=4, segments=10. */
  Array<float3, 125> verts;
  /** How many anchors (often abbreviated nv). */
  int anchors;
  /** How many segments (often abbreviated ns). */
  int segments;

  AdjVerts(const int anchors, const int segments) : anchors(anchors), segments(segments)
  {
    verts.reinitialize(v_total_verts(anchors, segments));
    verts.fill(float3(0.0f, 0.0f, 0.0f));
  }

  const int anchor_offset_to_outer_ring_vert(const int anchor, const int offset) const;

  const float3 &outer_ring_vert(const int anchor, const int offset) const
  {
    return verts[anchor_offset_to_outer_ring_vert(anchor, offset)];
  }

  float3 &mutable_outer_ring_vert(const int anchor, const int offset)
  {
    return verts[anchor_offset_to_outer_ring_vert(anchor, offset)];
  }

  const int ring_anchor_offset_to_vert(const int ring, const int anchor, const int offset) const;

  const float3 &vert(const int ring, const int anchor, const int offset) const
  {
    return verts[ring_anchor_offset_to_vert(ring, anchor, offset)];
  }

  float3 &mutable_vert(const int ring, const int anchor, const int offset)
  {
    return verts[ring_anchor_offset_to_vert(ring, anchor, offset)];
  }
};

const int AdjVerts::anchor_offset_to_outer_ring_vert(const int anchor, const int offset) const
{
  const int ring = v_num_rings(this->segments) - 1;
  return rao_to_vert(ring, anchor, offset, this->anchors, this->segments);
}

const int AdjVerts::ring_anchor_offset_to_vert(const int ring,
                                               const int anchor,
                                               const int offset) const
{
  return rao_to_vert(ring, anchor, offset, this->anchors, this->segments);
}

[[maybe_unused]] static void draw_adj(const AdjVerts &adjverts)
{
  constexpr uint life = draw::drw_debug_persistent_lifetime;
  for (const int i : adjverts.verts.index_range()) {
    float3 co = adjverts.verts[i];
    draw::drw_debug_point(co, 0.01f, {1, 0, 0, 1}, life);
  }
}

static float3 avg4(const float3 &v0, const float3 &v1, const float3 &v2, const float3 &v3)
{
  return 0.25f * (v0 + v1 + v2 + v3);
}

/* Gamma needed for smooth Catmull-Clark, Sabin modification. */
static float sabin_gamma(int n)
{
  /* Precalculated for common cases of n. */
  if (n < 3) {
    return 0.0f;
  }
  if (n == 3) {
    return 0.065247584f;
  }
  if (n == 4) {
    return 0.25f;
  }
  if (n == 5) {
    return 0.401983447f;
  }
  if (n == 6) {
    return 0.523423277f;
  }
  const double k = cos(math::numbers::pi / double(n));
  /* Need x, real root of x^3 + (4k^2 - 3)x - 2k = 0.
   * Answer calculated via Wolfram Alpha. */
  const double k2 = k * k;
  const double k4 = k2 * k2;
  const double k6 = k4 * k2;
  const double y = pow(
      math::numbers::sqrt3 * math::sqrt(64.0 * k6 - 144.0 * k4 + 135.0 * k2 - 27.0) + 9.0 * k,
      1.0 / 3.0);
  const double x = 0.480749856769136 * y - (0.231120424783545 * (12.0 * k2 - 9.0)) / y;
  return (k * x + 2.0 * k2 - 1.0) / (x * x * (k * x + 1.0));
}

/* Fill \a adjverts using recursive cubic subdivision, until reach the
 * base case where \a adjverts_2_segs is the answer.
 * adjverts.segments should be a power of 2, and >= 2.
 */
static void fill_adjverts(AdjVerts &adjverts,
                          const AdjVerts &adjverts_2_segs,
                          const profile::AnchorProfiles &profiles)
{
  BLI_assert(adjverts_2_segs.segments == 2 && adjverts.segments % 2 == 0);
  BLI_assert(adjverts.segments < 1000); /* TODO: stop too-deep recursion. */
  if (adjverts.segments < 2) {
    return;
  }
  if (adjverts.segments == 2) {
    std::copy(adjverts_2_segs.verts.begin(), adjverts_2_segs.verts.end(), adjverts.verts.begin());
    return;
  }
  int nhalf = adjverts.segments / 2;
  BLI_assert(nhalf % 2 == 0);
  AdjVerts adjverts_half = AdjVerts(adjverts.anchors, nhalf);
  fill_adjverts(adjverts_half, adjverts_2_segs, profiles);
  /* Do a step of cubic subdivision (Catmull-Clark) with special rules
   * at boundaries. See Levin 1999 paper "Filling an N-sided hole using combined
   * subdivision schemes".
   */
  const int n_boundary = adjverts_half.anchors;
  const int ns_in = adjverts_half.segments;
  const int ns_out = adjverts.segments;

  /* First adjust the boundary vertices of the input, storing in the output. */
  for (const int a : IndexRange(n_boundary)) {
    adjverts.mutable_outer_ring_vert(a, 0) = adjverts_half.outer_ring_vert(a, 0);
    for (const int o : IndexRange::from_begin_end(1, ns_in)) {
      float3 co = adjverts_half.outer_ring_vert(a, o);
      /* Smooth boundary rule for even verts. TODO: Custom profiles. */
      const float3 co1 = adjverts_half.outer_ring_vert(a, o - 1);
      const float3 co2 = adjverts_half.outer_ring_vert(a, o + 1);
      co = co - (1.0f / 6.0f) * (co1 + co2 - 2.0f * co);
      adjverts.mutable_outer_ring_vert(a, 2 * o) = co;
    }
  }
  /* Now set odd boundary verts, using input profiles. */
  for (const int a : IndexRange(n_boundary)) {
    const profile::Profile &profile = profiles[a];
    for (int o = 1; o < ns_out; o += 2) {
      float3 co = profile::get_profile_point(profile, o, ns_out);
      /* Smooth boundary rule for odd verts. TODO: Custom profiles. */
      const float3 co1 = adjverts.outer_ring_vert(a, o - 1);
      const float3 co2 = adjverts.outer_ring_vert(a, o + 1);
      co = co - (1.0f / 6.0f) * (co1 + co2 - 2.0f * co);
      adjverts.mutable_outer_ring_vert(a, o) = co;
    }
  }
  /* Copy adjusted boundary verts back into adjverts_half, prior to subdivision. */
  for (const int a : IndexRange(n_boundary)) {
    for (const int o : IndexRange(ns_in)) {
      adjverts_half.mutable_outer_ring_vert(a, o) = adjverts.outer_ring_vert(a, 2 * o);
    }
  }
  /* Now we do the internal vertices, using standard Catmull-Clark
   * and assuming all boundary vertices have valence 4. */

  /* The new face-center vertices.
   * Made as average of four corners of quads of the input mesh.
   */
  const int nrings_in = v_num_rings(ns_in);
  for (const int r : IndexRange::from_begin_end(1, nrings_in)) {
    const int side = v_anchor_div(r, n_boundary, ns_in);
    for (const int a : IndexRange(n_boundary)) {
      const int aprev = a == 0 ? n_boundary - 1 : a - 1;
      /* Corner faces have a different indexing pattern. */
      adjverts.mutable_vert(2 * r - 1, a, 0) = avg4(adjverts_half.vert(r, a, 0),
                                                    adjverts_half.vert(r, a, 1),
                                                    adjverts_half.vert(r - 1, r == 1 ? 0 : a, 0),
                                                    adjverts_half.vert(r, aprev, side - 1));
      for (int o = 1; o < side - 1; o++) {
        adjverts.mutable_vert(2 * r - 1, a, 2 * o) = avg4(adjverts_half.vert(r, a, o),
                                                          adjverts_half.vert(r, a, o + 1),
                                                          adjverts_half.vert(r - 1, a, o),
                                                          adjverts_half.vert(r - 1, a, o - 1));
      }
    }
  }

  /* The new cross-ring edge vertices.
   * Made as average of ends of cross-ring edges of the input mesh and the adjacent newly-made
   * face-center vertices of the output mesh left and right of that edge.
   */
  for (const int r : IndexRange::from_begin_end(1, nrings_in)) {
    const int side = v_anchor_div(r, n_boundary, ns_in);
    for (const int a : IndexRange(n_boundary)) {
      for (int o = 1; o < side; o++) {
        adjverts.mutable_vert(2 * r - 1, a, 2 * o - 1) = avg4(
            adjverts_half.vert(r, a, o),
            adjverts_half.vert(r - 1, a, o - 1),
            adjverts.vert(2 * r - 1, a, 2 * o - 2),
            adjverts.vert(2 * r - 1, a, 2 * o));
      }
    }
  }

  /* The new edge-ring edge vertices.
   * Made as the average of the ends of the edge-ring edges of input mesh and the
   * face-center vertices of the output mesh above and below that edge.
   */
  for (int r = 1; r < nrings_in - 1; r++) {
    const int side = v_anchor_div(r, n_boundary, ns_in);
    for (const int a : IndexRange(n_boundary)) {
      for (const int o : IndexRange(side)) {
        adjverts.mutable_vert(2 * r, a, 2 * o + 1) = avg4(adjverts_half.vert(r, a, o),
                                                          adjverts_half.vert(r, a, o + 1),
                                                          adjverts.vert(2 * r + 1, a, 2 * o + 2),
                                                          adjverts.vert(2 * r - 1, a, 2 * o));
      }
    }
  }

  /* The new vertex-vertex vertices (transformation of input interior vertices).
   * Made as a function of the average of the ends of the four new edges leading into the
   * vertex, and the average of the four new face vertices surrounding the vertex.
   */
  float gamma = sabin_gamma(4);
  float beta = -gamma;
  for (int r = 1; r < nrings_in - 1; r++) {
    const int iside = v_anchor_div(r, n_boundary, ns_in);
    const int oside = v_anchor_div(2 * r, n_boundary, ns_out);
    for (const int a : IndexRange(n_boundary)) {
      const int aprev = a == 0 ? n_boundary - 1 : a - 1;
      float3 ev_centroid = avg4(adjverts.vert(2 * r, aprev, oside - 1),
                                adjverts.vert(2 * r, a, 1),
                                adjverts.vert(2 * r + 1, aprev, oside + 1),
                                adjverts.vert(2 * r + 1, a, 1));
      float3 ef_centroid = avg4(adjverts.vert(2 * r - 1, a, 0),
                                adjverts.vert(2 * r + 1, a, 0),
                                adjverts.vert(2 * r + 1, a, 2),
                                adjverts.vert(2 * r + 1, aprev, oside));
      adjverts.mutable_vert(2 * r, a, 0) = ev_centroid + beta * ef_centroid +
                                           gamma * adjverts_half.vert(r, a, 0);
      for (int o = 1; o < iside; o++) {
        ev_centroid = avg4(adjverts.vert(2 * r, a, 2 * o - 1),
                           adjverts.vert(2 * r, a, 2 * o + 1),
                           adjverts.vert(2 * r + 1, a, 2 * o + 1),
                           adjverts.vert(2 * r - 1, a, 2 * o - 1));
        ef_centroid = avg4(adjverts.vert(2 * r - 1, a, 2 * o - 2),
                           adjverts.vert(2 * r - 1, a, 2 * o),
                           adjverts.vert(2 * r + 1, a, 2 * o),
                           adjverts.vert(2 * r + 1, a, 2 * o + 2));
        adjverts.mutable_vert(2 * r, a, 2 * o) = ev_centroid + beta * ef_centroid +
                                                 gamma * adjverts_half.vert(r, a, o);
      }
    }
  }

  /* The center vertex-vertex is like the ones calculated for the quads, but
   * the gamma is different because the valence is not necessarily 4.
   */
  gamma = sabin_gamma(n_boundary);
  beta = -gamma;
  /* Get centroids of new edge-verts and new face-verts adjacent to center. */
  float3 ev_centroid(0.0f, 0.0f, 0.0f);
  float3 ef_centroid(0.0f, 0.0f, 0.0f);
  for (const int a : IndexRange(n_boundary)) {
    ev_centroid = ev_centroid + adjverts.vert(1, a, 0);
    ef_centroid = ef_centroid + adjverts.vert(1, a, 1);
  }
  ev_centroid = ev_centroid / float(n_boundary);
  ef_centroid = ef_centroid / float(n_boundary);
  adjverts.mutable_vert(0, 0, 0) = ev_centroid + beta * ef_centroid +
                                   gamma * adjverts_half.vert(0, 0, 0);

  /* Final step: Copy the profile vertices to the output boundary. */
  for (const int a : IndexRange(n_boundary)) {
    const profile::Profile &profile = profiles[a];
    for (const int o : IndexRange(ns_out)) {
      adjverts.mutable_outer_ring_vert(a, o) = profile::get_profile_point(profile, o, ns_out);
    }
  }
}

/** Given a ruler whose values are \a div apart and start at 0, what is the index of the
 * tick on that ruler such that  \a value is between the tick of that index (inclusive) and
 * the next tick (exclusive), and also what is the remainder of value past that tick?
 * Assume value >= 0 and div > 0. */
static std::pair<int, float> tick_and_remainder(const float value, const float div)
{
  BLI_assert(value >= 0.0f && div > 0);
  const int tick = int(std::floor(value / div));
  return std::pair<int, float>(tick, value - tick * div);
}

/** Interpolate four coordinates by u along the co00 -> c01 direction and v along the c01 -> c10 direction. */
static float3 interp_bilinear_quad(const float3 &co00,
                                   const float3 &co01,
                                   const float3 &co11,
                                   const float3 &co10,
                                   const float u,
                                   const float v)
{
  return (1 - u) * (1 - v) * co00 +
         u * (1 - v) * co01 +
         u * v * co11 +
        (1 - u) * v * co10;
}

/** Interpolate given \a adjverts_in to make one with a different number of segments, in \a adjverts.
 * Assume they both have the same number of anchors, that adjverts_in has more segments,
 * and that adjverts_in has an even number of segments.
 */
static void interp_adj(AdjVerts &adjverts, const AdjVerts &adjverts_in, const profile::AnchorProfiles &profiles)
{
  const int na = adjverts_in.anchors;
  const int ns_in = adjverts_in.segments;
  const int ns_out = adjverts.segments;
  fmt::println("interp_adj, ns_in={} ns_out={}", ns_in, ns_out);
  BLI_assert(adjverts.anchors == na && (ns_in % 2) == 0 && ns_out < ns_in);
  const int num_rings_in = v_num_rings(ns_in);
  const int num_rings_out = v_num_rings(ns_out);
  fmt::println("num_rings_in={} num_rings_out={}", num_rings_in, num_rings_out);
  bool odd_out = (ns_out % 2) == 1;
  /* The radial_divs are the fractional spacing between rings.
   * For odd number of segments, there is a half space at the center. */
  const float radial_div_in = 2.0f / ns_in;
  const float radial_div_out = odd_out ? 1.0f / ((ns_out / 2) + 0.5f) : 2.0f / ns_out;
  fmt::println("div_in={}, div_out={}", radial_div_in, radial_div_out);
  for (int r_out = 0; r_out < num_rings_out - 1; r_out++) {
    if (r_out == 0 && !odd_out) {
      /* Just copy the center vertex. */
      adjverts.mutable_vert(0, 0, 0) = adjverts_in.vert(0, 0, 0);
      continue;
    }
    const float radial_frac_out = (num_rings_out - r_out - 1) * radial_div_out;
    auto [tick_in, remainder_in] = tick_and_remainder(radial_frac_out, radial_div_in);
    const int r_in = num_rings_in - 1 - tick_in;
    fmt::println("  r_out={} -> r_in={}, rem={}", r_out, r_in, remainder_in);
    for (const int a : IndexRange(na)) {
      const int iside = v_anchor_div(r_in, na, ns_in);
      const int aprev = a == 0 ? na - 1 : a - 1;
      float3 co00 = adjverts_in.vert(r_in, a, 0);
      float3 co01 = adjverts_in.vert(r_in, a, 1);
      float3 co11 = adjverts_in.vert(r_in - 1, r_in == 1 ? 0 : a, 0);
      float3 co10 = adjverts_in.vert(r_in, aprev, iside - 1);
      float3 co_interp = interp_bilinear_quad(co00, co01, co11, co10,
                                              remainder_in, remainder_in);
      adjverts.mutable_vert(r_out, a, 0) = co_interp;
    }
    for (const int a : IndexRange(na)) {
      const int iside = v_anchor_div(r_in, na, ns_in);
      const int oside = v_anchor_div(r_out, na, ns_out);
      fmt::println("  a={} iside={} oside={}", a, iside, oside);
      const int anext = a == na - 1 ? 0 : na + 1;
      const float ilen = math::length(adjverts_in.vert(r_in, anext, 0) - adjverts_in.vert(r_in, a, 0));
      const float idiv =  ilen / iside;
      const float olen = math::length(adjverts.vert(r_out, anext, 0) - adjverts.vert(r_out, a, 0));
      const float odiv =  olen / oside;
      for (int o_out = 1; o_out < oside; o_out++) {
        const float o_out_len = (o_out + remainder_in) * odiv;
        auto [o_in, rem] = tick_and_remainder(o_out_len, idiv);
        const float rem_frac = rem / idiv;
        fmt::println("  o_out={}, o_in={}, rem_frac={}", o_out, o_in, rem_frac);
      }
    }
  }
}

}  // namespace adj

/** Return a 4-tuple with the number of vertices, edges, faces, corners.  */
int4 MeshPattern::num_elements() const
{
  int4 ans;
  switch (this->kind) {
    case MeshKind::None:
      ans = int4(0, 0, 0, 0);
      break;
    case MeshKind::Line:
      ans = int4(num_anchors + num_segs - 1, num_segs, 0, 0);
      break;
    case MeshKind::TerminalPoly:
      ans = int4(
          num_anchors + num_segs - 2, num_anchors + num_segs - 1, 1, num_anchors + num_segs - 1);
      break;
    case MeshKind::Adj: {
      int totf = adj::f_total_faces(num_anchors, num_segs);
      ans = int4(adj::v_total_verts(num_anchors, num_segs),
                 adj::e_num_edges(num_anchors, num_segs),
                 totf,
                 4 * totf + (adj::odd(num_segs) ? num_anchors - 4 : 0));
      break;
    }
    case MeshKind::TriFan:
      ans = int4(num_segs + 2, num_segs + 2, num_segs, 3 * num_segs);
      break;
    case MeshKind::Cutoff:
      /* TODO */
      BLI_assert(false);
      ans = int4(0, 0, 0, 0);
  }
  return ans;
}

/** Return the index of the vertex in this MeshPattern for the given anchor vertex. */
int MeshPattern::anchor_vert(int anchor_index) const
{
  int ans;
  BLI_assert(0 <= anchor_index && anchor_index < this->num_anchors);
  switch (this->kind) {
    case MeshKind::None:
      ans = 0;
      break;
    case MeshKind::Line:
      ans = anchor_index == 0 ? 0 : this->num_segs;
      break;
    case MeshKind::TerminalPoly:
      ans = anchor_index == 0 ? 0 : this->num_segs + anchor_index;
      break;
    case MeshKind::Adj: {
      const int ring = adj::v_num_rings(this->num_segs) - 1;
      ans = adj::v_ringstart(ring, this->num_anchors, this->num_segs);
      ans += anchor_index * adj::v_anchor_div(ring, this->num_anchors, this->num_segs);
      break;
    }
    case MeshKind::TriFan:
      ans = anchor_index == 0 ? 0 : this->num_segs + anchor_index;
      break;
    case MeshKind::Cutoff:
      /* TODO */
      BLI_assert(false);
      ans = 0;
      break;
  }
  BLI_assert(0 <= ans && ans < this->num_elements()[0]);
  return ans;
}

namespace topology {
/** Functions for analyzing the topology around bevels and setting up main bevel data structures
 * before construction. */

/** if there is common member of both values1 and values2, return the first lexicographically
 * where that is so. Else return -1.
 */
static int find_in_both(Span<int> values1, Span<int> values2)
{
  for (const int v1 : values1) {
    if (values2.contains(v1)) {
      return v1;
    }
  }
  return -1;
}

/** If \a edge is in the mesh face \a face, return the position (amount the corners of the face)
 * where it is. Else return -1. */
static int find_edge_pos_in_face(int edge, int face, const Mesh &mesh)
{
  const Span<int> face_edges = mesh.corner_edges().slice(mesh.faces()[face]);
  for (int e_index : face_edges.index_range()) {
    if (edge == face_edges[e_index]) {
      return e_index;
    }
  }
  return -1;
}

/** Assuming edge1 and edge2 are both in face, are they ccw around a common vertex in that face?
 */
static int edges_ccw_in_face(const int edge1, const int edge2, const int face, const Mesh &mesh)
{
  int pos1 = find_edge_pos_in_face(edge1, face, mesh);
  int pos2 = find_edge_pos_in_face(edge2, face, mesh);
  BLI_assert(pos1 != -1 && pos2 != -1);
  return ((pos1 + 1) % mesh.faces()[face].size()) == pos2;
}

/** Are edges counter clockwise around vert when seen from normal side
 * of the majority of faces involved? */
static bool edges_are_ccw(Span<int> edges, const BevelState &bs)
{
  if (edges.size() <= 1) {
    return true;
  }
  int ccw_test_sum = 0;
  for (const int e_index : edges.index_range()) {
    int edge = edges[e_index];
    int edge_next = edges[(e_index + 1) % edges.size()];
    int face = find_in_both(bs.mesh_info.edge_faces()[edge], bs.mesh_info.edge_faces()[edge_next]);
    if (face != -1) {
      bool ccw_in_face = edges_ccw_in_face(edge, edge_next, face, bs.mesh_info.mesh);
      /* Note that if edges are ccw in the face then that means the are
       * going cw around the bevel vertex.
       */
      ccw_test_sum += ccw_in_face ? -1 : 1;
    }
  }
  return ccw_test_sum >= 0;
}

static MultiValueMap<int, int> build_edge_to_edge_map(Span<int> edges, const BevelState &bs)
{
  GroupedSpan<int> edge_faces = bs.mesh_info.edge_faces();
  MultiValueMap<int, int> face_to_edge;
  for (const int e : edges) {
    for (const int f : edge_faces[e]) {
      face_to_edge.add(f, e);
    }
  }
  MultiValueMap<int, int> map;
  for (const int e : edges) {
    for (const int f : edge_faces[e]) {
      map.add_multiple(e, face_to_edge.lookup(f));
    }
  }
  return map;
}

/** Order the \a edges as specified in #bevvert_order_edges in the general
 * case where not all edges are manifold.
 * Most of the time, try_all_manifold_order will succeed. This is a fallback.
 * Whatever permutation is applied to edges, do the same thing to weights.
 */
static void general_edge_order(MutableSpan<int> edges, const BevelState &bs)
{
  MultiValueMap<int, int> edge_to_edge = build_edge_to_edge_map(edges, bs);
  int n = edges.size();
  Map<int, int> edge_to_index;
  /* An edge should appear at most once in edges. Make an edge -> position map. */
  edge_to_index.reserve(n);
  for (const int i : edges.index_range()) {
    edge_to_index.add_new(edges[i], i);
  }
  const Mesh &mesh = bs.mesh_info.mesh;
  GroupedSpan<int> edge_faces = bs.mesh_info.edge_faces();
  Array<bool, 20> placed(n, false);
  Vector<int, 20> chain_cw(n);
  Vector<int, 20> chain_ccw(n);
  Vector<int, 20> ordered_edges;
  ordered_edges.reserve(n);
  while (ordered_edges.size() != n) {
    /* Loop invariant: edges in ordered_edges are placed.
     * The loop always adds at least first_e to ordered_edges, so will not loop forever.
     * Start by finding the first index #i of an edge that has not yet been placed.
     */
    int first_i = std::find_if(placed.begin(), placed.end(), [](bool b) { return !b; }) -
                  placed.begin();
    BLI_assert(first_i < n);
    int first_e = edges[first_i];
    placed[first_i] = true;
    chain_cw.resize(0);
    chain_ccw.resize(0);
    /* Grow chain in two directions from first_e: cw when dir == 0, ccw when dir == 1. */
    for (const int dir : {0, 1}) {
      /* Follow chain from first_e going in dir as far as possible. */
      int cur_e = first_e;
      while (cur_e != -1) {
        Span<int> neighbors = edge_to_edge.lookup(cur_e);
        int next_e = -1;
        for (const int e : neighbors) {
          if (placed[edge_to_index.lookup(e)]) {
            continue;
          }
          int common_face = find_in_both(edge_faces[cur_e], edge_faces[e]);
          if (common_face != -1) {
            /* Note: ccw in face means cw around common vertex. */
            bool correct_dir = dir == 0 ? edges_ccw_in_face(cur_e, e, common_face, mesh) :
                                          !edges_ccw_in_face(cur_e, e, common_face, mesh);
            if (correct_dir) {
              next_e = e;
              break;
            }
          }
        }
        if (next_e != -1) {
          if (dir == 0) {
            chain_cw.append(next_e);
          }
          else {
            chain_ccw.append(next_e);
          }
          placed[edge_to_index.lookup(next_e)] = true;
        }
        cur_e = next_e;
      }
    }
    /* Append to ordered_edges: rev(chain_cw), first_e, chain_ccw. */
    for (int i = chain_cw.size() - 1; i >= 0; i--) {
      ordered_edges.append(chain_cw[i]);
    }
    ordered_edges.append(first_e);
    ordered_edges.extend(chain_ccw);
  }
  BLI_assert(ordered_edges.size() == n);
  std::copy(ordered_edges.begin(), ordered_edges.end(), edges.begin());
}

/** See if the "all manifold" strategy for ordering edges applies, and if it does
 * order the edges CCW around the vertex.
 */
static bool try_all_manifold_order(MutableSpan<int> edges, const BevelState &bs)
{
  bool edges_all_manifold = std::all_of(edges.begin(), edges.end(), [&](const int e) {
    return bs.mesh_info.edge_faces()[e].size() == 2;
  });
  if (!edges_all_manifold) {
    return false;
  }
  int n = edges.size();
  bool found_order;
  if (n <= 3) {
    /* There is only one order, modulo possible reversal. */
    found_order = true;
  }
  else if (n <= 6) {
    GroupedSpan<int> edge_faces = bs.mesh_info.edge_faces();
    /* Use O(n^2) method that needs no extra map. */
    for (int i = 1; i < n - 1; i++) {
      /* Loop invariant: edges[0..i) form a chain through faces and are placed.
       * We only loop to i == n-2 because at that point, whether the last edge
       * forms part of the chain or not, there isn't a better order.
       */
      int prev_e = edges[i - 1];
      BLI_assert(edge_faces[prev_e].size() == 2);
      int f1 = edge_faces[prev_e][0];
      int f2 = edge_faces[prev_e][1];
      /* Find an edge another edge among the unplaced edges, adjacent to f1 or f2.
       * If find one, swap the edge in position i with the found one.
       */
      int next_e_index = -1;
      for (int j = i; j < n; j++) {
        int e = edges[j];
        BLI_assert(edge_faces[e].size() == 2);
        int g1 = edge_faces[e][0];
        int g2 = edge_faces[e][1];
        if (f1 == g1 || f1 == g2 || f2 == g1 || f2 == g2) {
          next_e_index = j;
          break;
        }
      }
      if (next_e_index == -1) {
        found_order = false;
        break;
      }
      if (i != next_e_index) {
        std::swap(edges[i], edges[next_e_index]);
      }
    }
    found_order = true;
  }
  else {
    /* For larger n, speed up the search by making a map from
     * edges to other edges that share a face.
     */
    MultiValueMap<int, int> edge_to_edge = build_edge_to_edge_map(edges, bs);
    /* Similar loop to above, but using edge_to_edge map. */
    for (int i = 1; i < n - 1; i++) {
      int prev_e = edges[i - 1];
      int next_e_index = -1;
      for (int j = i; j < n; j++) {
        int e = edges[j];
        if (edge_to_edge.lookup(prev_e).contains(e)) {
          next_e_index = j;
          break;
        }
      }
      if (next_e_index == -1) {
        found_order = false;
        break;
      }
      if (i != next_e_index) {
        std::swap(edges[i], edges[next_e_index]);
      }
    }
    found_order = true;
  }
  if (found_order) {
    /* Possibly reverse the order. */
    if (!edges_are_ccw(edges, bs)) {
      std::reverse(edges.begin(), edges.end());
    }
  }
  return found_order;
}

/** Order the \a bevedges so that, as much as possible,
 * the edges are in a sequence such that they faces between them.
 * Also, so that with viewed from the positive face normal side of those faces,
 * we generally move counterclockwise around the bevel vertex.
 * The \a faces argument will be filled in with indices of shared faces
 * between the corresponding bevedge and its successor (cyclicly) if there is one,
 * else -1.
 */
static void bevvert_order_edges(BevelState &bs, MutableSpan<int> bevedges, MutableSpan<int> faces)
{
  Array<int, 20> edges(bevedges.size());
  for (const int i : bevedges.index_range()) {
    edges[i] = bs.bevedge_mesh_edges()[bevedges[i]];
  }
  if (!try_all_manifold_order(edges.as_mutable_span(), bs)) {
    general_edge_order(edges.as_mutable_span(), bs);
  }
  /* If we are edge beveling, we need to rotate so that a beveled edge is first. */
  if (bs.params.affect_type == BevelAffect::Edges) {
    auto first_beveled_edge_iter = std::find_if(edges.begin(), edges.end(), [&](int e) {
      int be = bs.edge_bevedges()[e];
      return bs.bevedge_weights()[be] > 0.0f;
    });
    BLI_assert(first_beveled_edge_iter != edges.end());
    std::rotate(edges.begin(), first_beveled_edge_iter, edges.end());
  }
  /* Convert the ordered edges to bevedges and copy the new order to bevedges argument.
   * Also set the shared face between bevedges, if any.
   */
  GroupedSpan<int> edge_faces = bs.mesh_info.edge_faces();
  for (const int i : bevedges.index_range()) {
    const int e = edges[i];
    const int next_e = edges[(i == bevedges.size() - 1) ? 0 : i + 1];
    bevedges[i] = bs.edge_bevedges()[e];
    faces[i] = find_in_both(edge_faces[e], edge_faces[next_e]);
  }
}

/** How many edges attached to BevVert \a bv are beveled? */
static int bevvert_num_beveled_edges(int bv, const BevelState &bs)
{
  Span<int> bevedges = bs.bevvert_bevedges()[bv];
  return std::accumulate(bevedges.begin(), bevedges.end(), 0, [&](int sum, int be) {
    return sum + (bs.bevedge_is_beveled(be) > 0.0f ? 1 : 0);
  });
}

/** What MeshPattern will be used for BevVert \a bv ? */
static MeshPattern find_meshpattern(int bv, BevelState &bs)
{
  MeshPattern pat;
  pat.num_segs = bs.params.segments;
  int num_edges = bs.bevvert_bevedges()[bv].size();
  BLI_assert(num_edges > 1 && pat.num_segs >= 1);
  if (bs.params.affect_type == BevelAffect::Vertices) {
    if (num_edges < 3) {
      pat.kind = MeshKind::Line;
      /* Anchors are on each attached edge.
       * They will be connected by a line of pat.num_segs segments. */
      pat.num_anchors = 2;
    }
    else {
      /* Anchors are on each attached edge. */
      pat.kind = MeshKind::Adj;
      pat.num_anchors = num_edges;
    }
  }
  else {
    int num_beveled = bevvert_num_beveled_edges(bv, bs);
    BLI_assert(num_beveled > 0 && num_beveled <= num_edges);
    if (num_beveled == 1) {
      /* Terminal edge case. */
      if (num_edges == 2) {
        /* Need to add an artificial vertex on the unbeveled edge. */
        pat.kind = MeshKind::TriFan;
        pat.num_anchors = 3;
      }
      else if (num_edges == 3) {
        pat.kind = MeshKind::Line;
        pat.num_anchors = 2;
      }
      else if (num_edges == 4) {
        pat.kind = MeshKind::TriFan;
        pat.num_anchors = 3;
      }
      else {
        pat.kind = MeshKind::TerminalPoly;
        pat.num_anchors = num_edges - 1;
      }
    }
    else {
      /* Edge bevel with at least two beveled edges. */
      if (num_beveled == 2) {
        /* The two beveled edges butt up against each other, along a line. */
        pat.kind = MeshKind::Line;
        pat.num_anchors = 2;
      }
      else {
        /* Edge bevel with three or more beveled edges.
         * There are anchors at the point where the right side of one beveled
         * edge attaches to the left side of the next beveled edge. */
        pat.kind = MeshKind::Adj;
        pat.num_anchors = num_beveled;
        /* However if there is mitering, we put additional anchor points
         * between the beveled edges: one additional for "arc" miters and two
         * additional for "patch" miters. We can only have pacth miters on
         * outer miters (where there is an obtuse angle). */
        if (bs.params.miter_inner != BevelMiterType::Sharp) {
          pat.num_anchors += num_beveled;
        }
        if (bs.params.miter_outer == BevelMiterType::Patch) {
          /* TODO: this is wrong if there is no obtuse angle between two
           * successive beveled edges! */
          pat.num_anchors += 1;
        }
      }
    }
  }
  return pat;
}

/* Helper function to find the coordinates between two beveled edges in the adj pattern. */
static float3 adj_edge_anchor_co(int bv,
                                 int cur_anchor_pos,
                                 int next_anchor_pos,
                                 int num_in_plane,
                                 int be_in_plane_pos,
                                 int num_not_in_plane,
                                 int be_not_in_plane_pos,
                                 const BevelState &bs)
{
  const int be_cur = bs.bevvert_bevedges()[bv][cur_anchor_pos];
  const int be_next = bs.bevvert_bevedges()[bv][next_anchor_pos];
  const int be_cur_end = bs.bevedge_vert_end(bv, be_cur);
  const int be_next_end = bs.bevedge_vert_end(bv, be_next);
  BLI_assert(bs.bevedge_is_beveled(be_cur) && bs.bevedge_is_beveled(be_next));
  bool offset_edge_between = false;
  bool offset_meet_edges_between = false;
  float3 co;
  float sin_ratio = 1.0f;
  int eon_pos = -1;
  const float spec_r = bs.bevedge_widths()[be_cur][2 * be_cur_end + 1];
  const float spec_next_l = bs.bevedge_widths()[be_next][2 * be_next_end];
  int offset_bepos_in_plane = -1;
  int offset_meet_face = bs.face_next(bv, cur_anchor_pos);
  if (num_not_in_plane > 0) {
    if (bs.params.loop_slide && num_not_in_plane == 1 &&
        geom::good_slide(
            bv, cur_anchor_pos, next_anchor_pos, be_not_in_plane_pos, spec_r, spec_next_l, bs))
    {
      offset_edge_between = true;
      offset_bepos_in_plane = be_not_in_plane_pos;
      eon_pos = be_not_in_plane_pos;
    }
    else {
      offset_meet_face = -1;
      offset_meet_edges_between = true;
    }
  }
  else if (num_in_plane > 0) {
    if (bs.params.loop_slide && num_in_plane == 1 &&
        geom::good_slide(
            bv, cur_anchor_pos, next_anchor_pos, be_in_plane_pos, spec_r, spec_next_l, bs))
    {
      offset_edge_between = true;
      eon_pos = be_in_plane_pos;
    }
    else {
      offset_meet_edges_between = false;
    }
  }
  if (offset_edge_between) {
    if (!geom::try_offset_on_edge_between(bv,
                                          cur_anchor_pos,
                                          next_anchor_pos,
                                          eon_pos,
                                          spec_r,
                                          spec_next_l,
                                          &co,
                                          &sin_ratio,
                                          bs))
    {
      eon_pos = -1;
      offset_edge_between = false;
    }
  }
  if (!offset_edge_between) {
    co = geom::offset_meet(bv,
                           cur_anchor_pos,
                           next_anchor_pos,
                           spec_r,
                           spec_next_l,
                           offset_meet_face,
                           offset_meet_edges_between,
                           offset_bepos_in_plane,
                           bs);
  }
  /* TODO: stuff with eon and mitering. */
  return co;
}

/** Set the positions of the anchor vertices for bv, and the bevedge attachment points.  */
static void build_vmesh_skeleton(int bv,
                                 MutableSpan<float3> bv_newvert_positions,
                                 MutableSpan<int2> be_attach_verts,
                                 const BevelState &bs)
{
  fmt::println("build_vmesh_skelton bv={}", bv);
  bv_newvert_positions.fill(float3(0.0f, 0.0f, 0.0f));
  Span<int> bevedges = bs.bevvert_bevedges()[bv];
  const int num_edges = bevedges.size();
  BLI_assert(num_edges > 1);
  MeshPattern pat = bs.bevvert_meshpatterns()[bv];
  if (bs.params.affect_type == BevelAffect::Vertices) {
    /* Vertex bevel. */
    for (const int epos : bevedges.index_range()) {
      int be = bevedges[epos];
      int be_end = bs.bevedge_vert_end(bv, be);
      float spec = bs.bevedge_widths()[be][be_end * 2];
      float3 co = geom::slide_dist(bv, epos, spec, bs);
      int nv_index = pat.anchor_vert(epos);
      bv_newvert_positions[nv_index] = co;
      be_attach_verts[be][be_end] = nv_index;
    }
  }
  else {
    /* Edge bevel. */;
    /* First gather the positions of all the actual beveled edges and their width specs. */
    Vector<int, 20> bevel_pos;
    Vector<float2, 20> widths;
    for (const int epos : bevedges.index_range()) {
      const int be = bevedges[epos];
      if (bs.bevedge_is_beveled(be)) {
        bevel_pos.append(epos);
        const float4 be_bev_widths = bs.bevedge_widths()[be];
        const int end = bs.bevedge_vert_end(bv, be);
        widths.append(float2(be_bev_widths[2 * end], be_bev_widths[2 * end + 1]));
      }
      else {
        widths.append(float2(0.0f, 0.0f));
      }
    }
    int num_beveled = bevel_pos.size();
    BLI_assert(num_beveled > 0 && bevel_pos[0] == 0);
    const int be0 = bevedges[0];
    const int be0_end = bs.bevedge_vert_end(bv, be0);
    const int be1 = bevedges[1];
    const int be1_end = bs.bevedge_vert_end(bv, be1);
    const int anchor0 = pat.anchor_vert(0);
    const int anchor1 = pat.anchor_vert(1);
    /* All cases have be0 attached to anchor 0. */
    be_attach_verts[be0][be0_end] = anchor0;
    bool leg_side = ELEM(
        bs.params.offset_type, BevelOffsetType::Percent, BevelOffsetType::Absolute);
    if (num_beveled == 1) {
      /* Terminal edge cases. */
      if (num_edges == 2) {
        BLI_assert(pat.kind == MeshKind::TriFan && pat.num_anchors == 3);
        be_attach_verts[be1][be1_end] = pat.anchor_vert(2);
        bv_newvert_positions[anchor0] = geom::offset_bevedge(bv, 0, widths[0][0], true, bs);
        bv_newvert_positions[anchor1] = geom::offset_bevedge(bv, 0, widths[0][1], false, bs);
        bv_newvert_positions[pat.anchor_vert(2)] = geom::slide_dist(bv, 1, widths[0][0], bs);
      }
      else if (num_edges == 3) {
        BLI_assert(pat.kind == MeshKind::Line && pat.num_anchors == 2);
        be_attach_verts[be1][be1_end] = anchor1;
        const int be2 = bevedges[2];
        be_attach_verts[be2][bs.bevedge_vert_end(bv, be2)] = anchor0;
        if (leg_side) {
          bv_newvert_positions[anchor0] = geom::slide_dist(bv, 2, widths[0][0], bs);
          bv_newvert_positions[anchor1] = geom::slide_dist(bv, 1, widths[0][1], bs);
        }
        else {
          bv_newvert_positions[anchor0] = geom::offset_meet(
              bv, 2, 0, 0.0f, widths[0][0], bs.face_prev(bv, 0), false, -1, bs);
          bv_newvert_positions[anchor1] = geom::offset_meet(
              bv, 0, 1, widths[0][1], 0.0f, bs.face_next(bv, 0), false, -1, bs);
        }
      }
      else if (num_edges == 4) {
        BLI_assert(pat.kind == MeshKind::TriFan && pat.num_anchors == 3);
        be_attach_verts[be1][be1_end] = anchor1;
        const int be2 = bevedges[2];
        be_attach_verts[be2][bs.bevedge_vert_end(bv, be2)] = pat.anchor_vert(2);
        const int be3 = bevedges[3];
        be_attach_verts[be3][bs.bevedge_vert_end(bv, be3)] = anchor0;
        if (leg_side) {
          bv_newvert_positions[anchor0] = geom::slide_dist(bv, 3, widths[0][0], bs);
          bv_newvert_positions[anchor1] = geom::slide_dist(bv, 1, widths[0][1], bs);
        }
        else {
          bv_newvert_positions[anchor0] = geom::offset_meet(
              bv, 3, 0, 0.0f, widths[0][0], bs.face_prev(bv, 0), false, -1, bs);
          bv_newvert_positions[anchor1] = geom::offset_meet(
              bv, 0, 1, widths[0][1], 0.0f, bs.face_next(bv, 0), false, -1, bs);
        }
      }
      else {
        BLI_assert(pat.kind == MeshKind::TerminalPoly && pat.num_anchors == bevedges.size() - 1);
        for (int i = 1; i < num_edges; i++) {
          const int be = bevedges[i];
          be_attach_verts[be][bs.bevedge_vert_end(bv, be)] = pat.anchor_vert(i);
        }
        if (leg_side) {
          bv_newvert_positions[anchor0] = geom::slide_dist(bv, num_edges - 1, widths[0][0], bs);
          bv_newvert_positions[anchor1] = geom::slide_dist(bv, 1, widths[0][1], bs);
        }
        else {
          bv_newvert_positions[anchor0] = geom::offset_meet(
              bv, num_edges - 1, 0, 0.0f, widths[0][0], bs.face_prev(bv, 0), false, -1, bs);
          bv_newvert_positions[anchor1] = geom::offset_meet(
              bv, 0, 1, widths[0][1], 0.0f, bs.face_next(bv, 0), false, -1, bs);
        }
        for (int i = 2; i < num_edges - 1; i++) {
          bv_newvert_positions[anchor1 + i - 1] = geom::slide_dist(bv, i, widths[0][0], bs);
        }
      }
    }
    else if (num_beveled == 2) {
      BLI_assert(pat.kind == MeshKind::Line && pat.num_anchors == 2);
      bool seen_second_bevel = false;
      for (int i = 1; i < num_edges; i++) {
        seen_second_bevel = seen_second_bevel || i == bevel_pos[1];
        const int be = bevedges[i];
        be_attach_verts[be][bs.bevedge_vert_end(bv, be)] = seen_second_bevel ? anchor0 : anchor1;
        bv_newvert_positions[anchor0] = geom::offset_meet(
            bv, 1, 0, widths[1][1], widths[0][0], bs.face_prev(bv, 0), false, -1, bs);
        bv_newvert_positions[anchor1] = geom::offset_meet(
            bv, 0, 1, widths[0][1], widths[1][0], bs.face_next(bv, 0), false, -1, bs);
      }
    }
    else {
      BLI_assert(pat.kind == MeshKind::Adj && pat.num_anchors >= 3);
      /* TODO: miters. */
      BLI_assert(bs.params.miter_inner == BevelMiterType::Sharp &&
                 bs.params.miter_outer == BevelMiterType::Sharp);
      int cur_anchor = 0;
      int next_anchor = 1;
      /* We need to know the number of edges between the current beveled edge and the
       * next one that are "in plane" (coplanar with their two adjactent faces) or
       * "not in plane". */
      int num_in_plane = 0;
      int num_not_in_plane = 0;
      int be_in_plane_pos = -1;
      int be_not_in_plane_pos = -1;
      for (int i = 1; i < num_edges; i++) {
        const int be = bevedges[i];
        const int be_end = bs.bevedge_vert_end(bv, be);
        be_attach_verts[be][be_end] = pat.anchor_vert(next_anchor);
        if (i == bevel_pos[next_anchor]) {
          /* We have reached the next beveled edge. */
          const int anchor_pos = pat.anchor_vert(next_anchor);
          fmt::println(
              "cur_anchor={} next_anchor={} anchor_pos={}", cur_anchor, next_anchor, anchor_pos);
          bv_newvert_positions[anchor_pos] = adj_edge_anchor_co(bv,
                                                                bevel_pos[cur_anchor],
                                                                bevel_pos[next_anchor],
                                                                num_in_plane,
                                                                be_in_plane_pos,
                                                                num_not_in_plane,
                                                                be_not_in_plane_pos,
                                                                bs);
          print_float3(bv_newvert_positions[anchor_pos]);
          fmt::println("");
          cur_anchor = next_anchor;
          next_anchor = (cur_anchor + 1) % num_beveled;
          num_in_plane = 0;
          num_not_in_plane = 0;
          be_in_plane_pos = -1;
          be_not_in_plane_pos = -1;
        }
        else {
          if (geom::bevedge_on_plane(bv, i, bs)) {
            num_in_plane++;
            be_in_plane_pos = i;
          }
          else {
            num_not_in_plane++;
            be_not_in_plane_pos = i;
          }
        }
      }
      const int anchor_pos = pat.anchor_vert(next_anchor);
      fmt::println(
          "cur_anchor={} next_anchor={} anchor_pos={}", cur_anchor, next_anchor, anchor_pos);
      bv_newvert_positions[anchor_pos] = adj_edge_anchor_co(bv,
                                                            bevel_pos[cur_anchor],
                                                            bevel_pos[next_anchor],
                                                            num_in_plane,
                                                            be_in_plane_pos,
                                                            num_not_in_plane,
                                                            be_not_in_plane_pos,
                                                            bs);
      print_float3(bv_newvert_positions[anchor_pos]);
      fmt::println("");
    }
  }
}

static void build_internal_adj(const int bv,
                               const profile::AnchorProfiles &profiles,
                               MutableSpan<float3> bv_newvert_positions,
                               const BevelState &bs)
{
  const MeshPattern &pat = bs.bevvert_meshpatterns()[bv];
  const int ns = bs.params.segments;
  const int na = pat.num_anchors;
  BLI_assert(ns > 1 && pat.kind == MeshKind::Adj);
  const int ns_power_2 = bs.pro_spacing.segments_power_2;

  /* First construct an initial control mesh with 2 segments. */
  adj::AdjVerts adj2(na, 2);
  float3 center(0.0f, 0.0f, 0.0f);
  for (const int a : IndexRange(na)) {
    float3 pos = bv_newvert_positions[pat.anchor_vert(a)];
    adj2.mutable_outer_ring_vert(a, 0) = bv_newvert_positions[pat.anchor_vert(a)];
    pos = profile::get_profile_point(profiles[a], 1, 2);
    adj2.mutable_outer_ring_vert(a, 1) = profile::get_profile_point(profiles[a], 1, 2);
    center = center + adj2.vert(0, a, 0);
  }
  center = center / float(na);

  /* To place the center vertex, let:
   * 'negative_fullest' = the original vertex across the boundverts' center.
   * 'fullness' = fraction of the way from the boundvert's centroid
   *  to the original vertex (if positive) or to negative_fullest (if negative).
   */
  const float3 orig_v = bs.mesh_info.mesh.vert_positions()[bs.bevvert_mesh_verts()[bv]];
  const float fullness = bs.pro_spacing.fullness;
  const float3 center_dir = orig_v - center;
  if (math::length_squared(center_dir) > geom::bevel_epsilon_sq) {
    adj2.mutable_vert(0, 0, 0) = center + fullness * center_dir;
    /* TODO: custom profile is different here. */
    /* const float3 negative_fullest = center + (center - orig_v); */
  }
  else {
    adj2.mutable_vert(0, 0, 0) = center;
  }

  /* Make and fill adj mesh with #ns_power_2 segements. */
  adj::AdjVerts adj_sup_power_2(na, ns_power_2);
  adj::fill_adjverts(adj_sup_power_2, adj2, profiles);

  /* Interpolate the mesh to the needed number of eegments, if necessary. */
  adj::AdjVerts adj(na, ns);
  if (ns == ns_power_2) {
    BLI_assert(adj_sup_power_2.verts.size() == adj.verts.size());
    std::copy(adj_sup_power_2.verts.begin(), adj_sup_power_2.verts.end(), adj.verts.begin());
  }
  else {
    adj::interp_adj(adj, adj_sup_power_2, profiles);
  }

  /* Snap the vertices of adj to the superellipsoid. */


  BLI_assert(adj.verts.size() == bv_newvert_positions.size());
  std::copy(adj.verts.begin(), adj.verts.end(), bv_newvert_positions.begin());
}

static void build_internal_vmesh(const int bv,
                                 const profile::AnchorProfiles &profiles,
                                 MutableSpan<float3> bv_newvert_positions,
                                 const BevelState &bs)
{
  const MeshPattern &pat = bs.bevvert_meshpatterns()[bv];
  const int ns = bs.params.segments;
  if (ns == 1) {
    /* The newverts are all anchor points, so already set. */
    return;
  }
  if (pat.kind == MeshKind::Line) {
    BLI_assert(pat.num_anchors == 2 && pat.num_elements()[0] == ns + 1);
    const profile::Profile pro = profiles[0];
    for (int i = 1; i <= ns; i++) {
      bv_newvert_positions[i] = pro.prof_co[i];
    }
  }
  else if (pat.kind == MeshKind::Adj) {
    build_internal_adj(bv, profiles, bv_newvert_positions, bs);
  }
  else {
    /* TODO: implement me. */
    BLI_assert(false);
  }
}

}  // end of namespace topology

namespace spec {

/** Use the bs.params offset type and amount to set initial widths at the
 * \a bv end of the edge at its position \a edge_pos.
 * Fill in either the first two elements of \a r_widths or the last two,
 * depending on whether bv is first in the edge or second.
 */
static void edge_bevel_init_widths_from_params(const int bv,
                                               const int edge_pos,
                                               float4 &r_widths,
                                               const BevelState &bs)
{
  const int be = bs.bevvert_bevedges()[bv][edge_pos];
  float2 specs;
  if (!bs.bevedge_is_beveled(be)) {
    specs[0] = 0.0f;
    specs[1] = 0.0f;
  }
  else {
    float offset = bs.params.offset;
    switch (bs.params.offset_type) {
      case BevelOffsetType::Offset:
        specs[0] = offset;
        specs[1] = offset;
        break;
      case BevelOffsetType::Width:
      case BevelOffsetType::Depth: {
        float z = bs.params.offset_type == BevelOffsetType::Width ?
                      math::abs(2.0f * math::sin(geom::edge_face_angle(bv, edge_pos, bs) / 2.0f)) :
                      math::abs(math::cos(geom::edge_face_angle(bv, edge_pos, bs) / 2.0f));
        if (z < geom::bevel_epsilon) {
          /* Undefined behavior, so just do a tiny bevel. */
          specs[0] = 0.01f * offset;
        }
        else {
          specs[0] = offset / z;
        }
        specs[1] = specs[0];
        break;
      }
      case BevelOffsetType::Percent: {
        /* Offset needs to meet adjacent edges at percentage of their lengths.
         * Since the width isn't constant, we don't store a width at all, but
         * rather the distance along the adjacent edge that we need to go
         * at this end of the edge.
         */
        const int edge_pos_prev = bs.prev_edge_pos(bv, edge_pos);
        const int edge_pos_next = bs.next_edge_pos(bv, edge_pos);
        specs[0] = geom::edge_length(bv, edge_pos_prev, bs) * offset / 100.0f;
        specs[1] = geom::edge_length(bv, edge_pos_next, bs) * offset / 100.0f;
        break;
      }
      case BevelOffsetType::Absolute:
        /* Spec is the absolute amount along previous and next edges. */
        specs[0] = offset;
        specs[1] = offset;
        break;
      default:
        BLI_assert_unreachable();
    }
    if (bs.params.use_weights) {
      /* TODO: implement me.
       * Before doing this, change the API..
       */
    }
  }
  const bool near_end = bs.bevedge_vert_end(bv, be) == 0;
  if (near_end) {
    r_widths[0] = specs[0];
    r_widths[1] = specs[1];
  }
  else {
    r_widths[2] = specs[0];
    r_widths[3] = specs[1];
  }
}

/** Like the previous, but for vertex beveling. */
static void vertex_bevel_init_widths_from_params(const int bv,
                                                 const int edge_pos,
                                                 float4 &r_widths,
                                                 const BevelState &bs)
{
  const int be = bs.bevvert_bevedges()[bv][edge_pos];
  const bool near_end = bs.bevedge_vert_end(bv, be) == 0;
  float spec;
  float offset = bs.params.offset;
  switch (bs.params.offset_type) {
    case BevelOffsetType::Offset:
      spec = offset;
      break;
    case BevelOffsetType::Width:
    case BevelOffsetType::Depth: {

      float3 edir = geom::edge_dir(bv, edge_pos, bs);
      float3 vert_axis = geom::bevvert_axis(bv, bs);
      if (near_end) {
        edir = -edir;
      }
      float z = bs.params.offset_type == BevelOffsetType::Width ?
                    math::abs(2.0f * math::sin(angle_v3v3(vert_axis, edir))) :
                    math::abs(math::cos(angle_v3v3(vert_axis, edir)));
      if (z < geom::bevel_epsilon) {
        /* Undefined behavior, so do tiny bevel. */
        spec = 0.01f * offset;
      }
      else {
        spec = offset / z;
      }
      break;
    }
    case BevelOffsetType::Percent: {
      spec = geom::edge_length(bv, edge_pos, bs) * offset / 100.0f;
      break;
    }
    case BevelOffsetType::Absolute:
      spec = offset;
      break;
    default:
      BLI_assert_unreachable();
  }
  float weight = 1.0f;
  if (bs.params.use_weights) {
    /* TODO: change bevel weight API.
     * Current code takes in a weight name and converts that
     * to a CustomData layer index.
     * We should just take in the weight name.
     */
  }
  else if (bs.params.dvert != nullptr) {
    /* TODO: change bevel vertex group API.
     * Current code uses MOD_get_vgroup to get an MDeformVert pointer
     * from a name.
     * We should just take in the vgroup name.
     */
  }
  if (weight != 1.0f) {
    /* TODO: multiple all specs by weight. */
  }

  if (near_end) {
    r_widths[0] = spec;
    r_widths[1] = 0.0f;
  }
  else {
    r_widths[2] = spec;
    r_widths[3] = 0.0f;
  }
}

}  // end of namespace spec

/** Initialize the bevel_edges and bevel_verts IndexMasks, and also some fields
 * that will control the basic operaiton of the whole bevel.
 */
BevelState::BevelState(const Mesh &src_mesh,
                       const IndexMask &selection,
                       const BevelParameters &bevel_params)
    : mesh_info(src_mesh)
{
  /* Copy the parameters so we can safely modify them if necessary. */
  this->params = bevel_params;
  bool vertex_only = params.affect_type == BevelAffect::Vertices;

  /* There is an adjustment procedure that modifies bevel amounts to make the best compromise
   * when all specifications cannot be met simultaneously. That procedure cannot be applied
   * in all circumstances. */
  this->offset_adjust = !(vertex_only || params.offset_type == BevelOffsetType::Percent ||
                          params.offset_type == BevelOffsetType::Absolute);

  /* Disable the miters with the cutoff vertex mesh method, the combination isn't useful anyway.
   */
  if (params.vmesh_method == BevelVmeshMethod::Cutoff) {
    params.miter_outer = BevelMiterType::Sharp;
    params.miter_inner = BevelMiterType::Sharp;
  }

  /* Get the IndexMasks for bevel-involved edges and bevel-involved vertices,
   * and set the weights that modify the bevel amount for those elements
   * actually beveled.
   *
   * The selection parameter selectes what is beveled depending on the affect type.
   * We filter the selection to only select legal elements to bevel:
   * for edge bevels, the edge must be manifold (attached to exactly 2 faces);
   * for vertex bevels, the vertex must have at least two edges attached to it.
   *
   * Note: the legacy bmesh bevel code assumed that the caller made these things true.
   */
  IndexMask beveled_edges_mask;
  if (vertex_only) {
    bevverts_mask_ = index_mask::IndexMask::from_predicate(
        selection, GrainSize(4096), this->memory_, [&](const int v) {
          return this->mesh_info.vert_edges()[v].size() >= 2;
        });
  }
  else {
    beveled_edges_mask = index_mask::IndexMask::from_predicate(
        selection, GrainSize(4096), this->memory_, [&](const int e) {
          return this->mesh_info.edge_faces()[e].size() == 2;
        });
    Array<bool> bevel_involved_vert(src_mesh.verts_num, false);
    beveled_edges_mask.foreach_index([&](const int e) {
      int2 edge_vs = src_mesh.edges()[e];
      bevel_involved_vert[edge_vs[0]] = true;
      bevel_involved_vert[edge_vs[1]] = true;
    });
    bevverts_mask_ = IndexMask::from_bools(bevel_involved_vert, memory_);
  }
  Array<bool> bevel_involved_edge(src_mesh.edges_num, false);
  bevverts_mask_.foreach_index([&](const int v) {
    Span<int> edges = this->mesh_info.vert_edges()[v];
    for (const int e : edges) {
      bevel_involved_edge[e] = true;
    }
  });
  bevedges_mask_ = IndexMask::from_bools(bevel_involved_edge, memory_);

  /* Create the arrays for bevverts and maps between verts and bevverts. */
  this->bevverts_num = bevverts_mask_.size();
  bevvert_mesh_verts_ = Array<int>(this->bevverts_num);
  vert_bevverts_ = Array<int>(this->mesh_info.mesh.verts_num, -1);
  if (vertex_only) {
    bevvert_weights_ = Array<float>(this->bevverts_num, 1.0f);
  }
  bevverts_mask_.foreach_index([&](const int v, const int mask) {
    vert_bevverts_[v] = mask;
    bevvert_mesh_verts_[mask] = v;
    if (vertex_only) {
      bevvert_weights_[mask] = selection.contains(v);
    }
  });

  /* Create the arrays for the bevedges and maps between edges and bevedges. */
  this->bevedges_num = bevedges_mask_.size();
  bevedge_mesh_edges_ = Array<int>(this->bevedges_num);
  edge_bevedges_ = Array<int>(this->mesh_info.mesh.edges_num, -1);
  if (!vertex_only) {
    bevedge_weights_ = Array<float>(this->bevedges_num);
    bevedge_widths_ = Array<float4>(this->bevedges_num);
  }
  bevedges_mask_.foreach_index([&](const int e, const int mask) {
    edge_bevedges_[e] = mask;
    bevedge_mesh_edges_[mask] = e;
    if (vertex_only) {
      bevedge_weights_[mask] = 0.0f;
    }
    else {
      bevedge_weights_[mask] = selection.contains(e);
    }
  });

  /* Find the bevedges attached to each bevvert. */

  /* First fill in the counts of the number of attached edges, then convert to offsets. */
  bevvert_bevedges_offsets_ = Array<int>(this->bevverts_num + 1);
  threading::parallel_for(IndexRange(this->bevverts_num), 20000, [&](IndexRange range) {
    for (const int bv : range) {
      int v = bevvert_mesh_verts_[bv];
      bevvert_bevedges_offsets_[bv] = this->mesh_info.vert_edges()[v].size();
    }
  });
  offset_indices::accumulate_counts_to_offsets(bevvert_bevedges_offsets_);
  OffsetIndices<int> offsets(bevvert_bevedges_offsets_);
  bevvert_bevedges_indices_ = Array<int>(offsets.total_size());
  threading::parallel_for(IndexRange(this->bevverts_num), 20000, [&](IndexRange range) {
    for (const int bv : range) {
      const int v = bevvert_mesh_verts_[bv];
      Span<int> v_edges = this->mesh_info.vert_edges()[v];
      int offset_first = offsets[bv].first();
      for (const int i : IndexRange(offsets[bv].size())) {
        bevvert_bevedges_indices_[offset_first + i] = edge_bevedges_[v_edges[i]];
      }
    }
  });
  bevvert_bevedges_ = GroupedSpan<int>(offsets, bevvert_bevedges_indices_);
}

/** Return the bevedge poisition of the last bevedge attached to newvert, or -1 if none. */
int BevelState::last_attached_bevedge_pos(const int bv, const int newvert) const
{
  Span<int> edges = bevvert_bevedges_[bv];
  int ans = -1;
  for (const int epos : edges.index_range()) {
    const int be = edges[epos];
    const int be_end = bevedge_vert_end(bv, be);
    const int attach_vert = bevedge_attach_verts_[be][be_end];
    if (attach_vert == newvert) {
      ans = epos;
    }
  }
  return ans;
}

/** Initialize the part of the state related to the profile curves that will be used in the
 * non-custom shapes of multisegment bevels.
 */
void BevelState::initialize_profile_data()
{
  /* Convert the input profile shape parameter to the actual exponent of a superellipse. */
  const float psr = -std::log(2.0) /
                    std::log(math::sqrt(this->params.profile > 0 ? this->params.profile : 1e-20f));
  this->pro_super_r = psr;

  /* Snap some ranges of profile to particular shapes. */
  if (this->params.profile >= 0.950f) { /* r ~ 692, so pro_square_r is 1e4 */
    this->pro_super_r = profile::pro_square_r;
  }
  else if (abs(psr - profile::pro_circle_r) < 1e-4f) {
    this->pro_super_r = profile::pro_circle_r;
  }
  else if (abs(psr - profile::pro_line_r) < 1e-4f) {
    this->pro_super_r = profile::pro_line_r;
  }
  else if (abs(psr) < 1e-4f) {
    this->pro_super_r = profile::pro_square_in_r;
  }

  profile::set_profile_spacing(
      this, &this->pro_spacing, this->params.profile_type == BevelProfileType::Custom);

  if (this->params.segments > 1) {
    this->pro_spacing.fullness = profile::find_profile_fullness(this);
  }

  /* Get separate non-custom profile samples for the miter profiles if they are needed */
  /* TODO: check: bmesh code seems wrong here: the first check is ==, not != */
  if (this->params.profile_type != BevelProfileType::Custom &&
      (this->params.miter_inner != BevelMiterType::Sharp ||
       this->params.miter_outer != BevelMiterType::Sharp))
  {
    profile::set_profile_spacing(this, &this->pro_spacing_miter, false);
  }
}

/** Order the bevedges around each bevvert so that they form manifold caps
 * or pieces of manifold caps, as much as possible.
 * In the simple usual case, there is one manifold cap (a cyclic ordering of edges
 * where there is exactly one face between each successive pair of edges).
 * In that simple case, the ordering is counterclockwise when looking at it from
 * the positive normal side of the vertex and faces.
 *
 * As well as ordering the bevedges, set up the bevvert_bevfaces GroupedSpan,
 * which gives the faces shared between a bevvedge in bevvert_bevedges and its
 * immediately following (in cyclic order) bevedge.
 */
void BevelState::order_bevedges()
{
  /* We'll shared the bevvert_bevedges_offsets_ for bevvert_faces_, since the faces
   * are parallel to the bevedges. */
  OffsetIndices<int> offsets(bevvert_bevedges_offsets_);
  bevvert_faces_indices_ = Array<int>(offsets.total_size());
  threading::parallel_for(IndexRange(this->bevverts_num), 4096, [&](IndexRange range) {
    for (const int bv : range) {
      IndexRange r = OffsetIndices<int>(bevvert_bevedges_offsets_)[bv];
      MutableSpan<int> edges = bevvert_bevedges_indices_.as_mutable_span().slice(r);
      MutableSpan<int> faces = bevvert_faces_indices_.as_mutable_span().slice(r);
      topology::bevvert_order_edges(*this, edges, faces);
    }
  });
  bevvert_faces_ = GroupedSpan<int>(offsets, bevvert_faces_indices_);
}

/** Figure out the topology of the vertex meshes at each bevvert.
 * Allocate the containers for newverts, newedges, and newfaces.
 */
void BevelState::set_bevvert_mesh_topology()
{
  bevvert_meshpatterns_ = Array<MeshPattern>(this->bevverts_num);
  bevvert_newverts_offsets_ = Array<int>(this->bevverts_num + 1);
  bevvert_newedges_offsets_ = Array<int>(this->bevverts_num + 1);
  bevvert_newfaces_offsets_ = Array<int>(this->bevverts_num + 1);
  threading::parallel_for(IndexRange(this->bevverts_num), 20'000, [&](IndexRange range) {
    for (const int bv : range) {
      bevvert_meshpatterns_[bv] = topology::find_meshpattern(bv, *this);
      print_meshpattern(bevvert_meshpatterns_[bv]);
      auto [numv, nume, numf, numc] = bevvert_meshpatterns_[bv].num_elements();
      /* In this loop, accumulate counts. Will covert to offsets later. */
      bevvert_newverts_offsets_[bv] = numv;
      bevvert_newedges_offsets_[bv] = nume;
      bevvert_newfaces_offsets_[bv] = numf;
    }
  });
  offset_indices::accumulate_counts_to_offsets(bevvert_newverts_offsets_);
  offset_indices::accumulate_counts_to_offsets(bevvert_newedges_offsets_);
  offset_indices::accumulate_counts_to_offsets(bevvert_newfaces_offsets_);
  bevvert_newverts_ = OffsetIndices<int>(bevvert_newverts_offsets_);
  bevvert_newedges_ = OffsetIndices<int>(bevvert_newedges_offsets_);
  bevvert_newfaces_ = OffsetIndices<int>(bevvert_newfaces_offsets_);
  newvert_positions_ = Array<float3>(bevvert_newverts_.total_size());
}

/** Set the widths of each half of each end of the beveled and non-beveled edges. */
void BevelState::set_bevedge_widths()
{
  bevedge_widths_ = Array<float4>(this->bevedges_num);
  bool is_edge_bevel = this->params.affect_type == BevelAffect::Edges;
  threading::parallel_for(IndexRange(bevverts_num), 10'000, [&](IndexRange range) {
    for (const int bv : range) {
      for (const int epos : IndexRange(bevvert_bevedges_[bv].size())) {
        const int be = bevvert_bevedges_[bv][epos];
        if (is_edge_bevel) {
          spec::edge_bevel_init_widths_from_params(bv, epos, bevedge_widths_[be], *this);
        }
        else {
          spec::vertex_bevel_init_widths_from_params(bv, epos, bevedge_widths_[be], *this);
        }
      }
    }
  });
}

/** Build all the data needed for the vertex meshes. */
void BevelState::build_vertex_meshes()
{
  bevedge_attach_verts_ = Array<int2>(this->bevedges_num);

  threading::parallel_for(IndexRange(bevverts_num), 20'000, [&](IndexRange range) {
    for (const int bv : range) {
      topology::build_vmesh_skeleton(
          bv,
          newvert_positions_.as_mutable_span().slice(bevvert_newverts_[bv]),
          bevedge_attach_verts_.as_mutable_span(),
          *this);
      const int anchors_num = bevvert_meshpatterns_[bv].num_anchors;
      profile::AnchorProfiles profiles(anchors_num);
      profile::calculate_profiles(bv, profiles, *this);
      topology::build_internal_vmesh(
          bv, profiles, newvert_positions_.as_mutable_span().slice(bevvert_newverts_[bv]), *this);
    }
  });
}

std::optional<Mesh *> mesh_bevel(const Mesh &src_mesh,
                                 const IndexMask &selection,
                                 const BevelParameters &params,
                                 const bke::AttributeFilter & /*attribute_filter*/)
{
  if (params.offset <= 0) {
    return std::nullopt;
  }
  std::cout << "\n\nBEVEL, offset = " << params.offset << "\n\n";

  BevelState state(src_mesh, selection, params);
  state.initialize_profile_data();
  dump_bevel_state(state, "after initialization");

  state.order_bevedges();
  state.set_bevedge_widths();
  dump_bevel_state(state, "after construction and ordering bevedges");
  state.set_bevvert_mesh_topology();
  state.build_vertex_meshes();
  dump_bevel_state(state, "after build_vertex_meshes");
  // draw_all_bevverts(state);
  return std::nullopt;
}

}  // namespace blender::geometry
