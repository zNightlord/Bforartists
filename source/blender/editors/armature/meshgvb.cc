/* SPDX-FileCopyrightText: 2024 Bforartists. All rights reserved.
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "meshgvb.h"

#include <algorithm>  /* std::sort, std::min, std::max */
#include <cfloat>     /* FLT_MAX */

#include "MEM_guardedalloc.h"

#include "DNA_mesh_types.h"
#include "DNA_object_types.h"

/* Modern BLI C++ headers */
#include "BLI_array.hh"
#include "BLI_index_range.hh"
#include "BLI_math_base.hh"
#include "BLI_math_vector.hh"
#include "BLI_math_vector_types.hh"
#include "BLI_span.hh"
#include "BLI_task.hh"
#include "BLI_vector.hh"

/* BLI_Heap still has no .hh equivalent */
#include "BLI_heap.hh"

#include "BKE_deform.hh"
#include "BKE_mesh.hh"

#include "ED_mesh.hh"
#include "ED_object_vgroup.hh"



namespace blender::ed::armature::gvb {

/* -------------------------------------------------------------------- */
/** \name Voxel grid
 * \{ */

enum class VoxelState : uint8_t {
  Exterior = 0,
  Interior = 1,
  Boundary = 2,
};

struct Grid {
  int3 res;
  float3 origin;
  float3 vsize;   /* per-axis voxel dimensions */
  float diag;     /* bounding-box diagonal = D for Eq. 7 normalization */
  Array<VoxelState> state;

  int total() const
  {
    return res.x * res.y * res.z;
  }

  int idx(int x, int y, int z) const
  {
    return x + y * res.x + z * res.x * res.y;
  }

  int idx(int3 c) const
  {
    return c.x + c.y * res.x + c.z * res.x * res.y;
  }

  bool in_bounds(int3 c) const
  {
    return c.x >= 0 && c.x < res.x &&
           c.y >= 0 && c.y < res.y &&
           c.z >= 0 && c.z < res.z;
  }

  float3 voxel_center(int flat_idx) const
  {
    int3 c = {flat_idx % res.x,
              (flat_idx / res.x) % res.y,
              flat_idx / (res.x * res.y)};
    return origin + (float3(c) + float3(0.5f)) * vsize;
  }

  int world_to_idx(const float3 &p) const
  {
    int3 c = int3((p - origin) / vsize);
    c.x = std::clamp(c.x, 0, res.x - 1);
    c.y = std::clamp(c.y, 0, res.y - 1);
    c.z = std::clamp(c.z, 0, res.z - 1);
    return idx(c);
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 1 — Majority-voting voxelization (Section 4, Eqs. 5 & 6)
 * \{ */

/**
 * Intersect an axis-aligned ray at perpendicular coords (u, v) with a triangle.
 * Returns depth along `axis` (0=X 1=Y 2=Z), or FLT_MAX on miss.
 *
 * Solves the 2×2 system:
 *   tv0[b] + s*(tv1[b]-tv0[b]) + t*(tv2[b]-tv0[b]) = u
 *   tv0[c] + s*(tv1[c]-tv0[c]) + t*(tv2[c]-tv0[c]) = v
 * then depth = tv0[a] + s*(tv1[a]-tv0[a]) + t*(tv2[a]-tv0[a]).
 */
static float ray_tri(int axis,
                     float u,
                     float v,
                     const float3 &tv0,
                     const float3 &tv1,
                     const float3 &tv2)
{
  const int a = axis;
  const int b = (axis + 1) % 3;
  const int c = (axis + 2) % 3;

  const float e1b = tv1[b] - tv0[b], e1c = tv1[c] - tv0[c];
  const float e2b = tv2[b] - tv0[b], e2c = tv2[c] - tv0[c];
  const float det = e1b * e2c - e1c * e2b;

  if (std::abs(det) < 1e-10f) {
    return FLT_MAX;
  }

  const float rb = u - tv0[b];
  const float rc = v - tv0[c];
  const float s  = (rb * e2c - rc * e2b) / det;
  const float t  = (e1b * rc - e1c * rb) / det;

  if (s < 0.0f || t < 0.0f || s + t > 1.0f) {
    return FLT_MAX;
  }

  return tv0[a] + s * (tv1[a] - tv0[a]) + t * (tv2[a] - tv0[a]);
}

/**
 * Cast all rays perpendicular to `axis` and write votes into `votes`.
 *
 * For each (u, v) ray column:
 *   1. Gather all triangle intersection depths along the axis (broad-phase
 *      culled with per-triangle bounding boxes in the perpendicular plane).
 *   2. Sort depths and parity-fill from both ends.
 *   3. OR the two half-results into the vote bit for this axis (Eq. 5).
 */
static void vote_axis(const Grid &grid,
                      MutableSpan<uint8_t> votes,
                      Span<float3> tv0,
                      Span<float3> tv1,
                      Span<float3> tv2,
                      int axis)
{
  const int a        = axis;
  const int b        = (axis + 1) % 3;
  const int c        = (axis + 2) % 3;
  const uint8_t bit  = uint8_t(1u << axis);
  const int    ntris = int(tv0.size());

  /* Precompute per-triangle bounds in the perpendicular (b, c) plane */
  Array<float> tmin_b(ntris), tmax_b(ntris);
  Array<float> tmin_c(ntris), tmax_c(ntris);

  for (int t = 0; t < ntris; t++) {
    tmin_b[t] = std::min({tv0[t][b], tv1[t][b], tv2[t][b]});
    tmax_b[t] = std::max({tv0[t][b], tv1[t][b], tv2[t][b]});
    tmin_c[t] = std::min({tv0[t][c], tv1[t][c], tv2[t][c]});
    tmax_c[t] = std::max({tv0[t][c], tv1[t][c], tv2[t][c]});
  }

  Vector<float> hits;  /* reused each column, avoids per-column allocation */

  for (int iu = 0; iu < grid.res[b]; iu++) {
    const float u = grid.origin[b] + (iu + 0.5f) * grid.vsize[b];

    for (int iv = 0; iv < grid.res[c]; iv++) {
      const float v = grid.origin[c] + (iv + 0.5f) * grid.vsize[c];

      hits.clear();
      for (int t = 0; t < ntris; t++) {
        /* Broad-phase: bounding box cull in perpendicular plane */
        if (u < tmin_b[t] || u > tmax_b[t]) continue;
        if (v < tmin_c[t] || v > tmax_c[t]) continue;

        const float depth = ray_tri(axis, u, v, tv0[t], tv1[t], tv2[t]);
        if (depth != FLT_MAX) {
          hits.append(depth);
        }
      }

      if (hits.is_empty()) {
        continue;
      }

      std::sort(hits.begin(), hits.end());
      const int nhits = int(hits.size());

      /* Write vote bit for voxel at axis-index ia, perp-indices iu/iv */
      const auto set_vote = [&](int ia) {
        int3 coords;
        coords[a] = ia;
        coords[b] = iu;
        coords[c] = iv;
        votes[grid.idx(coords)] |= bit;
      };

      /* Forward pass: walk low→high, toggling parity at each intersection */
      {
        int cross = 0, hi = 0;
        for (int ia = 0; ia < grid.res[a]; ia++) {
          const float center = grid.origin[a] + (ia + 0.5f) * grid.vsize[a];
          while (hi < nhits && hits[hi] <= center) { cross++; hi++; }
          if (cross & 1) set_vote(ia);
        }
      }

      /* Backward pass: same logic from the opposite end (Eq. 5 OR) */
      {
        int cross = 0, hi = nhits - 1;
        for (int ia = grid.res[a] - 1; ia >= 0; ia--) {
          const float center = grid.origin[a] + (ia + 0.5f) * grid.vsize[a];
          while (hi >= 0 && hits[hi] >= center) { cross++; hi--; }
          if (cross & 1) set_vote(ia);
        }
      }
    }
  }
}

/**
 * Classify all voxels as Interior or Exterior via 3-axis majority voting (Eq. 6):
 * interior if at least 2 of 3 axis pairs vote interior.
 */
static void voxelize(Grid &grid,
                     Span<float3> tv0,
                     Span<float3> tv1,
                     Span<float3> tv2)
{
  Array<uint8_t> votes(grid.total(), 0);

  for (int axis = 0; axis < 3; axis++) {
    vote_axis(grid, votes, tv0, tv1, tv2, axis);
  }

  for (int i = 0; i < grid.total(); i++) {
    const uint8_t v     = votes[i];
    const int vote_count = (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1);
    grid.state[i]        = (vote_count >= 2) ? VoxelState::Interior : VoxelState::Exterior;
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 2 — Boundary voxel extraction (Section 4 post-processing)
 *
 * Uses the full Akenine-Möller SAT triangle/AABB overlap test to find all
 * non-exterior voxels touched by mesh geometry and mark them Boundary.
 * \{ */

/**
 * Full 13-axis SAT triangle/AABB overlap test (Akenine-Möller 2001).
 * Axes tested: 3 AABB face normals, 1 triangle normal, 9 edge×cardinal.
 *
 * The `separated` lambda projects the three translated triangle vertices
 * and the box half-extents onto a candidate separating axis, returning
 * true if a gap is found (= shapes are disjoint along this axis).
 */
static bool sat_tri_box(const float3 &center,
                         const float3 &half,
                         const float3 &tv0,
                         const float3 &tv1,
                         const float3 &tv2)
{
  /* Translate triangle into box-local space */
  const float3 v0 = tv0 - center;
  const float3 v1 = tv1 - center;
  const float3 v2 = tv2 - center;

  const float3 e0 = v1 - v0;
  const float3 e1 = v2 - v1;
  const float3 e2 = v0 - v2;

  const auto separated = [&](const float3 &ax) -> bool {
    const float p0  = math::dot(ax, v0);
    const float p1  = math::dot(ax, v1);
    const float p2  = math::dot(ax, v2);
    /* Box projection radius = half · |axis| (component-wise) */
    const float r   = math::dot(math::abs(ax), half);
    return (std::min({p0, p1, p2}) > r || std::max({p0, p1, p2}) < -r);
  };

  /* 1. Three AABB face axes (cardinal directions) */
  if (separated({1, 0, 0}) || separated({0, 1, 0}) || separated({0, 0, 1})) {
    return false;
  }

  /* 2. Triangle face normal */
  if (separated(math::cross(e0, e1))) {
    return false;
  }

  /* 3. Nine cross-product axes: triangle edge × cardinal axis */
  const float3 edges[3] = {e0, e1, e2};
  for (const float3 &e : edges) {
    /* e × X = (0, -e.z, e.y) */
    if (separated({0.0f, -e.z, e.y}))  return false;
    /* e × Y = (e.z, 0, -e.x) */
    if (separated({e.z,  0.0f, -e.x})) return false;
    /* e × Z = (-e.y, e.x, 0) */
    if (separated({-e.y, e.x, 0.0f}))  return false;
  }

  return true; /* no separating axis found */
}

/**
 * For each triangle, iterate over the voxel range its bounding box covers
 * and run an exact SAT test. Non-exterior voxels that pass are marked Boundary.
 */
static void mark_boundary(Grid &grid,
                           Span<float3> tv0,
                           Span<float3> tv1,
                           Span<float3> tv2)
{
  const float3 half = grid.vsize * float3(0.5f);
  const int ntris   = int(tv0.size());

  for (int t = 0; t < ntris; t++) {
    /* Compute the range of voxel indices the triangle can touch */
    int3 vmin, vmax;
    for (int d = 0; d < 3; d++) {
      const float lo = std::min({tv0[t][d], tv1[t][d], tv2[t][d]}) - half[d];
      const float hi = std::max({tv0[t][d], tv1[t][d], tv2[t][d]}) + half[d];
      vmin[d] = std::max(0, int(std::floor((lo - grid.origin[d]) / grid.vsize[d])));
      vmax[d] = std::min(grid.res[d] - 1,
                         int(std::floor((hi - grid.origin[d]) / grid.vsize[d])));
    }

    for (int iz = vmin.z; iz <= vmax.z; iz++) {
      for (int iy = vmin.y; iy <= vmax.y; iy++) {
        for (int ix = vmin.x; ix <= vmax.x; ix++) {
          const int vidx = grid.idx(ix, iy, iz);

          if (grid.state[vidx] == VoxelState::Exterior  ||
              grid.state[vidx] == VoxelState::Boundary)
          {
            continue;
          }

          const float3 center = grid.origin +
                                (float3(int3(ix, iy, iz)) + float3(0.5f)) * grid.vsize;

          if (sat_tri_box(center, half, tv0[t], tv1[t], tv2[t])) {
            grid.state[vidx] = VoxelState::Boundary;
          }
        }
      }
    }
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 3 — Per-bone Dijkstra geodesic distances (Section 5, Alg. 1)
 * \{ */

/* 6-connected neighbor offsets */
static constexpr int3 NBRS[6] = {
  {-1, 0, 0}, {1, 0, 0},
  { 0,-1, 0}, {0, 1, 0},
  { 0, 0,-1}, {0, 0, 1},
};

/**
 * Closest point on segment (a→b) to point p.
 */
static float3 closest_on_segment(const float3 &p, const float3 &a, const float3 &b)
{
  const float3 ab     = b - a;
  const float  ab_sq  = math::dot(ab, ab);
  if (ab_sq < 1e-10f) {
    return a;
  }
  const float t = std::clamp(math::dot(p - a, ab) / ab_sq, 0.0f, 1.0f);
  return a + t * ab;
}

/**
 * Run Dijkstra from bone seeds through all non-exterior voxels.
 * Seeds are voxels whose centre lies within (half voxel diagonal) of the bone.
 * Step cost = Euclidean distance between voxel centres (6-connected).
 * Writes into `dist` (pre-allocated, size = grid.total()).
 * FLT_MAX = unreachable; -1.0f = exterior (skipped).
 */
static void dijkstra(const Grid &grid,
                     MutableSpan<float> dist,
                     const float3 &bone_root,
                     const float3 &bone_tip)
{
  /* Initialize */
  for (int i = 0; i < grid.total(); i++) {
    dist[i] = (grid.state[i] != VoxelState::Exterior) ? FLT_MAX : -1.0f;
  }

  /* Pre-compute 6 Euclidean step costs */
  float step[6];
  for (int n = 0; n < 6; n++) {
    step[n] = math::length(float3(NBRS[n]) * grid.vsize);
  }

  /* Seed radius = half the voxel diagonal */
  const float seed_radius = math::length(grid.vsize) * 0.5f;

  Heap *heap = BLI_heap_new();

  for (int iz = 0; iz < grid.res.z; iz++) {
    for (int iy = 0; iy < grid.res.y; iy++) {
      for (int ix = 0; ix < grid.res.x; ix++) {
        const int vidx = grid.idx(ix, iy, iz);
        if (grid.state[vidx] == VoxelState::Exterior) continue;

        const float3 center   = grid.origin +
                                (float3(int3(ix, iy, iz)) + float3(0.5f)) * grid.vsize;
        const float3 closest  = closest_on_segment(center, bone_root, bone_tip);
        const float  seg_dist = math::distance(center, closest);

        if (seg_dist <= seed_radius) {
          dist[vidx] = 0.0f;
          BLI_heap_insert(heap, 0.0f, POINTER_FROM_INT(vidx));
        }
      }
    }
  }

  /* Main Dijkstra loop — Algorithm 1 in the paper */
  while (!BLI_heap_is_empty(heap)) {
    const float cur_dist = BLI_heap_top_value(heap);
    const int cur = POINTER_AS_INT(BLI_heap_pop_min(heap));

    if (cur_dist > dist[cur] + 1e-6f) continue; /* stale heap entry */

    const int cx =  cur % grid.res.x;
    const int cy = (cur / grid.res.x) % grid.res.y;
    const int cz =  cur / (grid.res.x * grid.res.y);

    for (int n = 0; n < 6; n++) {
      const int3 nc = int3(cx, cy, cz) + NBRS[n];
      if (!grid.in_bounds(nc)) continue;

      const int nidx = grid.idx(nc);
      if (grid.state[nidx] == VoxelState::Exterior) continue;

      const float new_dist = dist[cur] + step[n];
      if (new_dist < dist[nidx]) {
        dist[nidx] = new_dist;
        BLI_heap_insert(heap, new_dist, POINTER_FROM_INT(nidx));
      }
    }
  }

  BLI_heap_free(heap, nullptr);
}

/** \} */

static int nearest_non_exterior_voxel(const Grid &grid, int start_idx)
{
  /* If start voxel is non-exterior, use it directly */
  if (grid.state[start_idx] != VoxelState::Exterior) {
    return start_idx;
  }

  /* BFS outward until we find interior or boundary */
  int cx =  start_idx % grid.res.x;
  int cy = (start_idx / grid.res.x) % grid.res.y;
  int cz =  start_idx / (grid.res.x * grid.res.y);

  /* Search expanding cubic shells — 3 shells covers any surface gap */
  for (int r = 1; r <= 3; r++) {
    for (int dz = -r; dz <= r; dz++) {
      for (int dy = -r; dy <= r; dy++) {
        for (int dx = -r; dx <= r; dx++) {
          if (std::abs(dx) != r && std::abs(dy) != r && std::abs(dz) != r) {
            continue; /* only the shell, not the interior */
          }
          const int3 nc = int3(cx + dx, cy + dy, cz + dz);
          if (!grid.in_bounds(nc)) continue;
          const int nidx = grid.idx(nc);
          if (grid.state[nidx] != VoxelState::Exterior) {
            return nidx;
          }
        }
      }
    }
  }

  return start_idx; /* fallback — distances will be FLT_MAX, vertex skipped */
}

/* -------------------------------------------------------------------- */
/** \name Phase 4 — Weight computation (Section 6, Equations 7 & 8)
 * \{ */

/**
 * Assign vertex weights from per-bone geodesic distance fields.
 *
 * Eq. 7: d_j^i = (d_v^i + |p_vertex − p_voxel|) / D
 * Eq. 8: ω_j^i = (1 / ((1−α)·d + α·d²))²
 */
static void assign_weights(Object *ob,
                            Mesh *mesh,
                            Span<float3> world_verts,
                            int numbones,
                            bDeformGroup **dgrouplist,
                            bDeformGroup **dgroupflip,
                            const bool *selected,
                            Span<Array<float>> bone_dists,
                            const Grid &grid)
{
  using namespace bke;

  const float inv_D = 1.0f / grid.diag;
  const float eps   = 1e-6f;
  constexpr float alpha = 0.7f;

  const bool use_topo_mirror = (mesh->editflag & ME_EDIT_MIRROR_TOPO) != 0;

  /* --- Weight paint mask — mirrors heat_bone_weighting exactly --- */
  const bool use_vert_sel = (mesh->editflag & ME_EDIT_PAINT_VERT_SEL) != 0;
  const bool use_face_sel = (mesh->editflag & ME_EDIT_PAINT_FACE_SEL) != 0;

  Array<bool> mask;   /* empty = no mask, all vertices active */

  if ((ob->mode & OB_MODE_WEIGHT_PAINT) && (use_vert_sel || use_face_sel)) {
    mask = Array<bool>(mesh->verts_num, false);

    const AttributeAccessor attributes = mesh->attributes();
    const OffsetIndices<int> faces    = mesh->faces();
    const Span<int>      corner_verts = mesh->corner_verts();

    if (use_vert_sel) {
      const VArray select_vert = *attributes.lookup_or_default<bool>(
          ".select_vert", AttrDomain::Point, false);
      for (const int i : faces.index_range()) {
        for (const int vert : corner_verts.slice(faces[i])) {
          mask[vert] = select_vert[vert];
        }
      }
    }
    else if (use_face_sel) {
      const VArray select_poly = *attributes.lookup_or_default<bool>(
          ".select_poly", AttrDomain::Face, false);
      for (const int i : faces.index_range()) {
        if (select_poly[i]) {
          for (const int vert : corner_verts.slice(faces[i])) {
            mask[vert] = true;
          }
        }
      }
    }
  }

  Array<float> weights(numbones);

  for (const int vi : IndexRange(mesh->verts_num)) {

    /* Skip masked vertices in weight paint mode */
    if (!mask.is_empty() && !mask[vi]) {
      continue;
    }

    const float3 &p_vert     = world_verts[vi];
    const int     vidx_raw   = grid.world_to_idx(p_vert);
    const int     vidx        = nearest_non_exterior_voxel(grid, vidx_raw);
    const float3  p_vox       = grid.voxel_center(vidx);
    const float   sub_vox     = math::distance(p_vert, p_vox);

    float weight_sum = 0.0f;
    weights.fill(0.0f);

    for (int j = 0; j < numbones; j++) {
      if (!selected[j] || !dgrouplist[j] || bone_dists[j].is_empty()) continue;

      const float d_v = bone_dists[j][vidx];
      if (d_v < 0.0f || d_v == FLT_MAX) continue;

      const float d     = std::max((d_v + sub_vox) * inv_D, eps);
      const float inner = std::max((1.0f - alpha) * d + alpha * d * d, eps);
      const float w     = 1.0f / inner;
      weights[j]        = w * w;
      weight_sum       += weights[j];
    }

    if (weight_sum < 1e-10f) continue;

    const int vi_flip = dgroupflip ?
        mesh_get_x_mirror_vert(ob, nullptr, vi, use_topo_mirror) : -1;

    for (int j = 0; j < numbones; j++) {
      if (!dgrouplist[j] || weights[j] == 0.0f) continue;

      const float w = weights[j] / weight_sum;

      if (w > 0.0f) {
        ed::object::vgroup_vert_add(ob, dgrouplist[j], vi, w, WEIGHT_REPLACE);
      }
      else {
        ed::object::vgroup_vert_remove(ob, dgrouplist[j], vi);
      }

      if (dgroupflip && dgroupflip[j] && vi_flip >= 0) {
        if (w > 0.0f) {
          ed::object::vgroup_vert_add(ob, dgroupflip[j], vi_flip, w, WEIGHT_REPLACE);
        }
        else {
          ed::object::vgroup_vert_remove(ob, dgroupflip[j], vi_flip);
        }
      }
    }
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Triangle list builder (shared by Phases 1 and 2)
 * \{ */

struct TriList {
  Array<float3> v0, v1, v2;
  int size() const { return int(v0.size()); }
};

static TriList build_tris(Mesh *mesh, Span<float3> world_verts)
{
  const Span<int>         corner_verts = mesh->corner_verts();
  const OffsetIndices<int> faces       = mesh->faces();

  /* Count triangles (fan triangulation) */
  int ntris = 0;
  for (const int fi : faces.index_range()) {
    ntris += std::max(0, int(faces[fi].size()) - 2);
  }

  TriList tris;
  tris.v0 = Array<float3>(ntris);
  tris.v1 = Array<float3>(ntris);
  tris.v2 = Array<float3>(ntris);

  int ti = 0;
  for (const int fi : faces.index_range()) {
    const IndexRange face = faces[fi];
    if (face.size() < 3) continue;
    const int i0 = corner_verts[face[0]];
    for (int t = 1; t + 1 < int(face.size()); t++) {
      tris.v0[ti] = world_verts[i0];
      tris.v1[ti] = world_verts[corner_verts[face[t]]];
      tris.v2[ti] = world_verts[corner_verts[face[t + 1]]];
      ti++;
    }
  }

  return tris;
}

/** \} */

}  /* namespace blender::ed::armature::gvb */

/* -------------------------------------------------------------------- */
/** \name Public entry point
 * \{ */

void geodesic_voxel_bone_weighting(blender::Object *ob,
                                   blender::Mesh *mesh,
                                   float (*verts)[3],
                                   int numbones,
                                   blender::bDeformGroup **dgrouplist,
                                   blender::bDeformGroup **dgroupflip,
                                   float (*root)[3],
                                   float (*tip)[3],
                                   const bool *selected,
                                   const char **r_error_str)
{
  using namespace blender;
  using namespace blender::ed::armature::gvb;

  *r_error_str = nullptr;

  /* Reinterpret C arrays as typed spans — layout is identical */
  const Span<float3> verts_span(reinterpret_cast<const float3 *>(verts),
                                 mesh->verts_num);
  const Span<float3> root_span(reinterpret_cast<const float3 *>(root), numbones);
  const Span<float3> tip_span(reinterpret_cast<const float3 *>(tip),   numbones);

  /* --- Build world-space AABB --- */
  float3 bbox_min(FLT_MAX);
  float3 bbox_max(-FLT_MAX);

  for (const float3 &v : verts_span) {
    bbox_min = math::min(bbox_min, v);
    bbox_max = math::max(bbox_max, v);
  }

  /* Pad by one voxel on every side so surface voxels are never on the face */
  const int3    res    = int3(GVB_DEFAULT_RES_X, GVB_DEFAULT_RES_Y, GVB_DEFAULT_RES_Z);
  const float3  pad    = (bbox_max - bbox_min) / float3(res) + float3(1e-4f);
  bbox_min -= pad;
  bbox_max += pad;

  const float3 extent = bbox_max - bbox_min;

  /* --- Initialize grid --- */
  Grid grid;
  grid.res    = res;
  grid.origin = bbox_min;
  grid.vsize  = extent / float3(res);
  grid.diag   = math::length(extent); /* D in Eq. 7 */
  grid.state  = Array<VoxelState>(grid.total(), VoxelState::Exterior);

  /* Build triangle list once — shared by Phases 1 and 2 */
  const TriList tris = build_tris(mesh, verts_span);

  if (tris.size() == 0) {
    *r_error_str = "GVB: mesh produced no triangles after fan triangulation";
    return;
  }

  /* --- Phase 1: Majority-voting voxelization (Section 4, Eq. 5 & 6) --- */
  voxelize(grid, tris.v0, tris.v1, tris.v2);

  /* --- Phase 2: Boundary voxel extraction (Section 4 post-processing) --- */
  mark_boundary(grid, tris.v0, tris.v1, tris.v2);

  /* --- Phase 3: Per-bone Dijkstra, parallelized across bones (Section 5) --- */
  Array<Array<float>> bone_dists(numbones);

  for (int j = 0; j < numbones; j++) {
    if (!selected[j] || !dgrouplist[j]) continue;
    bone_dists[j] = Array<float>(grid.total(), 0.0f);
  }

  threading::parallel_for(IndexRange(numbones), 1, [&](IndexRange range) {
    for (const int j : range) {
      if (!selected[j] || !dgrouplist[j] || bone_dists[j].is_empty()) continue;
      dijkstra(grid, bone_dists[j], root_span[j], tip_span[j]);
    }
  });

  /* --- Phase 4: Weight assignment (Section 6, Eq. 7 & 8) --- */
  assign_weights(ob, mesh, verts_span, numbones,
                 dgrouplist, dgroupflip, selected,
                 bone_dists, grid);
}

/** \} */
