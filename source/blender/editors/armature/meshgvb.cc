/* SPDX-FileCopyrightText: 2024 Bforartists. All rights reserved.
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edarmature
 *
 * Four-phase pipeline (paper section numbers in comments):
 *   Phase 1 — Voxelization via 3-axis majority voting   (Section 4)
 *   Phase 2 — Boundary voxel extraction via SAT test    (Section 4 post-proc)
 *   Phase 3 — Per-bone Dijkstra geodesic distances      (Section 5, Algorithm 1)
 *   Phase 4 — Weight assignment via Equations 7 & 8     (Section 6)
 */

#include "meshgvb.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "MEM_guardedalloc.h"

#include "DNA_mesh_types.h"
#include "DNA_object_types.h"

#include "BLI_heap.h"
#include "BLI_math_geom.h"
#include "BLI_math_vector.h"
#include "BLI_task.h"

#include "BKE_deform.hh"
#include "BKE_mesh.hh"

#include "ED_object_vgroup.hh"

/* -------------------------------------------------------------------- */
/** \name Internal grid types
 * \{ */

enum GVBState : uint8_t {
  GVB_EXTERIOR = 0, /* outside mesh, never used for distance */
  GVB_INTERIOR = 1, /* inside mesh, traversed by Dijkstra */
  GVB_BOUNDARY = 2, /* on the mesh surface, receives weights */
};

/* Vote bits — one per axis pair (Eq. 5) */
#define GVB_VOTE_X (1u << 0)
#define GVB_VOTE_Y (1u << 1)
#define GVB_VOTE_Z (1u << 2)

struct GVBGrid {
  int res[3];          /* voxel counts along each axis */
  float origin[3];     /* world-space AABB minimum corner */
  float vsize[3];      /* per-axis voxel dimensions */
  float diag;          /* bbox diagonal = D in Eq. 7, for distance normalization */
  GVBState *state;     /* flat array [res[0]*res[1]*res[2]] */
};

static int gvb_idx(const GVBGrid &g, int x, int y, int z)
{
  return x + y * g.res[0] + z * g.res[0] * g.res[1];
}

static bool gvb_in_bounds(const GVBGrid &g, int x, int y, int z)
{
  return (x >= 0 && x < g.res[0]) &&
         (y >= 0 && y < g.res[1]) &&
         (z >= 0 && z < g.res[2]);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 1 — Majority-voting voxelization (Section 4, Eqs. 5 & 6)
 *
 * For each of the 3 axis pairs (±X, ±Y, ±Z) we cast parallel rays through
 * the mesh and use odd-even parity counting to vote each voxel inside or
 * outside. A voxel is classified interior if at least 2 of the 3 pairs
 * agree (Eq. 6). This is robust to non-manifold, non-watertight, and
 * intersecting geometry because the majority vote breaks tie cases.
 * \{ */

/**
 * Intersect an axis-aligned ray at (u, v) with a triangle (tv0..tv2).
 * Returns the coordinate along `axis` where the ray crosses the triangle,
 * or FLT_MAX if it misses.
 *
 * axis 0=X  → ray parallel to X, u=Y, v=Z
 * axis 1=Y  → ray parallel to Y, u=Z, v=X
 * axis 2=Z  → ray parallel to Z, u=X, v=Y
 *
 * Derivation: solve  tv0[b] + s*(tv1[b]-tv0[b]) + t*(tv2[b]-tv0[b]) = u
 *                    tv0[c] + s*(tv1[c]-tv0[c]) + t*(tv2[c]-tv0[c]) = v
 * for barycentric (s,t), then depth = tv0[a] + s*e1[a] + t*e2[a].
 */
static float gvb_ray_tri(int axis,
                          float u,
                          float v,
                          const float tv0[3],
                          const float tv1[3],
                          const float tv2[3])
{
  const int a  = axis;           /* parallel index */
  const int b  = (axis + 1) % 3;
  const int c  = (axis + 2) % 3;

  float e1b = tv1[b] - tv0[b], e1c = tv1[c] - tv0[c];
  float e2b = tv2[b] - tv0[b], e2c = tv2[c] - tv0[c];
  float det = e1b * e2c - e1c * e2b;

  if (fabsf(det) < 1e-10f) {
    return FLT_MAX; /* degenerate triangle or ray parallel to it */
  }

  float rb = u - tv0[b];
  float rc = v - tv0[c];
  float s  = (rb * e2c - rc * e2b) / det;
  float t  = (e1b * rc - e1c * rb) / det;

  if (s < 0.0f || t < 0.0f || s + t > 1.0f) {
    return FLT_MAX; /* ray misses triangle */
  }

  return tv0[a] + s * (tv1[a] - tv0[a]) + t * (tv2[a] - tv0[a]);
}

/**
 * Cast all rays perpendicular to `axis` and accumulate votes into `votes[]`.
 *
 * For each ray column at perpendicular coords (iu, iv):
 *   1. Find all triangle intersection depths along this axis, sort them.
 *   2. Walk forward (low→high): after an odd number of crossings, vote interior.
 *   3. Walk backward (high→low): same parity from the opposite direction.
 *   4. OR both passes into the vote bit for this axis (Eq. 5).
 *
 * Triangle bounding boxes in the perpendicular plane are used to skip
 * irrelevant triangles for each column (the main performance optimization).
 */
static void gvb_vote_axis(const GVBGrid &grid,
                           uint8_t *votes,
                           int ntris,
                           const float (*tv0)[3],
                           const float (*tv1)[3],
                           const float (*tv2)[3],
                           int axis)
{
  const int a  = axis;
  const int b  = (axis + 1) % 3;
  const int c  = (axis + 2) % 3;

  const uint8_t vote_bit = (uint8_t)(1u << axis);

  /* Precompute per-triangle bounding intervals in (b, c) directions */
  std::vector<float> tmin_b(ntris), tmax_b(ntris);
  std::vector<float> tmin_c(ntris), tmax_c(ntris);

  for (int t = 0; t < ntris; t++) {
    tmin_b[t] = min_fff(tv0[t][b], tv1[t][b], tv2[t][b]);
    tmax_b[t] = max_fff(tv0[t][b], tv1[t][b], tv2[t][b]);
    tmin_c[t] = min_fff(tv0[t][c], tv1[t][c], tv2[t][c]);
    tmax_c[t] = max_fff(tv0[t][c], tv1[t][c], tv2[t][c]);
  }

  std::vector<float> hits; /* reused per column */

  for (int iu = 0; iu < grid.res[b]; iu++) {
    float u = grid.origin[b] + (iu + 0.5f) * grid.vsize[b];

    for (int iv = 0; iv < grid.res[c]; iv++) {
      float v = grid.origin[c] + (iv + 0.5f) * grid.vsize[c];

      hits.clear();

      for (int t = 0; t < ntris; t++) {
        /* Broad-phase: skip triangles whose (b,c) bbox excludes this ray */
        if (u < tmin_b[t] || u > tmax_b[t]) continue;
        if (v < tmin_c[t] || v > tmax_c[t]) continue;

        float depth = gvb_ray_tri(axis, u, v, tv0[t], tv1[t], tv2[t]);
        if (depth != FLT_MAX) {
          hits.push_back(depth);
        }
      }

      if (hits.empty()) {
        continue;
      }

      std::sort(hits.begin(), hits.end());
      const int nhits = (int)hits.size();

      /* Helper: coords[a]=ia, coords[b]=iu, coords[c]=iv → flat index */
      auto set_vote = [&](int ia) {
        int coords[3];
        coords[a] = ia;
        coords[b] = iu;
        coords[c] = iv;
        votes[gvb_idx(grid, coords[0], coords[1], coords[2])] |= vote_bit;
      };

      /* Forward pass (+axis direction): odd crossing count = inside */
      {
        int cross = 0, hi = 0;
        for (int ia = 0; ia < grid.res[a]; ia++) {
          float center = grid.origin[a] + (ia + 0.5f) * grid.vsize[a];
          while (hi < nhits && hits[hi] <= center) { cross++; hi++; }
          if (cross & 1) set_vote(ia);
        }
      }

      /* Backward pass (-axis direction): same logic from the other end */
      {
        int cross = 0, hi = nhits - 1;
        for (int ia = grid.res[a] - 1; ia >= 0; ia--) {
          float center = grid.origin[a] + (ia + 0.5f) * grid.vsize[a];
          while (hi >= 0 && hits[hi] >= center) { cross++; hi--; }
          if (cross & 1) set_vote(ia);
        }
      }
    }
  }
}

static void gvb_voxelize(GVBGrid &grid,
                          int ntris,
                          const float (*tv0)[3],
                          const float (*tv1)[3],
                          const float (*tv2)[3])
{
  const int total = grid.res[0] * grid.res[1] * grid.res[2];

  /* Per-voxel vote accumulator: bits 0/1/2 = X/Y/Z pair voted interior */
  uint8_t *votes = MEM_calloc_arrayN<uint8_t>(total, "gvb_votes");

  /* Run voting for all 3 axis pairs */
  for (int axis = 0; axis < 3; axis++) {
    gvb_vote_axis(grid, votes, ntris, tv0, tv1, tv2, axis);
  }

  /* Eq. 6: interior if at least 2 of the 3 pairs agree */
  for (int i = 0; i < total; i++) {
    uint8_t v = votes[i];
    int count = (v & 1) + ((v >> 1) & 1) + ((v >> 2) & 1);
    grid.state[i] = (count >= 2) ? GVB_INTERIOR : GVB_EXTERIOR;
  }

  MEM_freeN(votes);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 2 — Boundary voxel extraction (Section 4 post-processing)
 *
 * After interior/exterior classification, identify boundary voxels:
 * non-exterior voxels that are actually touched by mesh geometry.
 * We use the Akenine-Möller separating axis theorem (SAT) triangle/box
 * overlap test to find these accurately (paper: octree + SAT).
 * \{ */

/**
 * Full Akenine-Möller SAT triangle-AABB overlap test.
 * Returns true if triangle (tv0,tv1,tv2) overlaps box at `center` ± `half`.
 *
 * Tests 13 separating axes:
 *   3  — AABB face normals (cardinal axes)
 *   1  — triangle face normal
 *   9  — cross products of 3 triangle edges × 3 cardinal axes
 */
static bool gvb_sat_tri_box(const float center[3],
                              const float half[3],
                              const float tv0[3],
                              const float tv1[3],
                              const float tv2[3])
{
  /* Translate triangle into box-local space */
  float v0[3], v1[3], v2[3];
  sub_v3_v3v3(v0, tv0, center);
  sub_v3_v3v3(v1, tv1, center);
  sub_v3_v3v3(v2, tv2, center);

  float e0[3], e1[3], e2[3];
  sub_v3_v3v3(e0, v1, v0);
  sub_v3_v3v3(e1, v2, v1);
  sub_v3_v3v3(e2, v0, v2);

  /* Project 3 triangle points and box onto a separating axis.
   * Return true (= SEPARATED) if the projections don't overlap. */
  auto separated = [&](float ax, float ay, float az) -> bool {
    float p0 = ax * v0[0] + ay * v0[1] + az * v0[2];
    float p1 = ax * v1[0] + ay * v1[1] + az * v1[2];
    float p2 = ax * v2[0] + ay * v2[1] + az * v2[2];
    float r  = half[0] * fabsf(ax) + half[1] * fabsf(ay) + half[2] * fabsf(az);
    return (min_fff(p0, p1, p2) > r || max_fff(p0, p1, p2) < -r);
  };

  /* 1. Three AABB face axes */
  if (separated(1, 0, 0) || separated(0, 1, 0) || separated(0, 0, 1)) return false;

  /* 2. Triangle face normal */
  float n[3]; cross_v3_v3v3(n, e0, e1);
  if (separated(n[0], n[1], n[2])) return false;

  /* 3. Nine edge × cardinal cross-product axes */
  const float edges[3][3] = {{e0[0], e0[1], e0[2]},
                               {e1[0], e1[1], e1[2]},
                               {e2[0], e2[1], e2[2]}};
  for (int ei = 0; ei < 3; ei++) {
    /* cross with X axis = (0, -e[2], e[1]) */
    if (separated(0.0f, -edges[ei][2], edges[ei][1])) return false;
    /* cross with Y axis = (e[2], 0, -e[0]) */
    if (separated(edges[ei][2], 0.0f, -edges[ei][0])) return false;
    /* cross with Z axis = (-e[1], e[0], 0) */
    if (separated(-edges[ei][1], edges[ei][0], 0.0f)) return false;
  }

  return true; /* no separating axis → overlapping */
}

/**
 * Mark all non-exterior voxels that contain mesh geometry as GVB_BOUNDARY.
 * Iterates each triangle, computes the range of voxels it can touch from
 * its bounding box, then does an exact SAT test per candidate voxel.
 */
static void gvb_mark_boundary(GVBGrid &grid,
                               int ntris,
                               const float (*tv0)[3],
                               const float (*tv1)[3],
                               const float (*tv2)[3])
{
  float half[3] = {grid.vsize[0] * 0.5f,
                   grid.vsize[1] * 0.5f,
                   grid.vsize[2] * 0.5f};

  for (int t = 0; t < ntris; t++) {
    /* Voxel index range the triangle's bbox can touch */
    int imin[3], imax[3];
    for (int d = 0; d < 3; d++) {
      float lo = min_fff(tv0[t][d], tv1[t][d], tv2[t][d]) - half[d];
      float hi = max_fff(tv0[t][d], tv1[t][d], tv2[t][d]) + half[d];
      imin[d] = max_ii(0, (int)floorf((lo - grid.origin[d]) / grid.vsize[d]));
      imax[d] = min_ii(grid.res[d] - 1,
                       (int)floorf((hi - grid.origin[d]) / grid.vsize[d]));
    }

    for (int iz = imin[2]; iz <= imax[2]; iz++) {
      for (int iy = imin[1]; iy <= imax[1]; iy++) {
        for (int ix = imin[0]; ix <= imax[0]; ix++) {

          int vidx = gvb_idx(grid, ix, iy, iz);
          if (grid.state[vidx] == GVB_EXTERIOR) continue;
          if (grid.state[vidx] == GVB_BOUNDARY) continue; /* already marked */

          float center[3] = {
            grid.origin[0] + (ix + 0.5f) * grid.vsize[0],
            grid.origin[1] + (iy + 0.5f) * grid.vsize[1],
            grid.origin[2] + (iz + 0.5f) * grid.vsize[2],
          };

          if (gvb_sat_tri_box(center, half, tv0[t], tv1[t], tv2[t])) {
            grid.state[vidx] = GVB_BOUNDARY;
          }
        }
      }
    }
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 3 — Per-bone Dijkstra geodesic distances (Section 5, Alg. 1)
 *
 * For each bone, seed the voxels it intersects (distance = 0) and run
 * Dijkstra through all non-exterior voxels using 6-connectivity with
 * Euclidean step costs between voxel centers (paper Section 8 discussion:
 * Euclidean outperformed Manhattan and gave better results for non-uniform
 * voxel grids).
 * \{ */

/* 6-connected neighbors */
static const int GVB_NBRS[6][3] = {
  {-1, 0, 0}, {1, 0, 0},
  { 0,-1, 0}, {0, 1, 0},
  { 0, 0,-1}, {0, 0, 1},
};

/**
 * Fill `dist[total]` with the geodesic distance from bone segment
 * (root→tip) to every non-exterior voxel.
 *
 * Seeds: voxels whose center is within (voxel diagonal * 0.5) of the
 * bone segment — enough to catch any voxel the segment passes through.
 *
 * dist[] is pre-allocated by the caller (size = res[0]*res[1]*res[2]).
 * Unreachable voxels keep the value FLT_MAX.
 */
static void gvb_dijkstra(const GVBGrid &grid,
                          float *dist,
                          const float root[3],
                          const float tip[3])
{
  const int total = grid.res[0] * grid.res[1] * grid.res[2];

  /* Initialize: non-exterior = FLT_MAX, exterior = -1 (never entered) */
  for (int i = 0; i < total; i++) {
    dist[i] = (grid.state[i] != GVB_EXTERIOR) ? FLT_MAX : -1.0f;
  }

  /* Pre-compute the 6 Euclidean step costs */
  float step[6];
  for (int n = 0; n < 6; n++) {
    float d[3] = {GVB_NBRS[n][0] * grid.vsize[0],
                  GVB_NBRS[n][1] * grid.vsize[1],
                  GVB_NBRS[n][2] * grid.vsize[2]};
    step[n] = len_v3(d);
  }

  BLI_Heap *heap = BLI_heap_new();

  /* Bone intersection radius: half voxel diagonal */
  float seed_radius = len_v3(grid.vsize) * 0.5f;

  /* Seed all voxels that the bone segment passes through */
  for (int iz = 0; iz < grid.res[2]; iz++) {
    for (int iy = 0; iy < grid.res[1]; iy++) {
      for (int ix = 0; ix < grid.res[0]; ix++) {
        int vidx = gvb_idx(grid, ix, iy, iz);
        if (grid.state[vidx] == GVB_EXTERIOR) continue;

        float center[3] = {
          grid.origin[0] + (ix + 0.5f) * grid.vsize[0],
          grid.origin[1] + (iy + 0.5f) * grid.vsize[1],
          grid.origin[2] + (iz + 0.5f) * grid.vsize[2],
        };

        float closest[3];
        closest_to_line_segment_v3(closest, center, root, tip);

        if (len_v3v3(center, closest) <= seed_radius) {
          dist[vidx] = 0.0f;
          BLI_heap_insert(heap, 0.0f, POINTER_FROM_INT(vidx));
        }
      }
    }
  }

  /* Dijkstra main loop — Algorithm 1 in the paper */
  while (!BLI_heap_is_empty(heap)) {
    float cur_dist;
    int cur = POINTER_AS_INT(BLI_heap_pop_min(heap, &cur_dist));

    /* Skip stale heap entries */
    if (cur_dist > dist[cur] + 1e-6f) continue;

    int cx =  cur % grid.res[0];
    int cy = (cur / grid.res[0]) % grid.res[1];
    int cz =  cur / (grid.res[0] * grid.res[1]);

    for (int n = 0; n < 6; n++) {
      int nx = cx + GVB_NBRS[n][0];
      int ny = cy + GVB_NBRS[n][1];
      int nz = cz + GVB_NBRS[n][2];

      if (!gvb_in_bounds(grid, nx, ny, nz)) continue;

      int nidx = gvb_idx(grid, nx, ny, nz);
      if (grid.state[nidx] == GVB_EXTERIOR) continue;

      float new_dist = dist[cur] + step[n];
      if (new_dist < dist[nidx]) {
        dist[nidx] = new_dist;
        BLI_heap_insert(heap, new_dist, POINTER_FROM_INT(nidx));
      }
    }
  }

  BLI_heap_free(heap, nullptr);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Phase 4 — Weight computation (Section 6, Equations 7 & 8)
 *
 * Eq. 7 normalizes raw geodesic distance and adds a sub-voxel correction
 *       to compensate for grid coarseness when multiple vertices land in
 *       the same voxel.
 *
 * Eq. 8 converts normalized distance to weight influence. The alpha
 *       parameter controls bind stiffness (paper default = 0.7).
 *
 * Weights are then normalized across all influencing bones before assignment.
 * \{ */

static int gvb_world_to_voxel(const GVBGrid &grid, const float p[3])
{
  int ix = (int)((p[0] - grid.origin[0]) / grid.vsize[0]);
  int iy = (int)((p[1] - grid.origin[1]) / grid.vsize[1]);
  int iz = (int)((p[2] - grid.origin[2]) / grid.vsize[2]);
  CLAMP(ix, 0, grid.res[0] - 1);
  CLAMP(iy, 0, grid.res[1] - 1);
  CLAMP(iz, 0, grid.res[2] - 1);
  return gvb_idx(grid, ix, iy, iz);
}

static void gvb_voxel_center(const GVBGrid &grid, int idx, float r[3])
{
  int ix =  idx % grid.res[0];
  int iy = (idx / grid.res[0]) % grid.res[1];
  int iz =  idx / (grid.res[0] * grid.res[1]);
  r[0] = grid.origin[0] + (ix + 0.5f) * grid.vsize[0];
  r[1] = grid.origin[1] + (iy + 0.5f) * grid.vsize[1];
  r[2] = grid.origin[2] + (iz + 0.5f) * grid.vsize[2];
}

static void gvb_assign_weights(Object *ob,
                                Mesh *mesh,
                                float (*world_verts)[3],
                                int numbones,
                                bDeformGroup **dgrouplist,
                                bDeformGroup **dgroupflip,
                                const bool *selected,
                                float **bone_dists,
                                const GVBGrid &grid,
                                float alpha)
{
  const float inv_D = 1.0f / grid.diag;  /* Eq. 7 normalization factor */
  const float eps   = 1e-6f;             /* minimum distance (avoids div/0 in Eq. 8) */

  blender::Array<float> weights(numbones);

  for (int vi = 0; vi < mesh->verts_num; vi++) {
    const float *p_vert = world_verts[vi];

    /* Find the voxel containing this vertex */
    int vidx = gvb_world_to_voxel(grid, p_vert);

    /* Sub-voxel correction: exact vertex-to-voxel-center offset (Eq. 7) */
    float p_vox[3];
    gvb_voxel_center(grid, vidx, p_vox);
    float sub_vox = len_v3v3(p_vert, p_vox);

    float weight_sum = 0.0f;
    weights.fill(0.0f);

    for (int j = 0; j < numbones; j++) {
      if (!selected[j] || !dgrouplist[j] || !bone_dists[j]) continue;

      float d_v = bone_dists[j][vidx];
      if (d_v < 0.0f || d_v == FLT_MAX) continue;

      /* Eq. 7: normalized distance incorporating sub-voxel correction */
      float d = (d_v + sub_vox) * inv_D;
      d = max_ff(d, eps);

      /* Eq. 8: falloff function with alpha blend-smoothness control */
      float inner = (1.0f - alpha) * d + alpha * (d * d);
      inner = max_ff(inner, eps);
      float w = 1.0f / inner;
      weights[j] = w * w;   /* squared to get the (...)^2 in Eq. 8 */

      weight_sum += weights[j];
    }

    if (weight_sum < 1e-10f) continue;

    /* Normalize and assign to deform groups */
    for (int j = 0; j < numbones; j++) {
      if (!dgrouplist[j] || weights[j] == 0.0f) continue;

      float w = weights[j] / weight_sum;
      blender::ed::object::vgroup_vert_add(ob, dgrouplist[j], vi, w, WEIGHT_REPLACE);

      /* Mirror vertex weight if flip group is present */
      if (dgroupflip && dgroupflip[j]) {
        int vi_flip = mesh_get_x_mirror_vert(
            ob, nullptr, vi, (mesh->editflag & ME_EDIT_MIRROR_TOPO) != 0);
        if (vi_flip != -1) {
          blender::ed::object::vgroup_vert_add(
              ob, dgroupflip[j], vi_flip, w, WEIGHT_REPLACE);
        }
      }
    }
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Triangle list builder (shared by Phases 1 and 2)
 * \{ */

struct GVBTriangles {
  int count;
  float (*v0)[3];
  float (*v1)[3];
  float (*v2)[3];
};

static GVBTriangles gvb_build_tris(const Mesh *mesh, float (*world_verts)[3])
{
  using namespace blender;
  const Span<int> corner_verts  = mesh->corner_verts();
  const OffsetIndices faces     = mesh->faces();

  /* Count triangles from fan triangulation */
  int ntris = 0;
  for (const int fi : faces.index_range()) {
    ntris += max_ii(0, (int)faces[fi].size() - 2);
  }

  GVBTriangles tris;
  tris.count = ntris;
  tris.v0 = MEM_calloc_arrayN<float[3]>(ntris, "gvb_tv0");
  tris.v1 = MEM_calloc_arrayN<float[3]>(ntris, "gvb_tv1");
  tris.v2 = MEM_calloc_arrayN<float[3]>(ntris, "gvb_tv2");

  int ti = 0;
  for (const int fi : faces.index_range()) {
    const IndexRange face = faces[fi];
    if (face.size() < 3) continue;
    int i0 = corner_verts[face[0]];
    for (int t = 1; t + 1 < (int)face.size(); t++) {
      int i1 = corner_verts[face[t]];
      int i2 = corner_verts[face[t + 1]];
      copy_v3_v3(tris.v0[ti], world_verts[i0]);
      copy_v3_v3(tris.v1[ti], world_verts[i1]);
      copy_v3_v3(tris.v2[ti], world_verts[i2]);
      ti++;
    }
  }
  return tris;
}

static void gvb_free_tris(GVBTriangles &tris)
{
  MEM_freeN(tris.v0);
  MEM_freeN(tris.v1);
  MEM_freeN(tris.v2);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Parallel Dijkstra task data
 * \{ */

struct GVBDijkstraTask {
  const GVBGrid *grid;
  float **bone_dists;
  float (*root)[3];
  float (*tip)[3];
  const bool *selected;
  bDeformGroup **dgrouplist;
};

static void gvb_dijkstra_task_cb(void *__restrict userdata,
                                  const int j,
                                  const TaskParallelTLS *__restrict /*tls*/)
{
  GVBDijkstraTask *td = static_cast<GVBDijkstraTask *>(userdata);
  if (!td->selected[j] || !td->dgrouplist[j] || !td->bone_dists[j]) return;
  gvb_dijkstra(*td->grid, td->bone_dists[j], td->root[j], td->tip[j]);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Main entry point
 * \{ */

void geodesic_voxel_bone_weighting(Object *ob,
                                   Mesh *mesh,
                                   float (*verts)[3],
                                   int numbones,
                                   bDeformGroup **dgrouplist,
                                   bDeformGroup **dgroupflip,
                                   float (*root)[3],
                                   float (*tip)[3],
                                   const bool *selected,
                                   float alpha,
                                   const int res[3],
                                   const char **error_str)
{
  *error_str = nullptr;

  /* --- Build world-space AABB --- */
  float bbox_min[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
  float bbox_max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  for (int i = 0; i < mesh->verts_num; i++) {
    minmax_v3v3_v3(bbox_min, bbox_max, verts[i]);
  }

  /* Pad by one voxel so surface voxels are never on the AABB face */
  float extent[3];
  sub_v3_v3v3(extent, bbox_max, bbox_min);
  for (int d = 0; d < 3; d++) {
    float pad = (extent[d] / res[d]) + 1e-4f;
    bbox_min[d] -= pad;
    bbox_max[d] += pad;
    extent[d]    = bbox_max[d] - bbox_min[d];
  }

  /* --- Initialize grid --- */
  GVBGrid grid = {};
  grid.res[0] = res[0]; grid.res[1] = res[1]; grid.res[2] = res[2];
  copy_v3_v3(grid.origin, bbox_min);
  for (int d = 0; d < 3; d++) {
    grid.vsize[d] = extent[d] / grid.res[d];
  }
  grid.diag = len_v3(extent); /* D in Eq. 7 */

  const int total = grid.res[0] * grid.res[1] * grid.res[2];
  grid.state = MEM_calloc_arrayN<GVBState>(total, "gvb_state");

  /* Triangulate mesh once — shared by Phases 1 and 2 */
  GVBTriangles tris = gvb_build_tris(mesh, verts);

  if (tris.count == 0) {
    *error_str = "GVB: mesh has no triangles after triangulation";
    MEM_freeN(grid.state);
    return;
  }

  /* ---- Phase 1: Majority-voting voxelization ---- */
  gvb_voxelize(grid, tris.count,
               (const float (*)[3])tris.v0,
               (const float (*)[3])tris.v1,
               (const float (*)[3])tris.v2);

  /* ---- Phase 2: Boundary voxel marking ---- */
  gvb_mark_boundary(grid, tris.count,
                    (const float (*)[3])tris.v0,
                    (const float (*)[3])tris.v1,
                    (const float (*)[3])tris.v2);

  gvb_free_tris(tris);

  /* ---- Phase 3: Per-bone Dijkstra (parallelized across bones) ---- */
  float **bone_dists = MEM_calloc_arrayN<float *>(numbones, "gvb_bone_dists");
  for (int j = 0; j < numbones; j++) {
    if (!selected[j] || !dgrouplist[j]) continue;
    bone_dists[j] = MEM_calloc_arrayN<float>(total, "gvb_dist_j");
  }

  GVBDijkstraTask td = {&grid, bone_dists, root, tip, selected, dgrouplist};

  TaskParallelSettings settings;
  BLI_parallel_range_settings_defaults(&settings);
  settings.use_threading = (numbones > 4); /* parallel only when worthwhile */
  BLI_task_parallel_range(0, numbones, &td, gvb_dijkstra_task_cb, &settings);

  /* ---- Phase 4: Weight computation and vertex group assignment ---- */
  gvb_assign_weights(ob, mesh, verts, numbones,
                     dgrouplist, dgroupflip, selected,
                     bone_dists, grid, alpha);

  /* --- Free per-bone distance fields --- */
  for (int j = 0; j < numbones; j++) {
    if (bone_dists[j]) MEM_freeN(bone_dists[j]);
  }
  MEM_freeN(bone_dists);
  MEM_freeN(grid.state);
}

/** \} */
