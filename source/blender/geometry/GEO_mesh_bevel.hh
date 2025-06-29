/* SPDX-FileCopyrightText: 2019 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once
#include <optional>

#include "BLI_index_mask.hh"

#include "BKE_attribute.hh"

struct MDeformVert;
struct Mesh;
struct CurveProfile;

namespace blender::geometry {

/** \warning Each of following enum values  are saved in files.
 * Currently, they also need to be in synch with values in bmesh_opdefines.h
 * and DNA_modifier_types.h. */
enum class BevelOffsetType {
  Offset = 0,
  Width = 1,
  Depth = 2,
  Percent = 3,
  Absolute = 4,
};

enum class BevelProfileType {
  Superellipse = 0,
  Custom = 1,
};

enum class BevelAffect {
  Vertices = 0,
  Edges = 1,
};

enum class BevelMiterType {
  Sharp = 0,
  Patch = 1,
  Arc = 2,
};

enum class BevelVmeshMethod {
  Adj = 0,
  Cutoff = 1,
};

enum class BevelFaceStrengthMode {
  None = 0,
  New = 1,
  Affected = 2,
  All = 3,
};

enum class BevelFaceStrengthValue {
  Weak = -16384,
  Medium = 0,
  Strong = 16384,
};

struct BevelParameters {
  /** Blender units to offset each side of a beveled edge. */
  float offset = 0.0;
  /** How offset is measured. */
  BevelOffsetType offset_type = BevelOffsetType::Offset;
  /** Profile type. */
  BevelProfileType profile_type = BevelProfileType::Superellipse;
  /** Number of segments in beveled edge or beveled vertex profile. */
  int segments = 1;
  /** Profile shape parameter.  */
  float profile = 0.5f;
  /** Bevel vertices only or edges. */
  BevelAffect affect_type = BevelAffect::Edges;
  /** Should we multiply edge width specs by bevel edge weight attribute?
   * Or, if vertex-only, should we multiply distance along arms by vertex bevel weight attribute?
   */
  bool use_weights = false;
  /** Should offsets be limited by collisions? */
  bool limit_offset = false;
  /** Vertex group array, maybe set if vertex only bevel to control which vertices are beveled and
   * by how much. */
  MDeformVert *dvert = nullptr;
  /** Vertex group index, maybe set if vertex only to say which group to get weights from.. */
  int vertex_group = -1;
  /** Material slot to use (if not -1) for material of newly created faces. */
  int mat = -1;
  /** Should bevel prefer to slide along edges rather than keep widths spec? */
  bool loop_slide = true;
  /** Should we propagate seam edge markings? */
  bool mark_seam = false;
  /** Should we propagate sharp edge markings? */
  bool mark_sharp = false;
  /** Should we harden normals? */
  bool harden_normals = false;
  /** Method for setting face strength if not None. */
  BevelFaceStrengthMode face_strength_mode = BevelFaceStrengthMode::None;
  /** What kind of miter pattern to use on reflex angles. */
  BevelMiterType miter_outer = BevelMiterType::Sharp;
  /** What kind of miter pattern to use on non-reflex angles. */
  BevelMiterType miter_inner = BevelMiterType::Sharp;
  /** Amount to spread when doing inside miter. */
  float spread = 0.1f;
  /** The custom profile input, if not null. */
  CurveProfile *custom_profile = nullptr;
  /** The method to use for vertex mesh creation */
  BevelVmeshMethod vmesh_method = BevelVmeshMethod::Adj;
  /** If not -1, the custom data offset for vertex weights. */
  int bweight_offset_vert = -1;
  /** If not -1, the custom data offset for edge weights. */
  int bweight_offset_edge = -1;
};

/**
 * \return #std::nullopt if the mesh is not changed (when every selected face is already a
 * triangle).
 */
std::optional<Mesh *> mesh_bevel(const Mesh &src_mesh,
                                 const IndexMask &selection,
                                 const BevelParameters &params,
                                 const bke::AttributeFilter &attribute_filter);

}  // namespace blender::geometry
