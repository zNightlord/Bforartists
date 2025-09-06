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
enum class BevelAffect {
  Vertices = 0,
  Edges = 1,
  Faces = 2,
};

enum class BevelMiterType {
  Sharp = 0,
  Patch = 1,
  Arc = 2,
};

struct BevelParameters {
  int segments = 1;
  /** Profile shape parameter.  */
  float shape = 0.5f;
  /** Bevel vertices only or edges. */
  BevelAffect affect_type = BevelAffect::Edges;
  /** The custom profile input, if not null. */
  CurveProfile *custom_profile = nullptr;
  /** Blender units to offset eacjh end of each edge.
   * A 4d Array of Arrays indexed by mesh edge id.
   * If affect_type is Edges or Faces, these are in order: source end (left, right), destination
   * end (left, right). If affect_type is Vertices, these are the amounts to move along each edge
   * from the vertex, and only the first and third values are used. */
  Array<Array<float>, 4> offsets;
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
