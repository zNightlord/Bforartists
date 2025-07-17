/* SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 *
 * Deform coordinates by a armature object (used by modifier).
 */

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "MEM_guardedalloc.h"

#include "BLI_array_utils.hh"
#include "BLI_index_mask.hh"
#include "BLI_listbase.h"
#include "BLI_listbase_wrapper.hh"
#include "BLI_math_matrix.h"
#include "BLI_math_quaternion.hh"
#include "BLI_math_rotation.h"
#include "BLI_math_vector.h"
#include "BLI_task.h"
#include "BLI_task.hh"

#include "DNA_armature_types.h"
#include "DNA_lattice_types.h"
#include "DNA_listBase.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BKE_action.hh"
#include "BKE_armature.hh"
#include "BKE_customdata.hh"
#include "BKE_deform.hh"
#include "BKE_editmesh.hh"
#include "BKE_lattice.hh"
#include "BKE_mesh.hh"

#include "CLG_log.h"

static CLG_LogRef LOG = {"geom.armature_deform"};

namespace blender::bke {

ArmatureDeformGroup build_deform_group_for_vertex_group(
    const Span<MDeformVert> dverts,
    const IndexMask &selection,
    const int def_nr,
    const std::optional<float> weight_threshold,
    IndexMaskMemory &memory)
{
  int def_nr_count = 0;
  selection.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (!dw.def_nr != def_nr) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      ++def_nr_count;
    }
  });

  Vector<int> indices;
  indices.reserve(def_nr_count);
  Array<float> weights(def_nr_count);

  selection.foreach_index(GrainSize(1024), [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    for (const MDeformWeight &dw : dweights) {
      if (!dw.def_nr == def_nr) {
        continue;
      }
      if (weight_threshold && dw.weight < *weight_threshold) {
        continue;
      }

      indices.append_unchecked(index);
      weights[indices.size() - 1] = dw.weight;
    }
  });

  if (discard_empty_groups && entries.is_empty()) {
    continue;
  }

  ArmatureDeformGroup group;
  const int def_nr = def_nr_range[def_nr_i];
  IndexMask mask = IndexMask::from_indices(all_indices.as_span().slice(entries), memory);
  Array<float> weights = all_weights.as_span().slice(entries);
  group_masks.add_new(def_nr, {std::move(mask), std::move(weights)});

  return group_masks;
}

void BKE_armature_deform_vectors(
    const Object &ob_arm,
    const blender::float4x4 &target_to_world,
    const bool use_envelope,
    const bool use_quaternion,
    const blender::IndexMask &selection,
    std::optional<ArmatureDeformVertexGroupParams> vertex_group_params,
    std::optional<blender::Span<float>> weights,
    blender::MutableSpan<blender::float3> vectors)
{
}

}  // namespace blender::bke

/** \} */
