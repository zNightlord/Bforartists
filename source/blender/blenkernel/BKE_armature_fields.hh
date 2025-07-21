/* SPDX-FileCopyrightText: 2001-2002 NaN Holding BV. All rights reserved.
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bke
 */

#include "BLI_listbase_wrapper.hh"
#include "BLI_math_matrix_types.hh"

#include "BKE_armature.hh"
#include "BKE_curves.hh"
#include "BKE_geometry_fields.hh"
#include "BKE_grease_pencil.hh"
#include "BKE_mesh.hh"

#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

struct bArmature;
struct Bone;
struct bPoseChannel;

namespace blender::bke {

namespace armature_deform {

struct VertexGroupData {
  ConstListBaseWrapper<bDeformGroup> vertex_groups;
  blender::Span<MDeformVert> dverts;
};

inline std::optional<VertexGroupData> vertex_groups_from_field_context(
    const bke::GeometryFieldContext &context)
{
  switch (context.type()) {
    case GeometryComponent::Type::Mesh:
      if (const Mesh *mesh = context.mesh()) {
        return VertexGroupData{{mesh->vertex_group_names}, mesh->deform_verts()};
      }
      break;
    case GeometryComponent::Type::GreasePencil: {
      const GreasePencil *grease_pencil = context.grease_pencil();
      const bke::greasepencil::Drawing *drawing = context.grease_pencil_layer_drawing();
      if (grease_pencil && drawing) {
        return VertexGroupData{{grease_pencil->vertex_group_names},
                               drawing->geometry.wrap().deform_verts()};
      }
      break;
    }
    default:
      break;
  }
  return std::nullopt;
}

inline std::optional<int> find_vertex_group_def_nr(
    const ConstListBaseWrapper<bDeformGroup> &vertex_groups, StringRef name)
{
  int def_nr = 0;
  for (const bDeformGroup *def_grp : vertex_groups) {
    if (def_grp->name == name) {
      return def_nr;
    }
    ++def_nr;
  }
  return std::nullopt;
}

inline void weights_from_deform_verts(const IndexMask &mask,
                                      const Span<MDeformVert> dverts,
                                      const int def_nr,
                                      MutableSpan<float> weights)
{
  constexpr GrainSize grain_size(1024);

  BLI_assert(dverts.size() >= mask.min_array_size());

  mask.foreach_index(grain_size, [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    weights[index] = 0.0f;
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr == def_nr) {
        weights[index] = dw.weight;
        break;
      }
    }
  });
}

inline void selection_from_deform_verts(const IndexMask &mask,
                                        const Span<MDeformVert> dverts,
                                        const int def_nr,
                                        const float threshold,
                                        MutableSpan<bool> selection)
{
  constexpr GrainSize grain_size(1024);

  BLI_assert(dverts.size() >= mask.min_array_size());

  mask.foreach_index(grain_size, [&](const int index) {
    const MDeformVert &dvert = dverts[index];
    const Span<MDeformWeight> dweights(dvert.dw, dvert.totweight);
    selection[index] = false;
    for (const MDeformWeight &dw : dweights) {
      if (dw.def_nr == def_nr && dw.weight > threshold) {
        selection[index] = true;
        break;
      }
    }
  });
}

/* XXX Getting vertex groups in field evaluation isn't possible because it relies on special vertex
 * group data arrays (MDeformVert) and the declared list of vertex group names (bDeformGroup) to
 * map to the pose channels of the armature.
 *
 * The only feasible way to get vertex group information is in fields, where the field context can
 * provide the necessary access to both the `dverts` array and the `vertex_group_names` list.
 *
 * Future implementation of sparse attributes should solve these issues. For now this should be
 * considered a proof-of-concept. */

class VertexGroupWeightInput final : public bke::GeometryFieldInput {
 private:
  const Bone *bone_;

 public:
  VertexGroupWeightInput(const Bone &bone)
      : bke::GeometryFieldInput(CPPType::get<float>(), "Vertex Group Weight Input"), bone_(&bone)
  {
  }

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const
  {
    const auto default_weights = VArray<float>::from_single(0.0f, mask.min_array_size());
    std::optional<VertexGroupData> vgroup_data = vertex_groups_from_field_context(context);
    if (!vgroup_data) {
      return default_weights;
    }
    const std::optional<int> def_nr = find_vertex_group_def_nr(vgroup_data->vertex_groups,
                                                               bone_->name);
    if (!def_nr) {
      return default_weights;
    }

    Array<float> weights(mask.min_array_size());
    weights_from_deform_verts(mask, vgroup_data->dverts, *def_nr, weights);
    return VArray<float>::from_container(std::move(weights));
  }

  uint64_t hash() const
  {
    return get_default_hash(bone_);
  }

  bool is_equal_to(const fn::FieldNode &other) const
  {
    if (const auto *other_field = dynamic_cast<const VertexGroupWeightInput *>(&other)) {
      return other_field->bone_ == bone_;
    }
    return false;
  }

  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
  {
    return AttrDomain::Point;
  }
};

class VertexGroupSelectionInput final : public bke::GeometryFieldInput {
 private:
  const Bone *bone_;
  float threshold_;

 public:
  VertexGroupSelectionInput(const Bone &bone, const float threshold)
      : bke::GeometryFieldInput(CPPType::get<bool>(), "Vertex Group Selection Input"),
        bone_(&bone),
        threshold_(threshold)
  {
  }

  GVArray get_varray_for_context(const bke::GeometryFieldContext &context,
                                 const IndexMask &mask) const
  {
    const auto default_selection = VArray<bool>::from_single(false, mask.min_array_size());
    std::optional<VertexGroupData> vgroup_data = vertex_groups_from_field_context(context);
    if (!vgroup_data) {
      return default_selection;
    }
    const std::optional<int> def_nr = find_vertex_group_def_nr(vgroup_data->vertex_groups,
                                                               bone_->name);
    if (!def_nr) {
      return default_selection;
    }

    Array<bool> selection(mask.min_array_size());
    selection_from_deform_verts(mask, vgroup_data->dverts, *def_nr, threshold_, selection);
    return VArray<bool>::from_container(std::move(selection));
  }

  uint64_t hash() const
  {
    return get_default_hash(bone_, threshold_);
  }

  bool is_equal_to(const fn::FieldNode &other) const
  {
    if (const auto *other_field = dynamic_cast<const VertexGroupSelectionInput *>(&other)) {
      return other_field->bone_ == bone_ && other_field->threshold_ == threshold_;
    }
    return false;
  }

  std::optional<AttrDomain> preferred_domain(const GeometryComponent & /*component*/) const
  {
    return AttrDomain::Point;
  }
};

class BoneEnvelopeMultiFunction : public mf::MultiFunction {
 private:
  float4x4 target_to_armature_;
  const Bone *bone_;
  float threshold_;

 public:
  BoneEnvelopeMultiFunction(const float4x4 &target_to_armature,
                            const Bone &bone,
                            const float threshold)
      : target_to_armature_(target_to_armature), bone_(&bone), threshold_(threshold)
  {
    static mf::Signature signature_;
    mf::SignatureBuilder builder{"Bone Envelope", signature_};
    builder.single_input<float3>("Position");
    builder.single_output<float>("Weight");
    builder.single_output<bool>("Selection");
    this->set_signature(&signature_);
  }

  void call(const IndexMask &mask, mf::Params params, mf::Context /*context*/) const override
  {
    const VArraySpan<float3> positions = params.readonly_single_input<float3>(0, "Position");
    MutableSpan<float> weights = params.uninitialized_single_output<float>(1, "Weight");
    MutableSpan<bool> selection = params.uninitialized_single_output<bool>(2, "Selection");
    mask.foreach_index([&](const int64_t i) {
      const float3 &pos = positions[i];
      const float3 arm_pos = math::transform_point(target_to_armature_, pos);
      weights[i] = distfactor_to_bone(
          arm_pos, bone_->head, bone_->tail, bone_->rad_head, bone_->rad_tail, bone_->dist);
      selection[i] = weights[i] > threshold_;
    });
  }
};

}  // namespace armature_deform

}  // namespace blender::bke
