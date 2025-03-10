/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spoutliner
 */

#pragma once

#include "tree_element.hh"

struct bArmature;
struct BoneCollection;

namespace blender::ed::outliner {

class TreeElementBoneCollectionBase final : public AbstractTreeElement {
  bArmature &armature_;

 public:
  TreeElementBoneCollectionBase(TreeElement &legacy_te, bArmature &armature);
  void expand(SpaceOutliner & /*soops*/) const override;
};

class TreeElementBoneCollection final : public AbstractTreeElement {
  bArmature &armature_;
  BoneCollection &bcoll_;

 public:
  TreeElementBoneCollection(TreeElement &legacy_te, bArmature &armature, BoneCollection &bcoll);
  void expand(SpaceOutliner & /*soops*/) const override;
};

}  // namespace blender::ed::outliner
