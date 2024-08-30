/* SPDX-License-Identifier: GPL-2.0-or-later
* Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * Structures to identify unique data blocks.
 * The keys are unique so we are able to match ids across frame updates.
 */

#pragma once

#include "DEG_depsgraph_query.hh"

// Copied from eevee_sync.hh

namespace blender::npr::strokegen
{
  /* -------------------------------------------------------------------- */
  /** \name ObjectKey
   *
   * Unique key to identify each object in the hash-map.
   * Note that we get a unique key for each object component.
   * \{ */
  class ObjectKey {
    /** Hash value of the key. */
    uint64_t hash_value_ = 0;
    /** Original Object or source object for duplis. */
    Object *ob_ = nullptr;
    /** Original Parent object for duplis. */
    Object *parent_ = nullptr;
    /** Dupli objects recursive unique identifier */
    int id_[MAX_DUPLI_RECUR];
    /** Used for particle system hair. */
    int sub_key_ = 0;

   public:
    ObjectKey() = default;

    ObjectKey(Object *ob, int sub_key = 0)
    {
      ob_ = DEG_get_original_object(ob);
      hash_value_ = BLI_ghashutil_ptrhash(ob_);

      if (DupliObject *dupli = DRW_object_get_dupli(ob)) {
        parent_ = DRW_object_get_dupli_parent(ob);
        hash_value_ = BLI_ghashutil_combine_hash(hash_value_, BLI_ghashutil_ptrhash(parent_));
        for (int i : IndexRange(MAX_DUPLI_RECUR)) {
          id_[i] = dupli->persistent_id[i];
          if (id_[i] == INT_MAX) {
            break;
          }
          hash_value_ = BLI_ghashutil_combine_hash(hash_value_, BLI_ghashutil_inthash(id_[i]));
        }
      }

      if (sub_key != 0) {
        sub_key_ = sub_key;
        hash_value_ = BLI_ghashutil_combine_hash(hash_value_, BLI_ghashutil_inthash(sub_key_));
      }
    }

    uint64_t hash() const
    {
      return hash_value_;
    }

    bool operator<(const ObjectKey &k) const
    {
      if (hash_value_ != k.hash_value_) {
        return hash_value_ < k.hash_value_;
      }
      if (ob_ != k.ob_) {
        return (ob_ < k.ob_);
      }
      if (parent_ != k.parent_) {
        return (parent_ < k.parent_);
      }
      if (sub_key_ != k.sub_key_) {
        return (sub_key_ < k.sub_key_);
      }
      if (parent_) {
        for (int i : IndexRange(MAX_DUPLI_RECUR)) {
          if (id_[i] < k.id_[i]) {
            return true;
          }
          if (id_[i] == INT_MAX) {
            break;
          }
        }
      }
      return false;
    }

    bool operator==(const ObjectKey &k) const
    {
      if (hash_value_ != k.hash_value_) {
        return false;
      }
      if (ob_ != k.ob_) {
        return false;
      }
      if (parent_ != k.parent_) {
        return false;
      }
      if (sub_key_ != k.sub_key_) {
        return false;
      }
      if (parent_) {
        for (int i : IndexRange(MAX_DUPLI_RECUR)) {
          if (id_[i] != k.id_[i]) {
            return false;
          }
          if (id_[i] == INT_MAX) {
            break;
          }
        }
      }
      return true;
    }
  };
  /** \} */



  /* -------------------------------------------------------------------- */
  /** \name Engine-Specific Draw-Data
   * \{ */
  struct BaseHandle {
    unsigned int recalc;
  };

  struct ObjectHandle : BaseHandle {
    ObjectKey object_key;
  };

  struct BnprDrawData : public DrawData
  {
    ObjectKey object_key;

    void reset_recalc_flag()
    {
      if (recalc != 0)
      {
        recalc = 0;
      }
    }
  };

  struct WorldHandle : public DrawData
  {
    void reset_recalc_flag()
    {
      if (recalc != 0)
      {
        recalc = 0;
      }
    }
  };

  struct SceneHandle : public DrawData
  {
    void reset_recalc_flag()
    {
      if (recalc != 0)
      {
        recalc = 0;
      }
    }
  };
  /** \} */
}

/** \} */
