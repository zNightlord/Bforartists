/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bke
 * \brief Low-level operations for the 'Sequence' data-block.
 */

#include "DNA_sequence_types.h"

struct Main;

Sequence *BKE_sequence_add(Main &bmain, const char *name);
bool BKE_sequence_can_be_removed(const Main &bmain, const Sequence &sequence);
