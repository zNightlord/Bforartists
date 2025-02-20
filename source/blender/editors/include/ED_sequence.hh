/* SPDX-FileCopyrightText: 2025 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup editors
 */

#pragma once

#include "BKE_sequence.hh"

enum class AddSequenceMode : int8_t { Blank = 0, CopySettings, FullCopy };

void ED_operatortypes_sequence();
