# SPDX-FileCopyrightText: 2025 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

set(INC
  PUBLIC .
)

set(INC_SYS
)

set(SRC
  mixbox.cpp
  mixbox.h
)

set(LIB
)


blender_add_lib(extern_mixbox "${SRC}" "${INC}" "" "${LIB}")
add_library(bf::extern::mixbox ALIAS extern_mixbox)
