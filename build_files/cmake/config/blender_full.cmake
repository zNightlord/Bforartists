# SPDX-FileCopyrightText: 2014-2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

# Custom config reproducing the printed CMake configuration for the
# Bforartists Android cross-build.
#
# Usage:
#   make android BUILD_CMAKE_ARGS="-C<path>/blender_android_full.cmake"
# or directly:
#   cmake -C<path>/blender_android_full.cmake ...

# ---------------------------------------------------------------------------
# Build Options
set(WITH_ALEMBIC             ON  CACHE BOOL "" FORCE)
set(WITH_BULLET               ON  CACHE BOOL "" FORCE)
set(WITH_CLANG                OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES               ON  CACHE BOOL "" FORCE)
set(WITH_FFTW3                ON  CACHE BOOL "" FORCE)
set(WITH_FREESTYLE            ON  CACHE BOOL "" FORCE)
set(WITH_GMP                  ON  CACHE BOOL "" FORCE)
set(WITH_HARU                 ON  CACHE BOOL "" FORCE)
set(WITH_IK_ITASC             ON  CACHE BOOL "" FORCE)
set(WITH_IK_SOLVER             ON  CACHE BOOL "" FORCE)
set(WITH_INPUT_NDOF           ON  CACHE BOOL "" FORCE)
set(WITH_INPUT_IME            ON  CACHE BOOL "" FORCE)
set(WITH_INTERNATIONAL        ON  CACHE BOOL "" FORCE)
set(WITH_MANIFOLD             ON  CACHE BOOL "" FORCE)
set(WITH_OPENIMAGEDENOISE     ON  CACHE BOOL "" FORCE)
set(WITH_OPENSUBDIV           ON  CACHE BOOL "" FORCE)
set(WITH_OPENVDB              ON  CACHE BOOL "" FORCE)
set(WITH_POTRACE              ON  CACHE BOOL "" FORCE)
set(WITH_PUGIXML              ON  CACHE BOOL "" FORCE)
set(WITH_QUADRIFLOW           ON  CACHE BOOL "" FORCE)
set(WITH_TBB                  ON  CACHE BOOL "" FORCE)
set(WITH_USD                  ON  CACHE BOOL "" FORCE)
set(WITH_MATERIALX            ON  CACHE BOOL "" FORCE)
set(WITH_XR_OPENXR            ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Compiler Options
set(WITH_BUILDINFO            ON  CACHE BOOL "" FORCE)
set(WITH_OPTIMIZED_BUILD_TOOLS ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# System Options
# WITH_INSTALL_PORTABLE was blank/unset in the printed config - left as
# platform default rather than forced here.
set(WITH_TBB_MALLOC_PROXY     ON  CACHE BOOL "" FORCE)
set(WITH_MEM_VALGRIND         OFF CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# GHOST Options
set(WITH_GHOST_DEBUG          OFF CACHE BOOL "" FORCE)
set(WITH_GHOST_SDL            OFF CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Profiling Options
set(WITH_TRACY                OFF CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Image Formats
set(WITH_IMAGE_CINEON         ON  CACHE BOOL "" FORCE)
set(WITH_IMAGE_OPENJPEG       ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Audio
set(WITH_AUDASPACE            ON  CACHE BOOL "" FORCE)
set(WITH_RUBBERBAND           ON  CACHE BOOL "" FORCE)
set(WITH_CODEC_FFMPEG         ON  CACHE BOOL "" FORCE)
set(WITH_CODEC_SNDFILE        ON  CACHE BOOL "" FORCE)
set(WITH_COREAUDIO            OFF CACHE BOOL "" FORCE)
set(WITH_JACK                 OFF CACHE BOOL "" FORCE)
set(WITH_OPENAL               ON  CACHE BOOL "" FORCE)
set(WITH_PULSEAUDIO           OFF CACHE BOOL "" FORCE)
set(WITH_SDL_AUDIO            OFF CACHE BOOL "" FORCE)
set(WITH_WASAPI               ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Python
set(WITH_PYTHON_INSTALL              ON  CACHE BOOL "" FORCE)
set(WITH_PYTHON_INSTALL_NUMPY        ON  CACHE BOOL "" FORCE)
set(WITH_PYTHON_INSTALL_ZSTANDARD    ON  CACHE BOOL "" FORCE)
set(WITH_PYTHON_MODULE               OFF CACHE BOOL "" FORCE)
set(WITH_PYTHON_SAFETY               OFF CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Modifiers
set(WITH_MOD_FLUID            ON  CACHE BOOL "" FORCE)
set(WITH_MOD_OCEANSIM         ON  CACHE BOOL "" FORCE)
set(WITH_MOD_REMESH           ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Rendering
set(WITH_HYDRA                ON  CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Rendering (Cycles)
set(WITH_CYCLES_OSL           ON  CACHE BOOL "" FORCE)
set(WITH_CYCLES_EMBREE        ON  CACHE BOOL "" FORCE)
set(WITH_CYCLES_PATH_GUIDING  ON  CACHE BOOL "" FORCE)
# Cycles CPU-only: every GPU compute backend forced OFF. None of
# CUDA/OptiX/HIP/HIPRT/oneAPI have an Android ARM equivalent anyway, and
# this file is meant to render exclusively on-CPU.
set(WITH_CYCLES_DEVICE_OPTIX  OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_DEVICE_CUDA   OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_CUDA_BINARIES OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_DEVICE_ONEAPI OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_ONEAPI_BINARIES OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_DEVICE_HIP    OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_HIP_BINARIES  OFF CACHE BOOL "" FORCE)
set(WITH_CYCLES_DEVICE_HIPRT  OFF CACHE BOOL "" FORCE)

# ---------------------------------------------------------------------------
# Tests
# All disabled: WITH_GTESTS is off (avoids the TEST_PYTHON_EXE configure
# error hit earlier), so anything gated behind it below is a no-op anyway.
set(WITH_GTESTS                       OFF CACHE BOOL "" FORCE)
set(CYCLES_TEST_DEVICES               "CPU" CACHE STRING "" FORCE)
set(WITH_CYCLES_TEST_OSL              OFF CACHE BOOL "" FORCE)
set(WITH_GPU_RENDER_TESTS             OFF CACHE BOOL "" FORCE)
set(WITH_GPU_BACKEND_TESTS            OFF CACHE BOOL "" FORCE)
set(WITH_GPU_DRAW_TESTS               OFF CACHE BOOL "" FORCE)
set(WITH_GPU_COMPOSITOR_TESTS         OFF CACHE BOOL "" FORCE)
set(WITH_GPU_MESH_PAINT_TESTS         OFF CACHE BOOL "" FORCE)
set(WITH_UI_TESTS                     OFF CACHE BOOL "" FORCE)
set(WITH_UI_TESTS_HEADLESS            OFF CACHE BOOL "" FORCE)
# Was ON in the printed config but meaningless with WITH_GTESTS OFF above.
set(WITH_UI_TESTS_VULKAN              OFF CACHE BOOL "" FORCE)
set(WITH_GPU_RENDER_TESTS_VULKAN      OFF CACHE BOOL "" FORCE)
set(WITH_LINUX_OFFICIAL_RELEASE_TESTS OFF CACHE BOOL "" FORCE)
set(WITH_IO_PLY              ON  CACHE BOOL "" FORCE)
set(WITH_IO_STL              ON  CACHE BOOL "" FORCE)
set(WITH_IO_WAVEFRONT_OBJ    ON  CACHE BOOL "" FORCE)
set(WITH_LIBMV               ON  CACHE BOOL "" FORCE)
set(WITH_LIBMV_SCHUR_SPECIALIZATIONS ON CACHE BOOL "" FORCE)
set(WITH_MANIFOLD            ON  CACHE BOOL "" FORCE)
set(WITH_MOD_FLUID           ON  CACHE BOOL "" FORCE)
set(WITH_MOD_OCEANSIM        ON  CACHE BOOL "" FORCE)
set(WITH_MOD_REMESH          ON  CACHE BOOL "" FORCE)
set(WITH_UV_SLIM             ON  CACHE BOOL "" FORCE)
set(WITH_NANOVDB             ON  CACHE BOOL "" FORCE)
set(WITH_OPENAL              ON  CACHE BOOL "" FORCE)
set(WITH_OPENIMAGEDENOISE    ON  CACHE BOOL "" FORCE)
set(WITH_OPENSUBDIV          ON  CACHE BOOL "" FORCE)
set(WITH_OPENVDB             ON  CACHE BOOL "" FORCE)
set(WITH_OPENVDB_BLOSC       ON  CACHE BOOL "" FORCE)
set(WITH_POTRACE             ON  CACHE BOOL "" FORCE)
set(WITH_PUGIXML             ON  CACHE BOOL "" FORCE)
set(WITH_PYTHON_INSTALL      ON  CACHE BOOL "" FORCE)
set(WITH_QUADRIFLOW          ON  CACHE BOOL "" FORCE)
set(WITH_RUBBERBAND          ON  CACHE BOOL "" FORCE)
set(WITH_SDL_AUDIO           OFF CACHE BOOL "" FORCE)
set(WITH_TBB                 ON  CACHE BOOL "" FORCE)
set(WITH_USD                 ON  CACHE BOOL "" FORCE)
set(WITH_MATERIALX           ON  CACHE BOOL "" FORCE)
set(WITH_HYDRA               ON  CACHE BOOL "" FORCE)
set(WITH_XR_OPENXR           ON  CACHE BOOL "" FORCE)

# platform dependent options
if(APPLE)
  set(WITH_COREAUDIO           ON  CACHE BOOL "" FORCE)
  set(WITH_CYCLES_DEVICE_METAL ON  CACHE BOOL "" FORCE)
  set(WITH_BLENDER_THUMBNAILER ON  CACHE BOOL "" FORCE)
endif()
if(WIN32)
  set(WITH_WASAPI              ON  CACHE BOOL "" FORCE)
endif()
if(UNIX AND NOT (APPLE OR ANDROID))
  set(WITH_JACK                ON  CACHE BOOL "" FORCE)
  set(WITH_DOC_MANPAGE         ON  CACHE BOOL "" FORCE)
  set(WITH_GHOST_XDND          ON  CACHE BOOL "" FORCE)
  set(WITH_PULSEAUDIO          ON  CACHE BOOL "" FORCE)
  set(WITH_X11_XINPUT          ON  CACHE BOOL "" FORCE)
endif()
if(NOT (APPLE OR ANDROID))
  # Can't use CMAKE_SYSTEM_PROCESSOR here as it's not set yet,
  # so fall back to checking the env for vcvarsall's VSCMD_ARG_TGT_ARCH
  if(NOT (WIN32 AND "$ENV{VSCMD_ARG_TGT_ARCH}" STREQUAL "arm64"))
    set(WITH_TBB_MALLOC_PROXY       ON  CACHE BOOL "" FORCE)
  endif()
endif()
