# File: CMakeLists.txt
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
#
# Fortran 2018 production core (ARCHITECTURE.md `src/`). Build + test:
#
#   cmake -S . -B build && cmake --build build && ctest --test-dir build --output-on-failure
#
cmake_minimum_required(VERSION 3.20)
# VERSION is a bare numeric triple matching the v0.1.0 release series.
project(CableDynCore VERSION 0.1.0 LANGUAGES C Fortran
        DESCRIPTION "CableDyn geometrically nonlinear cable/mooring solver — Fortran core")

option(CABLEDYN_ENABLE_COVERAGE "Instrument the GNU Fortran core for gcov/gcovr source coverage" OFF)
option(CABLEDYN_ENABLE_OPENMP "Enable OpenMP acceleration when the compiler provides it" ON)
option(CABLEDYN_STATIC_WINDOWS_RELEASE
       "Build a single-file non-OpenMP Windows release with Intel IFX/MSVC" OFF)

if(CABLEDYN_STATIC_WINDOWS_RELEASE)
  if(NOT WIN32 OR NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR NOT MSVC)
    message(FATAL_ERROR
      "CABLEDYN_STATIC_WINDOWS_RELEASE requires Windows, Intel IFX, and the MSVC C compiler")
  endif()
  # Match OpenFAST's redistributable-free Release configuration: static IFX,
  # static MSVC CRT, and no OpenMP runtime dependency (IFX has no static OpenMP
  # runtime on Windows).
  set(CABLEDYN_ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded" CACHE STRING "" FORCE)
  # Large finite-EI meshes create compiler-generated array temporaries in deeply
  # nested residual/assembly calls. Keep those temporaries off the process stack;
  # the released executable also receives a generous reserve below. The size
  # threshold keeps run-time-sized arrays and temporaries on the heap but leaves
  # compile-time-sized locals up to 1 MiB on the stack: with a bare /heap-arrays,
  # IFX 2025.3 heap-allocates the host frame of every procedure with internal
  # procedures and never frees it, leaking memory on each call.
  #
  # /Qipo: the element, hydrodynamic, and assembly kernels live in separate modules
  # and are called per element and per Newton iteration; without interprocedural
  # optimisation IFX cannot inline across those module boundaries. The core archive
  # then holds IPO objects, so it links only through an IPO-capable IFX link (the
  # driver below, or the OpenFAST Release|x64 project, which enables whole-program
  # optimisation).
  string(APPEND CMAKE_Fortran_FLAGS " /libs:static /heap-arrays:1024 /Qipo")
  # CMake does not pass the compile flags to an IFX executable link on Windows.
  string(APPEND CMAKE_EXE_LINKER_FLAGS " /Qipo")
  # LAPACK/BLAS: the netlib reference implementation compiled with this IFX
  # (release/build_static_windows.ps1 builds it), not MKL. Every CableDyn system is
  # banded with a small bandwidth, and MKL's sequential DGBTRF/DGBTRS/DTBSV spend
  # most of their time in per-column DGER/DAXPY dispatch on vectors of a few
  # entries: with MKL, 30-65% of the run time of the release benchmark decks was
  # spent inside MKL level-1/2 kernels, against 6-15% with the reference routines.
  set(CABLEDYN_STATIC_LAPACK_LIBRARIES "" CACHE STRING
      "Static reference LAPACK and BLAS libraries (in that order) built with the release IFX")
  if(NOT CABLEDYN_STATIC_LAPACK_LIBRARIES)
    message(FATAL_ERROR
      "CABLEDYN_STATIC_WINDOWS_RELEASE needs -DCABLEDYN_STATIC_LAPACK_LIBRARIES=<liblapack.lib>;<libblas.lib> "
      "(reference LAPACK built with IFX; release/build_static_windows.ps1 provides them)")
  endif()
  foreach(_cabledyn_static_lib IN LISTS CABLEDYN_STATIC_LAPACK_LIBRARIES)
    if(NOT EXISTS "${_cabledyn_static_lib}")
      message(FATAL_ERROR "CABLEDYN_STATIC_LAPACK_LIBRARIES entry not found: ${_cabledyn_static_lib}")
    endif()
  endforeach()
endif()

# src/cabledyn_mutex.c uses C11 <stdatomic.h>, which the MSVC C compiler provides only
# behind this flag, in every MSVC configuration (static release or not).
if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  string(APPEND CMAKE_C_FLAGS " /experimental:c11atomics")
endif()

# Default to an optimised build unless the caller (or a multi-config generator) sets one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  if(CABLEDYN_ENABLE_COVERAGE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
  else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

# The Debug tier is a hardening tier, not just an unoptimised variant: alongside the
# non-short-circuit .AND./.OR. behaviour the -O0 CI cell already exposes, enable runtime
# bounds / pointer / array-temporary checking so a latent out-of-bounds access (e.g. a
# substring formed before the guard that should protect it) fails loudly here rather than
# reading stray memory under Release. Release builds are unaffected.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=all -fbacktrace")
endif()

# Emit/consume all .mod files from one directory so targets see each other's modules.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

# Strict, standard-conforming flags for gfortran (ifx/ifort can be added later).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-std=f2018>
                      $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
                      $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
                      $<$<COMPILE_LANGUAGE:Fortran>:-fimplicit-none>
                      $<$<COMPILE_LANGUAGE:C>:-Wall>
                      $<$<COMPILE_LANGUAGE:C>:-Wextra>)
endif()
# GNU builds on Windows (MinGW-w64): every executable gets a 64 MiB stack reserve instead of
# the 2 MiB MinGW default. With OpenMP, gfortran places fixed-size local arrays on the stack
# (-fopenmp implies -frecursive), and the stack depth of the OpenBLAS kernels varies with the
# CPU, so a small reserve can overflow on one machine and pass on another. A reserve is address
# space, not committed memory. The library itself must run within a small stack, since a host
# such as python.exe sets the stack of the threads that call it: the stack_* tests check that.
# The static IFX release sets its own reserve on the driver.
set(CABLEDYN_GNU_WINDOWS_STACK OFF)
if(MINGW AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CABLEDYN_GNU_WINDOWS_STACK ON)
  set(CABLEDYN_EXE_STACK_RESERVE 67108864)
  set(CABLEDYN_SMALL_STACK_RESERVE 1048576)
  add_link_options(
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:-Wl,--stack,${CABLEDYN_EXE_STACK_RESERVE}>")
endif()

if(CABLEDYN_ENABLE_COVERAGE AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "CABLEDYN_ENABLE_COVERAGE requires GNU Fortran/gcov instrumentation")
endif()

if(CABLEDYN_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS Fortran)
endif()
if(OpenMP_Fortran_FOUND)
  message(STATUS "CableDyn: OpenMP Fortran enabled for element-parallel assembly")
else()
  message(STATUS "CableDyn: OpenMP Fortran not found; element-parallel assembly directives will be ignored")
endif()
# Linear-algebra backend for the static/dynamic solves (DGBSV wrappers for both
# the EI=0 and finite-EI paths). Prefer a system LAPACK (apt liblapack on CI);
# fall back to a conda OpenBLAS (which bundles BLAS+LAPACK) located via CONDA_PREFIX
# for local Windows builds where find_package(LAPACK) cannot locate it unaided. The
# static Windows release links the IFX-built reference LAPACK given above.
if(CABLEDYN_STATIC_WINDOWS_RELEASE)
  set(CABLEDYN_LAPACK ${CABLEDYN_STATIC_LAPACK_LIBRARIES})
  message(STATUS "CableDyn: static release LAPACK: ${CABLEDYN_STATIC_LAPACK_LIBRARIES}")
else()
  find_package(LAPACK QUIET)
endif()
if(CABLEDYN_STATIC_WINDOWS_RELEASE)
  # CABLEDYN_LAPACK is set above.
elseif(LAPACK_FOUND)
  set(CABLEDYN_LAPACK LAPACK::LAPACK)
else()
  find_library(CABLEDYN_OPENBLAS NAMES openblas
    HINTS ENV CONDA_PREFIX ENV CMAKE_PREFIX_PATH
    PATH_SUFFIXES Library/lib lib)
  if(NOT CABLEDYN_OPENBLAS)
    message(FATAL_ERROR
      "No LAPACK found and no OpenBLAS in CONDA_PREFIX/CMAKE_PREFIX_PATH. "
      "Install liblapack-dev (Linux) or activate the conda env (Windows).")
  endif()
  set(CABLEDYN_LAPACK ${CABLEDYN_OPENBLAS})
  message(STATUS "CableDyn: using OpenBLAS fallback for LAPACK: ${CABLEDYN_OPENBLAS}")
endif()
# A threaded OpenBLAS commits one work buffer per worker thread (about 4 GB on a 32-CPU
# machine) although every CableDyn solve is small; src/cabledyn_blas.c limits it to one
# thread at the first CableDyn_Create when the linked back end provides the setter.
include(CheckFunctionExists)
if(CABLEDYN_STATIC_WINDOWS_RELEASE)
  # Reference LAPACK has no thread pool.
  set(CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS FALSE)
else()
  if(LAPACK_FOUND)
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  else()
    set(CMAKE_REQUIRED_LIBRARIES ${CABLEDYN_OPENBLAS})
  endif()
  check_function_exists(openblas_set_num_threads CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS)
  unset(CMAKE_REQUIRED_LIBRARIES)
endif()

set(CABLEDYN_RUNTIME_DLLS)
if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Collect the runtime bin dirs to search, independent of HOW the toolchain was resolved.
  # OpenBLAS and the GNU compiler runtimes do not necessarily live in the same directory: in
  # the conda gfortran toolchain openblas.dll is under Library/bin, while libgfortran/
  # libquadmath/libgcc_s/libgomp/libwinpthread are under Library/mingw-w64/bin. Gather prefix
  # roots from every place LAPACK/OpenBLAS could have been found -- the OpenBLAS fallback lib,
  # an activated conda env (CONDA_PREFIX), an explicit -DCMAKE_PREFIX_PATH configure, and the
  # libraries find_package(LAPACK/BLAS) actually resolved -- then expand each root into the bin
  # dirs the DLLs may sit in and stage each DLL from the first directory that has it, so
  # build/runtime is self-contained.
  set(_cabledyn_dll_roots)
  if(CABLEDYN_OPENBLAS)
    get_filename_component(_cabledyn_lib_dir "${CABLEDYN_OPENBLAS}" DIRECTORY)
    get_filename_component(_cabledyn_root "${_cabledyn_lib_dir}" DIRECTORY)
    list(APPEND _cabledyn_dll_roots "${_cabledyn_root}")
  endif()
  if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND _cabledyn_dll_roots "$ENV{CONDA_PREFIX}")
  endif()
  foreach(_cabledyn_prefix IN LISTS CMAKE_PREFIX_PATH)
    list(APPEND _cabledyn_dll_roots "${_cabledyn_prefix}")
  endforeach()
  foreach(_cabledyn_lib IN LISTS LAPACK_LIBRARIES BLAS_LIBRARIES)
    if(EXISTS "${_cabledyn_lib}")
      get_filename_component(_cabledyn_lib_dir "${_cabledyn_lib}" DIRECTORY)
      get_filename_component(_cabledyn_root "${_cabledyn_lib_dir}" DIRECTORY)
      list(APPEND _cabledyn_dll_roots "${_cabledyn_root}")
    endif()
  endforeach()
  # Environment variables on Windows commonly use backslashes. Convert every
  # search root before it is propagated into install(FILES ...); otherwise the
  # generated cmake_install.cmake can interpret sequences such as \U as escapes.
  set(_cabledyn_normalized_roots)
  foreach(_cabledyn_root IN LISTS _cabledyn_dll_roots)
    file(TO_CMAKE_PATH "${_cabledyn_root}" _cabledyn_root_cmake)
    list(APPEND _cabledyn_normalized_roots "${_cabledyn_root_cmake}")
  endforeach()
  set(_cabledyn_dll_roots ${_cabledyn_normalized_roots})
  # Expand each root into the bin dirs the runtime DLLs may sit in. The four candidates cover a
  # root that is the env (Library/bin, Library/mingw-w64/bin), a root that is a conda Library
  # dir (bin, mingw-w64/bin), and a generic install prefix (bin). Nonexistent dirs are skipped
  # by the per-DLL EXISTS probe below.
  set(_cabledyn_dll_search_dirs)
  foreach(_cabledyn_root IN LISTS _cabledyn_dll_roots)
    list(APPEND _cabledyn_dll_search_dirs
         "${_cabledyn_root}/bin"
         "${_cabledyn_root}/mingw-w64/bin"
         "${_cabledyn_root}/Library/bin"
         "${_cabledyn_root}/Library/mingw-w64/bin")
  endforeach()
  if(_cabledyn_dll_search_dirs)
    list(REMOVE_DUPLICATES _cabledyn_dll_search_dirs)
  endif()
  foreach(_cabledyn_dll IN ITEMS
      openblas.dll
      libgomp-1.dll
      libgfortran-5.dll
      libquadmath-0.dll
      libgcc_s_seh-1.dll
      libwinpthread-1.dll)   # libgfortran-5.dll depends on it; without it the bundle is not self-contained
    set(_cabledyn_dll_found "")
    foreach(_cabledyn_dir IN LISTS _cabledyn_dll_search_dirs)
      if(EXISTS "${_cabledyn_dir}/${_cabledyn_dll}")
        set(_cabledyn_dll_found "${_cabledyn_dir}/${_cabledyn_dll}")
        break()
      endif()
    endforeach()
    if(_cabledyn_dll_found)
      list(APPEND CABLEDYN_RUNTIME_DLLS "${_cabledyn_dll_found}")
    elseif(_cabledyn_dll_search_dirs)
      # Surface an incomplete bundle rather than silently shipping a runtime/ that fails
      # outside conda with a missing GNU runtime DLL.
      message(WARNING "CableDyn: runtime DLL ${_cabledyn_dll} not found in any of "
                      "${_cabledyn_dll_search_dirs}; build/runtime may not be self-contained.")
    endif()
  endforeach()
endif()

# Run-time loading of OpenBLAS (Windows GNU builds that stage openblas.dll). A threaded
# OpenBLAS sizes its worker pool from OPENBLAS_NUM_THREADS while its DLL loads, and each
# worker commits a ~128 MB buffer at once (about 4 GB on a 32-CPU machine); Windows also
# maps delay-load imports in the background right after process start. With
# CABLEDYN_BLAS_RUNTIME_LOAD, src/cabledyn_blas.c defines the LAPACK entry points CableDyn
# calls and loads openblas.dll itself on first use with a one-thread pool, so neither the
# shared library nor any executable linking cabledyn_core imports it.
set(CABLEDYN_BLAS_RUNTIME_LOAD OFF)
if(WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS)
  foreach(_cabledyn_dll_src IN LISTS CABLEDYN_RUNTIME_DLLS)
    get_filename_component(_cabledyn_dll_name "${_cabledyn_dll_src}" NAME)
    if(_cabledyn_dll_name STREQUAL "openblas.dll")
      set(CABLEDYN_BLAS_RUNTIME_LOAD ON)
    endif()
  endforeach()
  if(CABLEDYN_BLAS_RUNTIME_LOAD)
    message(STATUS "CableDyn: openblas.dll is loaded at run time with a one-thread pool")
  endif()
endif()

# Stage the shared runtime DLL set exactly once, through a single custom command whose copies
# are declared OUTPUTs.  The build graph then serializes it, so parallel POST_BUILD steps cannot
# race on shared files.  Put Windows development/test executables and one DLL set in bin/, keeping
# the build root free of an OpenBLAS DLL that could shadow NumPy's own runtime in Python-driven
# tests.  A second self-contained deployment set remains in runtime/. Each packaged executable
# copies only its OWN exe into runtime/ and depends on this target.
if(CABLEDYN_RUNTIME_DLLS)
  set(_cabledyn_runtime_dir "${CMAKE_BINARY_DIR}/runtime")
  set(_cabledyn_build_bin "${CMAKE_BINARY_DIR}/bin")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${_cabledyn_build_bin}")
  set(_cabledyn_staged_dlls)
  set(_cabledyn_build_dlls)
  foreach(_cabledyn_dll_src IN LISTS CABLEDYN_RUNTIME_DLLS)
    get_filename_component(_cabledyn_dll_name "${_cabledyn_dll_src}" NAME)
    list(APPEND _cabledyn_staged_dlls "${_cabledyn_runtime_dir}/${_cabledyn_dll_name}")
    list(APPEND _cabledyn_build_dlls "${_cabledyn_build_bin}/${_cabledyn_dll_name}")
  endforeach()
  add_custom_command(
    OUTPUT ${_cabledyn_staged_dlls} ${_cabledyn_build_dlls}
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_cabledyn_runtime_dir}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${_cabledyn_build_bin}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CABLEDYN_RUNTIME_DLLS} "${_cabledyn_runtime_dir}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CABLEDYN_RUNTIME_DLLS} "${_cabledyn_build_bin}"
    DEPENDS ${CABLEDYN_RUNTIME_DLLS}
    COMMENT "Staging Windows runtime DLLs into ${_cabledyn_build_bin} and ${_cabledyn_runtime_dir}"
    VERBATIM)
  add_custom_target(cabledyn_stage_runtime_dll_set ALL
                    DEPENDS ${_cabledyn_staged_dlls} ${_cabledyn_build_dlls})
endif()

function(cabledyn_stage_runtime_dlls target_name)
  if(CABLEDYN_RUNTIME_DLLS)
    # The shared DLL set is staged once by cabledyn_stage_runtime_dll_set; depend on it so the
    # bundle is complete when this target builds. Copy only this target's OWN executable here --
    # a distinct output per target, so no cross-target race on shared files.
    add_dependencies(${target_name} cabledyn_stage_runtime_dll_set)
    add_custom_command(TARGET ${target_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory
              "${CMAKE_BINARY_DIR}/runtime"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              $<TARGET_FILE:${target_name}>
              "${CMAKE_BINARY_DIR}/runtime"
      VERBATIM)
  endif()
endfunction()

# --- the certified production core ---
# ONE compilation of the core sources (an OBJECT library), TWO link products:
#   cabledyn_core   (STATIC)  -- every in-repo executable/test and the OpenFAST-tree
#                                ExternalProject link against this, unchanged;
#   cabledyn_shared (SHARED)  -- the C-ABI product (CableDyn_CAPI's BIND(C) surface)
#                                consumed by the ctypes-based `cabledyn` Python package
#                                and any CFD-side C caller, the MoorDyn-C deployment
#                                pattern. PIC on the object library serves both.
add_library(cabledyn_obj OBJECT
  src/cabledyn_mutex.c
  src/cabledyn_blas.c
  src/cabledyn_crt_locale.c
  src/cabledyn_path.c
  src/cabledyn_fatal.c
  src/CableDyn_Precision.f90
  src/CableDyn_FatalReport.f90
  src/CableDyn_PathIO.f90
  src/CableDyn_Banner.f90
  src/CableDyn_Linalg.f90
  src/CableDyn_AD.f90
  src/CableDyn_SO3.f90
  src/CableDyn_HermiteArch.f90
  src/CableDyn_EndConnection.f90
  src/CableDyn_Vessel.f90
  src/CableDyn_HermiteCable.f90
  src/CableDyn_HermiteTorsion.f90
  src/CableDyn_HermiteCableStatic.f90
  src/CableDyn_HermiteCableDynamic.f90
  src/CableDyn_Conventions.f90
  src/CableDyn_RigidKinematics.f90
  src/CableDyn_Cosserat.f90
  src/CableDyn_CosseratAssemble.f90
  src/CableDyn_CosseratStatic.f90
  src/CableDyn_CosseratDynamic.f90
  src/CableDyn_CosseratEMC.f90
  src/CableDyn_FiniteEIModel.f90
  src/CableDyn_SeabedContact.f90
  src/CableDyn_Mesh.f90
  src/CableDyn_CableElem.f90
  src/CableDyn_Assemble.f90
  src/CableDyn_Loads.f90
  src/CableDyn_Bathymetry.f90
  src/CableDyn_Hydro.f90
  src/CableDyn_WaveSpectra.f90
  src/CableDyn_StreamWave.f90
  src/CableDyn_Static.f90
  src/CableDyn_Catenary.f90
  src/CableDyn_FiniteEIStatic.f90
  src/CableDyn_CoupledPlatform.f90
  src/CableDyn_Line.f90
  src/CableDyn_Damping.f90
  src/CableDyn_Viscoelastic.f90
  src/CableDyn_Syrope.f90
  src/CableDyn_Dynamic.f90
  src/CableDyn_Model.f90
  src/CableDyn_Modal.f90
  src/CableDyn_System.f90
  src/CableDyn_OpenFAST_Types.f90
  src/CableDyn_OpenFAST.f90
  src/CableDyn_OpenFAST_Mesh.f90
  src/CableDyn_OpenFAST_FMF.f90
  src/CableDyn_OpenFAST_HermiteFMF.f90
  src/CableDyn_CAPI.f90
  src/CableDyn_Driver.f90
  src/CableDyn_RangeOutput.f90
  src/CableDyn_DeckDriver.f90
  src/CableDyn_OpenFAST_Aggregate.f90
)
set_target_properties(cabledyn_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(cabledyn_obj PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
target_include_directories(cabledyn_obj PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}")
target_link_libraries(cabledyn_obj PUBLIC ${CABLEDYN_LAPACK})
if(CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS)
  target_compile_definitions(cabledyn_obj PRIVATE CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS=1)
endif()
if(CABLEDYN_BLAS_RUNTIME_LOAD)
  target_compile_definitions(cabledyn_obj PRIVATE CABLEDYN_BLAS_RUNTIME_LOAD=1)
endif()
if(NOT WIN32)
  # src/cabledyn_mutex.c and src/cabledyn_blas.c use POSIX threads off Windows.
  find_package(Threads REQUIRED)
  target_link_libraries(cabledyn_obj PUBLIC Threads::Threads)
endif()
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cabledyn_obj PUBLIC OpenMP::OpenMP_Fortran)
endif()
if(CABLEDYN_ENABLE_COVERAGE)
  target_compile_options(cabledyn_obj PRIVATE -O0 -g --coverage)
  target_link_options(cabledyn_obj PUBLIC --coverage)
endif()

add_library(cabledyn_core STATIC $<TARGET_OBJECTS:cabledyn_obj>)
target_link_libraries(cabledyn_core PUBLIC cabledyn_obj)
if(MINGW)
  # Pull src/cabledyn_crt_locale.c into every executable that links the archive, so
  # its load-time constructor installs the libgfortran locale guard before any I/O.
  target_link_options(cabledyn_core INTERFACE "LINKER:--undefined=cabledyn_crt_locale_guard")
endif()

add_library(cabledyn_shared SHARED $<TARGET_OBJECTS:cabledyn_obj>)
target_link_libraries(cabledyn_shared PUBLIC cabledyn_obj)
set_target_properties(cabledyn_shared PROPERTIES
  OUTPUT_NAME cabledyn
  VERSION ${PROJECT_VERSION}
  SOVERSION 1)
# The shared library exports exactly the BIND(C) entry points declared in
# src/CableDyn_CAPI.h. Fortran module procedures and the registry mutex stay
# internal so they cannot become accidental ABI or collide with host symbols.
if(WIN32)
  target_sources(cabledyn_shared PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.def")
elseif(APPLE)
  target_link_options(cabledyn_shared PRIVATE
    "LINKER:-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.exp")
else()
  target_link_options(cabledyn_shared PRIVATE
    "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.map")
endif()
set_property(TARGET cabledyn_shared APPEND PROPERTY LINK_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.def"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.exp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cabledyn.map")

# Install the supported C-ABI product as a conventional runtime/development pair.
# The static core remains an internal/OpenFAST integration artifact; external
# consumers receive the versioned shared library and its frozen public header.
include(GNUInstallDirs)
install(TARGETS cabledyn_shared
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES src/CableDyn_CAPI.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# --- unit tests (CTest assertion programs) ---
#
# The module unit tests cross-validate against independent reference values, so they carry
# `fortran;unit;reference`; capability labels (`static`, `dynamic`) select a subset.
# The driver test is an end-to-end `integration` test (file -> solve -> CSV), not
# a reference-parity check, so it carries `fortran;integration;static` instead.
# Useful selections:
#   ctest --test-dir build -L fortran -LE slow  # fast Fortran developer loop
#   ctest --test-dir build -L reference    # the reference-parity unit checks
#   ctest --test-dir build -L integration  # the end-to-end driver check
#   ctest --test-dir build -L slow --parallel 4     # long validation cells
#   ctest --test-dir build -L 'static|dynamic'  # the solver checks (regex OR)
enable_testing()
add_executable(test_ad tests/test_ad.f90)
target_link_libraries(test_ad PRIVATE cabledyn_core)
add_test(NAME ad COMMAND test_ad)
set_tests_properties(ad PROPERTIES LABELS "fortran;unit")

add_executable(test_linalg tests/test_linalg.f90)
target_link_libraries(test_linalg PRIVATE cabledyn_core)
add_test(NAME linalg COMMAND test_linalg)
set_tests_properties(linalg PROPERTIES LABELS "fortran;unit")

# CD_All_Finite on signalling NaN and infinities with floating-point trapping on in every
# build type: a floating-point compare inside the finite check would stop the program.
add_executable(test_all_finite_traps tests/test_all_finite_traps.f90)
target_link_libraries(test_all_finite_traps PRIVATE cabledyn_core)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(test_all_finite_traps PRIVATE -ffpe-trap=invalid,zero,overflow)
endif()
add_test(NAME all_finite_traps COMMAND test_all_finite_traps)
set_tests_properties(all_finite_traps PROPERTIES LABELS "fortran;unit")

add_executable(test_so3 tests/test_so3.f90)
target_link_libraries(test_so3 PRIVATE cabledyn_core)
add_test(NAME so3 COMMAND test_so3)
set_tests_properties(so3 PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_conventions tests/test_conventions.f90)
target_link_libraries(test_conventions PRIVATE cabledyn_core)
add_test(NAME conventions COMMAND test_conventions)
set_tests_properties(conventions PROPERTIES LABELS "fortran;unit")

add_executable(test_rigid_kinematics tests/test_rigid_kinematics.f90)
target_link_libraries(test_rigid_kinematics PRIVATE cabledyn_core)
add_test(NAME rigid_kinematics COMMAND test_rigid_kinematics)
set_tests_properties(rigid_kinematics PROPERTIES LABELS "fortran;unit;openfast")

add_executable(test_cosserat_strains tests/test_cosserat_strains.f90)
target_link_libraries(test_cosserat_strains PRIVATE cabledyn_core)
add_test(NAME cosserat_strains COMMAND test_cosserat_strains)
set_tests_properties(cosserat_strains PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cosserat_force tests/test_cosserat_force.f90)
target_link_libraries(test_cosserat_force PRIVATE cabledyn_core)
add_test(NAME cosserat_force COMMAND test_cosserat_force)
set_tests_properties(cosserat_force PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cosserat_assemble tests/test_cosserat_assemble.f90)
target_link_libraries(test_cosserat_assemble PRIVATE cabledyn_core)
add_test(NAME cosserat_assemble COMMAND test_cosserat_assemble)
set_tests_properties(cosserat_assemble PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cosserat_mass tests/test_cosserat_mass.f90)
target_link_libraries(test_cosserat_mass PRIVATE cabledyn_core)
add_test(NAME cosserat_mass COMMAND test_cosserat_mass)
set_tests_properties(cosserat_mass PROPERTIES LABELS "fortran;unit;reference;dynamic")

add_executable(test_cosserat_dynamic tests/test_cosserat_dynamic.f90)
target_link_libraries(test_cosserat_dynamic PRIVATE cabledyn_core)
add_test(NAME cosserat_dynamic COMMAND test_cosserat_dynamic)
set_tests_properties(cosserat_dynamic PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_finite_ei_large_rotation tests/test_finite_ei_large_rotation.f90)
target_link_libraries(test_finite_ei_large_rotation PRIVATE cabledyn_core)
add_test(NAME finite_ei_large_rotation COMMAND test_finite_ei_large_rotation)
set_tests_properties(finite_ei_large_rotation PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_precession tests/test_finite_ei_precession.f90)
target_link_libraries(test_finite_ei_precession PRIVATE cabledyn_core)
add_test(NAME finite_ei_precession COMMAND test_finite_ei_precession)
set_tests_properties(finite_ei_precession PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_gyroscopic tests/test_finite_ei_gyroscopic.f90)
target_link_libraries(test_finite_ei_gyroscopic PRIVATE cabledyn_core)
add_test(NAME finite_ei_gyroscopic COMMAND test_finite_ei_gyroscopic)
set_tests_properties(finite_ei_gyroscopic PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_mult_tangent tests/test_finite_ei_mult_tangent.f90)
target_link_libraries(test_finite_ei_mult_tangent PRIVATE cabledyn_core)
add_test(NAME finite_ei_mult_tangent COMMAND test_finite_ei_mult_tangent)
set_tests_properties(finite_ei_mult_tangent PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_mult_step tests/test_finite_ei_mult_step.f90)
target_link_libraries(test_finite_ei_mult_step PRIVATE cabledyn_core)
add_test(NAME finite_ei_mult_step COMMAND test_finite_ei_mult_step)
set_tests_properties(finite_ei_mult_step PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_gyro_primitives tests/test_finite_ei_gyro_primitives.f90)
target_link_libraries(test_finite_ei_gyro_primitives PRIVATE cabledyn_core)
add_test(NAME finite_ei_gyro_primitives COMMAND test_finite_ei_gyro_primitives)
set_tests_properties(finite_ei_gyro_primitives PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_gyro_tangent tests/test_finite_ei_gyro_tangent.f90)
target_link_libraries(test_finite_ei_gyro_tangent PRIVATE cabledyn_core)
add_test(NAME finite_ei_gyro_tangent COMMAND test_finite_ei_gyro_tangent)
set_tests_properties(finite_ei_gyro_tangent PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_gyro_dynamics tests/test_finite_ei_gyro_dynamics.f90)
target_link_libraries(test_finite_ei_gyro_dynamics PRIVATE cabledyn_core)
add_test(NAME finite_ei_gyro_dynamics COMMAND test_finite_ei_gyro_dynamics)
set_tests_properties(finite_ei_gyro_dynamics PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_emc tests/test_emc.f90)
target_link_libraries(test_emc PRIVATE cabledyn_core)
add_test(NAME emc COMMAND test_emc)
set_tests_properties(emc PROPERTIES LABELS "fortran;dynamic;finite_ei")

add_executable(test_finite_ei_model tests/test_finite_ei_model.f90)
target_link_libraries(test_finite_ei_model PRIVATE cabledyn_core)
add_test(NAME finite_ei_model COMMAND test_finite_ei_model)
set_tests_properties(finite_ei_model PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_hermite_arch tests/test_hermite_arch.f90)
target_link_libraries(test_hermite_arch PRIVATE cabledyn_core)
add_test(NAME hermite_arch COMMAND test_hermite_arch)
set_tests_properties(hermite_arch PROPERTIES LABELS "fortran;unit;static")

add_executable(test_hermite_cable tests/test_hermite_cable.f90)
target_link_libraries(test_hermite_cable PRIVATE cabledyn_core)
add_test(NAME hermite_cable COMMAND test_hermite_cable)
set_tests_properties(hermite_cable PROPERTIES LABELS "fortran;unit;static")

add_executable(test_hermite_torsion_kernel tests/test_hermite_torsion_kernel.f90)
target_link_libraries(test_hermite_torsion_kernel PRIVATE cabledyn_core)
add_test(NAME hermite_torsion_kernel
         COMMAND test_hermite_torsion_kernel ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/torsion_oracle_reference.txt)
set_tests_properties(hermite_torsion_kernel PROPERTIES LABELS "fortran;unit;static")

add_executable(test_hermite_torsion_static tests/test_hermite_torsion_static.f90)
target_link_libraries(test_hermite_torsion_static PRIVATE cabledyn_core)
add_test(NAME hermite_torsion_static
         COMMAND test_hermite_torsion_static ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/torsion_analytic_refs.txt)
set_tests_properties(hermite_torsion_static PROPERTIES LABELS "fortran;unit;static;validation")

add_executable(test_torsion_validation tests/test_torsion_validation.f90)
target_link_libraries(test_torsion_validation PRIVATE cabledyn_core)
add_test(NAME torsion_validation
         COMMAND test_torsion_validation ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/torsion_analytic_refs.txt)
set_tests_properties(torsion_validation PROPERTIES LABELS "fortran;unit;static;validation")
add_executable(test_hermite_torsion_dynamic tests/test_hermite_torsion_dynamic.f90)
target_link_libraries(test_hermite_torsion_dynamic PRIVATE cabledyn_core)
add_test(NAME hermite_torsion_dynamic COMMAND test_hermite_torsion_dynamic)
set_tests_properties(hermite_torsion_dynamic PROPERTIES LABELS "fortran;unit;dynamic;validation")

add_executable(test_hermite_cable_static tests/test_hermite_cable_static.f90)
target_link_libraries(test_hermite_cable_static PRIVATE cabledyn_core)
add_test(NAME hermite_cable_static COMMAND test_hermite_cable_static)
set_tests_properties(hermite_cable_static PROPERTIES LABELS "fortran;unit;static")

# Finite-EI mesh-resolution diagnosis: h*kappa peak / curvature-jump / contact-island metrics
# and the power-of-two refinement recommendation -- the preflight trigger for the adaptive mesh
# policy. Synthetic controlled profiles + fail-closed validation.
add_executable(test_hermite_resolution tests/test_hermite_resolution.f90)
target_link_libraries(test_hermite_resolution PRIVATE cabledyn_core)
add_test(NAME hermite_resolution COMMAND test_hermite_resolution)
# Finite-EI static physical-branch safeguards: fold/loop and kink audit, compression
# screens, equilibrium stability and the EI continuation from the catenary seed.
add_executable(test_hermite_branch_audit tests/test_hermite_branch_audit.f90)
target_link_libraries(test_hermite_branch_audit PRIVATE cabledyn_core)
add_test(NAME hermite_branch_audit COMMAND test_hermite_branch_audit)
set_tests_properties(hermite_branch_audit PROPERTIES LABELS "fortran;hermite;static;unit")
set_tests_properties(hermite_resolution PROPERTIES LABELS "fortran;unit;static")

# Metric-gated mesh economics: solve coarse, auto-refine to the minimum safe mesh only where the
# resolution metric demands it (refine-when-fragile, no over-refinement, cap, equivalence, fail-closed).
add_executable(test_hermite_automesh tests/test_hermite_automesh.f90)
target_link_libraries(test_hermite_automesh PRIVATE cabledyn_core)
add_test(NAME hermite_automesh COMMAND test_hermite_automesh)
set_tests_properties(hermite_automesh PROPERTIES LABELS "fortran;unit;static")

# Mesh-sequenced static entry: coarse-first continuation + Hermite prolongation + per-level
# polish -- equivalence vs the direct solve on the committed 800 m lazy-wave rig, determinism,
# finest-only constraint enforcement, fail-closed validation. Needs the committed seed file.
add_executable(test_hermite_sequenced tests/test_hermite_sequenced.f90)
target_link_libraries(test_hermite_sequenced PRIVATE cabledyn_core)
add_test(NAME hermite_sequenced COMMAND test_hermite_sequenced WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(hermite_sequenced PROPERTIES LABELS "fortran;unit;static")

add_executable(test_hermite_cable_dynamic tests/test_hermite_cable_dynamic.f90)
target_link_libraries(test_hermite_cable_dynamic PRIVATE cabledyn_core)
add_test(NAME hermite_cable_dynamic COMMAND test_hermite_cable_dynamic)
set_tests_properties(hermite_cable_dynamic PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_hermite_drag_dynamics tests/test_hermite_drag_dynamics.f90)
target_link_libraries(test_hermite_drag_dynamics PRIVATE cabledyn_core)
add_test(NAME hermite_drag_dynamics COMMAND test_hermite_drag_dynamics)
set_tests_properties(hermite_drag_dynamics PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_hermite_axial_damping tests/test_hermite_axial_damping.f90)
target_link_libraries(test_hermite_axial_damping PRIVATE cabledyn_core)
add_test(NAME hermite_axial_damping COMMAND test_hermite_axial_damping)
set_tests_properties(hermite_axial_damping PROPERTIES LABELS "fortran;unit;dynamic;damping")

add_executable(test_hermite_added_mass tests/test_hermite_added_mass.f90)
target_link_libraries(test_hermite_added_mass PRIVATE cabledyn_core)
add_test(NAME hermite_added_mass COMMAND test_hermite_added_mass)
set_tests_properties(hermite_added_mass PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_hermite_wave_forcing tests/test_hermite_wave_forcing.f90)
target_link_libraries(test_hermite_wave_forcing PRIVATE cabledyn_core)
add_test(NAME hermite_wave_forcing COMMAND test_hermite_wave_forcing)
set_tests_properties(hermite_wave_forcing PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_cable_mesh tests/test_cable_mesh.f90)
target_link_libraries(test_cable_mesh PRIVATE cabledyn_core)
add_test(NAME cable_mesh COMMAND test_cable_mesh)
set_tests_properties(cable_mesh PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cable_elem tests/test_cable_elem.f90)
target_link_libraries(test_cable_elem PRIVATE cabledyn_core)
add_test(NAME cable_elem COMMAND test_cable_elem)
set_tests_properties(cable_elem PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cable_assemble tests/test_cable_assemble.f90)
target_link_libraries(test_cable_assemble PRIVATE cabledyn_core)
add_test(NAME cable_assemble COMMAND test_cable_assemble)
set_tests_properties(cable_assemble PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_cable_loads tests/test_cable_loads.f90)
target_link_libraries(test_cable_loads PRIVATE cabledyn_core)
add_test(NAME cable_loads COMMAND test_cable_loads)
set_tests_properties(cable_loads PROPERTIES LABELS "fortran;unit;reference")

add_executable(test_bathymetry tests/test_bathymetry.f90)
target_link_libraries(test_bathymetry PRIVATE cabledyn_core)
add_test(NAME bathymetry COMMAND test_bathymetry)
set_tests_properties(bathymetry PROPERTIES LABELS "fortran;unit;seabed")

add_executable(test_hydro tests/test_hydro.f90)
target_link_libraries(test_hydro PRIVATE cabledyn_core)
add_test(NAME hydro COMMAND test_hydro)
set_tests_properties(hydro PROPERTIES LABELS "fortran;unit;hydro")

add_executable(test_cable_static tests/test_cable_static.f90)
target_link_libraries(test_cable_static PRIVATE cabledyn_core)
add_test(NAME cable_static COMMAND test_cable_static)
set_tests_properties(cable_static PROPERTIES LABELS "fortran;unit;reference;static")

add_executable(test_catenary tests/test_catenary.f90)
target_link_libraries(test_catenary PRIVATE cabledyn_core)
add_test(NAME catenary COMMAND test_catenary)
set_tests_properties(catenary PROPERTIES LABELS "fortran;unit;reference;static")

add_executable(test_robust_init tests/test_robust_init.f90)
target_link_libraries(test_robust_init PRIVATE cabledyn_core)
add_test(NAME robust_init COMMAND test_robust_init)
set_tests_properties(robust_init PROPERTIES LABELS "fortran;integration;static")

add_executable(test_line tests/test_line.f90)
target_link_libraries(test_line PRIVATE cabledyn_core)
add_test(NAME line COMMAND test_line)
set_tests_properties(line PROPERTIES LABELS "fortran;unit;reference;static")

add_executable(test_l2_chain_parity tests/test_l2_chain_parity.f90)
target_link_libraries(test_l2_chain_parity PRIVATE cabledyn_core)
add_test(NAME l2_chain_parity COMMAND test_l2_chain_parity)
set_tests_properties(l2_chain_parity PROPERTIES LABELS "fortran;l2;validation;static")

add_executable(test_l2_volturnus_mooring tests/test_l2_volturnus_mooring.f90)
target_link_libraries(test_l2_volturnus_mooring PRIVATE cabledyn_core)
add_test(NAME l2_volturnus_mooring COMMAND test_l2_volturnus_mooring)
set_tests_properties(l2_volturnus_mooring PROPERTIES LABELS "fortran;l2;validation;static;external_ref")

add_executable(test_l2_lozon80_mooring tests/test_l2_lozon80_mooring.f90)
target_link_libraries(test_l2_lozon80_mooring PRIVATE cabledyn_core)
add_test(NAME l2_lozon80_mooring COMMAND test_l2_lozon80_mooring)
set_tests_properties(l2_lozon80_mooring PROPERTIES LABELS "fortran;l2;validation;static;external_ref")

add_executable(test_l2_lozon200_semitaut tests/test_l2_lozon200_semitaut.f90)
target_link_libraries(test_l2_lozon200_semitaut PRIVATE cabledyn_core)
add_test(NAME l2_lozon200_semitaut COMMAND test_l2_lozon200_semitaut)
set_tests_properties(l2_lozon200_semitaut PROPERTIES LABELS "fortran;l2;validation;static;external_ref")

add_executable(test_l2_lozon200_semitaut_deck tests/test_l2_lozon200_semitaut_deck.f90)
target_link_libraries(test_l2_lozon200_semitaut_deck PRIVATE cabledyn_core)
add_test(NAME l2_lozon200_semitaut_deck COMMAND test_l2_lozon200_semitaut_deck)
set_tests_properties(l2_lozon200_semitaut_deck PROPERTIES LABELS "fortran;l2;validation;static;external_ref")

add_executable(test_l3_orcaflex_static tests/test_l3_orcaflex_static.f90)
target_link_libraries(test_l3_orcaflex_static PRIVATE cabledyn_core)
add_test(NAME l3_orcaflex_static COMMAND test_l3_orcaflex_static)
set_tests_properties(l3_orcaflex_static PROPERTIES LABELS "fortran;l3;validation;static;orcaflex")

add_executable(test_l3_orcaflex_current_static tests/test_l3_orcaflex_current_static.f90)
target_link_libraries(test_l3_orcaflex_current_static PRIVATE cabledyn_core)
add_test(NAME l3_orcaflex_current_static COMMAND test_l3_orcaflex_current_static)
set_tests_properties(l3_orcaflex_current_static PROPERTIES LABELS "fortran;l3;validation;static;orcaflex;hydro")

# L3-4 finite-EI touchdown on the PRODUCTION cubic-Hermite bending-cable path (tensions via
# EA(|m|-1), max curvature resampled onto OrcaFlex's segment-midpoint grid).
add_executable(test_l3_hermite_touchdown tests/test_l3_hermite_touchdown.f90)
target_link_libraries(test_l3_hermite_touchdown PRIVATE cabledyn_core)
add_test(NAME l3_finite_ei_touchdown COMMAND test_l3_hermite_touchdown)
set_tests_properties(l3_finite_ei_touchdown PROPERTIES LABELS "fortran;l3;validation;static;orcaflex")

# Same case on the secondary Cosserat (rotation-DOF) path, kept as regression coverage.
add_executable(test_l3_orcaflex_finite_ei_static tests/test_l3_orcaflex_finite_ei_static.f90)
target_link_libraries(test_l3_orcaflex_finite_ei_static PRIVATE cabledyn_core)
add_test(NAME cosserat_orcaflex_finite_ei_static COMMAND test_l3_orcaflex_finite_ei_static)
set_tests_properties(cosserat_orcaflex_finite_ei_static PROPERTIES LABELS "fortran;cosserat;static;orcaflex")

# Finite-EI lazy-wave power-cable curvature across the three Lozon (2025) depths (80/200/800 m,
# VolturnUS-S) vs OrcaFlex. The committed suspended-span seeds live next to the test executable;
# the solve's buoyancy continuation makes them portable.
configure_file(tests/data/gomex80_lazywave_seed.xyz
               ${CMAKE_BINARY_DIR}/gomex80_lazywave_seed.xyz COPYONLY)
configure_file(tests/data/gomaine200_lazywave_seed.xyz
               ${CMAKE_BINARY_DIR}/gomaine200_lazywave_seed.xyz COPYONLY)
configure_file(tests/data/humboldt_lazywave_seed.xyz
               ${CMAKE_BINARY_DIR}/humboldt_lazywave_seed.xyz COPYONLY)
configure_file(tests/data/lozon_coupled_harmonics.dat
               ${CMAKE_BINARY_DIR}/lozon_coupled_harmonics.dat COPYONLY)
add_executable(test_lazywave_dynamics tests/test_lazywave_dynamics.f90)
target_link_libraries(test_lazywave_dynamics PRIVATE cabledyn_core)
add_test(NAME lazywave_dynamics COMMAND test_lazywave_dynamics WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(lazywave_dynamics PROPERTIES LABELS "fortran;dynamic;slow")

add_executable(test_hermite_axial_drag_arbiter tests/test_hermite_axial_drag_arbiter.f90)
target_link_libraries(test_hermite_axial_drag_arbiter PRIVATE cabledyn_core)
add_test(NAME hermite_axial_drag_arbiter COMMAND test_hermite_axial_drag_arbiter)
set_tests_properties(hermite_axial_drag_arbiter PROPERTIES LABELS "fortran;validation;l3;dynamic;arbiter")

add_executable(test_l3_longduration_mooring tests/test_l3_longduration_mooring.f90)
target_link_libraries(test_l3_longduration_mooring PRIVATE cabledyn_core)
add_test(NAME l3_longduration_mooring COMMAND test_l3_longduration_mooring)
set_tests_properties(l3_longduration_mooring PROPERTIES LABELS "fortran;validation;l3;l5;dynamic;slow")

add_executable(test_l3_lazywave_dynamic tests/test_l3_lazywave_dynamic.f90)
target_link_libraries(test_l3_lazywave_dynamic PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_dynamic COMMAND test_l3_lazywave_dynamic WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_dynamic PROPERTIES LABELS "fortran;validation;l3;dynamic;slow")
# The same parity gate on the default force blend (CD_L3_FORCE_BLEND=1) with the same Newton
# tolerance; it writes its own curvature-profile file, so the two can run concurrently.
add_test(NAME l3_lazywave_dynamic_force_blend COMMAND test_l3_lazywave_dynamic WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_dynamic_force_blend PROPERTIES LABELS "fortran;validation;l3;dynamic;slow"
                     ENVIRONMENT "CD_L3_FORCE_BLEND=1")

# L4-1c: the lazy-wave cables under a realistic COUPLED FOWT hang-off motion (the shared
# N-harmonic surge+heave table), scored vs OrcaFlex at all three depths.
add_executable(test_l4_lazywave_coupled tests/test_l4_lazywave_coupled.f90)
target_link_libraries(test_l4_lazywave_coupled PRIVATE cabledyn_core)
add_test(NAME l4_lazywave_coupled COMMAND test_l4_lazywave_coupled WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l4_lazywave_coupled PROPERTIES LABELS "fortran;validation;l4;dynamic;slow")

add_executable(test_l3_lazywave_curvature tests/test_l3_lazywave_curvature.f90)
target_link_libraries(test_l3_lazywave_curvature PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_curvature COMMAND test_l3_lazywave_curvature)
set_tests_properties(l3_lazywave_curvature PROPERTIES
                     LABELS "fortran;l3;validation;static;orcaflex;external_ref;slow"
                     TIMEOUT 480)

# L3-5 companion: Humboldt 800 m sag-bend curvature MESH CONVERGENCE (480 -> 960 elements by
# mesh-sequencing continuation; documents that the peak is mesh-converged inside the OrcaFlex band).
add_executable(test_l3_lazywave_refinement tests/test_l3_lazywave_refinement.f90)
target_link_libraries(test_l3_lazywave_refinement PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_refinement COMMAND test_l3_lazywave_refinement WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_refinement PROPERTIES
                     LABELS "fortran;l3;validation;static;orcaflex;slow"
                     TIMEOUT 900)

# L3-5 closure cell: the SAME 800 m cable solved as the FULL INSTALLED cable (fairlead + seabed
# anchor + grounded run) vs the OrcaFlex installed-cable reference -- the fourth corner of the
# pinned-span-vs-installed protocol attribution, measured on both protocols. (The executable
# name avoids the substring "install": Windows' UAC installer-detection heuristic refuses to
# launch unelevated executables whose file name matches install/setup/update/patch.)
add_executable(test_l3_lazywave_fullcable tests/test_l3_lazywave_fullcable.f90)
target_link_libraries(test_l3_lazywave_fullcable PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_fullcable COMMAND test_l3_lazywave_fullcable WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_fullcable PROPERTIES
                     LABELS "fortran;l3;validation;static;orcaflex;slow"
                     TIMEOUT 900)

# STATIC-INITIALIZATION STRESS MATRIX: the robustness matrix for the sequenced
# static entry -- 122 parametric cases across depth / slackness / buoyancy ratio / EA-EI /
# seabed stiffness / span ratio / buoyant fraction, each solved from a trivial bowed seed
# with the production settings and gated on convergence + a turn-angle kink detector +
# penalty-consistent penetration (the committed reference gates certify specific cases;
# this certifies the initialization CLASS).
add_executable(test_static_init_sweep tests/test_static_init_sweep.f90)
target_link_libraries(test_static_init_sweep PRIVATE cabledyn_core)
add_test(NAME static_init_sweep COMMAND test_static_init_sweep WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
# The hang guard allows for a checked Debug build, about ten times slower than Release.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(_static_init_sweep_timeout 5400)
else()
  set(_static_init_sweep_timeout 900)
endif()
set_tests_properties(static_init_sweep PROPERTIES
                     LABELS "fortran;hermite;static;robustness;slow"
                     TIMEOUT ${_static_init_sweep_timeout})

# L3-6e: OUT-OF-PLANE (3D) dynamic lazy-wave -- the committed dynamic gates are planar;
# this one frees the y DOFs and drives the 80 m hang-off with an elliptical heave+sway
# whirl vs the OrcaFlex 3D reference.
add_executable(test_l3_lazywave_3d_dynamic tests/test_l3_lazywave_3d_dynamic.f90)
target_link_libraries(test_l3_lazywave_3d_dynamic PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_3d_dynamic COMMAND test_l3_lazywave_3d_dynamic
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_3d_dynamic PROPERTIES
                     LABELS "fortran;l3;validation;dynamic;orcaflex;slow"
                     TIMEOUT 900)

# L3-6d: MOVING-TOUCHDOWN dynamic lazy-wave -- the unpinned successor to the pinned L3-6
# protocol (Lozon 80 m installed cable on the live penalty seabed, heave-driven TDP
# migration) vs the OrcaFlex installed-cable dynamic reference.
add_executable(test_l3_lazywave_touchdown_dynamic tests/test_l3_lazywave_touchdown_dynamic.f90)
target_link_libraries(test_l3_lazywave_touchdown_dynamic PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_touchdown_dynamic COMMAND test_l3_lazywave_touchdown_dynamic
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_touchdown_dynamic PROPERTIES
                     LABELS "fortran;l3;validation;dynamic;orcaflex;slow"
                     TIMEOUT 900)

# Regression for the pinned-end tangent mode observed in 200 m coupled OpenFAST runs.
# A section-aligned installed lazy wave receives a smooth lateral hang-off translation
# and is then held for a 30 s march.  The motion must enter the suspended span without
# a mesh-local curvature spike or collapse of |dr/ds| at the moving endpoint.
add_executable(test_pinned_touchdown_stability tests/test_pinned_touchdown_stability.f90)
target_link_libraries(test_pinned_touchdown_stability PRIVATE cabledyn_core)
add_test(NAME pinned_touchdown_stability COMMAND test_pinned_touchdown_stability
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(pinned_touchdown_stability PROPERTIES
                     LABELS "fortran;hermite;dynamic;contact;regression;slow"
                     TIMEOUT 900)

# Component (irregular) wave path: the general component-table kinematics kernel and its
# wiring into the Hermite dynamic model (the deterministic component-matched irregular-sea
# machinery behind L3-6c).
add_executable(test_hermite_irregular_waves tests/test_hermite_irregular_waves.f90)
target_link_libraries(test_hermite_irregular_waves PRIVATE cabledyn_core)
add_test(NAME hermite_irregular_waves COMMAND test_hermite_irregular_waves)
set_tests_properties(hermite_irregular_waves PROPERTIES LABELS "fortran;unit;dynamic")

# L3-6c: irregular-sea (JONSWAP) dynamic lazy-wave vs OrcaFlex -- the component-matched
# deterministic protocol (CableDyn's seeded JONSWAP realisation, given to both codes).
add_executable(test_l3_lazywave_irregular tests/test_l3_lazywave_irregular.f90)
target_link_libraries(test_l3_lazywave_irregular PRIVATE cabledyn_core)
add_test(NAME l3_lazywave_irregular COMMAND test_l3_lazywave_irregular
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(l3_lazywave_irregular PROPERTIES
                     LABELS "fortran;l3;validation;dynamic;orcaflex;slow"
                     TIMEOUT 900)

add_executable(test_l4_platform_mooring tests/test_l4_platform_mooring.f90)
target_link_libraries(test_l4_platform_mooring PRIVATE cabledyn_core)
add_test(NAME l4_platform_mooring COMMAND test_l4_platform_mooring)
set_tests_properties(l4_platform_mooring PROPERTIES LABELS "fortran;l4;validation;static")

add_executable(test_damping tests/test_damping.f90)
target_link_libraries(test_damping PRIVATE cabledyn_core)
add_test(NAME damping COMMAND test_damping)
set_tests_properties(damping PROPERTIES LABELS "fortran;unit;reference;dynamic")

add_executable(test_viscoelastic tests/test_viscoelastic.f90)
target_link_libraries(test_viscoelastic PRIVATE cabledyn_core)
add_test(NAME viscoelastic COMMAND test_viscoelastic)
set_tests_properties(viscoelastic PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_syrope tests/test_syrope.f90)
target_link_libraries(test_syrope PRIVATE cabledyn_core)
add_test(NAME syrope COMMAND test_syrope)
set_tests_properties(syrope PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_core_release_regressions tests/test_core_release_regressions.f90)
target_link_libraries(test_core_release_regressions PRIVATE cabledyn_core)
add_test(NAME core_release_regressions COMMAND test_core_release_regressions)
set_tests_properties(core_release_regressions PROPERTIES LABELS "fortran;unit;dynamic;static")

add_executable(test_core_solver_regressions tests/test_core_solver_regressions.f90)
target_link_libraries(test_core_solver_regressions PRIVATE cabledyn_core)
add_test(NAME core_solver_regressions COMMAND test_core_solver_regressions)
set_tests_properties(core_solver_regressions PROPERTIES LABELS "fortran;unit;dynamic;static")

add_executable(test_aggregate_step_atomicity tests/test_aggregate_step_atomicity.f90)
target_link_libraries(test_aggregate_step_atomicity PRIVATE cabledyn_core)
add_test(NAME aggregate_step_atomicity COMMAND test_aggregate_step_atomicity
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(aggregate_step_atomicity PROPERTIES LABELS "fortran;unit;dynamic;capi")

add_executable(test_core_input_guards tests/test_core_input_guards.f90)
target_link_libraries(test_core_input_guards PRIVATE cabledyn_core)
add_test(NAME core_input_guards COMMAND test_core_input_guards WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(core_input_guards PROPERTIES LABELS "fortran;unit;static")

add_executable(test_prescribed_motion tests/test_prescribed_motion.f90)
target_link_libraries(test_prescribed_motion PRIVATE cabledyn_core)
add_test(NAME prescribed_motion COMMAND test_prescribed_motion)
set_tests_properties(prescribed_motion PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_model tests/test_model.f90)
target_link_libraries(test_model PRIVATE cabledyn_core)
add_test(NAME model COMMAND test_model)
set_tests_properties(model PROPERTIES LABELS "fortran;unit;dynamic")

add_executable(test_system tests/test_system.f90)
target_link_libraries(test_system PRIVATE cabledyn_core)
add_test(NAME system COMMAND test_system)
set_tests_properties(system PROPERTIES LABELS "fortran;unit;dynamic")

# A massless Free junction between lines carrying large Morison added mass stays stable at
# every dt (attached end-node added-mass self block treated implicitly by the point step).
add_executable(test_point_added_mass tests/test_point_added_mass.f90)
target_link_libraries(test_point_added_mass PRIVATE cabledyn_core)
add_test(NAME point_added_mass COMMAND test_point_added_mass)
set_tests_properties(point_added_mass PROPERTIES LABELS "fortran;unit;dynamic")

# Deck-level: a massless Connect joining a rope and a buoyancy section (Can = 1) marches
# 60 s under waves at a typical dtM with finite, bounded loads.
add_executable(test_deck_massless_connect tests/test_deck_massless_connect.f90)
target_link_libraries(test_deck_massless_connect PRIVATE cabledyn_core)
add_test(NAME deck_massless_connect COMMAND test_deck_massless_connect WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(deck_massless_connect PROPERTIES LABELS "fortran;driver;integration;dynamic")

add_executable(test_openfast_shell tests/test_openfast_shell.f90)
target_link_libraries(test_openfast_shell PRIVATE cabledyn_core)
add_test(NAME openfast_shell COMMAND test_openfast_shell)
set_tests_properties(openfast_shell PROPERTIES LABELS "fortran;unit;dynamic;openfast")

# The shell step keeps the point store on the committed step: prescribed points follow a
# non-dynamic step, and a failed dynamic-point step restores the step-entry points.
add_executable(test_fast_point_store tests/test_fast_point_store.f90)
target_link_libraries(test_fast_point_store PRIVATE cabledyn_core)
add_test(NAME fast_point_store COMMAND test_fast_point_store)
set_tests_properties(fast_point_store PROPERTIES LABELS "fortran;unit;dynamic;openfast")

# FMF lifecycle over the Hermite dynamic power cable: shell transparency (bit-identical to
# the direct step under the same prescribed kinematics), endpoint reaction balance and
# tension consistency at rest, dynamic fairlead-load sanity, fail-closed battery. Needs the
# committed 80 m seed file.
add_executable(test_hermite_fmf tests/test_hermite_fmf.f90)
target_link_libraries(test_hermite_fmf PRIVATE cabledyn_core)
add_test(NAME hermite_fmf COMMAND test_hermite_fmf WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(hermite_fmf PROPERTIES LABELS "fortran;unit;dynamic;openfast")

add_executable(test_deck_hermite_cable tests/test_deck_hermite_cable.f90)
target_link_libraries(test_deck_hermite_cable PRIVATE cabledyn_core)
add_test(NAME deck_hermite_cable COMMAND test_deck_hermite_cable WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(deck_hermite_cable PROPERTIES LABELS "fortran;unit;dynamic;openfast")
# OPTION cable_statics: the continuation and sequenced finite-EI static routes agree on the
# shipped Lozon power cables.
add_executable(test_cable_statics_modes tests/test_cable_statics_modes.f90)
target_link_libraries(test_cable_statics_modes PRIVATE cabledyn_core)
add_test(NAME cable_statics_modes COMMAND test_cable_statics_modes ${CMAKE_CURRENT_SOURCE_DIR}/examples)
set_tests_properties(cable_statics_modes PROPERTIES LABELS "fortran;driver;hermite;static")

# 6-DOF prescribed vessel motion (vesselMotion / vesselRAO): rigid-body kinematics, RAO
# phase convention and spectrum, pure-rotation bit-identity with motionFile rows, and the
# clamped hang-off moment under vessel pitch.
add_executable(test_vessel_motion tests/test_vessel_motion.f90)
target_link_libraries(test_vessel_motion PRIVATE cabledyn_core)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/vessel_motion)
add_test(NAME vessel_motion COMMAND test_vessel_motion WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/vessel_motion)
set_tests_properties(vessel_motion PROPERTIES LABELS "fortran;unit;dynamic;validation")

# Discrete ATTACHMENTS (buoyancy modules, clumps) on a cubic-Hermite cable: convergence to
# the smeared section as the module pitch shrinks, resolved local curvature between
# modules, stable dynamics, and the deck grammar.
add_executable(test_line_attachments tests/test_line_attachments.f90)
target_link_libraries(test_line_attachments PRIVATE cabledyn_core)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/line_attachments)
add_test(NAME line_attachments COMMAND test_line_attachments WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/line_attachments)
set_tests_properties(line_attachments PROPERTIES LABELS "fortran;static;dynamic;validation")

# OpenFAST CompMooring=5 aggregate: one deck carrying BOTH EI=0 mooring lines AND finite-EI
# Hermite cables composed behind the same FMF moving-point facade. Proves pure-mooring
# equivalence to CD_FMF (bit-for-bit), mixed-deck composition + column routing, pure-cable
# build, coupled finite-EI touchdown/contact, and mixed-deck-with-dynamic-point fail-closed.
add_executable(test_openfast_aggregate tests/test_openfast_aggregate.f90)
target_link_libraries(test_openfast_aggregate PRIVATE cabledyn_core)
add_test(NAME openfast_aggregate COMMAND test_openfast_aggregate WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(openfast_aggregate PROPERTIES LABELS "fortran;unit;dynamic;openfast")
add_test(NAME openfast_aggregate_touchdown COMMAND test_openfast_aggregate --touchdown-only
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(openfast_aggregate_touchdown PROPERTIES
                     LABELS "fortran;unit;dynamic;openfast;contact")
# both runs of the executable write and re-read the same scratch decks in the build root
set_tests_properties(openfast_aggregate openfast_aggregate_touchdown PROPERTIES
                     RESOURCE_LOCK openfast_aggregate_files)

# Maintained IEA-15 MW mixed deck: the finite-EI power cable has a genuine multi-node
# grounded run and marches through the OpenFAST aggregate at 40 Hz with host current.
add_executable(test_openfast_touchdown tests/test_openfast_touchdown.f90)
target_link_libraries(test_openfast_touchdown PRIVATE cabledyn_core)
add_test(NAME openfast_touchdown
         COMMAND test_openfast_touchdown ${CMAKE_SOURCE_DIR}/examples/iea15mw_umaine_mixed_cabledyn.dat)
set_tests_properties(openfast_touchdown PROPERTIES
                     LABELS "fortran;integration;dynamic;openfast;contact"
                     # Fix the OpenBLAS dispatch so the trajectory retains coverage of
                     # the round-off-limited Newton state found on older x86 runners.
                     # Other BLAS implementations ignore this environment variable.
                     ENVIRONMENT "OPENBLAS_CORETYPE=Sandybridge")

# The same deck initialises without a static compression/mesh note: the touchdown
# contact kink is element-local, the smoothed element-mean axial force is tensile.
add_executable(test_flagship_compression_note tests/test_flagship_compression_note.f90)
target_link_libraries(test_flagship_compression_note PRIVATE cabledyn_core)
add_test(NAME flagship_compression_note
         COMMAND test_flagship_compression_note ${CMAKE_SOURCE_DIR}/examples/iea15mw_umaine_mixed_cabledyn.dat)
set_tests_properties(flagship_compression_note PROPERTIES LABELS "fortran;integration;static;contact")

add_executable(test_hfmf_mirror tests/test_hfmf_mirror.f90)
target_link_libraries(test_hfmf_mirror PRIVATE cabledyn_core)
add_test(NAME hfmf_mirror COMMAND test_hfmf_mirror WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(hfmf_mirror PROPERTIES LABELS "fortran;unit;dynamic;openfast")

add_executable(test_c_api tests/test_c_api.f90)
target_link_libraries(test_c_api PRIVATE cabledyn_core)
add_test(NAME c_api COMMAND test_c_api)
set_tests_properties(c_api PROPERTIES LABELS "fortran;unit;dynamic;capi")

# Compile and run through the installed-style C header against the SHARED ABI
# product. This catches missing Windows exports, C/Fortran prototype drift, ABI
# version drift, and loader failures that the Fortran-native lifecycle gate cannot.
add_executable(test_c_api_shared tests/test_c_api_shared.c)
target_include_directories(test_c_api_shared PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(test_c_api_shared PRIVATE cabledyn_shared)
add_test(NAME c_api_shared
         COMMAND test_c_api_shared ${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_volturnus_mooring.dat
                 ${CMAKE_CURRENT_SOURCE_DIR}/examples/dynamic_chain_waves.dat)
set_tests_properties(c_api_shared PROPERTIES LABELS "c;unit;capi;shared")

find_package(Threads QUIET)
if(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
  add_executable(test_c_api_concurrency tests/test_c_api_concurrency.c)
  target_include_directories(test_c_api_concurrency PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
  target_link_libraries(test_c_api_concurrency PRIVATE cabledyn_core Threads::Threads)
  add_test(NAME c_api_concurrency COMMAND test_c_api_concurrency)
  set_tests_properties(c_api_concurrency PROPERTIES LABELS "c;unit;capi;threading")
else()
  message(STATUS "CableDyn: pthread C API concurrency stress gate unavailable on this toolchain; "
                 "the portable C API lifecycle gate remains enabled")
endif()

# Portable (Win32 or POSIX threads) stack bound of the shared library: each deck is
# initialised, stepped, and evaluated on a thread with a 1 MiB stack (a host program, not
# CableDyn, chooses the stack of the threads that call the library). The decks are a
# 2048-element finite-EI cable (12294 DOFs), an EI=0 + finite-EI aggregate, and an EI=0
# mooring; a stack need that grew with the model would overflow here.
add_executable(test_c_api_stack tests/test_c_api_stack.c)
target_include_directories(test_c_api_stack PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(test_c_api_stack PRIVATE $<TARGET_LINKER_FILE:cabledyn_shared>)
add_dependencies(test_c_api_stack cabledyn_shared)
if(NOT WIN32)
  target_link_libraries(test_c_api_stack PRIVATE Threads::Threads m)
endif()
add_test(NAME stack_c_api
         COMMAND test_c_api_stack 1024
                 0.05 ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/stack_humboldt800_2048_elements.dat
                 0.025 ${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_umaine_mixed_cabledyn.dat
                 0.1 ${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_volturnus_mooring.dat)
set_tests_properties(stack_c_api PROPERTIES LABELS "c;unit;capi;stack" TIMEOUT 900)

# Portable (Win32 or POSIX threads) gate: concurrent InitDeck of one deck and of decks
# sharing a Syrope OWC table in the user's locale, concurrent Step on distinct handles,
# the missing-file diagnostic, and the BLAS thread/memory policy of CableDyn_Create.
add_executable(test_c_api_threads tests/test_c_api_threads.c)
target_include_directories(test_c_api_threads PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
# Link the import library alone: the usage requirements of cabledyn_shared carry the
# OpenBLAS import library, whose re-exported C runtime symbols (printf) would otherwise
# bind this C executable to openblas.dll at load time and defeat the memory check.
target_link_libraries(test_c_api_threads PRIVATE $<TARGET_LINKER_FILE:cabledyn_shared>)
add_dependencies(test_c_api_threads cabledyn_shared)
if(WIN32)
  target_link_libraries(test_c_api_threads PRIVATE psapi)
else()
  target_link_libraries(test_c_api_threads PRIVATE Threads::Threads m)
endif()
add_test(NAME c_api_threads
         COMMAND test_c_api_threads ${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat
                 ${CMAKE_CURRENT_SOURCE_DIR}/examples/syrope_polyester_mooring.dat
                 $<IF:$<BOOL:${CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS}>,1,0>)
# TIMEOUT: a process hung in its exit path cannot be terminated and keeps its output
# pipes open; without a timeout CTest waits for it forever.
set_tests_properties(c_api_threads PROPERTIES LABELS "c;unit;capi;threading" TIMEOUT 600)
if(CABLEDYN_BLAS_RUNTIME_LOAD)
  # openblas.dll must not be mapped before the first CableDyn call.
  target_compile_definitions(test_c_api_threads PRIVATE CABLEDYN_BLAS_RUNTIME_LOAD=1)
endif()
# The same gate with a 32-thread OpenMP setting, which a pthreads OpenBLAS loaded outside
# the one-thread policy would size its pool from: the memory bound then fails on any
# machine, not only on one with 32 logical CPUs.
add_test(NAME c_api_threads_omp32
         COMMAND test_c_api_threads ${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat
                 ${CMAKE_CURRENT_SOURCE_DIR}/examples/syrope_polyester_mooring.dat
                 $<IF:$<BOOL:${CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS}>,1,0>)
set_tests_properties(c_api_threads_omp32 PROPERTIES LABELS "c;unit;capi;threading"
                     ENVIRONMENT "OMP_NUM_THREADS=32" TIMEOUT 600)
# Repeated concurrent phases in one process, under the in-process watchdog.
add_test(NAME c_api_threads_stress
         COMMAND test_c_api_threads ${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat
                 ${CMAKE_CURRENT_SOURCE_DIR}/examples/syrope_polyester_mooring.dat
                 $<IF:$<BOOL:${CABLEDYN_HAVE_OPENBLAS_SET_NUM_THREADS}>,1,0>)
set_tests_properties(c_api_threads_stress PROPERTIES LABELS "c;unit;capi;threading" TIMEOUT 600
                     ENVIRONMENT "CABLEDYN_THREAD_TEST_LOOPS=5;CABLEDYN_THREAD_TEST_TIMEOUT=540")
if(WIN32)
  # Host threads that ran OpenMP teams, then an immediate process exit: MinGW-w64
  # libgomp could hang in its exit handler when a pool thread was still ending.
  add_test(NAME c_api_threads_exit
           COMMAND test_c_api_threads --exit-stress ${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat 1000)
  set_tests_properties(c_api_threads_exit PROPERTIES LABELS "c;unit;capi;threading" TIMEOUT 600
                       ENVIRONMENT "CABLEDYN_THREAD_TEST_TIMEOUT=540")
endif()
# Every LAPACK routine the Fortran sources call must be forwarded by src/cabledyn_blas.c;
# one that is not makes the shared library or the driver import openblas.dll at load time.
if(CABLEDYN_BLAS_RUNTIME_LOAD AND CMAKE_OBJDUMP)
  add_test(NAME blas_runtime_imports
           COMMAND ${CMAKE_COMMAND} -DOBJDUMP=${CMAKE_OBJDUMP}
                   "-DBINARIES=$<TARGET_FILE:cabledyn_shared>|$<TARGET_FILE:cabledyn>"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_blas_runtime_imports.cmake)
  set_tests_properties(blas_runtime_imports PROPERTIES LABELS "c;unit;capi")
endif()
# Without openblas.dll, the driver and the C API must fail at initialisation with the
# load diagnostic (library, location tried, load error, requesting entry point).
if(CABLEDYN_BLAS_RUNTIME_LOAD)
  add_executable(test_c_api_blas_missing tests/test_c_api_blas_missing.c)
  target_include_directories(test_c_api_blas_missing PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
  target_link_libraries(test_c_api_blas_missing PRIVATE $<TARGET_LINKER_FILE:cabledyn_shared>)
  add_dependencies(test_c_api_blas_missing cabledyn_shared)
  add_executable(test_blas_missing_solve tests/test_blas_missing_solve.f90)
  target_link_libraries(test_blas_missing_solve PRIVATE cabledyn_core)
  add_test(NAME blas_missing_runtime
           COMMAND ${CMAKE_COMMAND} "-DBIN_DIR=$<TARGET_FILE_DIR:cabledyn_shared>"
                   "-DWORK_DIR=${CMAKE_BINARY_DIR}/blas_missing_runtime"
                   "-DDRIVER=$<TARGET_FILE:cabledyn>" "-DCAPI_TEST=$<TARGET_FILE:test_c_api_blas_missing>"
                   "-DSOLVE_TEST=$<TARGET_FILE:test_blas_missing_solve>"
                   "-DDECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_blas_missing_runtime.cmake)
  set_tests_properties(blas_missing_runtime PROPERTIES LABELS "c;unit;capi;integration")
endif()

# A raw line with no external load converges at very small dt (round-off force scale).
add_executable(test_c_api_raw_small_dt tests/test_c_api_raw_small_dt.c)
target_include_directories(test_c_api_raw_small_dt PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(test_c_api_raw_small_dt PRIVATE cabledyn_core)
if(NOT WIN32)
  target_link_libraries(test_c_api_raw_small_dt PRIVATE m)
endif()
add_test(NAME c_api_raw_small_dt COMMAND test_c_api_raw_small_dt)
set_tests_properties(c_api_raw_small_dt PROPERTIES LABELS "c;unit;capi;dynamic")

# Long coupled marches must not grow the process memory: catches per-step leaks in the
# C API deck route (EI=0 chain) and the OpenFAST aggregate route (EI=0 moorings plus a
# finite-EI Hermite cable, with snapshot, fluid fields, outputs, and channels).
add_library(cabledyn_test_process_memory STATIC tests/process_memory.c)
if(WIN32)
  target_link_libraries(cabledyn_test_process_memory PUBLIC psapi)
endif()
add_executable(test_memory_growth tests/test_memory_growth.c)
target_include_directories(test_memory_growth PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(test_memory_growth PRIVATE cabledyn_core cabledyn_test_process_memory)
if(NOT WIN32)
  target_link_libraries(test_memory_growth PRIVATE m)
endif()
# Growth bounds [MiB]: observed 0.0 MiB over 20000 C-API steps and 0.11 MiB over 1000
# aggregate steps; 2 MiB leaves room for allocator page granularity.
add_test(NAME memory_growth
         COMMAND test_memory_growth 2 ${CMAKE_CURRENT_SOURCE_DIR}/examples/spread_3line_chain.dat 20000)
set_tests_properties(memory_growth PROPERTIES LABELS "c;capi;dynamic;memory" SKIP_RETURN_CODE 77)
add_executable(test_memory_growth_agg tests/test_memory_growth_agg.f90)
target_link_libraries(test_memory_growth_agg PRIVATE cabledyn_core cabledyn_test_process_memory)
add_test(NAME memory_growth_aggregate
         COMMAND test_memory_growth_agg ${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_umaine_mixed_cabledyn.dat
                 200.0 1000 2)
set_tests_properties(memory_growth_aggregate PROPERTIES LABELS "fortran;openfast;dynamic;memory"
                     SKIP_RETURN_CODE 77)

# Zero heap allocations in the per-step aggregate advance after warm-up, counted through
# GNU ld --wrap malloc/calloc/realloc wrappers (GNU toolchains only).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
  add_executable(test_step_allocations tests/test_step_allocations.f90 tests/allocation_counter.c)
  target_link_libraries(test_step_allocations PRIVATE cabledyn_core)
  target_link_options(test_step_allocations PRIVATE "-Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc")
  add_test(NAME step_allocations
           COMMAND test_step_allocations ${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_umaine_mixed_cabledyn.dat
                   200.0 200)
  set_tests_properties(step_allocations PROPERTIES LABELS "fortran;openfast;dynamic;memory")
  # The same gate with the range graphs of every line accumulated after each step.
  add_test(NAME step_allocations_ranges
           COMMAND test_step_allocations ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/range_mixed_alloc.dat
                   200.0 200 ranges)
  set_tests_properties(step_allocations_ranges PROPERTIES LABELS "fortran;openfast;dynamic;memory")
  # The same gate on a finite-EI cable with condensed torsion (a turning parent and a changing twist).
  add_test(NAME step_allocations_torsion COMMAND test_step_allocations torsion)
  set_tests_properties(step_allocations_torsion PROPERTIES LABELS "fortran;dynamic;memory")
endif()

add_executable(test_snap_load tests/test_snap_load.f90)
target_link_libraries(test_snap_load PRIVATE cabledyn_core)
add_test(NAME snap_load COMMAND test_snap_load)
set_tests_properties(snap_load PROPERTIES LABELS "fortran;l5;validation;dynamic")

add_executable(test_l2_2_dynamic_parity tests/test_l2_2_dynamic_parity.f90)
target_link_libraries(test_l2_2_dynamic_parity PRIVATE cabledyn_core)
add_test(NAME l2_2_dynamic_parity_wd0050_a1 COMMAND test_l2_2_dynamic_parity WD0050/A1)
add_test(NAME l2_2_dynamic_parity_wd0050_a2 COMMAND test_l2_2_dynamic_parity WD0050/A2)
add_test(NAME l2_2_dynamic_parity_wd0200_a1 COMMAND test_l2_2_dynamic_parity WD0200/A1)
add_test(NAME l2_2_dynamic_parity_wd0200_a2 COMMAND test_l2_2_dynamic_parity WD0200/A2)
add_test(NAME l2_2_dynamic_parity_wd0600_a1 COMMAND test_l2_2_dynamic_parity WD0600/A1)
add_test(NAME l2_2_dynamic_parity_wd0600_a2 COMMAND test_l2_2_dynamic_parity WD0600/A2)
set_tests_properties(l2_2_dynamic_parity_wd0050_a1 l2_2_dynamic_parity_wd0050_a2
                     l2_2_dynamic_parity_wd0200_a1 l2_2_dynamic_parity_wd0200_a2
                     l2_2_dynamic_parity_wd0600_a1 l2_2_dynamic_parity_wd0600_a2
                     PROPERTIES
                       TIMEOUT 600
                       WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                       ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1")
set_tests_properties(l2_2_dynamic_parity_wd0050_a1 l2_2_dynamic_parity_wd0050_a2
                     l2_2_dynamic_parity_wd0200_a1 l2_2_dynamic_parity_wd0200_a2
                     l2_2_dynamic_parity_wd0600_a1 l2_2_dynamic_parity_wd0600_a2
                     PROPERTIES LABELS "fortran;l2;validation;dynamic")
# Non-vacuity of the L2-2 OrcaFlex summary gate: the same WD0050/A1 cell with half the normal
# drag must be rejected by the summary comparison (the run prints the failing statistic).
add_test(NAME l2_2_dynamic_parity_perturbed_rejected COMMAND test_l2_2_dynamic_parity WD0050/A1 cdn=0.5)
set_tests_properties(l2_2_dynamic_parity_perturbed_rejected PROPERTIES
                       TIMEOUT 600
                       WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                       ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1"
                       PASS_REGULAR_EXPRESSION "MISMATCH \\[WD0050/A1:fairlead:orcaflex-summary-h1-"
                       LABELS "fortran;l2;validation;dynamic")

# Analytical-catenary arbiter self-validation: the canonical closed-form grounded-catenary
# static reference (submerged weight supplied explicitly, no buoyancy convention baked in) that
# backs the static parity gates. Runs the tool's self-check, which exits non-zero if the
# equations/constants regress past its gate -- so a stale arbiter fails loudly in the validation
# path instead of silently passing. Registered only when a Python interpreter is found (the tool
# needs numpy/scipy in the validation env).
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  # A Python interpreter selected from a Windows conda environment is not necessarily
  # launched through ``conda activate``. In that case extension modules cannot find the
  # DLLs in <prefix>/Library/bin (NumPy/OpenBLAS is the common failure) even though Python
  # itself starts normally. Keep that interpreter runtime ahead of the inherited PATH for
  # Python-driven tests only. Do not add CableDyn's deployment directory here: its OpenBLAS
  # DLL may otherwise pre-empt the BLAS runtime already loaded by NumPy.
  set(_cabledyn_python_test_environment "")
  set(_cabledyn_python_runtime_path "")
  if(WIN32)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import sys; print(sys.prefix)"
      RESULT_VARIABLE _cabledyn_python_prefix_rc
      OUTPUT_VARIABLE _cabledyn_python_prefix
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(_cabledyn_python_prefix_rc EQUAL 0 AND
       EXISTS "${_cabledyn_python_prefix}/Library/bin")
      set(_cabledyn_python_runtime_path
          "${_cabledyn_python_prefix}/Library/bin;$ENV{PATH}")
      string(REPLACE ";" "\\;" _cabledyn_python_runtime_path_escaped
             "${_cabledyn_python_runtime_path}")
      set(_cabledyn_python_test_environment
          "PATH=${_cabledyn_python_runtime_path_escaped}")
    endif()
  endif()

  # Public examples are executable input specifications. Keep their OPTIONS
  # rows in the documented OpenFAST value/key/comment/unit/choice layout.
  add_test(NAME example_deck_style
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_example_deck_style.py)
  set_tests_properties(example_deck_style PROPERTIES
    LABELS "python;documentation;examples"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Every OpenMP parallel region prepares its threads for the abnormal-end report.
  add_test(NAME omp_fatal_init
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_omp_fatal_init.py)
  set_tests_properties(omp_fatal_init PROPERTIES
    LABELS "python;unit;stack"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  add_test(NAME documentation_sources
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_documentation.py)
  set_tests_properties(documentation_sources PROPERTIES
    LABELS "python;documentation"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Audit the exact Git archive that would form a source release.
  add_test(NAME release_tree_audit
           COMMAND ${Python3_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/release/audit_release_tree.py
                   --tree ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(release_tree_audit PROPERTIES
    LABELS "python;documentation;release"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Solver-neutral cross-code comparison contract. This stdlib-only gate keeps
  # MoorDyn/OrcaFlex channel spelling at the adapter edge and rejects their
  # internal POINT/ROD discretization as CableDyn regression objects.
  add_test(NAME fem_observable_compare
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_fem_observable_compare.py)
  set_tests_properties(fem_observable_compare PROPERTIES
    LABELS "python;validation;cross-code;fem-observable"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  add_test(NAME validation_records
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_validation_records.py)
  set_tests_properties(validation_records PROPERTIES
    LABELS "python;validation;physical;performance"
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Finite-EI cable statics robustness matrix: sampled robustness decks (admissible,
  # inadmissible, tensile-rejected, metamorphic and at-rest twins) through the driver.
  # The per-deck bound is on CPU time, so the verdict does not depend on machine load;
  # PROCESSORS reserves the four driver processes the matrix runs at once, and the
  # wall-clock TIMEOUT only stops a hang. The Release bound of 150 s (Debug 1500 s) scales every
  # per-deck bound by 2.5 so that slower hosted CI runners pass; it still catches a pathological
  # slowdown.
  add_test(NAME cable_statics_matrix
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cable_statics_matrix.py
                   --exe $<TARGET_FILE:cabledyn>
                   --cases ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/cable_statics_matrix.json
                   --output ${CMAKE_BINARY_DIR}/cable_statics_matrix
                   --workers 4
                   --cpu-limit $<IF:$<CONFIG:Debug>,1500,150>)
  set_tests_properties(cable_statics_matrix PROPERTIES
    LABELS "python;hermite;static;robustness;slow"
    PROCESSORS 4
    TIMEOUT 1800
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Finite-EI laid cable in a current held by seabed friction (0.1-1.5 m/s x mu 0.1-1.0,
  # isotropic and anisotropic, three headings): converges, balances globally within the
  # Coulomb capacity and stays at rest when marched from the static state.
  add_test(NAME laid_cable_friction_current
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_laid_cable_friction_current.py
                   --exe $<TARGET_FILE:cabledyn>
                   --output ${CMAKE_BINARY_DIR}/laid_cable_friction_current
                   --workers 4)
  set_tests_properties(laid_cable_friction_current PROPERTIES
    LABELS "python;hermite;static;robustness;seabed"
    PROCESSORS 4
    TIMEOUT 900
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Finite-EI generalised-alpha blends: force-blend dtM convergence under prescribed surge.
  add_test(NAME hermite_force_blend_dt
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hermite_force_blend_dt.py
                   --exe $<TARGET_FILE:cabledyn>
                   --deck ${CMAKE_CURRENT_SOURCE_DIR}/examples/lozon_gomex80_power_cable_motion.dat
                   --output ${CMAKE_BINARY_DIR}/hermite_force_blend_dt)
  set_tests_properties(hermite_force_blend_dt PROPERTIES
    LABELS "python;hermite;dynamic;validation"
    TIMEOUT 600
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Finite-EI dynamic tensile monitor on the element-mean axial force: a tensile lazy wave
  # passes on a coarse mesh, a genuine compression is rejected on coarse and fine meshes.
  add_test(NAME hermite_tensile_monitor
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hermite_tensile_monitor.py
                   --exe $<TARGET_FILE:cabledyn>
                   --deck ${CMAKE_CURRENT_SOURCE_DIR}/examples/lozon_gomex80_power_cable_motion.dat
                   --output ${CMAKE_BINARY_DIR}/hermite_tensile_monitor)
  set_tests_properties(hermite_tensile_monitor PROPERTIES
    LABELS "python;hermite;dynamic;validation"
    TIMEOUT 600
    ENVIRONMENT "${_cabledyn_python_test_environment}")

  # Range graphs (LINES flag r) and touchdown channels on every standalone route: envelopes
  # against the node channels, the analytic catenary TDP, quasi-static tracking under a
  # moving fairlead, and agreement with cabledyn.touchdown_history (needs numpy).
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy"
    RESULT_VARIABLE _range_tdp_numpy_rc
    OUTPUT_QUIET ERROR_QUIET)
  if(_range_tdp_numpy_rc EQUAL 0)
    add_test(NAME range_graph_tdp
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_range_tdp.py
                     --exe $<TARGET_FILE:cabledyn>
                     --examples ${CMAKE_CURRENT_SOURCE_DIR}/examples
                     --output ${CMAKE_BINARY_DIR}/range_graph_tdp)
    set_tests_properties(range_graph_tdp PROPERTIES
      LABELS "python;integration;static;dynamic;validation"
      TIMEOUT 600
      ENVIRONMENT "${_cabledyn_python_test_environment}")
  else()
    message(STATUS "range_graph_tdp: skipped (numpy not importable in ${Python3_EXECUTABLE})")
  endif()

  # The arbiter imports numpy + scipy at load time; only register the test when they are
  # importable, so an image with python3 but without those packages does not fail the
  # validation suite with ModuleNotFoundError before any CableDyn gate runs.
  if(_cabledyn_python_runtime_path)
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E env "PATH=${_cabledyn_python_runtime_path}"
              "${Python3_EXECUTABLE}" -c "import numpy, scipy"
      RESULT_VARIABLE _arbiter_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  else()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import numpy, scipy"
      RESULT_VARIABLE _arbiter_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  endif()
  if(_arbiter_deps_rc EQUAL 0)
    add_test(NAME analytical_catenary_arbiter
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validation/scripts/analytical_catenary.py)
    set_tests_properties(analytical_catenary_arbiter PROPERTIES
      LABELS "python;validation;l2;static;arbiter"
      ENVIRONMENT "${_cabledyn_python_test_environment}")
  else()
    message(STATUS "analytical_catenary_arbiter: skipped (numpy/scipy not importable in ${Python3_EXECUTABLE})")
  endif()

  # cabledyn PYTHON PACKAGE gate: pytest over the ctypes binding (python/cabledyn) against
  # the built shared library -- the MoorDyn-C-style Python deployment surface. The suite
  # pins the binding to the SAME physics the Fortran gates certify (static fairlead pull,
  # fixed-point hold, deck-writer round-trip through the production parser), so a
  # prototype drift between CableDyn_CAPI and the binding fails as a physics mismatch,
  # not just a crash.
  if(_cabledyn_python_runtime_path)
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E env "PATH=${_cabledyn_python_runtime_path}"
              "${Python3_EXECUTABLE}" -c "import numpy, pytest"
      RESULT_VARIABLE _pybind_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  else()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import numpy, pytest"
      RESULT_VARIABLE _pybind_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  endif()
  if(_pybind_deps_rc EQUAL 0)
    add_test(NAME python_package
             COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/python/tests -q)
    if(WIN32 AND CABLEDYN_RUNTIME_DLLS)
      # Python 3.8+ deliberately ignores PATH for dependent-DLL resolution.
      # Exercise the self-contained deployment set whose shared library and GNU/OpenBLAS
      # runtimes are co-located, rather than the bare link product in the build root.
      set(_cabledyn_python_library
          "${CMAKE_BINARY_DIR}/runtime/$<TARGET_FILE_NAME:cabledyn_shared>")
    else()
      set(_cabledyn_python_library "$<TARGET_FILE:cabledyn_shared>")
    endif()
    set_tests_properties(python_package PROPERTIES
      LABELS "fortran;python;capi;integration"
      ENVIRONMENT
        "${_cabledyn_python_test_environment};CABLEDYN_LIBRARY=${_cabledyn_python_library};\
CABLEDYN_TEST_DRIVER=$<TARGET_FILE:cabledyn>;PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python")
  else()
    message(STATUS "python_package: skipped (numpy/pytest not importable in ${Python3_EXECUTABLE})")
  endif()
endif()

add_executable(test_cable_dynamic tests/test_cable_dynamic.f90)
target_link_libraries(test_cable_dynamic PRIVATE cabledyn_core)
add_test(NAME cable_dynamic COMMAND test_cable_dynamic)
set_tests_properties(cable_dynamic PROPERTIES LABELS "fortran;unit;reference;dynamic")

add_executable(test_cable_driver tests/test_cable_driver.f90)
target_link_libraries(test_cable_driver PRIVATE cabledyn_core)
add_test(NAME cable_driver COMMAND test_cable_driver)
set_tests_properties(cable_driver PROPERTIES LABELS "fortran;integration;static")

add_executable(test_line_driver tests/test_line_driver.f90)
target_link_libraries(test_line_driver PRIVATE cabledyn_core)
add_test(NAME line_driver COMMAND test_line_driver)
set_tests_properties(line_driver PROPERTIES LABELS "fortran;integration;static")

# Standalone analytical L1 gates scored on src/ against closed-form solutions (not
# the reference); `ctest -L l1` selects them. All ten L1 cases are gated on the two
# PRODUCTION paths: the EI=0 cable path carries L1-1 axial patch, L1-5 taut-string
# transverse vibration, L1-7 static catenary and L1-9 seabed reaction; the cubic-Hermite
# bending-cable path carries L1-2 bending/axial patch, L1-3 linear cantilever, L1-4
# large-deflection elastica, L1-6 Euler-beam vibration, L1-8 energy stability and
# L1-10 convergence order.
add_executable(test_l1_catenary tests/test_l1_catenary.f90)
target_link_libraries(test_l1_catenary PRIVATE cabledyn_core)
add_test(NAME l1_catenary COMMAND test_l1_catenary)
set_tests_properties(l1_catenary PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_axial_patch tests/test_l1_axial_patch.f90)
target_link_libraries(test_l1_axial_patch PRIVATE cabledyn_core)
add_test(NAME l1_axial_patch COMMAND test_l1_axial_patch)
set_tests_properties(l1_axial_patch PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_seabed tests/test_l1_seabed.f90)
target_link_libraries(test_l1_seabed PRIVATE cabledyn_core)
add_test(NAME l1_seabed COMMAND test_l1_seabed)
set_tests_properties(l1_seabed PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_string tests/test_l1_string.f90)
target_link_libraries(test_l1_string PRIVATE cabledyn_core)
add_test(NAME l1_string COMMAND test_l1_string)
set_tests_properties(l1_string PROPERTIES LABELS "fortran;l1;validation;dynamic")

# Finite-EI L1 gates on the PRODUCTION cubic-Hermite bending-cable path.
add_executable(test_l1_hermite_bending_patch tests/test_l1_hermite_bending_patch.f90)
target_link_libraries(test_l1_hermite_bending_patch PRIVATE cabledyn_core)
add_test(NAME l1_bending_patch COMMAND test_l1_hermite_bending_patch)
set_tests_properties(l1_bending_patch PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_hermite_cantilever tests/test_l1_hermite_cantilever.f90)
target_link_libraries(test_l1_hermite_cantilever PRIVATE cabledyn_core)
add_test(NAME l1_cantilever COMMAND test_l1_hermite_cantilever)
set_tests_properties(l1_cantilever PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_hermite_convergence tests/test_l1_hermite_convergence.f90)
target_link_libraries(test_l1_hermite_convergence PRIVATE cabledyn_core)
add_test(NAME l1_convergence_order COMMAND test_l1_hermite_convergence)
set_tests_properties(l1_convergence_order PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_hermite_elastica tests/test_l1_hermite_elastica.f90)
target_link_libraries(test_l1_hermite_elastica PRIVATE cabledyn_core)
add_test(NAME l1_large_deflection COMMAND test_l1_hermite_elastica)
set_tests_properties(l1_large_deflection PROPERTIES LABELS "fortran;l1;validation;static")

add_executable(test_l1_hermite_beam_vibration tests/test_l1_hermite_beam_vibration.f90)
target_link_libraries(test_l1_hermite_beam_vibration PRIVATE cabledyn_core)
add_test(NAME l1_euler_beam COMMAND test_l1_hermite_beam_vibration)
set_tests_properties(l1_euler_beam PROPERTIES LABELS "fortran;l1;validation;dynamic")

add_executable(test_l1_hermite_energy_stability tests/test_l1_hermite_energy_stability.f90)
target_link_libraries(test_l1_hermite_energy_stability PRIVATE cabledyn_core)
add_test(NAME l1_energy_stability COMMAND test_l1_hermite_energy_stability)
set_tests_properties(l1_energy_stability PROPERTIES LABELS "fortran;l1;validation;dynamic")

# The same six cases on the secondary Cosserat (rotation-DOF) path, kept as regression
# coverage for that element; `ctest -L cosserat` selects them.
add_executable(test_l1_bending_patch tests/test_l1_bending_patch.f90)
target_link_libraries(test_l1_bending_patch PRIVATE cabledyn_core)
add_test(NAME cosserat_bending_patch COMMAND test_l1_bending_patch)
set_tests_properties(cosserat_bending_patch PROPERTIES LABELS "fortran;cosserat;static")

# Linear isotropic rotational end connection: the spring is a configuration-dependent
# generalized load, so its force, consistent tangent, symmetry and the exact
# stretch-freedom property f.d == 0 are gated directly rather than inferred from a
# system-level result.
add_executable(test_endconn_spring tests/test_endconn_spring.f90)
target_link_libraries(test_endconn_spring PRIVATE cabledyn_core)
add_test(NAME endconn_spring COMMAND test_endconn_spring)
set_tests_properties(endconn_spring PROPERTIES LABELS "fortran;unit;endconnection;static")

# The same connection inside the finite-EI static solve: an analytic rigid-body rotation
# against the spring (a case that only has an isolated equilibrium BECAUSE the connection
# is there), plus the system-level stretch-freedom and pinned-limit statements.
add_executable(test_endconn_static tests/test_endconn_static.f90)
target_link_libraries(test_endconn_static PRIVATE cabledyn_core)
add_test(NAME endconn_static COMMAND test_endconn_static)
set_tests_properties(endconn_static PROPERTIES LABELS "fortran;endconnection;static;integration")

add_executable(test_l1_cantilever tests/test_l1_cantilever.f90)
target_link_libraries(test_l1_cantilever PRIVATE cabledyn_core)
add_test(NAME cosserat_cantilever COMMAND test_l1_cantilever)
set_tests_properties(cosserat_cantilever PROPERTIES LABELS "fortran;cosserat;static")

add_executable(test_l1_convergence_order tests/test_l1_convergence_order.f90)
target_link_libraries(test_l1_convergence_order PRIVATE cabledyn_core)
add_test(NAME cosserat_convergence_order COMMAND test_l1_convergence_order)
set_tests_properties(cosserat_convergence_order PROPERTIES LABELS "fortran;cosserat;static")

add_executable(test_l1_large_deflection tests/test_l1_large_deflection.f90)
target_link_libraries(test_l1_large_deflection PRIVATE cabledyn_core)
add_test(NAME cosserat_large_deflection COMMAND test_l1_large_deflection)
set_tests_properties(cosserat_large_deflection PROPERTIES LABELS "fortran;cosserat;static")

add_executable(test_l1_euler_beam tests/test_l1_euler_beam.f90)
target_link_libraries(test_l1_euler_beam PRIVATE cabledyn_core)
add_test(NAME cosserat_euler_beam COMMAND test_l1_euler_beam)
set_tests_properties(cosserat_euler_beam PROPERTIES LABELS "fortran;cosserat;dynamic")

add_executable(test_l1_energy_stability tests/test_l1_energy_stability.f90)
target_link_libraries(test_l1_energy_stability PRIVATE cabledyn_core)
add_test(NAME cosserat_energy_stability COMMAND test_l1_energy_stability)
set_tests_properties(cosserat_energy_stability PROPERTIES LABELS "fortran;cosserat;dynamic")

add_executable(test_l3_wave_smoke tests/test_l3_wave_smoke.f90)
target_link_libraries(test_l3_wave_smoke PRIVATE cabledyn_core)
add_test(NAME l3_wave_smoke COMMAND test_l3_wave_smoke)
set_tests_properties(l3_wave_smoke PROPERTIES LABELS "fortran;l3;dynamic;hydro")

# --- the standalone deck driver (cabledyn <deck.dat> <out_root>), released as CableDyn_driver ---
add_executable(cabledyn app/cabledyn.f90)
target_link_libraries(cabledyn PRIVATE cabledyn_core)
set(CABLEDYN_DRIVER_UTF8_MANIFEST OFF)
if(WIN32 AND (MSVC OR CMAKE_Fortran_SIMULATE_ID STREQUAL "MSVC"))
  # MSVC-style linkers embed a .manifest source, merged with the linker's own manifest.
  # This one makes UTF-8 the active code page, so the Fortran runtime's narrow file calls
  # open a working folder or file whose name is outside the system ANSI code page. The
  # MinGW runtime of the GNU build reaches such names through src/cabledyn_path.c.
  target_sources(cabledyn PRIVATE app/utf8_code_page.manifest)
  set(CABLEDYN_DRIVER_UTF8_MANIFEST ON)
endif()
if(CABLEDYN_DRIVER_UTF8_MANIFEST)
  # With an MSVC-style linker the driver and the shared library (OUTPUT_NAME cabledyn)
  # would write cabledyn.pdb and cabledyn.ilk to one directory; the product name keeps
  # them apart.
  set_target_properties(cabledyn PROPERTIES OUTPUT_NAME "CableDyn_driver")
endif()
if(CABLEDYN_STATIC_WINDOWS_RELEASE)
  # A 952-element production lazy-wave deck needs more than the PE default 1 MB
  # stack before IFX can enter its first dynamic step. Reserving address space is
  # cheap on x64 and matches the large-model contract of the standalone product.
  target_link_options(cabledyn PRIVATE "LINKER:/STACK:268435456")
endif()
install(TARGETS cabledyn RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(WIN32 AND CABLEDYN_RUNTIME_DLLS)
  # A conventional GNU/conda Windows install must be usable outside the build
  # environment, not only from the separately assembled release bundle.
  install(FILES ${CABLEDYN_RUNTIME_DLLS} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

add_executable(test_deck_driver tests/test_deck_driver.f90)
target_link_libraries(test_deck_driver PRIVATE cabledyn_core)
add_test(NAME deck_driver COMMAND test_deck_driver)
set_tests_properties(deck_driver PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_executable(test_stream_wave tests/test_stream_wave.f90)
target_link_libraries(test_stream_wave PRIVATE cabledyn_core)
add_test(NAME stream_wave COMMAND test_stream_wave)
set_tests_properties(stream_wave PROPERTIES LABELS "fortran;unit;driver;dynamic;validation;orcaflex")

add_executable(test_modal tests/test_modal.f90 tests/blas_config_probe.c)
target_link_libraries(test_modal PRIVATE cabledyn_core)
add_test(NAME modal COMMAND test_modal)
set_tests_properties(modal PROPERTIES LABELS "fortran;unit;driver;static;validation")

# Stack bound (GNU Windows builds): the modal test and the driver relinked with a 1 MiB stack
# reserve, half the MinGW default, instead of the 64 MiB reserve of the other executables. A
# large local array (as the dense matrices test_modal once kept on the stack) or a library
# frame that grows with the model overflows here on every machine. The driver runs a
# 2048-element finite-EI cable. Linux and macOS executables take their main-thread stack
# from the environment (ulimit -s), so these links are Windows-only; test_c_api_stack covers
# the library on every platform.
if(CABLEDYN_GNU_WINDOWS_STACK)
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/stack_small")
  add_executable(test_modal_small_stack tests/test_modal.f90 tests/blas_config_probe.c)
  target_link_libraries(test_modal_small_stack PRIVATE cabledyn_core)
  add_executable(cabledyn_small_stack app/cabledyn.f90)
  target_link_libraries(cabledyn_small_stack PRIVATE cabledyn_core)
  foreach(_cabledyn_small IN ITEMS test_modal_small_stack cabledyn_small_stack)
    # Appended after the 64 MiB directory-wide option; the last --stack wins (checked below).
    target_link_options(${_cabledyn_small} PRIVATE "-Wl,--stack,${CABLEDYN_SMALL_STACK_RESERVE}")
  endforeach()
  add_test(NAME stack_modal COMMAND test_modal_small_stack
           WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/stack_small")
  set_tests_properties(stack_modal PROPERTIES LABELS "fortran;unit;driver;static;stack")
  add_test(NAME stack_driver_large_deck
           COMMAND cabledyn_small_stack
                   "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/stack_humboldt800_2048_elements.dat"
                   "${CMAKE_BINARY_DIR}/stack_small/stack_driver_large_deck")
  set_tests_properties(stack_driver_large_deck PROPERTIES
                       LABELS "fortran;driver;integration;dynamic;stack" TIMEOUT 900)
  # The reserve each link actually received, read from the PE headers: <binary>=<bytes>.
  string(JOIN "|" _cabledyn_stack_checks
         "$<TARGET_FILE:cabledyn>=${CABLEDYN_EXE_STACK_RESERVE}"
         "$<TARGET_FILE:test_modal>=${CABLEDYN_EXE_STACK_RESERVE}"
         "$<TARGET_FILE:test_c_api_stack>=${CABLEDYN_EXE_STACK_RESERVE}"
         "$<TARGET_FILE:test_modal_small_stack>=${CABLEDYN_SMALL_STACK_RESERVE}"
         "$<TARGET_FILE:cabledyn_small_stack>=${CABLEDYN_SMALL_STACK_RESERVE}")
  add_test(NAME stack_reserve
           COMMAND ${CMAKE_COMMAND} "-DCHECKS=${_cabledyn_stack_checks}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_stack_reserve.cmake)
  set_tests_properties(stack_reserve PROPERTIES LABELS "unit;driver;stack")
  # The check itself must see a reserve that differs from the expected one.
  add_test(NAME stack_reserve_mismatch
           COMMAND ${CMAKE_COMMAND}
                   "-DCHECKS=$<TARGET_FILE:test_modal_small_stack>=${CABLEDYN_EXE_STACK_RESERVE}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_stack_reserve.cmake)
  set_tests_properties(stack_reserve_mismatch PROPERTIES LABELS "unit;driver;stack"
                       PASS_REGULAR_EXPRESSION
                       "stack reserve ${CABLEDYN_SMALL_STACK_RESERVE} bytes, expected ${CABLEDYN_EXE_STACK_RESERVE}")
endif()

add_executable(test_seabed_friction_aniso tests/test_seabed_friction_aniso.f90)
target_link_libraries(test_seabed_friction_aniso PRIVATE cabledyn_core)
add_test(NAME seabed_friction_aniso COMMAND test_seabed_friction_aniso)
set_tests_properties(seabed_friction_aniso PROPERTIES LABELS "fortran;unit;driver;dynamic;seabed")

add_executable(test_wave_spectra tests/test_wave_spectra.f90)
target_link_libraries(test_wave_spectra PRIVATE cabledyn_core)
add_test(NAME wave_spectra COMMAND test_wave_spectra)
set_tests_properties(wave_spectra PROPERTIES LABELS "fortran;unit;driver;dynamic;validation")

add_executable(test_dynamics_regressions tests/test_dynamics_regressions.f90)
target_link_libraries(test_dynamics_regressions PRIVATE cabledyn_core)
add_test(NAME dynamics_regressions COMMAND test_dynamics_regressions)
set_tests_properties(dynamics_regressions PROPERTIES LABELS "fortran;driver;integration;dynamic")

add_executable(test_deck_reader_io tests/test_deck_reader_io.f90)
target_link_libraries(test_deck_reader_io PRIVATE cabledyn_core)
add_test(NAME deck_reader_io COMMAND test_deck_reader_io)
set_tests_properties(deck_reader_io PROPERTIES LABELS "fortran;driver;integration;static")

add_executable(test_body_rod_physics tests/test_body_rod_physics.f90)
target_link_libraries(test_body_rod_physics PRIVATE cabledyn_core)
add_test(NAME body_rod_physics COMMAND test_body_rod_physics)
set_tests_properties(body_rod_physics PROPERTIES LABELS "fortran;driver;integration;dynamic")
add_executable(test_monolithic_gates tests/test_monolithic_gates.f90)
target_link_libraries(test_monolithic_gates PRIVATE cabledyn_core)
add_test(NAME monolithic_gates COMMAND test_monolithic_gates)
set_tests_properties(monolithic_gates PROPERTIES LABELS "fortran;driver;integration;dynamic")

# MoorDyn-C WaveKin 3/7 and Currents 1 vs MoorDyn-C v2.7.1 on shared inputs (the reference rows
# come from validation/scripts/wavekin_mdc_probe.cpp)
file(GLOB CD_WAVEKIN_MDC_FILES RELATIVE ${CMAKE_SOURCE_DIR}/tests/data/wavekin_mdc
     ${CMAKE_SOURCE_DIR}/tests/data/wavekin_mdc/*)
foreach(f IN LISTS CD_WAVEKIN_MDC_FILES)
  configure_file(tests/data/wavekin_mdc/${f} ${CMAKE_BINARY_DIR}/wavekin_mdc/${f} COPYONLY)
endforeach()
add_executable(test_wavekin_modes tests/test_wavekin_modes.f90)
target_link_libraries(test_wavekin_modes PRIVATE cabledyn_core)
add_test(NAME wavekin_modes COMMAND test_wavekin_modes WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(wavekin_modes PROPERTIES LABELS "fortran;driver;validation;dynamic")

add_executable(test_clamped_cable_ends tests/test_clamped_cable_ends.f90)
target_link_libraries(test_clamped_cable_ends PRIVATE cabledyn_core)
add_test(NAME clamped_cable_ends COMMAND test_clamped_cable_ends)
set_tests_properties(clamped_cable_ends PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_executable(test_torsion_deck tests/test_torsion_deck.f90)
target_link_libraries(test_torsion_deck PRIVATE cabledyn_core)
add_test(NAME torsion_deck COMMAND test_torsion_deck)
set_tests_properties(torsion_deck PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_executable(test_torsion_orcaflex tests/test_torsion_orcaflex.f90)
target_link_libraries(test_torsion_orcaflex PRIVATE cabledyn_core)
add_test(NAME torsion_orcaflex
         COMMAND test_torsion_orcaflex ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/torsion_orcaflex_refs.txt)
set_tests_properties(torsion_orcaflex PROPERTIES LABELS "fortran;driver;integration;static;validation;l3")

add_executable(test_rod_hydro tests/test_rod_hydro.f90)
target_link_libraries(test_rod_hydro PRIVATE cabledyn_core)
add_test(NAME rod_hydro COMMAND test_rod_hydro)
set_tests_properties(rod_hydro PROPERTIES LABELS "fortran;unit;dynamic;static")

add_executable(test_body_hydro tests/test_body_hydro.f90)
target_link_libraries(test_body_hydro PRIVATE cabledyn_core)
add_test(NAME body_hydro COMMAND test_body_hydro)
set_tests_properties(body_hydro PROPERTIES LABELS "fortran;unit;dynamic;static")

add_executable(test_shared_anchor tests/test_shared_anchor.f90)
target_link_libraries(test_shared_anchor PRIVATE cabledyn_core)
add_test(NAME shared_anchor COMMAND test_shared_anchor)
set_tests_properties(shared_anchor PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_executable(test_bodies_rods tests/test_bodies_rods.f90)
target_link_libraries(test_bodies_rods PRIVATE cabledyn_core)
add_test(NAME bodies_rods COMMAND test_bodies_rods)
set_tests_properties(bodies_rods PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_test(
  NAME driver_nonconvergence
  COMMAND ${CMAKE_COMMAND}
    -DDRIVER=$<TARGET_FILE:cabledyn>
    -DDECK=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonconvergent_motion.dat
    -DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/driver_nonconvergence
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_nonconvergence.cmake)
set_tests_properties(driver_nonconvergence PROPERTIES
  LABELS "fortran;driver;integration;dynamic;failure")

# libgfortran's per-statement setlocale restore reads freed UCRT memory on MinGW-w64
# unless src/cabledyn_crt_locale.c redirects its import. The unit check requires the
# redirect there and releases the freed pages mid-WRITE, which faults deterministically
# without the guard; the repeat run requires a clean exit from every process.
add_executable(test_crt_locale_guard tests/test_crt_locale_guard.f90 tests/crt_locale_guard_heap.c)
target_link_libraries(test_crt_locale_guard PRIVATE cabledyn_core)
add_test(NAME crt_locale_guard COMMAND test_crt_locale_guard $<IF:$<BOOL:${MINGW}>,1,0>)
set_tests_properties(crt_locale_guard PROPERTIES LABELS "fortran;unit")
# An abnormal end of the driver process is never silent: a stack overflow (on a 1 MiB stack
# in Windows GNU builds), an access violation and, off Windows, SIGTERM each end the process
# with a report on stderr of the cause and the simulated time reached; the march records the
# time of every committed step for that report.
add_executable(test_fatal_report tests/test_fatal_report.f90 tests/fatal_report_faults.c)
target_link_libraries(test_fatal_report PRIVATE cabledyn_core)
if(CABLEDYN_GNU_WINDOWS_STACK)
  target_link_options(test_fatal_report PRIVATE "-Wl,--stack,${CABLEDYN_SMALL_STACK_RESERVE}")
endif()
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/fatal_report")
add_test(NAME fatal_report_march
         COMMAND test_fatal_report march ${CMAKE_CURRENT_SOURCE_DIR}/examples/dynamic_chain_held.dat
                 ${CMAKE_BINARY_DIR}/fatal_report/march)
set_tests_properties(fatal_report_march PROPERTIES LABELS "fortran;unit;driver;dynamic")
# The cause names the fault as each platform delivers it (macOS may deliver SIGBUS).
foreach(_fatal IN ITEMS "overflow|during|(stack overflow|segmentation fault|bus error)"
                        "overflow|before|(stack overflow|segmentation fault|bus error)"
                        "null|during|(access violation|segmentation fault|bus error)"
                        "term|during|a termination request")
  string(REPLACE "|" ";" _fatal_case "${_fatal}")
  list(GET _fatal_case 0 _fatal_mode)
  list(GET _fatal_case 1 _fatal_when)
  list(SUBLIST _fatal_case 2 -1 _fatal_cause)
  string(JOIN "|" _fatal_cause ${_fatal_cause})
  add_test(NAME fatal_report_${_fatal_mode}_${_fatal_when}
           COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:test_fatal_report> -DMODE=${_fatal_mode}
                   -DWHEN=${_fatal_when} "-DCAUSE=${_fatal_cause}"
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_fatal_report.cmake)
  set_tests_properties(fatal_report_${_fatal_mode}_${_fatal_when} PROPERTIES
                       LABELS "fortran;unit;driver;stack" TIMEOUT 120)
endforeach()
# A stack overflow on an OpenMP worker thread with a 1 MiB stack is reported too: every
# thread that enters a parallel region gets room for the handlers (tests/test_omp_fatal_init.py
# checks that every region of the library and the driver asks for it).
add_test(NAME fatal_report_omp_overflow_during
         COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:test_fatal_report> -DMODE=omp_overflow
                 -DWHEN=during "-DCAUSE=(stack overflow|segmentation fault|bus error)"
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_fatal_report.cmake)
set_tests_properties(fatal_report_omp_overflow_during PROPERTIES
                     LABELS "fortran;unit;driver;stack" TIMEOUT 120
                     ENVIRONMENT "OMP_NUM_THREADS=2;OMP_STACKSIZE=1M")
# POSIX: a hardware fault keeps its own siginfo and context through the report.
add_test(NAME fatal_report_fault_context COMMAND test_fatal_report context)
set_tests_properties(fatal_report_fault_context PROPERTIES
                     LABELS "fortran;unit;driver;stack" TIMEOUT 120
                     FAIL_REGULAR_EXPRESSION "FAIL:")
# POSIX: the report keeps every previous disposition -- an ignored signal (with or without
# SA_SIGINFO) stays ignored, a default one still ends the process, a handler still runs.
add_test(NAME fatal_report_altstack COMMAND test_fatal_report altstack)
set_tests_properties(fatal_report_altstack PROPERTIES
                     LABELS "fortran;unit;driver;stack" TIMEOUT 120
                     FAIL_REGULAR_EXPRESSION "FAIL:")
add_test(NAME fatal_report_dispositions COMMAND test_fatal_report dispositions)
set_tests_properties(fatal_report_dispositions PROPERTIES
                     LABELS "fortran;unit;driver" TIMEOUT 120
                     FAIL_REGULAR_EXPRESSION "FAIL:")
add_test(
  NAME driver_repeat_exit
  COMMAND ${CMAKE_COMMAND}
    -DDRIVER=$<TARGET_FILE:cabledyn>
    -DDECK=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/three_line_composite_tmax0.dat
    -DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/driver_repeat_exit
    -DRUNS=200
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_repeat_exit.cmake)
set_tests_properties(driver_repeat_exit PROPERTIES
  LABELS "fortran;driver;integration;static")

# Command-line contract: unknown options, an output root whose .out would overwrite
# the deck, and an unwritable output directory all fail before any solve.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonconvergent_motion.dat
               ${CMAKE_CURRENT_BINARY_DIR}/driver_cli_self.out COPYONLY)
# Each case checks the exit code (1) and the diagnostic together.
add_test(NAME driver_cli_unknown_option
  COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:cabledyn> -DEXPECT_RC=1
    "-DARGS=--no-such-option|${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonconvergent_motion.dat|cli_out"
    "-DMUST_MATCH=unknown option \"--no-such-option\""
    -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/driver_cli_unknown_option
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_output.cmake)
add_test(NAME driver_cli_deck_overwrite
  COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:cabledyn> -DEXPECT_RC=1
    "-DARGS=${CMAKE_CURRENT_BINARY_DIR}/driver_cli_self.out|${CMAKE_CURRENT_BINARY_DIR}/driver_cli_self"
    "-DMUST_MATCH=would overwrite the input deck"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_output.cmake)
add_test(NAME driver_cli_missing_output_dir
  COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:cabledyn> -DEXPECT_RC=1
    "-DARGS=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/nonconvergent_motion.dat|${CMAKE_CURRENT_BINARY_DIR}/driver_cli_no_such_dir/run"
    "-DMUST_MATCH=cannot write output files"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_output.cmake)
set_tests_properties(driver_cli_unknown_option driver_cli_deck_overwrite driver_cli_missing_output_dir
  PROPERTIES LABELS "fortran;driver;cli")
# The mixed EI=0/finite-EI standalone route writes no modal files: OPTION nModes fails closed.
add_test(NAME driver_cli_mixed_nmodes
  COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:cabledyn> -DEXPECT_RC=1
    "-DARGS=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/mixed_nmodes.dat|${CMAKE_CURRENT_BINARY_DIR}/driver_cli_mixed_nmodes"
    "-DMUST_MATCH=nModes is not supported for mixed EI=0/finite-EI decks"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_output.cmake)
set_tests_properties(driver_cli_mixed_nmodes PROPERTIES LABELS "fortran;driver;cli")
# Non-ASCII, reserved-device, input-clobbering, and concurrently used output paths.
add_test(
  NAME driver_cli_paths
  COMMAND ${CMAKE_COMMAND}
    -DDRIVER=$<TARGET_FILE:cabledyn>
    -DSTATIC_DECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/chain_catenary_r3_100m.dat
    -DDYNAMIC_DECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/dynamic_chain_held.dat
    -DUTF8_MANIFEST=${CABLEDYN_DRIVER_UTF8_MANIFEST}
    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/driver_cli_paths
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_paths.cmake)
set_tests_properties(driver_cli_paths PROPERTIES LABELS "fortran;driver;cli")

# Advisories that do not stop a run: the start-up shock note of a prescribed motion that
# starts moving, and the no-seabed warning of an OpenFAST MooringFile run standalone.
add_test(NAME driver_advisories
  COMMAND ${CMAKE_COMMAND}
    -DDRIVER=$<TARGET_FILE:cabledyn>
    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/driver_advisories
    -DNO_SEABED_DECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/iea15mw_umaine_openfast_cabledyn.dat
    -DSEABED_DECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/chain_catenary_shallow_30m.dat
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_advisories.cmake)
set_tests_properties(driver_advisories PROPERTIES LABELS "fortran;driver;integration")
# A clamped finite-EI cable with a buoyant section on a Rigid6 body: the body static solve
# holds the objects at their deck poses, warns once, and reports no balanced clamped end loads.
add_test(NAME driver_clamped_buoyant_held
  COMMAND ${CMAKE_COMMAND} -DDRIVER=$<TARGET_FILE:cabledyn> -DEXPECT_RC=0
    "-DARGS=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/clamped_buoyant_cable_held.dat|held"
    "-DONCE=held at their deck poses"
    "-DMUST_NOT=balanced in"
    -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/driver_clamped_buoyant_held
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_driver_output.cmake)
set_tests_properties(driver_clamped_buoyant_held PROPERTIES LABELS "fortran;driver;integration;static")

if(WIN32 AND CABLEDYN_RUNTIME_DLLS)
  add_test(
    NAME windows_install_runtime
    COMMAND ${CMAKE_COMMAND}
      -DBUILD_DIR=${CMAKE_BINARY_DIR}
      -DINSTALL_PREFIX=${CMAKE_BINARY_DIR}/install-smoke
      -DEXE_NAME=$<TARGET_FILE_NAME:cabledyn>
      -DDECK=${CMAKE_CURRENT_SOURCE_DIR}/examples/chain_catenary_r3_100m.dat
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_windows_install_runtime.cmake)
  set_tests_properties(windows_install_runtime PROPERTIES
    LABELS "fortran;driver;integration;installation;windows")
endif()

# Maintained reference standalone decks are executable specifications. These
# smoke gates catch parser drift, missing auxiliary motion files, and loss of
# the intended lazy-wave equilibrium branch.
foreach(_lozon_case IN ITEMS
    lozon_gomex80_power_cable
    lozon_gomaine200_power_cable
    lozon_humboldt800_power_cable
    lozon_gomex80_mooring
    lozon_gomaine200_mooring
    lozon_humboldt800_mooring
    lozon_gomex80_power_cable_motion)
  add_test(
    NAME example_${_lozon_case}
    COMMAND cabledyn
      "${CMAKE_CURRENT_SOURCE_DIR}/examples/${_lozon_case}.dat"
      "${CMAKE_CURRENT_BINARY_DIR}/example_${_lozon_case}")
  set_tests_properties(example_${_lozon_case} PROPERTIES
    LABELS "fortran;driver;integration;example;lozon")
endforeach()

# Torsion example: a lazy-wave cable twisted two turns at the hang-off through the motionFile
# roll column (condensed torsion, quasi-static).
add_test(
  NAME example_torsion_lazy_wave_hangoff_twist
  COMMAND cabledyn
    "${CMAKE_CURRENT_SOURCE_DIR}/examples/torsion_lazy_wave_hangoff_twist.dat"
    "${CMAKE_CURRENT_BINARY_DIR}/example_torsion_lazy_wave_hangoff_twist")
set_tests_properties(example_torsion_lazy_wave_hangoff_twist PROPERTIES
  LABELS "fortran;driver;integration;example;dynamic")

# The teaching/reference deck enumerates all active defaults and documents every
# mutually exclusive alternative.  Keep it runnable so the canonical OPTIONS
# example cannot drift away from the production parser.
add_test(
  NAME example_cabledyn_options_reference
  COMMAND cabledyn
    "${CMAKE_CURRENT_SOURCE_DIR}/examples/cabledyn_options_reference.dat"
    "${CMAKE_CURRENT_BINARY_DIR}/example_cabledyn_options_reference")
set_tests_properties(example_cabledyn_options_reference PROPERTIES
  LABELS "fortran;driver;integration;example;options")

if(Python3_Interpreter_FOUND)
  add_test(
    NAME lozon_example_equilibrium
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_lozon_example_equilibrium.py
      --exe $<TARGET_FILE:cabledyn>
      --source ${CMAKE_CURRENT_SOURCE_DIR}
      --output ${CMAKE_CURRENT_BINARY_DIR}/lozon_example_equilibrium)
  set_tests_properties(lozon_example_equilibrium PROPERTIES
    LABELS "python;driver;integration;example;lozon;validation")
  # Static-initialization robustness matrix: a deterministic sample of EI=0 mooring regimes
  # (catenary, semi-taut, taut, suspended, composite, slack/vertical, anchors off the seabed,
  # sloped seabed, fine mesh, soft seabed, current, one/two segments, Free/Connect points)
  # must converge to a balanced tension-only equilibrium that matches the elastic catenary.
  # Needs numpy + scipy in the interpreter (the same probe as the catenary arbiter).
  if(_cabledyn_python_runtime_path)
    execute_process(
      COMMAND "${CMAKE_COMMAND}" -E env "PATH=${_cabledyn_python_runtime_path}"
              "${Python3_EXECUTABLE}" -c "import numpy, scipy"
      RESULT_VARIABLE _statics_matrix_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  else()
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import numpy, scipy"
      RESULT_VARIABLE _statics_matrix_deps_rc
      OUTPUT_QUIET ERROR_QUIET)
  endif()
  if(_statics_matrix_deps_rc EQUAL 0)
    add_test(
      NAME statics_robustness_matrix
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_statics_robustness.py
        --exe $<TARGET_FILE:cabledyn>
        --output ${CMAKE_CURRENT_BINARY_DIR}/statics_robustness_matrix)
    # PROCESSORS reserves the four driver processes the matrix runs at once (--jobs 4); the
    # TIMEOUT only stops a hang.
    set_tests_properties(statics_robustness_matrix PROPERTIES
      LABELS "python;driver;integration;static;validation;slow"
      PROCESSORS 4
      TIMEOUT 3600
      ENVIRONMENT "${_cabledyn_python_test_environment}")
  else()
    message(STATUS "statics_robustness_matrix: skipped (numpy/scipy not importable in ${Python3_EXECUTABLE})")
  endif()
endif()

add_executable(test_driver_channels tests/test_driver_channels.f90)
target_link_libraries(test_driver_channels PRIVATE cabledyn_core)
add_test(NAME driver_channels COMMAND test_driver_channels)
set_tests_properties(driver_channels PROPERTIES LABELS "fortran;driver;integration;static;dynamic")

add_executable(test_l3_orcaflex_wave_dynamic tests/test_l3_orcaflex_wave_dynamic.f90)
target_link_libraries(test_l3_orcaflex_wave_dynamic PRIVATE cabledyn_core)
add_test(NAME l3_orcaflex_wave_dynamic COMMAND test_l3_orcaflex_wave_dynamic)
set_tests_properties(l3_orcaflex_wave_dynamic PROPERTIES
  LABELS "fortran;l3;validation;dynamic;hydro;orcaflex" TIMEOUT 180)

# PORTABLE-COUNTER performance gate on the production EI=0 aggregate path: 60 s of the
# module cadence on the committed VolturnUS-S deck, gating COUNTS (Newton/step bands,
# zero substep fallbacks, finite loads) rather than machine-bound wall-clock seconds.
add_executable(test_of_mooring_counters tests/test_of_mooring_counters.f90)
target_link_libraries(test_of_mooring_counters PRIVATE cabledyn_core)
add_test(NAME openfast_mooring_counters
         COMMAND test_of_mooring_counters ${CMAKE_SOURCE_DIR}/examples/iea15mw_volturnus_mooring.dat
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(openfast_mooring_counters PROPERTIES
                     LABELS "fortran;openfast;perf;dynamic")

# PORTABLE-COUNTER mesh-scaling gate on the standalone EI=0 dynamic route: the committed
# held chain re-meshed to 82 and 820 segments, gating COUNTS (zero temporal subdivisions,
# Newton/step bands that do not grow with the mesh, no dense n_dof**2 mass workspace)
# rather than wall-clock seconds. A solver that cannot recognise the fine-mesh residual
# round-off floor re-solves every step by subdivision and trips it.
add_executable(test_ei0_mesh_scaling_counters tests/test_ei0_mesh_scaling_counters.f90)
target_link_libraries(test_ei0_mesh_scaling_counters PRIVATE cabledyn_core)
add_test(NAME ei0_mesh_scaling_counters
         COMMAND test_ei0_mesh_scaling_counters ${CMAKE_SOURCE_DIR}/examples/dynamic_chain_held.dat
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(ei0_mesh_scaling_counters PROPERTIES
                     LABELS "fortran;perf;dynamic" TIMEOUT 300)

if(CABLEDYN_RUNTIME_DLLS)
  foreach(_cabledyn_target IN ITEMS
      cabledyn
      cabledyn_shared)
    if(TARGET ${_cabledyn_target})
      cabledyn_stage_runtime_dlls(${_cabledyn_target})
    endif()
  endforeach()
endif()

# Static-dynamic consistency for the end connection: an equilibrium of the static solve
# must persist under a dynamic march, with a control run that MUST drift so the gate is
# demonstrably sensitive to the contribution being present.
add_executable(test_endconn_dynamic tests/test_endconn_dynamic.f90)
target_link_libraries(test_endconn_dynamic PRIVATE cabledyn_core)
add_test(NAME endconn_dynamic COMMAND test_endconn_dynamic)
set_tests_properties(endconn_dynamic PROPERTIES LABELS "fortran;endconnection;dynamic;integration")

# These regression programs deliberately assert bit identity and exact-zero limits. GNU's
# -Wcompare-reals warning is useful in the solver core but not for those stated test contracts.
# GCC 15 also diagnoses allocatable descriptors nested in a default-initialized aggregate as
# uninitialized even though Fortran initializes every allocatable component as unallocated.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  foreach(_cabledyn_exact_real_test IN ITEMS
      test_hermite_cable
      test_hermite_drag_dynamics
      test_hermite_axial_damping
      test_endconn_spring)
    target_compile_options(${_cabledyn_exact_real_test} PRIVATE -Wno-compare-reals)
  endforeach()
  target_compile_options(test_openfast_aggregate PRIVATE
                         -Wno-compare-reals
                         -Wno-uninitialized
                         -Wno-maybe-uninitialized)
  target_compile_options(test_torsion_deck PRIVATE -Wno-uninitialized -Wno-maybe-uninitialized)
endif()

if(CABLEDYN_ENABLE_COVERAGE)
  find_program(CABLEDYN_GCOVR gcovr)
  if(CABLEDYN_GCOVR)
    set(CABLEDYN_FORTRAN_COVERAGE_TARGETS
      test_ad
      test_linalg
      test_so3
      test_conventions
      test_rigid_kinematics
      test_emc
      test_hermite_arch
      test_cosserat_strains
      test_cosserat_force
      test_cosserat_assemble
      test_cosserat_mass
      test_cosserat_dynamic
      test_finite_ei_large_rotation
      test_finite_ei_precession
      test_finite_ei_gyroscopic
      test_finite_ei_mult_tangent
      test_finite_ei_mult_step
      test_finite_ei_gyro_primitives
      test_finite_ei_gyro_tangent
      test_finite_ei_gyro_dynamics
      test_finite_ei_model
      test_cable_mesh
      test_cable_elem
      test_cable_assemble
      test_cable_loads
      test_bathymetry
      test_hydro
      test_cable_static
      test_catenary
      test_robust_init
      test_line
      test_hermite_cable
      test_hermite_torsion_kernel
      test_hermite_torsion_static
      test_hermite_torsion_dynamic
      test_hermite_cable_static
      test_hermite_sequenced
      test_hermite_resolution
      test_hermite_automesh
      test_hermite_cable_dynamic
      test_hermite_drag_dynamics
      test_hermite_axial_damping
      test_hermite_added_mass
      test_hermite_wave_forcing
      test_lazywave_dynamics
      test_l3_lazywave_dynamic
      test_l3_longduration_mooring
      test_hermite_axial_drag_arbiter
      test_l2_chain_parity
      test_l2_volturnus_mooring
      test_l2_lozon80_mooring
      test_l2_lozon200_semitaut
      test_l2_lozon200_semitaut_deck
      test_l3_orcaflex_static
      test_l3_orcaflex_current_static
      test_l3_orcaflex_finite_ei_static
      test_l3_lazywave_curvature
      test_l3_lazywave_refinement
      test_l3_lazywave_fullcable
      test_l4_platform_mooring
      test_damping
      test_viscoelastic
      test_syrope
      test_prescribed_motion
      test_model
      test_system
      test_openfast_shell
      test_hermite_fmf
      test_openfast_aggregate
      test_of_mooring_counters
      test_ei0_mesh_scaling_counters
      test_c_api
      test_cable_dynamic
      test_cable_driver
      test_line_driver
      test_l1_catenary
      test_l1_axial_patch
      test_l1_seabed
      test_l1_string
      test_l1_hermite_bending_patch
      test_l1_hermite_cantilever
      test_l1_hermite_convergence
      test_l1_hermite_elastica
      test_l1_hermite_beam_vibration
      test_l1_hermite_energy_stability
      test_l3_hermite_touchdown
      test_l1_bending_patch
      test_l1_cantilever
      test_l1_convergence_order
      test_l1_large_deflection
      test_l1_euler_beam
      test_l1_energy_stability
      test_l3_wave_smoke
      test_deck_driver
      test_l3_orcaflex_wave_dynamic
      test_snap_load
      test_endconn_spring
      test_endconn_static
      test_endconn_dynamic
    )
    add_custom_target(fortran_coverage
      COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
              -L fortran -LE slow --output-on-failure
      COMMAND ${CABLEDYN_GCOVR}
              --root ${CMAKE_SOURCE_DIR}
              --filter ${CMAKE_SOURCE_DIR}/src/
              --exclude ${CMAKE_SOURCE_DIR}/tests/
              --exclude ${CMAKE_SOURCE_DIR}/app/
              --gcov-ignore-parse-errors negative_hits.warn_once_per_file
              --txt
              --xml-pretty
              --xml ${CMAKE_BINARY_DIR}/fortran-coverage.xml
              --html-details ${CMAKE_BINARY_DIR}/fortran-coverage.html
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Running non-slow Fortran tests and reporting src/ coverage"
      DEPENDS ${CABLEDYN_FORTRAN_COVERAGE_TARGETS}
      VERBATIM
    )
  else()
    message(WARNING
      "CABLEDYN_ENABLE_COVERAGE is ON, but gcovr was not found. "
      "Install gcovr in the active environment to enable the fortran_coverage target.")
  endif()
endif()

# Dynamic steps at the residual round-off floor (a held stiff cable at rest, a strongly
# damped line on a refined mesh), and every shipped dynamic example must complete its march.
foreach(_floor_case IN ITEMS held_start_gomex80_power_cable held_start_humboldt800_power_cable
    ve_polyester_waves_fine_mesh umaine_power_cable_alone_static umaine_power_cable_alone_dynamic)
  add_test(NAME ${_floor_case}
    COMMAND cabledyn "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/${_floor_case}.dat"
      "${CMAKE_CURRENT_BINARY_DIR}/${_floor_case}")
  set_tests_properties(${_floor_case} PROPERTIES LABELS "fortran;driver;integration;dynamic")
endforeach()
foreach(_dyn_example IN ITEMS dynamic_chain_held dynamic_chain_current dynamic_chain_waves
    ve_polyester_dynamic_waves syrope_polyester_mooring clump_weight_free_point connect_weighted_point)
  add_test(NAME example_${_dyn_example}
    COMMAND cabledyn "${CMAKE_CURRENT_SOURCE_DIR}/examples/${_dyn_example}.dat"
      "${CMAKE_CURRENT_BINARY_DIR}/example_${_dyn_example}")
  set_tests_properties(example_${_dyn_example} PROPERTIES LABELS "fortran;driver;integration;example;dynamic")
endforeach()
# Feature examples (vessel motion and RAOs, discrete buoyancy, range graphs and touchdown, wave
# spectra, stream waves, anisotropic friction, modal analysis, clamped ends on bodies, MoorDyn-C
# kinematics files, bodies and rods, line failure): every shipped deck must run to completion.
foreach(_feature_example IN ITEMS
    lazy_wave_vessel_motion lazy_wave_vessel_rao lazy_wave_buoyancy_modules lazy_wave_buoyancy_smeared
    chain_range_tdp chain_torsethaugen_spread chain_two_train_sea stream_wave_shallow_chain
    airy_wave_shallow_chain laid_cable_cross_current chain_modes lazy_wave_modes buoy_clamped_cable
    moordynC_wavekin/chain_wavekin7_currents1 spar_pinned_rods mixed_body_rods_points
    als_volturnus_line_break_time als_volturnus_line_break_tension)
  get_filename_component(_feature_name "${_feature_example}" NAME)
  add_test(NAME example_${_feature_name}
    COMMAND cabledyn "${CMAKE_CURRENT_SOURCE_DIR}/examples/${_feature_example}.dat"
      "${CMAKE_CURRENT_BINARY_DIR}/example_${_feature_name}")
  set_tests_properties(example_${_feature_name} PROPERTIES LABELS "fortran;driver;integration;example;dynamic")
endforeach()

# Composite grounded moorings must initialize at every mesh refinement of the shipped decks.
add_executable(test_static_mesh_refinement tests/test_static_mesh_refinement.f90)
target_link_libraries(test_static_mesh_refinement PRIVATE cabledyn_core)
add_test(NAME static_mesh_refinement
  COMMAND test_static_mesh_refinement "${CMAKE_CURRENT_SOURCE_DIR}/examples")
set_tests_properties(static_mesh_refinement PROPERTIES LABELS "fortran;driver;integration;static")

add_executable(test_hermite_strain_bound tests/test_hermite_strain_bound.f90)
target_link_libraries(test_hermite_strain_bound PRIVATE cabledyn_core)
add_test(NAME hermite_strain_bound COMMAND test_hermite_strain_bound)
set_tests_properties(hermite_strain_bound PROPERTIES LABELS "fortran;unit;finite_ei")

# --- test time limits ---------------------------------------------------------------------
# Tests that run for minutes in an optimised build carry the `slow` label, so the fast
# loop (`ctest -L fortran -LE slow`) leaves them out. CTest applies no timeout of its own,
# so every test without an explicit TIMEOUT gets a default that stops a hang; the long
# tests get a larger one. Both leave room for the Debug (-fcheck=all) build.
foreach(_slow_test IN ITEMS example_als_volturnus_line_break_time example_als_volturnus_line_break_tension
                            static_init_sweep)
  if(TEST ${_slow_test})
    set_property(TEST ${_slow_test} APPEND PROPERTY LABELS slow)
  endif()
endforeach()
foreach(_long_test IN ITEMS l3_lazywave_dynamic l3_lazywave_dynamic_force_blend l4_lazywave_coupled
                            l3_longduration_mooring lazywave_dynamics static_init_sweep body_rod_physics
                            example_als_volturnus_line_break_time example_als_volturnus_line_break_tension
                            lozon_example_equilibrium memory_growth memory_growth_aggregate driver_repeat_exit)
  if(TEST ${_long_test})
    get_test_property(${_long_test} TIMEOUT _long_timeout)
    if(NOT _long_timeout)
      set_tests_properties(${_long_test} PROPERTIES TIMEOUT 5400)
    endif()
  endif()
endforeach()
get_property(_all_tests DIRECTORY PROPERTY TESTS)
foreach(_any_test IN LISTS _all_tests)
  get_test_property(${_any_test} TIMEOUT _any_timeout)
  if(NOT _any_timeout)
    set_tests_properties(${_any_test} PROPERTIES TIMEOUT 1800)
  endif()
endforeach()
