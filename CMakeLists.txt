cmake_minimum_required(VERSION 3.28)

# The version lives in include/myb/core.h only.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/myb/core.h" _myb_version_lines
  REGEX "^#define MYB_VERSION_(MAJOR|MINOR|PATCH|TWEAK) [0-9]+$")
foreach(_line IN LISTS _myb_version_lines)
  string(REGEX MATCH "MYB_VERSION_(MAJOR|MINOR|PATCH|TWEAK) ([0-9]+)" _ "${_line}")
  set(_myb_version_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
endforeach()

project(myb
  VERSION "${_myb_version_MAJOR}.${_myb_version_MINOR}.${_myb_version_PATCH}.${_myb_version_TWEAK}"
  DESCRIPTION "Motion of bodies in space: N-body, orbits, transfers - rigorous physics in C23"
  LANGUAGES C)

# C++ is used only to check that the public headers compile as C++17.
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
  enable_language(CXX)
endif()

set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(MybFlags)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# ---- Vendored correctly rounded functions (CORE-MATH, MIT) ------------------------------
file(GLOB _myb_coremath_sources CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/core-math/binary64/*/*.c")
add_library(myb_coremath OBJECT ${_myb_coremath_sources})
target_include_directories(myb_coremath PRIVATE include src third_party/core-math)
# Never part of the API of a binary: the functions are reached through myb/cmath.h only. The
# objects go into both the static and the shared library, always next to the functions of the
# core they call: no import or export attributes (MYB_STATIC).
set_target_properties(myb_coremath PROPERTIES C_VISIBILITY_PRESET hidden)
target_compile_definitions(myb_coremath PRIVATE MYB_STATIC)
myb_apply_flags(myb_coremath THIRD_PARTY)

# ---- The library ----------------------------------------------------------------------
if(WIN32)
  set(_myb_platform_sources src/platform/platform_win32.c)
else()
  set(_myb_platform_sources src/platform/platform_posix.c)
endif()

set(_myb_sources
  src/attitude/attitude.c
  src/attitude/free.c
  src/body/body.c
  src/core/alloc.c
  src/core/check.c
  src/core/pool.c
  src/core/status.c
  src/cr3bp/cr3bp.c
  src/cr3bp/torus.c
  src/dyn/model.c
  src/dyn/nbody.c
  src/dyn/system.c
  src/dyn/whfast.c
  src/env/env.c
  src/env/terms.c
  src/env/weather.c
  src/estim/filter.c
  src/estim/fit.c
  src/estim/light.c
  src/estim/measure.c
  src/frames/bpck.c
  src/frames/cip.c
  src/frames/equinox.c
  src/frames/eop.c
  src/frames/geodetic.c
  src/frames/graph.c
  src/frames/site.c
  src/frames/orientation.c
  src/frames/orbital.c
  src/frames/rotation.c
  src/frames/tkframe.c
  src/force/atmosphere.c
  src/force/nongrav.c
  src/force/force.c
  src/force/magnetic.c
  src/force/charging.c
  src/force/radiation.c
  src/force/surface.c
  src/force/qlaw.c
  src/force/thrust.c
  src/frames/wgccre.c
  src/geom/geometry.c
  src/geom/view.c
  src/grav/ellipsoidal.c
  src/grav/figure.c
  src/grav/mascon.c
  src/grav/polyhedron.c
  src/grav/sh.c
  src/grav/tides.c
  src/integ/dop853.c
  src/integ/event.c
  src/integ/gauss_jackson.c
  src/integ/ias15.c
  src/integ/rkn12.c
  src/integ/symplectic.c
  src/io/ccsds.c
  src/io/daf.c
  src/io/kpool.c
  src/kepler/anomaly.c
  src/kin/ck.c
  src/kin/sclk.c
  src/kin/segset.c
  src/kin/spk.c
  src/kin/sp3.c
  src/kin/spk_segment.c
  src/kin/spk_writer.c
  src/kepler/elements.c
  src/kepler/propagate.c
  src/kepler/iod.c
  src/kepler/mean.c
  src/kepler/regular.c
  src/kepler/relative.c
  src/kepler/stumpff.c
  src/io/numformat.c
  src/io/numparse.c
  src/io/table.c
  src/io/text.c
  src/num/dd.c
  src/num/elliptic.c
  src/num/erf.c
  src/num/linalg.c
  src/num/poisson.c
  src/num/root.c
  src/num/series.c
  src/observe/observe.c
  src/orbit/events.c
  src/orbit/orbit.c
  src/num/sum.c
  src/rel/rel.c
  src/cmath/cmath.c
  src/time/instant.c
  src/transfer/indirect.c
  src/transfer/lambert.c
  src/transfer/nlp.c
  src/transfer/sims_flanagan.c
  src/transfer/transfer.c
  src/uncert/budget.c
  src/uncert/conjunction.c
  src/uncert/uncert.c
  src/time/timesys.c
  ${_myb_platform_sources})

add_library(myb STATIC ${_myb_sources} $<TARGET_OBJECTS:myb_coremath>)
add_library(myb::myb ALIAS myb)
target_include_directories(myb
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
         "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE src third_party/core-math)
target_compile_definitions(myb PUBLIC MYB_STATIC)
myb_apply_flags(myb)

# ---- The shared library (optional) ------------------------------------------------------
# The same sources compiled once more with hidden visibility: only the public API (MYB_API)
# and the few support functions the libraries of ext/ need leave the binary. The static
# library stays the one the tests link, since they reach into the internals.
# myb_symbol_versions(<shared target> <map>): ELF symbol versions (cmake/myb.map explains).
# ELF only: PE and Mach-O have no such thing.
function(myb_symbol_versions target map)
  if(UNIX AND NOT APPLE)
    target_link_options(${target} PRIVATE "LINKER:--version-script=${map}")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${map}")
  endif()
endfunction()

option(MYB_BUILD_SHARED "Also build the shared library (libmyb.so, libmyb.dll)" OFF)
if(MYB_BUILD_SHARED)
  set_target_properties(myb_coremath PROPERTIES POSITION_INDEPENDENT_CODE ON)
  add_library(myb_shared SHARED ${_myb_sources} $<TARGET_OBJECTS:myb_coremath>)
  add_library(myb::myb_shared ALIAS myb_shared)
  set_target_properties(myb_shared PROPERTIES
    OUTPUT_NAME myb
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    C_VISIBILITY_PRESET hidden)
  target_include_directories(myb_shared
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    PRIVATE src third_party/core-math)
  target_compile_definitions(myb_shared PRIVATE MYB_BUILDING)
  myb_apply_flags(myb_shared)
  myb_symbol_versions(myb_shared "${CMAKE_CURRENT_SOURCE_DIR}/cmake/myb.map")
endif()

# Stack frames reported per function (Clang -fstack-usage); tests/tools/check_stack.cmake keeps
# them within a budget. Off under the sanitizers, whose redzones enlarge every frame.
if(MYB_SANITIZE STREQUAL "")
  set(_myb_stack_default ON)
else()
  set(_myb_stack_default OFF)
endif()
option(MYB_STACK_USAGE "Report stack frames (-fstack-usage) and check their budget" ${_myb_stack_default})
if(MYB_STACK_USAGE)
  target_compile_options(myb PRIVATE -fstack-usage)
  target_compile_options(myb_coremath PRIVATE -fstack-usage)
endif()
if(NOT WIN32)
  # src/platform/platform_posix.c: the worker threads of myb_pool.
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries(myb PUBLIC Threads::Threads)
  # glibc keeps the floating-point environment functions (fegetround, feraiseexcept, used by
  # CORE-MATH) in libm; tests/tools/check_imports.cmake still forbids its inexact functions.
  target_link_libraries(myb PUBLIC m)
  if(MYB_BUILD_SHARED)
    target_link_libraries(myb_shared PRIVATE Threads::Threads m)
  endif()
endif()

# ---- Public headers stay valid C17 and C++17 -------------------------------------------
add_library(myb_headers_c17 OBJECT tests/compat/headers_c17.c)
set_target_properties(myb_headers_c17 PROPERTIES C_STANDARD 17)
target_link_libraries(myb_headers_c17 PRIVATE myb)
target_compile_options(myb_headers_c17 PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(CMAKE_CXX_COMPILER)
  add_library(myb_headers_cxx17 OBJECT tests/compat/headers_cxx17.cpp)
  set_target_properties(myb_headers_cxx17 PROPERTIES CXX_STANDARD 17 CXX_EXTENSIONS OFF)
  target_link_libraries(myb_headers_cxx17 PRIVATE myb)
  target_compile_options(myb_headers_cxx17 PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

# ---- Optional modules of particular bodies (ext/) -------------------------------------
option(MYB_BUILD_EXT "Build the body-specific models of ext/ (NRLMSISE-00, SGP4, troposphere, ionosphere)" ON)
if(MYB_BUILD_EXT)
  add_subdirectory(ext)
endif()

# ---- Bindings -------------------------------------------------------------------------
add_subdirectory(bindings)

# ---- Examples -------------------------------------------------------------------------
option(MYB_BUILD_EXAMPLES "Build the programs of examples/" OFF)
if(MYB_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ---- Installation -----------------------------------------------------------------------
# Headers, the libraries, the data kernels and a CMake package (find_package(myb): targets
# myb::myb, myb::myb_shared, myb::sgp4, ...; myb_DATA_DIR) with a pkg-config file.
option(MYB_INSTALL "Generate the install rules" ${PROJECT_IS_TOP_LEVEL})
if(MYB_INSTALL)
  include(MybInstall)
endif()

# ---- Tests ----------------------------------------------------------------------------
option(MYB_BUILD_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
if(MYB_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
