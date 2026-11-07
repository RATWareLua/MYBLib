cmake_minimum_required(VERSION 3.28)

# The version lives in include/myb/core.h only.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/myb/core.h" _myb_version_lines
  REGEX "^#define MYB_VERSION_(MAJOR|MINOR|PATCH) [0-9]+$")
foreach(_line IN LISTS _myb_version_lines)
  string(REGEX MATCH "MYB_VERSION_(MAJOR|MINOR|PATCH) ([0-9]+)" _ "${_line}")
  set(_myb_version_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
endforeach()

project(myb
  VERSION "${_myb_version_MAJOR}.${_myb_version_MINOR}.${_myb_version_PATCH}"
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

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# ---- Vendored correctly rounded functions (CORE-MATH, MIT) ------------------------------
file(GLOB _myb_coremath_sources CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/core-math/binary64/*/*.c")
add_library(myb_coremath OBJECT ${_myb_coremath_sources})
target_include_directories(myb_coremath PRIVATE src third_party/core-math)
myb_apply_flags(myb_coremath THIRD_PARTY)

# ---- The library ----------------------------------------------------------------------
if(WIN32)
  set(_myb_platform_sources src/platform/platform_win32.c)
else()
  set(_myb_platform_sources src/platform/platform_posix.c)
endif()

add_library(myb STATIC
  src/body/body.c
  src/core/alloc.c
  src/core/check.c
  src/core/status.c
  src/io/kpool.c
  src/kepler/anomaly.c
  src/kepler/elements.c
  src/kepler/propagate.c
  src/kepler/stumpff.c
  src/io/numparse.c
  src/num/dd.c
  src/num/root.c
  src/num/series.c
  src/num/sum.c
  src/cmath/cmath.c
  src/time/instant.c
  src/transfer/lambert.c
  src/time/timesys.c
  ${_myb_platform_sources}
  $<TARGET_OBJECTS:myb_coremath>)
add_library(myb::myb ALIAS myb)
target_include_directories(myb
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  PRIVATE src third_party/core-math)
target_compile_definitions(myb PUBLIC MYB_STATIC)
myb_apply_flags(myb)

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

# ---- Tests ----------------------------------------------------------------------------
option(MYB_BUILD_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
if(MYB_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
