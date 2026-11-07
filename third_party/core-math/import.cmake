# import.cmake — imports CORE-MATH binary64 functions into third_party/core-math.
#
# Usage:  cmake -DUPSTREAM=<path to a core-math checkout> -P third_party/core-math/import.cmake
#
# For every function listed below it copies src/binary64/<f>/<f>.c and the local headers that
# file includes ("x.h"), applying exactly these mechanical changes (README.myb.md):
#   1. '#include "myb_coremath_prelude.h"' inserted before the first preprocessor line of the
#      .c file (so that it is never inside a conditional block);
#   2. __builtin_fma        -> MYB_CM_FMA    (exact fma, see the prelude);
#   3. printf( / exit(      -> MYB_CM_PRINTF( / MYB_CM_EXIT(   (no output; an "impossible"
#      worst case terminates through MYB_CHECK instead of exit(1)).
# Running it again on the same upstream revision must leave `git status` clean.
cmake_minimum_required(VERSION 3.21)

if(NOT UPSTREAM OR NOT EXISTS "${UPSTREAM}/src/binary64")
  message(FATAL_ERROR "usage: cmake -DUPSTREAM=<core-math checkout> -P import.cmake")
endif()

set(functions
  sin cos tan atan atan2 asin acos
  exp expm1 log log1p
  sinh cosh tanh asinh acosh atanh
  cbrt hypot)

set(target_root "${CMAKE_CURRENT_LIST_DIR}/binary64")

function(transform content out_var)
  # CMake regular expressions have no \b: the character before or after is matched explicitly.
  string(REGEX REPLACE "__builtin_fma([^a-zA-Z_0-9])" "MYB_CM_FMA\\1" content "${content}")
  string(REGEX REPLACE "([^a-zA-Z_0-9])printf([ ]*\\()" "\\1MYB_CM_PRINTF\\2" content "${content}")
  string(REGEX REPLACE "([^a-zA-Z_0-9])exit([ ]*\\()" "\\1MYB_CM_EXIT\\2" content "${content}")
  set(${out_var} "${content}" PARENT_SCOPE)
endfunction()

foreach(f IN LISTS functions)
  set(source_dir "${UPSTREAM}/src/binary64/${f}")
  set(target_dir "${target_root}/${f}")
  file(MAKE_DIRECTORY "${target_dir}")

  file(READ "${source_dir}/${f}.c" content)
  # The first line starting with '#' follows the licence comment block.
  string(FIND "${content}" "\n#" first_directive)
  if(first_directive LESS 0)
    message(FATAL_ERROR "${f}.c has no preprocessor line")
  endif()
  math(EXPR first_directive "${first_directive} + 1")
  string(SUBSTRING "${content}" 0 ${first_directive} head)
  string(SUBSTRING "${content}" ${first_directive} -1 tail)
  set(content "${head}#include \"myb_coremath_prelude.h\" /* myb: see third_party/core-math/README.myb.md */\n${tail}")
  transform("${content}" content)
  file(WRITE "${target_dir}/${f}.c" "${content}")

  # Local headers included by the function (dint.h, tint.h, ...).
  string(REGEX MATCHALL "#include \"[a-zA-Z_0-9]+\\.h\"" locals "${tail}")
  foreach(line IN LISTS locals)
    string(REGEX REPLACE "#include \"([a-zA-Z_0-9]+\\.h)\"" "\\1" header "${line}")
    if(header STREQUAL "myb_coremath_prelude.h")
      continue()
    endif()
    file(READ "${source_dir}/${header}" header_content)
    transform("${header_content}" header_content)
    file(WRITE "${target_dir}/${header}" "${header_content}")
  endforeach()

  # Every fma must be routed; no output or exit may remain.
  file(GLOB written "${target_dir}/*")
  foreach(file IN LISTS written)
    file(READ "${file}" check)
    if(check MATCHES "__builtin_fma[^a-zA-Z_0-9]" OR check MATCHES "[^a-zA-Z_0-9]printf[ ]*\\("
       OR check MATCHES "[^a-zA-Z_0-9]exit[ ]*\\(")
      message(FATAL_ERROR "untransformed call left in ${file}")
    endif()
  endforeach()
  message(STATUS "imported ${f}")
endforeach()
