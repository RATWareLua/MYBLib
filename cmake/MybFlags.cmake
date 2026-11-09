# MybFlags.cmake — compiler gate and flags for every target of myb.
#
# myb is built with Clang only (GNU-style driver), on Windows (mingw-w64 target) and Linux.
# The flags below are part of the numerical contract of the library, not a matter of taste:
#   -ffp-contract=off   no silent a*b+c -> fma fusion (otherwise results depend on -march)
#   -fno-fast-math      IEEE-754 semantics: no reassociation, NaN/Inf honoured, signed zeros
#   -fno-math-errno     sqrt becomes one instruction on every target (target defaults differ)
#   -march=x86-64-v2    one instruction-set baseline for every build (SSE4.2, POPCNT, no FMA)
# They are appended AFTER user flags (CMAKE_C_FLAGS, Gentoo make.conf CFLAGS), so a user
# -march=native or -O3 cannot re-enable contraction. -ffast-math from the user is rejected
# at compile time by src/core/fpguard.h.
#
#   myb_apply_flags(<target> [THIRD_PARTY])
#     THIRD_PARTY — vendored code: same numerical flags, no project warnings.
#
# Cache options:
#   MYB_WARNINGS_AS_ERRORS  ON/OFF   -Werror for project code (development and CI)
#   MYB_HARDENING           ON/OFF   stack protector, CET, auto-var-init, RELRO
#   MYB_X64_BASELINE        x86-64-v2 (default) | x86-64-v3 — ISA baseline on x86-64. It is
#                           not optional: below v2 (no SSE4.1) CORE-MATH's roundeven becomes a
#                           C runtime call, absent from mingw-w64. Results are bit-identical
#                           for v2 and v3 (tests/repro): contraction is off, fma is myb_fma.
#   MYB_SANITIZE            "" | address | undefined | leak
#       address   — ASan + UBSan (LeakSanitizer included on Linux), abort on first error; in
#                   Debug also the checks of pointer pairs (pointer-compare/-subtract)
#       undefined — UBSan: with its runtime (report and stack, abort on the first error) where
#                   the runtime exists; in trap mode on the Windows target, which has none
#       leak      — standalone LeakSanitizer (Linux): leaks without ASan's instrumentation
#       thread    — ThreadSanitizer (Linux): data races and lock-order inversions of the
#                   thread team (myb_pool) and of everything run on it
include_guard(GLOBAL)
include(CheckCCompilerFlag)

option(MYB_WARNINGS_AS_ERRORS "Treat warnings as errors in project code" OFF)
option(MYB_HARDENING "Enable hardening flags" ON)
set(MYB_X64_BASELINE "x86-64-v2" CACHE STRING "x86-64 ISA baseline: x86-64-v2 or x86-64-v3")
set_property(CACHE MYB_X64_BASELINE PROPERTY STRINGS x86-64-v2 x86-64-v3)
set(MYB_SANITIZE "" CACHE STRING "Sanitizer: empty, address, undefined, leak or thread")
set_property(CACHE MYB_SANITIZE PROPERTY STRINGS "" address undefined leak thread)

# ---- Compiler gate ----------------------------------------------------------------------
set(MYB_MIN_CLANG_VERSION 19)
get_property(_myb_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
foreach(_lang IN LISTS _myb_languages)
  if(NOT _lang MATCHES "^(C|CXX)$")
    continue()
  endif()
  if(NOT CMAKE_${_lang}_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
      "myb is built with Clang only; ${_lang} compiler is '${CMAKE_${_lang}_COMPILER_ID}' "
      "(${CMAKE_${_lang}_COMPILER}). Use a preset from CMakePresets.json or a toolchain file "
      "from cmake/toolchains/.")
  endif()
  if(CMAKE_${_lang}_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    message(FATAL_ERROR "clang-cl is not supported: use the GNU-style 'clang' driver "
      "(cmake/toolchains/clang-mingw-w64.cmake on Windows).")
  endif()
  if(CMAKE_${_lang}_COMPILER_VERSION VERSION_LESS MYB_MIN_CLANG_VERSION)
    message(FATAL_ERROR "myb needs Clang >= ${MYB_MIN_CLANG_VERSION} (C23 constexpr, #embed); "
      "found ${CMAKE_${_lang}_COMPILER_VERSION}.")
  endif()
endforeach()

# ---- Target architecture ------------------------------------------------------------------
include(CheckCSourceCompiles)
check_c_source_compiles("
#if !defined(__x86_64__)
#error not x86-64
#endif
int main(void) { return 0; }" MYB_ARCH_X86_64)

# Keeps only the flags the compiler accepts without a warning (the check also links).
function(_myb_supported_flags out_var)
  set(CMAKE_REQUIRED_FLAGS "-Werror")
  set(result "")
  foreach(flag IN LISTS ARGN)
    string(MAKE_C_IDENTIFIER "MYB_HAS_FLAG_${flag}" var)
    check_c_compiler_flag("${flag}" ${var})
    if(${var})
      list(APPEND result "${flag}")
    endif()
  endforeach()
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

# ---- Numerical contract (all targets, appended last) -----------------------------------
set(_myb_fp_flags -ffp-contract=off -fno-fast-math -fno-math-errno)
if(MYB_ARCH_X86_64)
  if(NOT MYB_X64_BASELINE MATCHES "^x86-64-v[23]$")
    message(FATAL_ERROR "MYB_X64_BASELINE must be x86-64-v2 or x86-64-v3, not "
      "'${MYB_X64_BASELINE}' (AVX-512 levels are outside the project's ISA ceiling).")
  endif()
  list(PREPEND _myb_fp_flags "-march=${MYB_X64_BASELINE}")
endif()
# __FILE__ in MYB_CHECK messages and debug info do not depend on the checkout location.
set(_myb_repro_flags
  "-ffile-prefix-map=${CMAKE_SOURCE_DIR}/="
  "-ffile-prefix-map=${CMAKE_BINARY_DIR}/=build/")

# ---- Warnings (project code only) -------------------------------------------------------
_myb_supported_flags(_myb_warnings
  -Wall -Wextra -Wpedantic
  -Wshadow -Wconversion -Wsign-conversion -Wdouble-promotion
  -Wformat=2 -Werror=format-security
  -Wimplicit-fallthrough -Wnull-dereference -Wcast-qual -Wundef
  -Wstrict-prototypes -Wmissing-prototypes -Wmissing-variable-declarations -Wvla
  -Wcomma -Wunreachable-code
  -Werror=implicit -Werror=incompatible-pointer-types -Werror=int-conversion)

# ---- Hardening ------------------------------------------------------------------------
set(_myb_hardening_candidates -fstack-protector-strong -fstack-clash-protection)
if(MYB_ARCH_X86_64)
  list(APPEND _myb_hardening_candidates -fcf-protection=full)
endif()
_myb_supported_flags(_myb_hardening ${_myb_hardening_candidates})
_myb_supported_flags(_myb_autoinit_debug -ftrivial-auto-var-init=pattern)
_myb_supported_flags(_myb_autoinit_release -ftrivial-auto-var-init=zero)

# ---- Sanitizers -------------------------------------------------------------------------
set(_myb_san_compile "")
set(_myb_san_link "")
if(MYB_SANITIZE STREQUAL "address")
  set(_myb_san_compile -fsanitize=address,undefined -fno-sanitize-recover=all
    -fno-omit-frame-pointer)
  set(_myb_san_link -fsanitize=address,undefined)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Relations between pointers into different objects (run time:
    # ASAN_OPTIONS=detect_invalid_pointer_pairs=2). Only without optimisation: at -O2 Clang 21
    # instruments comparisons the optimiser synthesises (a pointer against -9) and reports
    # pairs that do not exist in the source.
    list(APPEND _myb_san_compile -fsanitize=pointer-compare,pointer-subtract)
    list(APPEND _myb_san_link -fsanitize=pointer-compare,pointer-subtract)
  endif()
elseif(MYB_SANITIZE STREQUAL "undefined" AND WIN32)
  # Trap mode: UB becomes an immediate crash at the faulting instruction, no runtime needed.
  set(_myb_san_compile -fsanitize=undefined -fsanitize-trap=undefined -fno-omit-frame-pointer)
  set(_myb_san_link -fsanitize=undefined -fsanitize-trap=undefined)
elseif(MYB_SANITIZE STREQUAL "undefined")
  # With the runtime: a report with the stack, then abort (no recovery).
  set(_myb_san_compile -fsanitize=undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
  set(_myb_san_link -fsanitize=undefined)
elseif(MYB_SANITIZE STREQUAL "leak")
  set(_myb_san_compile -fsanitize=leak -fno-omit-frame-pointer)
  set(_myb_san_link -fsanitize=leak)
elseif(MYB_SANITIZE STREQUAL "thread")
  set(_myb_san_compile -fsanitize=thread -fno-omit-frame-pointer)
  set(_myb_san_link -fsanitize=thread)
elseif(NOT MYB_SANITIZE STREQUAL "")
  message(FATAL_ERROR "MYB_SANITIZE: unknown value '${MYB_SANITIZE}'")
endif()

function(myb_apply_flags target)
  cmake_parse_arguments(_arg "THIRD_PARTY" "" "" ${ARGN})
  if(_arg_THIRD_PARTY)
    target_compile_options(${target} PRIVATE -w)
  else()
    target_compile_options(${target} PRIVATE ${_myb_warnings})
    if(MYB_WARNINGS_AS_ERRORS)
      target_compile_options(${target} PRIVATE -Werror)
    endif()
  endif()

  if(MYB_SANITIZE STREQUAL "")
    if(MYB_HARDENING)
      target_compile_options(${target} PRIVATE ${_myb_hardening}
        "$<IF:$<CONFIG:Debug>,${_myb_autoinit_debug},${_myb_autoinit_release}>"
        "$<$<AND:$<PLATFORM_ID:Linux>,$<NOT:$<CONFIG:Debug>>>:-U_FORTIFY_SOURCE;-D_FORTIFY_SOURCE=3>")
      if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        get_target_property(_type ${target} TYPE)
        if(_type STREQUAL "EXECUTABLE" OR _type STREQUAL "SHARED_LIBRARY")
          target_link_options(${target} PRIVATE -Wl,-z,noexecstack -Wl,-z,relro -Wl,-z,now)
        endif()
      endif()
    endif()
  else()
    target_compile_options(${target} PRIVATE ${_myb_san_compile} -U_FORTIFY_SOURCE)
    target_link_options(${target} PRIVATE ${_myb_san_link})
  endif()

  # Last on the command line: nothing configured above or by the user may override these.
  target_compile_options(${target} PRIVATE ${_myb_repro_flags} ${_myb_fp_flags})
endfunction()
