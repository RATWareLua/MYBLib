# run_fuzzers.cmake — runs the libFuzzer harnesses of a build with MYB_FUZZ (preset linux-fuzz)
# for a time budget each and keeps what they find:
#
#   cmake -DBUILD_DIR=build/linux-fuzz [-DSECONDS=300] [-DTARGETS=spk;sp3] [-DFORK=0]
#         [-DMAX_LEN=0] [-DTIMEOUT=20] [-DRSS_MB=2048] -P tests/fuzz/run_fuzzers.cmake
#
# For every harness fuzz_<name> (all whose seed list <build>/tests/fuzz/fuzz_<name>.seeds the
# configuration wrote, or those of TARGETS):
#   <build>/fuzz/corpus/<name>/   the working corpus: made from the seeds on the first run,
#                                 grown by every run after (libFuzzer adds what reaches new code)
#   <build>/fuzz/crashes/<name>/  crash-*, timeout-*, oom-*, leak-*: the inputs that failed
#   <build>/fuzz/<name>.log       the output of the last run
# FORK > 0 runs that many processes per harness (-fork=N); MAX_LEN 0 lets libFuzzer take the
# length of the largest seed. Exit status: not 0 if a harness left a new failing input — copy
# it into tests/fuzz/corpus/<name>/ (regress-*) once the reader is fixed, so that every build
# replays it.
cmake_minimum_required(VERSION 3.28)

if(NOT BUILD_DIR)
  message(FATAL_ERROR "usage: cmake -DBUILD_DIR=<build tree of MYB_FUZZ> [-DSECONDS=300] "
    "[-DTARGETS=a;b] -P run_fuzzers.cmake")
endif()
get_filename_component(BUILD_DIR "${BUILD_DIR}" ABSOLUTE)
foreach(_default SECONDS:300 FORK:0 MAX_LEN:0 TIMEOUT:20 RSS_MB:2048)
  string(REPLACE ":" ";" _pair "${_default}")
  list(GET _pair 0 _name)
  list(GET _pair 1 _value)
  if(NOT DEFINED ${_name})
    set(${_name} ${_value})
  endif()
endforeach()

if(NOT EXISTS "${BUILD_DIR}/CMakeCache.txt")
  message(FATAL_ERROR "${BUILD_DIR} is not a build tree")
endif()
file(STRINGS "${BUILD_DIR}/CMakeCache.txt" _cache REGEX "^(MYB_FUZZ|MYB_HAS_LIBFUZZER):")
if(NOT "MYB_FUZZ:BOOL=ON" IN_LIST _cache OR NOT "MYB_HAS_LIBFUZZER:INTERNAL=1" IN_LIST _cache)
  message(FATAL_ERROR "${BUILD_DIR} is no libFuzzer build (MYB_FUZZ with Clang on Linux, "
    "preset linux-fuzz); the harnesses of other builds replay files (tests/fuzz/replay.cmake)")
endif()

set(_fuzz_dir "${BUILD_DIR}/tests/fuzz")
if(NOT TARGETS)
  file(GLOB _lists "${_fuzz_dir}/fuzz_*.seeds")
  foreach(_list IN LISTS _lists)
    get_filename_component(_stem "${_list}" NAME_WE)
    string(REGEX REPLACE "^fuzz_" "" _target "${_stem}")
    list(APPEND TARGETS ${_target})
  endforeach()
endif()
if(NOT TARGETS)
  message(FATAL_ERROR "no harnesses in ${_fuzz_dir} (configure and build with MYB_FUZZ)")
endif()

set(_found "")
foreach(_target IN LISTS TARGETS)
  set(_program "${_fuzz_dir}/fuzz_${_target}")
  if(NOT EXISTS "${_program}" AND EXISTS "${_program}.exe")
    set(_program "${_program}.exe")
  endif()
  if(NOT EXISTS "${_program}")
    message(WARNING "fuzz_${_target} is not built: skipped")
    continue()
  endif()
  set(_corpus "${BUILD_DIR}/fuzz/corpus/${_target}")
  set(_crashes "${BUILD_DIR}/fuzz/crashes/${_target}")
  file(MAKE_DIRECTORY "${_corpus}" "${_crashes}")
  # the seeds go into the working corpus once (libFuzzer takes directories, not lists)
  file(GLOB _present "${_corpus}/*")
  if(NOT _present)
    file(STRINGS "${_fuzz_dir}/fuzz_${_target}.seeds" _seeds)
    set(_k 0)
    foreach(_seed IN LISTS _seeds)
      if(EXISTS "${_seed}")
        get_filename_component(_base "${_seed}" NAME)
        file(COPY_FILE "${_seed}" "${_corpus}/seed-${_k}-${_base}")
        math(EXPR _k "${_k} + 1")
      endif()
    endforeach()
  endif()
  file(GLOB _before "${_crashes}/*")
  set(_args -max_total_time=${SECONDS} -timeout=${TIMEOUT} -rss_limit_mb=${RSS_MB}
            -print_final_stats=1 "-artifact_prefix=${_crashes}/")
  if(MAX_LEN GREATER 0)
    list(APPEND _args -max_len=${MAX_LEN})
  endif()
  if(FORK GREATER 0)
    list(APPEND _args -fork=${FORK} -ignore_crashes=1)
  endif()
  message(STATUS "fuzz_${_target}: ${SECONDS} s")
  execute_process(COMMAND "${_program}" ${_args} "${_corpus}"
    RESULT_VARIABLE _rc OUTPUT_FILE "${BUILD_DIR}/fuzz/${_target}.log"
    ERROR_FILE "${BUILD_DIR}/fuzz/${_target}.log")
  file(GLOB _after "${_crashes}/*")
  set(_new "")
  foreach(_file IN LISTS _after)
    if(NOT _file IN_LIST _before)
      list(APPEND _new "${_file}")
    endif()
  endforeach()
  file(STRINGS "${BUILD_DIR}/fuzz/${_target}.log" _stats REGEX "^stat::number_of_executed_units")
  message(STATUS "fuzz_${_target}: exit ${_rc}; ${_stats}")
  foreach(_file IN LISTS _new)
    message(STATUS "  NEW ${_file}")
  endforeach()
  list(APPEND _found ${_new})
endforeach()

if(_found)
  list(LENGTH _found _count)
  message(FATAL_ERROR "${_count} new failing input(s) under ${BUILD_DIR}/fuzz/crashes; "
    "reproduce one with: <build>/tests/fuzz/fuzz_<name> <file>")
endif()
message(STATUS "no new failing inputs")
