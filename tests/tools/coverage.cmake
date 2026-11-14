# coverage.cmake — source-based coverage of the library (src/) from a build with MYB_COVERAGE
# (presets linux-coverage, windows-coverage):
#
#   cmake -DBUILD_DIR=build/linux-coverage [-DRUN_TESTS=OFF] [-DJOBS=2] [-DTEST_ARGS=-LE;slow]
#         [-DLEAST=15] -P tests/tools/coverage.cmake
#   cmake -DSUMMARY=build/linux-coverage/coverage/summary.json -P tests/tools/coverage.cmake
#
# 1. runs the tests of the build tree (unless RUN_TESTS=OFF: the raw profiles of an earlier run
#    are taken) with LLVM_PROFILE_FILE=<build>/profiles/%p-%m.profraw;
# 2. merges the raw profiles with llvm-profdata into <build>/coverage/myb.profdata;
# 3. reports with llvm-cov over the instrumented programs (<build>/coverage-objects.txt, written
#    by tests/CMakeLists.txt), the files of src/ only: <build>/coverage/report.txt (per file and
#    total: regions, functions, lines, branches), <build>/coverage/html/index.html and
#    <build>/coverage/summary.json; prints the totals and the LEAST files least covered by
#    lines. With SUMMARY only the last step, from a summary written before.
# The tools are looked for next to the compiler of the build, then on PATH (LLVM_PROFDATA and
# LLVM_COV name them otherwise). Exit status: not 0 if the tests failed or a step did.
cmake_minimum_required(VERSION 3.28)

if(NOT DEFINED LEAST)
  set(LEAST 15)
endif()

# "12.3 % (covered/count)" of two integers.
function(_percent covered count out)
  if(count GREATER 0)
    math(EXPR tenths "${covered} * 1000 / ${count}")
  else()
    set(tenths 0)
  endif()
  math(EXPR whole "${tenths} / 10")
  math(EXPR fraction "${tenths} % 10")
  set(${out} "${whole}.${fraction} % (${covered}/${count})" PARENT_SCOPE)
endfunction()

# Prints the totals of a summary of llvm-cov export and the files least covered by lines.
function(print_summary json_file source)
  file(READ "${json_file}" json)
  set(line "")
  foreach(kind lines regions functions branches)
    string(JSON covered ERROR_VARIABLE error GET "${json}" data 0 totals ${kind} covered)
    string(JSON count ERROR_VARIABLE error GET "${json}" data 0 totals ${kind} count)
    if(NOT error)
      _percent(${covered} ${count} text)
      string(APPEND line "\n    ${kind} ${text}")
    endif()
  endforeach()
  message(STATUS "coverage of src/ (total):${line}")
  string(JSON files LENGTH "${json}" data 0 files)
  set(rows "")
  if(files GREATER 0)
    math(EXPR last "${files} - 1")
    foreach(i RANGE ${last})
      string(JSON name GET "${json}" data 0 files ${i} filename)
      string(JSON covered GET "${json}" data 0 files ${i} summary lines covered)
      string(JSON count GET "${json}" data 0 files ${i} summary lines count)
      if(count GREATER 0)
        math(EXPR tenths "${covered} * 1000 / ${count}")
        # zero-padded to four digits so that the rows sort by the percentage
        string(LENGTH "${tenths}" digits)
        math(EXPR pad "4 - ${digits}")
        string(REPEAT "0" ${pad} zeros)
        if(IS_ABSOLUTE "${name}" AND source)
          file(RELATIVE_PATH relative "${source}" "${name}")
          if(NOT relative MATCHES "^\\.\\.")
            set(name "${relative}")
          endif()
        endif()
        list(APPEND rows "${zeros}${tenths}|${name}|${covered}|${count}")
      endif()
    endforeach()
  endif()
  list(SORT rows)
  message(STATUS "coverage: the ${LEAST} files of src/ least covered by lines:")
  set(shown 0)
  foreach(row IN LISTS rows)
    if(shown LESS LEAST)
      string(REPLACE "|" ";" parts "${row}")
      list(GET parts 1 name)
      list(GET parts 2 covered)
      list(GET parts 3 count)
      _percent(${covered} ${count} text)
      math(EXPR missed "${count} - ${covered}")
      message(STATUS "    ${text}  ${name}  (${missed} lines not run)")
      math(EXPR shown "${shown} + 1")
    endif()
  endforeach()
endfunction()

if(SUMMARY)
  print_summary("${SUMMARY}" "")
  return()
endif()

if(NOT BUILD_DIR)
  message(FATAL_ERROR "usage: cmake -DBUILD_DIR=<build tree of MYB_COVERAGE> -P coverage.cmake")
endif()
get_filename_component(BUILD_DIR "${BUILD_DIR}" ABSOLUTE)
if(NOT EXISTS "${BUILD_DIR}/CMakeCache.txt")
  message(FATAL_ERROR "${BUILD_DIR} is not a build tree (no CMakeCache.txt)")
endif()
if(NOT DEFINED RUN_TESTS)
  set(RUN_TESTS ON)
endif()
if(NOT JOBS)
  set(JOBS 2)
endif()

# What the build tree knows: the compiler, the sources, whether it is instrumented.
file(STRINGS "${BUILD_DIR}/CMakeCache.txt" _cache
  REGEX "^(CMAKE_C_COMPILER|CMAKE_HOME_DIRECTORY|MYB_COVERAGE|CMAKE_CTEST_COMMAND):")
foreach(_line IN LISTS _cache)
  string(REGEX MATCH "^([A-Z_]+):[A-Z]+=(.*)$" _ "${_line}")
  set(_cache_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
endforeach()
if(NOT _cache_MYB_COVERAGE)
  message(FATAL_ERROR "${BUILD_DIR} was not configured with MYB_COVERAGE=ON")
endif()
set(_source "${_cache_CMAKE_HOME_DIRECTORY}")
get_filename_component(_compiler_dir "${_cache_CMAKE_C_COMPILER}" DIRECTORY)
if(NOT LLVM_PROFDATA)
  find_program(LLVM_PROFDATA NAMES llvm-profdata HINTS "${_compiler_dir}")
endif()
if(NOT LLVM_COV)
  find_program(LLVM_COV NAMES llvm-cov HINTS "${_compiler_dir}")
endif()
if(NOT LLVM_PROFDATA OR NOT LLVM_COV)
  message(FATAL_ERROR "llvm-profdata or llvm-cov not found (LLVM_PROFDATA, LLVM_COV)")
endif()
set(_ctest "${_cache_CMAKE_CTEST_COMMAND}")
if(NOT _ctest)
  set(_ctest "${CMAKE_CTEST_COMMAND}")
endif()

set(_profiles "${BUILD_DIR}/profiles")
set(_out "${BUILD_DIR}/coverage")
set(_failed 0)

# ---- 1. The tests --------------------------------------------------------------------------
if(RUN_TESTS)
  file(REMOVE_RECURSE "${_profiles}")
  file(MAKE_DIRECTORY "${_profiles}")
  message(STATUS "coverage: running the tests (${JOBS} jobs)")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env "LLVM_PROFILE_FILE=${_profiles}/%p-%m.profraw"
            "${_ctest}" --test-dir "${BUILD_DIR}" --parallel ${JOBS} --output-on-failure
            ${TEST_ARGS}
    RESULT_VARIABLE _tests_rc)
  if(NOT _tests_rc EQUAL 0)
    message(WARNING "coverage: some tests failed (${_tests_rc}); the report covers what ran")
    set(_failed 1)
  endif()
endif()

# ---- 2. Merge ----------------------------------------------------------------------------------
file(GLOB _raw "${_profiles}/*.profraw")
list(LENGTH _raw _raw_count)
if(_raw_count EQUAL 0)
  message(FATAL_ERROR "coverage: no raw profiles in ${_profiles}")
endif()
file(MAKE_DIRECTORY "${_out}")
list(JOIN _raw "\n" _raw_text)
file(WRITE "${_out}/profiles.txt" "${_raw_text}\n")
execute_process(
  COMMAND "${LLVM_PROFDATA}" merge -sparse -f "${_out}/profiles.txt" -o "${_out}/myb.profdata"
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "coverage: llvm-profdata merge failed (${_rc})")
endif()
message(STATUS "coverage: ${_raw_count} raw profiles merged into ${_out}/myb.profdata")

# ---- 3. Report ---------------------------------------------------------------------------------
if(NOT EXISTS "${BUILD_DIR}/coverage-objects.txt")
  message(FATAL_ERROR "coverage: ${BUILD_DIR}/coverage-objects.txt is missing (configure again)")
endif()
file(STRINGS "${BUILD_DIR}/coverage-objects.txt" _listed)
set(_objects "")
foreach(_program IN LISTS _listed)
  if(EXISTS "${_program}")
    if(_objects)
      list(APPEND _objects -object "${_program}")
    else()
      list(APPEND _objects "${_program}")
    endif()
  endif()
endforeach()
if(NOT _objects)
  message(FATAL_ERROR "coverage: none of the programs of coverage-objects.txt is built")
endif()
# src/ only: the paths are relative to the source tree where -ffile-prefix-map made them so,
# absolute otherwise; the directories left out are matched as components of the path
set(_ignore "(^|[/\\\\])(tests|third_party|ext|bindings|examples|bench|include|build)[/\\\\]")
set(_common "-instr-profile=${_out}/myb.profdata" "-ignore-filename-regex=${_ignore}")

execute_process(
  COMMAND "${LLVM_COV}" report ${_objects} ${_common} -show-region-summary -show-branch-summary
  WORKING_DIRECTORY "${_source}" OUTPUT_FILE "${_out}/report.txt" RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "coverage: llvm-cov report failed (${_rc})")
endif()
execute_process(
  COMMAND "${LLVM_COV}" show ${_objects} ${_common} -format=html -show-branches=count
          "-output-dir=${_out}/html"
  WORKING_DIRECTORY "${_source}" RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(WARNING "coverage: llvm-cov show (HTML) failed (${_rc})")
  set(_failed 1)
endif()
execute_process(
  COMMAND "${LLVM_COV}" export ${_objects} ${_common} -summary-only -format=text
  WORKING_DIRECTORY "${_source}" OUTPUT_FILE "${_out}/summary.json" RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "coverage: llvm-cov export failed (${_rc})")
endif()

print_summary("${_out}/summary.json" "${_source}")
message(STATUS "coverage: ${_out}/report.txt, ${_out}/html/index.html")
if(_failed)
  message(FATAL_ERROR "coverage: done, but a step failed (see above)")
endif()
