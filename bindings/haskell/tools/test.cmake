# SPDX-License-Identifier: BSD-3-Clause
# The test of the Haskell package (test_haskell_binding), run by ctest:
#   cmake -DCABAL=<cabal> -DGHC=<ghc> -DCLANG=<clang> -DSOURCE=<repository> -DBINARY=<directory>
#         -DLIBRARIES=<a|b|...> [-DRUNTIME_DIRS=<a|b|...>] [-DDATA_DIR=...] [-DTEST_DATA_DIR=...]
#         -P test.cmake
#
# Everything cabal writes - its configuration, its store, its logs, the build - stays under
# BINARY: CABAL_DIR points there, with a configuration that names no repository (the package
# needs the boot packages of GHC only). The libraries to link (LIBRARIES: libmyb.a or the
# import library of libmyb.dll, and those of ext/) are copied into BINARY/lib, the one
# directory given to the linker, so that it cannot take the other kind from a directory that
# holds both; RUNTIME_DIRS are put in front of PATH (LD_LIBRARY_PATH) for the shared ones.
#
# Three steps: the lower layer must agree with the headers (myb-generate --check), the test
# suite of the package runs, and the example of the README (examples/LowOrbit.hs).
foreach(_var CABAL GHC CLANG SOURCE BINARY LIBRARIES)
  if(NOT DEFINED ${_var})
    message(FATAL_ERROR "test.cmake: ${_var} is not given")
  endif()
endforeach()

set(_package "${SOURCE}/bindings/haskell")
file(MAKE_DIRECTORY "${BINARY}/cabal" "${BINARY}/lib" "${BINARY}/output")
# the files the test writes (an SPK file read back) stay in the build tree
file(TO_NATIVE_PATH "${BINARY}/output" _output)
set(ENV{MYB_TEST_OUTPUT_DIR} "${_output}")
if(NOT EXISTS "${BINARY}/cabal/config")
  file(WRITE "${BINARY}/cabal/config"
    "-- Written by bindings/haskell/tools/test.cmake: no repository, the package needs the\n"
    "-- boot packages of GHC only.\n")
endif()
file(TO_NATIVE_PATH "${BINARY}/cabal" _cabal_dir)
set(ENV{CABAL_DIR} "${_cabal_dir}")

# the libraries to link, alone in their directory
file(GLOB _old "${BINARY}/lib/*")
if(_old)
  file(REMOVE ${_old})
endif()
string(REPLACE "|" ";" _libraries "${LIBRARIES}")
foreach(_library IN LISTS _libraries)
  file(COPY "${_library}" DESTINATION "${BINARY}/lib")
endforeach()

if(DEFINED RUNTIME_DIRS AND NOT RUNTIME_DIRS STREQUAL "")
  string(REPLACE "|" ";" _runtime "${RUNTIME_DIRS}")
  if(CMAKE_HOST_WIN32)
    set(_native "")
    foreach(_dir IN LISTS _runtime)
      file(TO_NATIVE_PATH "${_dir}" _dir_native)
      list(APPEND _native "${_dir_native}")
    endforeach()
    list(JOIN _native "\;" _joined)
    set(ENV{PATH} "${_joined};$ENV{PATH}")
  else()
    list(JOIN _runtime ":" _joined)
    set(ENV{LD_LIBRARY_PATH} "${_joined}:$ENV{LD_LIBRARY_PATH}")
    set(ENV{DYLD_LIBRARY_PATH} "${_joined}:$ENV{DYLD_LIBRARY_PATH}")
  endif()
endif()
if(DEFINED DATA_DIR)
  set(ENV{MYB_DATA_DIR} "${DATA_DIR}")
endif()
if(DEFINED TEST_DATA_DIR)
  set(ENV{MYB_TEST_DATA_DIR} "${TEST_DATA_DIR}")
endif()
set(ENV{MYB_CLANG} "${CLANG}")

# the same flags for both steps, so that cabal configures the package once
file(TO_NATIVE_PATH "${BINARY}/dist" _builddir)
file(TO_NATIVE_PATH "${BINARY}/lib" _libdir)
set(_flags --builddir=${_builddir} -w ${GHC} -j1 --extra-lib-dirs=${_libdir}
  "--ghc-options=+RTS -M2g -RTS")

execute_process(
  COMMAND "${CABAL}" run ${_flags} myb-generate -- --check --clang "${CLANG}" "${SOURCE}"
  WORKING_DIRECTORY "${_package}"
  RESULT_VARIABLE _generated)
if(NOT _generated EQUAL 0)
  message(FATAL_ERROR "the lower layer does not agree with the headers (exit ${_generated}): "
    "run myb-generate as bindings/haskell/README.md says")
endif()

execute_process(
  COMMAND "${CABAL}" test ${_flags} binding --test-show-details=direct
  WORKING_DIRECTORY "${_package}"
  RESULT_VARIABLE _tested)
if(NOT _tested EQUAL 0)
  message(FATAL_ERROR "the test suite of the Haskell package failed (exit ${_tested})")
endif()

# the example of the README
execute_process(
  COMMAND "${CABAL}" run ${_flags} myb-example
  WORKING_DIRECTORY "${_package}"
  RESULT_VARIABLE _example)
if(NOT _example EQUAL 0)
  message(FATAL_ERROR "the example of the Haskell package failed (exit ${_example})")
endif()
