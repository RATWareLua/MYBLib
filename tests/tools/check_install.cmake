# check_install.cmake — the installed library is usable by an outside project: installs the
# build tree into a scratch prefix, configures examples/ against it with find_package(myb),
# builds the examples and runs them (each returns 0 when its own check holds).
#
# Usage: cmake -DBUILD_DIR=<build tree> -DEXAMPLES_DIR=<examples> -DWORK_DIR=<scratch>
#              -DGENERATOR=<generator> [-DTOOLCHAIN=<file>] [-DMAKE_PROGRAM=<make>]
#              [-DCONFIG=<configuration>] [-DSHARED=ON]
#              [-DLIBDIR=<lib>] [-DINCLUDEDIR=<include>] [-DDATADIR=<share>]
#              -P check_install.cmake
# LIBDIR, INCLUDEDIR and DATADIR are the install directories of the build relative to the
# prefix (CMAKE_INSTALL_LIBDIR is lib64 on some systems).
foreach(_required BUILD_DIR EXAMPLES_DIR WORK_DIR GENERATOR)
  if(NOT ${_required})
    message(FATAL_ERROR "check_install.cmake: ${_required} is not set")
  endif()
endforeach()
if(NOT LIBDIR)
  set(LIBDIR lib)
endif()
if(NOT INCLUDEDIR)
  set(INCLUDEDIR include)
endif()
if(NOT DATADIR)
  set(DATADIR share)
endif()

function(run what)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${what} failed (${rc}):\n${out}\n${err}")
  endif()
  set(run_output "${out}" PARENT_SCOPE)
endfunction()

file(REMOVE_RECURSE "${WORK_DIR}")
set(prefix "${WORK_DIR}/prefix")
set(config_args "")
if(CONFIG)
  set(config_args --config "${CONFIG}")
endif()
run("install" "${CMAKE_COMMAND}" --install "${BUILD_DIR}" --prefix "${prefix}" ${config_args})

foreach(_file ${INCLUDEDIR}/myb/myb.h ${INCLUDEDIR}/myb/myb.hpp ${INCLUDEDIR}/myb_sgp4.h
        ${LIBDIR}/cmake/myb/mybConfig.cmake ${LIBDIR}/pkgconfig/myb.pc
        ${LIBDIR}/pkgconfig/myb_sgp4.pc ${LIBDIR}/pkgconfig/myb_nrlmsise00.pc
        ${LIBDIR}/pkgconfig/myb_troposphere.pc ${LIBDIR}/pkgconfig/myb_ionosphere.pc
        ${DATADIR}/myb/data/time/time-scales.tk ${DATADIR}/licenses/myb/LICENSE)
  if(NOT EXISTS "${prefix}/${_file}")
    message(FATAL_ERROR "not installed: ${_file}")
  endif()
endforeach()
if(EXISTS "${prefix}/${DATADIR}/myb/data/download")
  message(FATAL_ERROR "the downloaded data sets must not be installed")
endif()

set(configure_args -S "${EXAMPLES_DIR}" -B "${WORK_DIR}/build" -G "${GENERATOR}"
  "-DCMAKE_PREFIX_PATH=${prefix}" -DCMAKE_BUILD_TYPE=Release)
if(TOOLCHAIN)
  list(APPEND configure_args "-DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN}")
endif()
if(MAKE_PROGRAM)
  list(APPEND configure_args "-DCMAKE_MAKE_PROGRAM=${MAKE_PROGRAM}")
endif()
if(SHARED)
  list(APPEND configure_args -DMYB_EXAMPLES_SHARED=ON)
endif()
run("configuring the examples" "${CMAKE_COMMAND}" ${configure_args})
run("building the examples" "${CMAKE_COMMAND}" --build "${WORK_DIR}/build" --config Release)

# The shared library is found next to the program on Windows (PATH); elsewhere the build tree
# carries its location (RPATH).
if(CMAKE_HOST_WIN32)
  set(ENV{PATH} "${prefix}/bin;$ENV{PATH}")
endif()
foreach(_example low_orbit time_scales)
  file(GLOB_RECURSE programs "${WORK_DIR}/build/myb_example_${_example}"
    "${WORK_DIR}/build/myb_example_${_example}.exe")
  list(LENGTH programs count)
  if(NOT count EQUAL 1)
    message(FATAL_ERROR "example ${_example}: expected one program, found '${programs}'")
  endif()
  run("example ${_example}" "${programs}")
  message(STATUS "${_example}:\n${run_output}")
endforeach()
# the C++ example, if the examples found a C++ compiler
file(GLOB_RECURSE programs "${WORK_DIR}/build/myb_example_low_orbit_cpp"
  "${WORK_DIR}/build/myb_example_low_orbit_cpp.exe")
list(LENGTH programs count)
if(count EQUAL 1)
  run("example low_orbit_cpp" "${programs}")
  message(STATUS "low_orbit_cpp:\n${run_output}")
elseif(REQUIRE_CXX)
  message(FATAL_ERROR "the C++ example was not built: '${programs}'")
endif()
message(STATUS "installed library: examples built and passed")
