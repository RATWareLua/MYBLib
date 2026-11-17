# check_gcc.cmake — a program built with GCC uses the shared library built by Clang. The headers
# are written for any compiler, and GCC reads some of them more strictly than Clang: in C23 it
# refuses a standard attribute after a GNU one, which Clang takes (the visibility of the
# declarations and [[nodiscard]], include/myb/core.h). Checked:
#   - the public headers, those of ext/ among them, compile as C17 and C23 under gcc and as
#     C++17, C++20 and C++23 under g++, with -Wall -Wextra -Wpedantic -Werror, both with the
#     attributes of the shared library and with MYB_STATIC;
#   - examples/low_orbit.c built by gcc as C23 runs against libmyb.so;
#   - the C++ face (bindings/cpp) built by g++, with exceptions off, passes its own test — the
#     orbit of the example through the classes bit for bit that of the C calls.
#
# Usage: cmake -DGCC=<gcc> -DGXX=<g++> -DSOURCE_DIR=<repository> -DLIB_DIR=<dir of libmyb.so>
#              -DWORK_DIR=<scratch> -P check_gcc.cmake
foreach(_required GCC GXX SOURCE_DIR LIB_DIR WORK_DIR)
  if(NOT ${_required})
    message(FATAL_ERROR "check_gcc.cmake: ${_required} is not set")
  endif()
endforeach()

function(run what)
  execute_process(COMMAND ${ARGN} WORKING_DIRECTORY "${WORK_DIR}"
    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${what} failed (${rc}):\n${out}\n${err}")
  endif()
  set(run_output "${out}" PARENT_SCOPE)
endfunction()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

# GCC itself, not another compiler under its name (Apple's gcc is Clang)
run("${GCC} -dM -E" "${GCC}" -dM -E -x c /dev/null)
if(run_output MATCHES "__clang__" OR NOT run_output MATCHES "__GNUC__")
  message(FATAL_ERROR "${GCC} is not GCC")
endif()

set(_include "-I${SOURCE_DIR}/include" "-I${SOURCE_DIR}/bindings/cpp/include")
foreach(_ext ionosphere jb2008 nrlmsise00 sgp4 troposphere)
  list(APPEND _include "-I${SOURCE_DIR}/ext/${_ext}")
endforeach()
set(_warnings -Wall -Wextra -Wpedantic -Werror)
set(_data "-DMYB_TEST_DATA_DIR=\"${SOURCE_DIR}/tests/data\""
  "-DMYB_DATA_DIR=\"${SOURCE_DIR}/data\"")
# the example calls sqrt itself: libm, which a C program on Linux names
set(_link "-L${LIB_DIR}" "-Wl,-rpath,${LIB_DIR}" -lmyb -lm)

foreach(_static "" "-DMYB_STATIC")
  run("the headers as C17 ${_static}" "${GCC}" -std=c17 ${_static} ${_warnings} ${_include}
    -fsyntax-only "${SOURCE_DIR}/tests/compat/headers_c17.c")
  run("the headers as C23 ${_static}" "${GCC}" -std=c2x ${_static} ${_warnings} ${_include}
    -fsyntax-only "${SOURCE_DIR}/tests/compat/headers_c23.c")
  foreach(_standard c++17 c++20 c++23)
    run("the headers as ${_standard} ${_static}" "${GXX}" -std=${_standard} ${_static}
      ${_warnings} ${_include} -fsyntax-only "${SOURCE_DIR}/tests/compat/headers_cxx17.cpp")
  endforeach()
endforeach()

run("examples/low_orbit.c" "${GCC}" -std=c2x ${_warnings} ${_include}
  "${SOURCE_DIR}/examples/low_orbit.c" -o low_orbit ${_link})
run("low_orbit" "${WORK_DIR}/low_orbit")
message(STATUS "examples/low_orbit.c built by gcc:\n${run_output}")

run("the C++ face" "${GXX}" -std=c++17 ${_warnings} -fno-exceptions ${_include} ${_data}
  "${SOURCE_DIR}/tests/bindings/test_cpp.cc" -o test_cpp ${_link})
run("test_cpp" "${WORK_DIR}/test_cpp")
if(NOT run_output MATCHES "checks, 0 failures")
  message(FATAL_ERROR "the C++ face built by g++:\n${run_output}")
endif()
message(STATUS "bindings/cpp built by g++:\n${run_output}")
