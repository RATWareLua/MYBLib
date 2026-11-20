# check_msvc.cmake — a program built with Microsoft's compiler uses the shared library built by
# Clang for the GNU ABI. The DLL imports only KERNEL32 and msvcrt and passes no object of a C
# runtime across its interface (memory is released by the library that took it), so the two
# toolchains meet at the C ABI of Windows x64 alone. Checked:
#   - the public headers compile as C17 and as C++17 and C++20 under cl /W4 /WX;
#   - examples/low_orbit.c built by cl runs against libmyb.dll;
#   - the C++ face (bindings/cpp) built by cl, with exceptions off, passes its own test — the
#     orbit of the example through the classes bit for bit that of the C calls.
# The import library of the GNU toolchain (libmyb.dll.a) is read by Microsoft's linker as it
# is; it only has to carry the name that linker looks for, myb.lib.
#
# Usage: cmake -DVCVARS=<vcvars64.bat> -DSOURCE_DIR=<repository> -DLIB_DIR=<dir of libmyb.dll>
#              -DWORK_DIR=<scratch> -P check_msvc.cmake
cmake_minimum_required(VERSION 3.28)
foreach(_required VCVARS SOURCE_DIR LIB_DIR WORK_DIR)
  if(NOT ${_required})
    message(FATAL_ERROR "check_msvc.cmake: ${_required} is not set")
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

# The environment of the compiler: what vcvars64.bat leaves behind, taken from `set`. Through a
# file of two lines: cmd does not read the quotes of a command line passed to it with escapes,
# and vcvars64.bat may end with a nonzero code over a missing helper (vswhere) while setting
# everything — the variables themselves are checked below.
file(TO_NATIVE_PATH "${VCVARS}" _vcvars)
file(WRITE "${WORK_DIR}/environment.cmd" "@call \"${_vcvars}\" >nul 2>nul\r\n@set\r\n")
file(TO_NATIVE_PATH "${WORK_DIR}/environment.cmd" _environment_cmd)
execute_process(COMMAND cmd /c "${_environment_cmd}"
  RESULT_VARIABLE rc OUTPUT_VARIABLE environment ERROR_QUIET)
# The names of the environment of Windows have no case (a cmd started from PowerShell writes
# Path): each is matched in any case.
foreach(_name INCLUDE LIB LIBPATH PATH)
  set(_any_case "")
  string(LENGTH "${_name}" _length)
  math(EXPR _last "${_length} - 1")
  foreach(_i RANGE ${_last})
    string(SUBSTRING "${_name}" ${_i} 1 _upper)
    string(TOLOWER "${_upper}" _lower)
    string(APPEND _any_case "[${_upper}${_lower}]")
  endforeach()
  if(NOT environment MATCHES "(^|\n)${_any_case}=([^\r\n]*)")
    message(FATAL_ERROR "vcvars64.bat left no ${_name}")
  endif()
  set(ENV{${_name}} "${CMAKE_MATCH_2}")
endforeach()
file(TO_NATIVE_PATH "${LIB_DIR}" _lib_dir)
set(ENV{PATH} "${_lib_dir};$ENV{PATH}")
find_program(_cl cl NO_CACHE NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH
  NO_PACKAGE_ROOT_PATH)
if(NOT _cl)
  message(FATAL_ERROR "cl.exe is not on the PATH of vcvars64.bat")
endif()

configure_file("${LIB_DIR}/libmyb.dll.a" "${WORK_DIR}/myb.lib" COPYONLY)

set(_include "/I${SOURCE_DIR}/include" "/I${SOURCE_DIR}/bindings/cpp/include")
set(_data "/DMYB_TEST_DATA_DIR=\"${SOURCE_DIR}/tests/data\"" "/DMYB_DATA_DIR=\"${SOURCE_DIR}/data\"")
set(_c /nologo /W4 /WX /utf-8 ${_include})
set(_cxx /nologo /W4 /WX /utf-8 /permissive- /Zc:__cplusplus ${_include})

run("the headers as C17" "${_cl}" ${_c} /std:c17 /c "${SOURCE_DIR}/tests/compat/headers_c17.c"
  /Foheaders_c17.obj)
foreach(_standard c++17 c++20)
  run("the headers as ${_standard}" "${_cl}" ${_cxx} /std:${_standard} /c
    "${SOURCE_DIR}/tests/compat/headers_cxx17.cpp" /Foheaders_${_standard}.obj)
endforeach()

run("examples/low_orbit.c" "${_cl}" ${_c} /std:c17 "${SOURCE_DIR}/examples/low_orbit.c" myb.lib
  /Fe:low_orbit.exe)
run("low_orbit.exe" "${WORK_DIR}/low_orbit.exe")
message(STATUS "examples/low_orbit.c built by cl:\n${run_output}")

# no /EHsc: the header handles no exceptions, and with them off the standard library it uses
# must not need them either (cl warns, C4530, and /WX makes that an error)
run("the C++ face" "${_cl}" ${_cxx} /std:c++17 ${_data} "${SOURCE_DIR}/tests/bindings/test_cpp.cc"
  myb.lib /Fe:test_cpp.exe)
run("test_cpp.exe" "${WORK_DIR}/test_cpp.exe")
if(NOT run_output MATCHES "checks, 0 failures")
  message(FATAL_ERROR "the C++ face built by cl:\n${run_output}")
endif()
message(STATUS "bindings/cpp built by cl:\n${run_output}")
