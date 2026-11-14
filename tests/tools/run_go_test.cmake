# run_go_test.cmake — runs the tests of the Go package (bindings/go) against the build tree.
#
# cgo compiles the C of the package with the project's Clang (for the target and the sysroot
# of the build: on Windows the mingw-w64 one, whose C runtime the library is built for) and
# links with the linker of the build; the library comes from the pkg-config file of the build
# tree (myb-static.pc). Everything Go writes — its module and build caches, its work files,
# the temporary files of the tests — stays inside WORK_DIR; no toolchain is downloaded and no
# telemetry is kept (nothing goes into the user's profile).
#
# Usage: cmake -DGO=<go> -DPKG_CONFIG=<pkg-config> -DPC_DIR=<dir of myb-static.pc>
#              -DCC=<C compiler> [-DCC_TARGET=<triple>] [-DSYSROOT=<dir>] [-DLDFLAGS=<flags>]
#              -DWORK_DIR=<dir> -DSOURCE_DIR=<bindings/go> -DREFERENCE=<low_orbit_reference>
#              -DDATA_DIR=<data> -DTEST_DATA_DIR=<tests/data> [-DPACKAGES=<go packages>]
#              -P run_go_test.cmake
foreach(_var GO PKG_CONFIG PC_DIR CC WORK_DIR SOURCE_DIR REFERENCE)
  if(NOT ${_var})
    message(FATAL_ERROR "run_go_test.cmake: ${_var} is not given")
  endif()
endforeach()

file(MAKE_DIRECTORY "${WORK_DIR}/tmp")

# the compiler as cgo takes it: a command with its arguments, the path quoted (spaces)
set(_cc "\"${CC}\"")
if(CC_TARGET)
  string(APPEND _cc " --target=${CC_TARGET}")
endif()
if(SYSROOT)
  string(APPEND _cc " --sysroot=${SYSROOT}")
endif()

set(ENV{CGO_ENABLED} "1")
set(ENV{CC} "${_cc}")
# the C of the package computes nothing; it is still compiled under the numerical contract of
# the library
# The build cache of Go knows the files of the package, not the headers it includes from
# elsewhere: after a structure of the library changes, a glue object compiled against the old
# one would be taken from the cache (and the library would read past its end). A hash of the
# public headers among the flags makes their contents part of the key of every cgo compile.
file(GLOB _headers "${SOURCE_DIR}/../../include/myb/*.h" "${SOURCE_DIR}/../../ext/*/myb_*.h")
list(SORT _headers)
set(_contents "")
foreach(_header IN LISTS _headers)
  file(READ "${_header}" _text)
  string(APPEND _contents "${_text}")
endforeach()
string(SHA256 _headers_hash "${_contents}")
string(SUBSTRING "${_headers_hash}" 0 16 _headers_hash)
set(ENV{CGO_CFLAGS} "-O2 -g -ffp-contract=off -fno-fast-math -DMYBGO_HEADERS=0x${_headers_hash}")
set(ENV{CGO_LDFLAGS} "${LDFLAGS}")
set(ENV{PKG_CONFIG} "${PKG_CONFIG}")
set(ENV{PKG_CONFIG_LIBDIR} "${PC_DIR}") # the build tree only: not an installed myb-static.pc
set(ENV{PKG_CONFIG_PATH} "")
set(ENV{GOTOOLCHAIN} "local")
set(ENV{GOTELEMETRY} "off")
set(ENV{GOENV} "off")
set(ENV{GOFLAGS} "-mod=mod")
set(ENV{GOPROXY} "off")
set(ENV{GOPATH} "${WORK_DIR}/path")
set(ENV{GOCACHE} "${WORK_DIR}/cache")
set(ENV{GOMODCACHE} "${WORK_DIR}/mod")
set(ENV{GOTMPDIR} "${WORK_DIR}/tmp")
set(ENV{TMPDIR} "${WORK_DIR}/tmp")
set(ENV{TMP} "${WORK_DIR}/tmp")
set(ENV{TEMP} "${WORK_DIR}/tmp")
set(ENV{MYB_GO_REFERENCE} "${REFERENCE}")
set(ENV{MYB_DATA_DIR} "${DATA_DIR}")
set(ENV{MYB_TEST_DATA_DIR} "${TEST_DATA_DIR}")

if(NOT PACKAGES)
  set(PACKAGES ./...)
endif()
execute_process(COMMAND "${GO}" test -count=1 -v ${PACKAGES}
  WORKING_DIRECTORY "${SOURCE_DIR}" RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "go test failed (${rc})")
endif()
