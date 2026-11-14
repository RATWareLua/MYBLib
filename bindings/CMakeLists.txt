# Bindings: the library in other languages. Each is a thin face over the C interface — no
# arithmetic of its own, so a result has the bits of the same calls made in C.
#
#   cpp/     a C++17 header, myb/myb.hpp: target myb::cpp
#   python/  a package over the shared library (cffi)
#   rust/    the crates myb-sys (declarations) and myb (ownership, outcomes)
#   haskell/ a cabal package: Myb.C (generated declarations, hsc2hs) and Myb (ownership,
#            outcomes); built and tested by cabal (tests/CMakeLists.txt: test_haskell_binding)
#   go/      a module of plain cgo over the static library (the coarse layer: one call, one task)

# ---- C++ ---------------------------------------------------------------------------------------
# Header-only: ownership (RAII), outcomes (myb::result<T>), lambdas for terms and events.
add_library(myb_cpp INTERFACE)
add_library(myb::cpp ALIAS myb_cpp)
set_target_properties(myb_cpp PROPERTIES EXPORT_NAME cpp)
target_include_directories(myb_cpp INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(myb_cpp INTERFACE myb)
if(CMAKE_CXX_COMPILER)
  target_compile_features(myb_cpp INTERFACE cxx_std_17)
endif()

# ---- Go ----------------------------------------------------------------------------------------
# cgo finds the libraries through pkg-config, the modules myb-static and myb_sgp4-static: the
# static libraries whatever else is there (myb.pc links the shared one when it is installed
# too) and the system libraries they need.
#   ${PROJECT_BINARY_DIR}/bindings/go/pkgconfig/*-static.pc   the build tree (test_go_binding)
#   <prefix>/<libdir>/pkgconfig/*-static.pc                    an installation
# Each archive is named by its path: that takes lib<name>.a even where lib<name>.dll.a or
# lib<name>.so lies next to it, on every linker (Apple's has no -Bstatic, and pkg-config drops
# the second of two -Wl,-Bstatic when one module requires another). A path with blanks is
# written with each blank escaped, as pkg-config and cgo split words.
set(_myb_go_system_libs "")
if(NOT WIN32)
  set(_myb_go_system_libs " -lm -pthread")
endif()
function(_myb_go_escape out_var text)
  string(REPLACE " " "\\ " escaped "${text}")
  set(${out_var} "${escaped}" PARENT_SCOPE)
endfunction()
# _myb_go_libs(<out> <dir> <name>): how a static library <dir>/lib<name>.a is linked.
function(_myb_go_libs out_var dir name)
  set(${out_var} "${dir}/lib${name}.a" PARENT_SCOPE)
endfunction()
# _myb_go_archive_dir(<out> <target>): the directory of the static library of a target in the
# build tree (its archive directory, or the binary directory it was defined in; with a
# configuration's subdirectory under a multi-configuration generator).
function(_myb_go_archive_dir out_var target)
  get_target_property(dir ${target} ARCHIVE_OUTPUT_DIRECTORY)
  if(NOT dir)
    get_target_property(dir ${target} BINARY_DIR)
  endif()
  get_property(multi GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(multi)
    string(APPEND dir "/$<CONFIG>")
  endif()
  _myb_go_escape(dir "${dir}")
  set(${out_var} "${dir}" PARENT_SCOPE)
endfunction()

_myb_go_escape(_myb_go_include "${PROJECT_SOURCE_DIR}/include")
_myb_go_archive_dir(_myb_go_dir myb)
_myb_go_libs(_myb_go_link "${_myb_go_dir}" myb)
file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/bindings/go/pkgconfig/myb-static.pc" CONTENT
"Name: myb-static
Description: ${PROJECT_DESCRIPTION} (the static library of the build tree)
Version: ${PROJECT_VERSION}
Cflags: -I${_myb_go_include} -DMYB_STATIC
Libs: ${_myb_go_link}${_myb_go_system_libs}
")
if(TARGET myb_sgp4)
  _myb_go_escape(_myb_go_sgp4_include "${PROJECT_SOURCE_DIR}/ext/sgp4")
  _myb_go_archive_dir(_myb_go_sgp4_dir myb_sgp4)
  _myb_go_libs(_myb_go_sgp4_link "${_myb_go_sgp4_dir}" myb_sgp4)
  file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/bindings/go/pkgconfig/myb_sgp4-static.pc" CONTENT
"Name: myb_sgp4-static
Description: SGP4/SDP4 propagation of two-line element sets for myb (the static library of the build tree)
Version: ${PROJECT_VERSION}
Requires: myb-static
Cflags: -I${_myb_go_sgp4_include}
Libs: ${_myb_go_sgp4_link}
")
endif()

# MYB_INSTALL is defined after this directory (at the top level), with the default given here.
if((DEFINED MYB_INSTALL AND MYB_INSTALL) OR (NOT DEFINED MYB_INSTALL AND PROJECT_IS_TOP_LEVEL))
  # The prefix from the place of the file itself, as myb.pc has it (cmake/MybInstall.cmake).
  if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    set(_myb_go_prefix "${CMAKE_INSTALL_PREFIX}")
    set(_myb_go_libdir "${CMAKE_INSTALL_LIBDIR}")
  else()
    string(REGEX REPLACE "[^/]+" ".." _myb_go_up "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
    set(_myb_go_prefix "\${pcfiledir}/${_myb_go_up}")
    set(_myb_go_libdir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
  endif()
  if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
    set(_myb_go_includedir "${CMAKE_INSTALL_INCLUDEDIR}")
  else()
    set(_myb_go_includedir "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
  endif()
  _myb_go_libs(_myb_go_link "\${libdir}" myb)
  file(WRITE "${PROJECT_BINARY_DIR}/bindings/go/myb-static.pc"
"prefix=${_myb_go_prefix}
libdir=${_myb_go_libdir}
includedir=${_myb_go_includedir}

Name: myb-static
Description: ${PROJECT_DESCRIPTION} (the static library, whatever else is installed)
URL: https://github.com/RATWareLua/MYBLib
Version: ${PROJECT_VERSION}
Cflags: -I\${includedir} -DMYB_STATIC
Libs: ${_myb_go_link}${_myb_go_system_libs}
")
  install(FILES "${PROJECT_BINARY_DIR}/bindings/go/myb-static.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  if(TARGET myb_sgp4)
    _myb_go_libs(_myb_go_sgp4_link "\${libdir}" myb_sgp4)
    file(WRITE "${PROJECT_BINARY_DIR}/bindings/go/myb_sgp4-static.pc"
"prefix=${_myb_go_prefix}
libdir=${_myb_go_libdir}
includedir=${_myb_go_includedir}

Name: myb_sgp4-static
Description: SGP4/SDP4 propagation of two-line element sets for myb (the static library)
URL: https://github.com/RATWareLua/MYBLib
Version: ${PROJECT_VERSION}
Requires: myb-static = ${PROJECT_VERSION}
Cflags: -I\${includedir}
Libs: ${_myb_go_sgp4_link}
")
    install(FILES "${PROJECT_BINARY_DIR}/bindings/go/myb_sgp4-static.pc"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  endif()
endif()
