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
# cgo finds the library through pkg-config, the module myb-static: the static library whatever
# else is there (myb.pc links the shared one when it is installed too) and the system libraries
# it needs.
#   ${PROJECT_BINARY_DIR}/bindings/go/pkgconfig/myb-static.pc   the build tree (test_go_binding)
#   <prefix>/<libdir>/pkgconfig/myb-static.pc                    an installation
# -Bstatic takes libmyb.a even where libmyb.dll.a or libmyb.so lies next to it.
set(_myb_go_system_libs "")
if(NOT WIN32)
  set(_myb_go_system_libs " -lm -pthread")
endif()
string(REPLACE " " "\\ " _myb_go_include "${PROJECT_SOURCE_DIR}/include") # as cgo splits words
file(GENERATE OUTPUT "${PROJECT_BINARY_DIR}/bindings/go/pkgconfig/myb-static.pc" CONTENT
"Name: myb-static
Description: ${PROJECT_DESCRIPTION} (the static library of the build tree)
Version: ${PROJECT_VERSION}
Cflags: -I${_myb_go_include} -DMYB_STATIC
Libs: -L$<TARGET_FILE_DIR:myb> -Wl,-Bstatic -lmyb -Wl,-Bdynamic${_myb_go_system_libs}
")

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
  file(WRITE "${PROJECT_BINARY_DIR}/bindings/go/myb-static.pc"
"prefix=${_myb_go_prefix}
libdir=${_myb_go_libdir}
includedir=${_myb_go_includedir}

Name: myb-static
Description: ${PROJECT_DESCRIPTION} (the static library, whatever else is installed)
URL: https://github.com/RATWareLua/MYBLib
Version: ${PROJECT_VERSION}
Cflags: -I\${includedir} -DMYB_STATIC
Libs: -L\${libdir} -Wl,-Bstatic -lmyb -Wl,-Bdynamic${_myb_go_system_libs}
")
  install(FILES "${PROJECT_BINARY_DIR}/bindings/go/myb-static.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endif()
