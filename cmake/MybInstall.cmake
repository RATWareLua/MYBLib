# MybInstall.cmake — install rules of myb: libraries, headers, data kernels, licences, the
# CMake package (find_package(myb)) and a pkg-config file.
#
#   <prefix>/include/myb/*.h, include/myb_<ext>.h; include/myb/myb.hpp (C++, target myb::cpp)
#   <prefix>/lib/libmyb.a [libmyb.so, libmyb.dll.a], libmyb_<ext>.a ...; bin/libmyb.dll
#   <prefix>/lib/cmake/myb/mybConfig.cmake ...     targets myb::myb, myb::myb_shared,
#                                                  myb::<ext>, myb::<ext>_shared; myb_DATA_DIR
#   <prefix>/lib/pkgconfig/myb.pc, myb_<ext>.pc
#   <prefix>/share/myb/data/...                    the kernels of data/ (not data/download)
#   <prefix>/share/licenses/myb/...
include_guard(GLOBAL)
include(CMakePackageConfigHelpers)

set(_myb_targets myb myb_cpp)
if(MYB_BUILD_SHARED)
  list(APPEND _myb_targets myb_shared)
endif()
get_property(_myb_ext_modules GLOBAL PROPERTY MYB_EXT_MODULES)
foreach(_ext IN LISTS _myb_ext_modules)
  list(APPEND _myb_targets myb_${_ext})
  if(MYB_BUILD_SHARED)
    list(APPEND _myb_targets myb_${_ext}_shared)
  endif()
  install(FILES "${PROJECT_SOURCE_DIR}/ext/${_ext}/myb_${_ext}.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(FILES "${PROJECT_SOURCE_DIR}/ext/${_ext}/NOTICE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/myb" RENAME "NOTICE.${_ext}")
endforeach()

install(TARGETS ${_myb_targets} EXPORT mybTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/myb"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" FILES_MATCHING PATTERN "*.h")
install(FILES "${PROJECT_SOURCE_DIR}/bindings/cpp/include/myb/myb.hpp"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/myb")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/myb/data"
  PATTERN "download" EXCLUDE)
install(FILES "${PROJECT_SOURCE_DIR}/LICENSE" "${PROJECT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/myb")

# ---- CMake package ------------------------------------------------------------------------
set(_myb_cmake_dir "${CMAKE_INSTALL_LIBDIR}/cmake/myb")
set(MYB_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}/myb/data")
install(EXPORT mybTargets NAMESPACE myb:: DESTINATION "${_myb_cmake_dir}")
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/mybConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/mybConfig.cmake"
  INSTALL_DESTINATION "${_myb_cmake_dir}"
  PATH_VARS MYB_INSTALL_DATADIR)
# The API is stable within a major version.
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/mybConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES "${PROJECT_BINARY_DIR}/mybConfig.cmake" "${PROJECT_BINARY_DIR}/mybConfigVersion.cmake"
  DESTINATION "${_myb_cmake_dir}")

# ---- pkg-config -----------------------------------------------------------------------------
# With the shared library installed `-lmyb` finds it; a static-only installation needs
# MYB_STATIC on Windows (no dllimport) and the system libraries of the static archive.
function(_myb_pc_path out_var base dir)
  if(IS_ABSOLUTE "${dir}")
    set(${out_var} "${dir}" PARENT_SCOPE)
  else()
    set(${out_var} "${base}/${dir}" PARENT_SCOPE)
  endif()
endfunction()
# The prefix is found from the place of the file itself (${pcfiledir}), so an installation
# moved or installed with another --prefix still answers correctly.
if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
  set(MYB_PC_PREFIX "${CMAKE_INSTALL_PREFIX}")
else()
  # one ".." for every directory between the prefix and the file
  string(REGEX REPLACE "[^/]+" ".." _myb_pc_up "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  set(MYB_PC_PREFIX "\${pcfiledir}/${_myb_pc_up}")
endif()
_myb_pc_path(MYB_PC_LIBDIR [=[${exec_prefix}]=] "${CMAKE_INSTALL_LIBDIR}")
_myb_pc_path(MYB_PC_INCLUDEDIR [=[${prefix}]=] "${CMAKE_INSTALL_INCLUDEDIR}")
_myb_pc_path(MYB_PC_DATADIR [=[${prefix}]=] "${CMAKE_INSTALL_DATADIR}/myb/data")
set(_myb_system_libs "")
if(NOT WIN32)
  set(_myb_system_libs "-lm -pthread")
endif()
if(MYB_BUILD_SHARED)
  set(MYB_PC_CFLAGS "")
  set(MYB_PC_LIBS "")
  set(MYB_PC_LIBS_PRIVATE "${_myb_system_libs}")
else()
  set(MYB_PC_CFLAGS "-DMYB_STATIC")
  set(MYB_PC_LIBS "${_myb_system_libs}")
  set(MYB_PC_LIBS_PRIVATE "")
endif()
configure_file("${PROJECT_SOURCE_DIR}/cmake/myb.pc.in" "${PROJECT_BINARY_DIR}/myb.pc" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/myb.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
# One file for each library of ext/: `pkg-config --libs myb_sgp4` gives -lmyb_sgp4 and, through
# Requires, what myb itself needs.
set(_myb_ext_text_nrlmsise00 "NRLMSISE-00 model of the Earth's upper atmosphere for myb")
set(_myb_ext_text_sgp4 "SGP4/SDP4 propagation of two-line element sets for myb")
set(_myb_ext_text_troposphere "Tropospheric delay of radio and optical signals (IERS Conventions 2010) for myb")
set(_myb_ext_text_ionosphere "Ionospheric total electron content maps (IONEX) for myb")
set(_myb_ext_text_jb2008 "JB2008 model of the Earth's thermosphere, implemented from the papers, for myb")
foreach(MYB_PC_EXT IN LISTS _myb_ext_modules)
  set(MYB_PC_EXT_DESCRIPTION "${_myb_ext_text_${MYB_PC_EXT}}")
  if(NOT MYB_PC_EXT_DESCRIPTION)
    set(MYB_PC_EXT_DESCRIPTION "Optional module ${MYB_PC_EXT} of myb")
  endif()
  configure_file("${PROJECT_SOURCE_DIR}/cmake/myb_ext.pc.in"
    "${PROJECT_BINARY_DIR}/myb_${MYB_PC_EXT}.pc" @ONLY)
  install(FILES "${PROJECT_BINARY_DIR}/myb_${MYB_PC_EXT}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
endforeach()
