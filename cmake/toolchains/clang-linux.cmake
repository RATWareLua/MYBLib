# Clang on Linux (Gentoo and others). The system libc and linker are used.
# lld is taken if it is installed; otherwise the default linker of the clang driver.
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)

find_program(_myb_lld NAMES ld.lld)
if(_myb_lld)
  foreach(_kind EXE SHARED MODULE)
    set(CMAKE_${_kind}_LINKER_FLAGS_INIT "-fuse-ld=lld")
  endforeach()
endif()
