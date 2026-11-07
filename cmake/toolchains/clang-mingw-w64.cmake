# Clang (LLVM for Windows) targeting x86_64-w64-windows-gnu.
#
# The compiler is clang; headers, C runtime and libgcc come from a mingw-w64 sysroot
# (MSYS2 MINGW64 by default: C:/msys64/mingw64). Visual Studio is not needed.
# Linker: lld. Make tool: mingw32-make from the same sysroot (or Ninja, if installed).
#
# Override the sysroot with the environment variable MYB_MINGW_ROOT or -DMYB_MINGW_ROOT=...
if(DEFINED ENV{MYB_MINGW_ROOT} AND NOT DEFINED MYB_MINGW_ROOT)
  set(MYB_MINGW_ROOT "$ENV{MYB_MINGW_ROOT}")
endif()
if(NOT DEFINED MYB_MINGW_ROOT)
  set(MYB_MINGW_ROOT "C:/msys64/mingw64")
endif()
if(NOT EXISTS "${MYB_MINGW_ROOT}/include/stdio.h")
  message(FATAL_ERROR "mingw-w64 sysroot not found at '${MYB_MINGW_ROOT}'. Install MSYS2 "
    "(MINGW64 environment) or set MYB_MINGW_ROOT.")
endif()

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_C_COMPILER_TARGET x86_64-w64-windows-gnu)
set(CMAKE_CXX_COMPILER_TARGET x86_64-w64-windows-gnu)
set(CMAKE_SYSROOT "${MYB_MINGW_ROOT}")

foreach(_kind EXE SHARED MODULE)
  set(CMAKE_${_kind}_LINKER_FLAGS_INIT "-fuse-ld=lld")
endforeach()

if(NOT CMAKE_MAKE_PROGRAM AND CMAKE_GENERATOR STREQUAL "MinGW Makefiles")
  find_program(CMAKE_MAKE_PROGRAM mingw32-make HINTS "${MYB_MINGW_ROOT}/bin" NO_DEFAULT_PATH)
endif()

# Libraries for tests (MPFR, GMP) live in the sysroot itself; programs (llvm-nm) on the host.
list(APPEND CMAKE_PREFIX_PATH "${MYB_MINGW_ROOT}")
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
