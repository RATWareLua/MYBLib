# Clang of LLVM on macOS (Apple silicon or Intel). The Clang of Xcode is older than the one the
# library needs (C23 constexpr, #embed: Clang >= 19), so the compiler is that of an LLVM
# installed beside it — Homebrew's llvm, found through MYB_LLVM_ROOT (environment or cache) or
# the usual Homebrew prefixes. The system's SDK, linker (ld64) and C library are used.
set(MYB_LLVM_ROOT "$ENV{MYB_LLVM_ROOT}" CACHE PATH "Root of the LLVM installation (bin/clang)")
if(NOT MYB_LLVM_ROOT)
  foreach(_myb_candidate /opt/homebrew/opt/llvm /usr/local/opt/llvm)
    if(NOT MYB_LLVM_ROOT AND EXISTS "${_myb_candidate}/bin/clang")
      set(MYB_LLVM_ROOT "${_myb_candidate}" CACHE PATH "Root of the LLVM installation (bin/clang)" FORCE)
    endif()
  endforeach()
endif()
if(MYB_LLVM_ROOT)
  set(CMAKE_C_COMPILER "${MYB_LLVM_ROOT}/bin/clang")
  set(CMAKE_CXX_COMPILER "${MYB_LLVM_ROOT}/bin/clang++")
else()
  set(CMAKE_C_COMPILER clang)
  set(CMAKE_CXX_COMPILER clang++)
endif()

# The SDK of the system: the headers and libraries an LLVM outside Xcode does not carry.
if(NOT CMAKE_OSX_SYSROOT)
  execute_process(COMMAND xcrun --show-sdk-path OUTPUT_VARIABLE _myb_sdk
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_myb_sdk)
    set(CMAKE_OSX_SYSROOT "${_myb_sdk}" CACHE PATH "macOS SDK")
  endif()
endif()

# Homebrew's MPFR and GMP for the tests against MPFR.
foreach(_myb_prefix /opt/homebrew /usr/local)
  if(EXISTS "${_myb_prefix}/include/mpfr.h")
    list(APPEND CMAKE_PREFIX_PATH "${_myb_prefix}")
  endif()
endforeach()
