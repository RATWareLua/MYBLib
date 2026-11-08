# Tests of myb. Every test is an executable: 0 = pass, 1 = fail, 77 = skipped.

function(myb_add_test name)
  cmake_parse_arguments(_t "" "" "SOURCES;LIBRARIES;INCLUDES" ${ARGN})
  add_executable(${name} ${_t_SOURCES})
  target_link_libraries(${name} PRIVATE myb ${_t_LIBRARIES})
  target_include_directories(${name} PRIVATE
    "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/support")
  if(_t_INCLUDES)
    target_include_directories(${name} SYSTEM PRIVATE ${_t_INCLUDES})
  endif()
  target_compile_definitions(${name} PRIVATE
    "MYB_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
    "MYB_DATA_DIR=\"${PROJECT_SOURCE_DIR}/data\""
    # C99-conformant printf (%a, %lld) from mingw-w64 instead of msvcrt's.
    $<$<BOOL:${WIN32}>:__USE_MINGW_ANSI_STDIO=1>)
  myb_apply_flags(${name})
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES SKIP_RETURN_CODE 77)
endfunction()

# Files compiled with -mfma: the hardware reference for myb_fma, and the probe proving that
# -ffp-contract=off survives a target that has FMA (as with a user's -march=native).
set(_myb_fma_option "")
if(MYB_ARCH_X86_64)
  set(_myb_fma_option -mfma)
endif()
set_source_files_properties(unit/fma_hw.c unit/fpenv_probe.c
  PROPERTIES COMPILE_OPTIONS "${_myb_fma_option}")

myb_add_test(test_core SOURCES unit/test_core.c)
myb_add_test(test_fpenv SOURCES unit/test_fpenv.c unit/fpenv_probe.c)
myb_add_test(test_fma SOURCES unit/test_fma.c unit/fma_hw.c)
myb_add_test(test_num SOURCES unit/test_num.c)
myb_add_test(test_parse SOURCES unit/test_parse.c)
myb_add_test(test_kpool SOURCES unit/test_kpool.c)
myb_add_test(test_time SOURCES unit/test_time.c)
myb_add_test(test_body SOURCES unit/test_body.c)
myb_add_test(test_kepler SOURCES unit/test_kepler.c)
myb_add_test(test_lambert SOURCES unit/test_lambert.c)
myb_add_test(test_dop853 SOURCES unit/test_dop853.c)
myb_add_test(test_ias15 SOURCES unit/test_ias15.c)
myb_add_test(test_symplectic SOURCES unit/test_symplectic.c)
myb_add_test(test_gauss_jackson SOURCES unit/test_gauss_jackson.c)
myb_add_test(test_rkn12 SOURCES unit/test_rkn12.c)
myb_add_test(test_repro SOURCES repro/test_repro.c)

# IAU SOFA reference through ERFA (BSD-3-Clause), compiled for the tests only.
file(GLOB _erfa_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa/*.c")
add_library(erfa_reference STATIC ${_erfa_sources})
target_include_directories(erfa_reference SYSTEM PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa")
myb_apply_flags(erfa_reference THIRD_PARTY)
myb_add_test(test_erfa_time SOURCES erfa/test_erfa_time.c LIBRARIES erfa_reference)

# MYB_CHECK must print one line and terminate the process.
add_executable(check_death unit/check_death.c)
target_link_libraries(check_death PRIVATE myb)
target_include_directories(check_death PRIVATE "${PROJECT_SOURCE_DIR}/src")
myb_apply_flags(check_death)
add_test(NAME test_check_death
  COMMAND "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:check_death>"
          "-DEXPECT=check_death\\.c:[0-9]+: main: CHECK failed: 1 \\+ 1 == 3"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/expect_death.cmake")

# The library must not import libm functions whose results differ between C runtimes.
get_filename_component(_myb_compiler_dir "${CMAKE_C_COMPILER}" DIRECTORY)
find_program(MYB_LLVM_NM NAMES llvm-nm HINTS "${_myb_compiler_dir}")
if(MYB_LLVM_NM)
  add_test(NAME test_no_libm_imports
    COMMAND "${CMAKE_COMMAND}" "-DNM=${MYB_LLVM_NM}" "-DLIB=$<TARGET_FILE:myb>"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_imports.cmake")
else()
  message(WARNING "llvm-nm not found: test_no_libm_imports is not registered")
endif()

# Correct rounding and error bounds against MPFR (GNU MPFR, arbitrary precision).
# Static archives first: a test must not depend on DLLs being on PATH.
find_path(MYB_MPFR_INCLUDE_DIR mpfr.h)
find_library(MYB_MPFR_LIBRARY NAMES libmpfr.a mpfr)
find_library(MYB_GMP_LIBRARY NAMES libgmp.a gmp)
if(MYB_MPFR_INCLUDE_DIR AND MYB_MPFR_LIBRARY AND MYB_GMP_LIBRARY)
  set(_mpfr_libs "${MYB_MPFR_LIBRARY}" "${MYB_GMP_LIBRARY}")
  if(WIN32 AND MYB_MPFR_LIBRARY MATCHES "\\.a$")
    # MSYS2 builds MPFR thread-safe on top of winpthreads.
    find_library(MYB_WINPTHREAD_LIBRARY NAMES libwinpthread.a winpthread)
    if(MYB_WINPTHREAD_LIBRARY)
      list(APPEND _mpfr_libs "${MYB_WINPTHREAD_LIBRARY}")
    endif()
  endif()
  myb_add_test(test_mpfr_cmath SOURCES mpfr/test_mpfr_cmath.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_num SOURCES mpfr/test_mpfr_num.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_parse SOURCES mpfr/test_mpfr_parse.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_kepler SOURCES mpfr/test_mpfr_kepler.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_radau SOURCES mpfr/test_mpfr_radau.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_yoshida SOURCES mpfr/test_mpfr_yoshida.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_gauss_jackson SOURCES mpfr/test_mpfr_gauss_jackson.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_rkn12 SOURCES mpfr/test_mpfr_rkn12.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
else()
  message(WARNING "MPFR/GMP not found: correct-rounding tests are NOT built "
    "(MYB_MPFR_INCLUDE_DIR=${MYB_MPFR_INCLUDE_DIR}, MYB_MPFR_LIBRARY=${MYB_MPFR_LIBRARY})")
endif()
