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
myb_add_test(test_linalg SOURCES unit/test_linalg.c)
myb_add_test(test_elliptic SOURCES unit/test_elliptic.c)
myb_add_test(test_attitude SOURCES unit/test_attitude.c)
myb_add_test(test_regular SOURCES unit/test_regular.c)
myb_add_test(test_thrust SOURCES unit/test_thrust.c)
myb_add_test(test_uncert SOURCES unit/test_uncert.c)
myb_add_test(test_cr3bp SOURCES unit/test_cr3bp.c)
myb_add_test(test_transfer SOURCES unit/test_transfer.c)
myb_add_test(test_spacecraft SOURCES unit/test_spacecraft.c)
myb_add_test(test_parse SOURCES unit/test_parse.c)
myb_add_test(test_format SOURCES unit/test_format.c)
myb_add_test(test_table SOURCES unit/test_table.c)
myb_add_test(test_ccsds SOURCES unit/test_ccsds.c)
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
myb_add_test(test_pool SOURCES unit/test_pool.c)
myb_add_test(test_nbody SOURCES unit/test_nbody.c)
myb_add_test(test_whfast SOURCES unit/test_whfast.c)
myb_add_test(test_spk SOURCES unit/test_spk.c)
myb_add_test(test_spk_writer SOURCES unit/test_spk_writer.c)
myb_add_test(test_pck_writer SOURCES unit/test_pck_writer.c)
myb_add_test(test_sp3 SOURCES unit/test_sp3.c)
myb_add_test(test_testpo SOURCES unit/test_testpo.c)
myb_add_test(test_eop SOURCES unit/test_eop.c)
myb_add_test(test_wgccre SOURCES unit/test_wgccre.c)
myb_add_test(test_earth_frames SOURCES unit/test_earth_frames.c)
myb_add_test(test_tkframe SOURCES unit/test_tkframe.c)
myb_add_test(test_rel SOURCES unit/test_rel.c)
myb_add_test(test_grav SOURCES unit/test_grav.c)
myb_add_test(test_force SOURCES unit/test_force.c)
myb_add_test(test_ellipsoidal SOURCES unit/test_ellipsoidal.c)
myb_add_test(test_lageos SOURCES unit/test_lageos.c)
myb_add_test(test_de440_rel SOURCES unit/test_de440_rel.c)
myb_add_test(test_moon_libration SOURCES unit/test_moon_libration.c)
myb_add_test(test_pn_binary SOURCES unit/test_pn_binary.c)
myb_add_test(test_pn_higher SOURCES unit/test_pn_higher.c)
myb_add_test(test_observe SOURCES unit/test_observe.c LIBRARIES erfa_reference)
myb_add_test(test_charge SOURCES unit/test_charge.c)
myb_add_test(test_plasma SOURCES unit/test_plasma.c)
myb_add_test(test_charging SOURCES unit/test_charging.c)
myb_add_test(test_gradient SOURCES unit/test_gradient.c)
myb_add_test(test_env SOURCES unit/test_env.c)
myb_add_test(test_orbit SOURCES unit/test_orbit.c)
myb_add_test(test_estim SOURCES unit/test_estim.c)
myb_add_test(test_events SOURCES unit/test_events.c)
myb_add_test(test_geometry SOURCES unit/test_geometry.c)
myb_add_test(test_nongrav SOURCES unit/test_nongrav.c)
myb_add_test(test_surface SOURCES unit/test_surface.c)
myb_add_test(test_shadow SOURCES unit/test_shadow.c)
myb_add_test(test_repro SOURCES repro/test_repro.c LIBRARIES myb_sgp4)

# IAU SOFA reference through ERFA (BSD-3-Clause), compiled for the tests only.
file(GLOB _erfa_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa/*.c")
add_library(erfa_reference STATIC ${_erfa_sources})
target_include_directories(erfa_reference SYSTEM PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa")
myb_apply_flags(erfa_reference THIRD_PARTY)
myb_add_test(test_erfa_time SOURCES erfa/test_erfa_time.c LIBRARIES erfa_reference)
myb_add_test(test_erfa_frames SOURCES erfa/test_erfa_frames.c LIBRARIES erfa_reference)
myb_add_test(test_erfa_equinox SOURCES erfa/test_erfa_equinox.c LIBRARIES erfa_reference)
myb_add_test(test_tides SOURCES unit/test_tides.c LIBRARIES erfa_reference)

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

# The stack budget of the library: 4 KB per frame when optimised (the deepest chains, e.g. a
# variational step through a force model, stay within a few tens of KB), 8 KB at -O0.
if(MYB_STACK_USAGE)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(_myb_stack_budget 8192)
  else()
    set(_myb_stack_budget 4096)
  endif()
  add_test(NAME test_stack_budget
    COMMAND "${CMAKE_COMMAND}"
            "-DDIRS=${PROJECT_BINARY_DIR}/CMakeFiles/myb.dir;${PROJECT_BINARY_DIR}/CMakeFiles/myb_coremath.dir"
            "-DBUDGET=${_myb_stack_budget}" -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_stack.cmake")
endif()

# The installed library serves an outside project: examples/ configured with
# find_package(myb) against a scratch prefix, built and run. Not under the sanitizers, whose
# runtime the outside project would have to link as well.
if(MYB_INSTALL AND MYB_SANITIZE STREQUAL "")
  set(_myb_install_args
    "-DBUILD_DIR=${PROJECT_BINARY_DIR}" "-DEXAMPLES_DIR=${PROJECT_SOURCE_DIR}/examples"
    "-DGENERATOR=${CMAKE_GENERATOR}" "-DTOOLCHAIN=${CMAKE_TOOLCHAIN_FILE}"
    "-DMAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}" "-DCONFIG=$<CONFIG>"
    "-DLIBDIR=${CMAKE_INSTALL_LIBDIR}" "-DINCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}"
    "-DDATADIR=${CMAKE_INSTALL_DATADIR}")
  if(CMAKE_CXX_COMPILER)
    list(APPEND _myb_install_args -DREQUIRE_CXX=ON) # the C++ example must build and pass too
  endif()
  add_test(NAME test_install
    COMMAND "${CMAKE_COMMAND}" ${_myb_install_args}
            "-DWORK_DIR=${PROJECT_BINARY_DIR}/install-test"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_install.cmake")
  set_tests_properties(test_install PROPERTIES RESOURCE_LOCK myb_install)
  if(MYB_BUILD_SHARED)
    add_test(NAME test_install_shared
      COMMAND "${CMAKE_COMMAND}" ${_myb_install_args} -DSHARED=ON
              "-DWORK_DIR=${PROJECT_BINARY_DIR}/install-test-shared"
              -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_install.cmake")
    set_tests_properties(test_install_shared PROPERTIES RESOURCE_LOCK myb_install)
  endif()
endif()

# The C++ face of the library (bindings/cpp): the classes against the same calls made in C.
if(CMAKE_CXX_COMPILER)
  add_executable(test_cpp_binding bindings/test_cpp.cc)
  set_target_properties(test_cpp_binding PROPERTIES
    CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  target_link_libraries(test_cpp_binding PRIVATE myb_cpp)
  target_compile_definitions(test_cpp_binding PRIVATE
    "MYB_TEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\""
    "MYB_DATA_DIR=\"${PROJECT_SOURCE_DIR}/data\"")
  myb_apply_flags(test_cpp_binding)
  # the header handles no exceptions: it compiles with them switched off
  target_compile_options(test_cpp_binding PRIVATE -fno-exceptions)
  if(WIN32)
    # the C++ runtime inside the test: it must not look for DLLs on PATH
    target_link_options(test_cpp_binding PRIVATE -static)
    target_link_libraries(test_cpp_binding PRIVATE winpthread) # for the static libstdc++
  endif()
  add_test(NAME test_cpp_binding COMMAND test_cpp_binding)
endif()

# The Python package (bindings/python) over the shared library: needs an interpreter with cffi
# (MYB_PYTHON names one; otherwise the first of those found that has it). Not with a
# sanitizer: its runtime would have to be loaded into the interpreter first.
if(MYB_BUILD_SHARED AND MYB_SANITIZE STREQUAL "")
  set(MYB_PYTHON "" CACHE FILEPATH "Python interpreter (with cffi) for the test of bindings/python")
  set(_myb_python_candidates "${MYB_PYTHON}")
  if(NOT MYB_PYTHON)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    # the interpreters of PATH as well: a toolchain's own Python often has no cffi
    find_program(MYB_PYTHON_ON_PATH NAMES python3 python NAMES_PER_DIR
      NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
      NO_PACKAGE_ROOT_PATH)
    mark_as_advanced(MYB_PYTHON_ON_PATH)
    set(_myb_python_candidates "${Python3_EXECUTABLE}" "${MYB_PYTHON_ON_PATH}")
  endif()
  set(_myb_python_ok FALSE)
  foreach(_myb_python IN LISTS _myb_python_candidates)
    if(_myb_python AND EXISTS "${_myb_python}" AND NOT _myb_python_ok)
      execute_process(COMMAND "${_myb_python}" -c "import cffi"
        RESULT_VARIABLE _myb_cffi_rc OUTPUT_QUIET ERROR_QUIET)
      if(_myb_cffi_rc EQUAL 0)
        set(_myb_python_ok TRUE)
        set(Python3_EXECUTABLE "${_myb_python}")
      endif()
    endif()
  endforeach()
  if(_myb_python_ok)
    add_test(NAME test_python_binding
      COMMAND "${CMAKE_COMMAND}" -E env
        "PYTHONPATH=${PROJECT_SOURCE_DIR}/bindings/python"
        "MYB_LIBRARY=$<TARGET_FILE:myb_shared>"
        "MYB_LIBRARY_PATH=$<TARGET_FILE_DIR:myb_sgp4_shared>"
        "MYB_DATA_DIR=${PROJECT_SOURCE_DIR}/data"
        "MYB_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/data"
        "MYB_CLANG=${CMAKE_C_COMPILER}"
        "${Python3_EXECUTABLE}" -m unittest discover -s
        "${PROJECT_SOURCE_DIR}/bindings/python/tests")
  else()
    message(STATUS "no Python with cffi found (set MYB_PYTHON): test_python_binding is not "
      "registered")
  endif()
endif()

# The Rust crates (bindings/rust) over the static library: needs cargo. On Windows the library
# is built for the GNU ABI, so the Rust toolchain of that ABI is asked for. Not with a
# sanitizer: Rust would have to link its runtime.
if(MYB_SANITIZE STREQUAL "")
  find_program(MYB_CARGO cargo)
  mark_as_advanced(MYB_CARGO)
  set(_myb_cargo_ok FALSE)
  set(_myb_cargo_toolchain "")
  if(MYB_CARGO)
    if(MINGW)
      set(_myb_cargo_toolchain "+stable-x86_64-pc-windows-gnu")
    endif()
    execute_process(COMMAND "${MYB_CARGO}" ${_myb_cargo_toolchain} --version
      RESULT_VARIABLE _myb_cargo_rc OUTPUT_QUIET ERROR_QUIET)
    if(_myb_cargo_rc EQUAL 0)
      set(_myb_cargo_ok TRUE)
    endif()
  endif()
  if(_myb_cargo_ok)
    # the directories of the libraries, separated as PATH is
    if(CMAKE_HOST_WIN32)
      set(_myb_path_separator "\;")
    else()
      set(_myb_path_separator ":")
    endif()
    add_test(NAME test_rust_binding
      COMMAND "${CMAKE_COMMAND}" -E env
        "MYB_LIB_DIR=$<TARGET_FILE_DIR:myb>${_myb_path_separator}$<TARGET_FILE_DIR:myb_sgp4>"
        "CARGO_TARGET_DIR=${PROJECT_BINARY_DIR}/rust"
        "MYB_DATA_DIR=${PROJECT_SOURCE_DIR}/data"
        "MYB_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/data"
        "${MYB_CARGO}" ${_myb_cargo_toolchain} test --offline --locked
        --manifest-path "${PROJECT_SOURCE_DIR}/bindings/rust/Cargo.toml" --features myb/sgp4)
    set_tests_properties(test_rust_binding PROPERTIES TIMEOUT 900)
  else()
    message(STATUS "cargo not found: test_rust_binding is not registered")
  endif()
endif()

# The shared libraries export the public API and nothing else.
if(MYB_BUILD_SHARED)
  find_program(MYB_LLVM_READOBJ NAMES llvm-readobj HINTS "${_myb_compiler_dir}")
  if(MYB_LLVM_NM AND (MYB_LLVM_READOBJ OR NOT WIN32))
    add_test(NAME test_exports
      COMMAND "${CMAKE_COMMAND}" "-DLIB=$<TARGET_FILE:myb_shared>" "-DPREFIX=^myb_"
              "-DNM=${MYB_LLVM_NM}" "-DREADOBJ=${MYB_LLVM_READOBJ}" -DMINIMUM=400
              "-DMAP=${PROJECT_SOURCE_DIR}/cmake/myb.map"
              -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_exports.cmake")
    get_property(_myb_ext_modules GLOBAL PROPERTY MYB_EXT_MODULES)
    foreach(_ext IN LISTS _myb_ext_modules)
      add_test(NAME test_exports_${_ext}
        COMMAND "${CMAKE_COMMAND}" "-DLIB=$<TARGET_FILE:myb_${_ext}_shared>" "-DPREFIX=^myb_"
                "-DNM=${MYB_LLVM_NM}" "-DREADOBJ=${MYB_LLVM_READOBJ}" -DMINIMUM=1
                "-DMAP=${PROJECT_SOURCE_DIR}/ext/${_ext}/myb_${_ext}.map"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tools/check_exports.cmake")
    endforeach()
  else()
    message(WARNING "llvm-nm or llvm-readobj not found: test_exports is not registered")
  endif()
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
  myb_add_test(test_mpfr_wh SOURCES mpfr/test_mpfr_wh.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_rkn12 SOURCES mpfr/test_mpfr_rkn12.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_frames SOURCES mpfr/test_mpfr_frames.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
  myb_add_test(test_mpfr_grav SOURCES mpfr/test_mpfr_grav.c
    LIBRARIES ${_mpfr_libs} INCLUDES "${MYB_MPFR_INCLUDE_DIR}")
else()
  message(WARNING "MPFR/GMP not found: correct-rounding tests are NOT built "
    "(MYB_MPFR_INCLUDE_DIR=${MYB_MPFR_INCLUDE_DIR}, MYB_MPFR_LIBRARY=${MYB_MPFR_LIBRARY})")
endif()

# ---- Optional modules (ext/) ----------------------------------------------------------------
if(TARGET myb_nrlmsise00)
  myb_add_test(test_nrlmsise00 SOURCES ext/test_nrlmsise00.c LIBRARIES myb_nrlmsise00)
endif()
if(TARGET myb_sgp4)
  myb_add_test(test_sgp4 SOURCES ext/test_sgp4.c LIBRARIES myb_sgp4)
endif()
if(TARGET myb_ionosphere)
  myb_add_test(test_ionosphere SOURCES ext/test_ionosphere.c LIBRARIES myb_ionosphere)
endif()
if(TARGET myb_troposphere)
  myb_add_test(test_troposphere SOURCES ext/test_troposphere.c LIBRARIES myb_troposphere)
endif()
