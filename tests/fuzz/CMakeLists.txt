# Fuzz harnesses of the readers of myb: one executable per reader, fuzz_<name>.
#
# With MYB_FUZZ and libFuzzer (Clang on Linux, preset linux-fuzz) a harness is a libFuzzer
# binary; tests/fuzz/run_fuzzers.cmake runs them for a time budget and keeps what they find.
# Everywhere else it links the replay main (replay.c). Either way every build registers one
# test per harness that runs its seeds — the files of tests/data and data/ named below and the
# hand-made cases and regressions of tests/fuzz/corpus/<name>/ — and, with the replay main,
# MYB_FUZZ_REPLAY_VARIANTS damaged variants of each (truncations and small damages).
#
# The seed lists are also written to <build>/tests/fuzz/fuzz_<name>.seeds for the scripts.

set(MYB_FUZZ_REPLAY_VARIANTS 24 CACHE STRING
  "Damaged variants of each seed that a fuzz test runs with the replay main")

set(_fuzz_corpus "${CMAKE_CURRENT_SOURCE_DIR}/corpus")
set(_fuzz_data "${PROJECT_SOURCE_DIR}/data")
set(_fuzz_test_data "${PROJECT_SOURCE_DIR}/tests/data")

# myb_add_fuzz(<name> [LIBRARIES <targets>...] SEEDS <files or directories>...)
#   The harness fuzz_<name>.c; its own corpus directory corpus/<name> is always a seed when it
#   exists. A directory stands for the files in it.
function(myb_add_fuzz name)
  cmake_parse_arguments(_f "" "" "LIBRARIES;SEEDS" ${ARGN})
  set(target fuzz_${name})
  foreach(library IN LISTS _f_LIBRARIES)
    if(NOT TARGET ${library})
      message(STATUS "${target} is not built: ${library} is not (MYB_BUILD_EXT)")
      return()
    endif()
  endforeach()
  if(MYB_FUZZ_LIBFUZZER)
    add_executable(${target} ${target}.c)
    target_link_options(${target} PRIVATE -fsanitize=fuzzer)
  else()
    add_executable(${target} ${target}.c replay.c)
  endif()
  target_link_libraries(${target} PRIVATE myb ${_f_LIBRARIES})
  target_include_directories(${target} PRIVATE
    "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
  # the kernels a harness builds its fixed models from (#embed <time/utc-tai.tk>)
  target_compile_options(${target} PRIVATE "--embed-dir=${_fuzz_data}")
  myb_apply_flags(${target})

  set(seeds "")
  set(sources ${_f_SEEDS})
  if(IS_DIRECTORY "${_fuzz_corpus}/${name}")
    list(APPEND sources "${_fuzz_corpus}/${name}")
  endif()
  foreach(seed IN LISTS sources)
    if(IS_DIRECTORY "${seed}")
      file(GLOB files CONFIGURE_DEPENDS LIST_DIRECTORIES false "${seed}/*")
      list(SORT files)
      list(APPEND seeds ${files})
    elseif(EXISTS "${seed}")
      list(APPEND seeds "${seed}")
    else()
      message(FATAL_ERROR "seed of fuzz_${name} not found: ${seed}")
    endif()
  endforeach()
  if(NOT seeds)
    message(FATAL_ERROR "fuzz_${name} has no seeds")
  endif()
  list(JOIN seeds "\n" text)
  file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${target}.seeds" CONTENT "${text}\n" @ONLY
    NEWLINE_STYLE UNIX)

  if(MYB_FUZZ_LIBFUZZER)
    # libFuzzer runs the files it is given once each
    add_test(NAME ${target} COMMAND ${target} ${seeds})
  else()
    add_test(NAME ${target} COMMAND ${target} -mutate=${MYB_FUZZ_REPLAY_VARIANTS} ${seeds})
  endif()
  set_tests_properties(${target} PROPERTIES LABELS fuzz TIMEOUT 900)
  set_property(GLOBAL APPEND PROPERTY MYB_TEST_PROGRAMS ${target})
endfunction()

myb_add_fuzz(number)
myb_add_fuzz(table SEEDS "${_fuzz_test_data}/weather/sw-excerpt.csv")
myb_add_fuzz(kpool SEEDS
  "${_fuzz_data}/time/time-scales.tk" "${_fuzz_data}/time/utc-tai.tk"
  "${_fuzz_data}/constants/codata-2022.tk" "${_fuzz_data}/frames/naif-inertial.tk"
  "${_fuzz_data}/frames/eop-tides.tk" "${_fuzz_test_data}/ck/test.tsc")
myb_add_fuzz(kernel_time SEEDS
  "${_fuzz_data}/time/utc-tai.tk" "${_fuzz_data}/time/time-scales.tk"
  "${_fuzz_data}/time/tdb-tt-fb1990.tk")
myb_add_fuzz(kernel_frames SEEDS
  "${_fuzz_data}/frames/naif-inertial.tk" "${_fuzz_data}/constants/codata-2022.tk"
  "${_fuzz_test_data}/wgccre/pck00011.tpc")
myb_add_fuzz(kernel_series SEEDS
  "${_fuzz_data}/frames/iau1980.tk" "${_fuzz_data}/frames/eop-tides.tk"
  "${_fuzz_data}/frames/iers2010-site.tk" "${_fuzz_data}/grav/iers2010-tides.tk"
  "${_fuzz_data}/frames/iau2006-cip.tk")
myb_add_fuzz(sclk SEEDS "${_fuzz_test_data}/ck/test.tsc")
myb_add_fuzz(spk LIBRARIES myb_sgp4 SEEDS
  "${_fuzz_test_data}/spk/t01.bsp" "${_fuzz_test_data}/spk/t02_03.bsp"
  "${_fuzz_test_data}/spk/t05.bsp" "${_fuzz_test_data}/spk/t08_09_12_13.bsp"
  "${_fuzz_test_data}/spk/t10.bsp" "${_fuzz_test_data}/spk/t14.bsp"
  "${_fuzz_test_data}/spk/t15.bsp" "${_fuzz_test_data}/spk/t17.bsp"
  "${_fuzz_test_data}/spk/t18_19.bsp" "${_fuzz_test_data}/spk/t20.bsp"
  "${_fuzz_test_data}/spk/chain_a.bsp" "${_fuzz_test_data}/spk/chain_b.bsp"
  "${_fuzz_test_data}/frames/graph.bsp")
myb_add_fuzz(pck SEEDS "${_fuzz_test_data}/frames/graph.bpc")
myb_add_fuzz(ck SEEDS "${_fuzz_test_data}/ck/test.bc" "${_fuzz_test_data}/ck/t5.bc")
myb_add_fuzz(icgem SEEDS
  "${_fuzz_test_data}/grav/EGM2008_n20.gfc" "${_fuzz_test_data}/grav/synthetic2.gfc")
myb_add_fuzz(ocean_tides SEEDS "${_fuzz_test_data}/grav/fes2004_excerpt.dat")
myb_add_fuzz(magnetic)
myb_add_fuzz(oem)
myb_add_fuzz(opm)
myb_add_fuzz(omm LIBRARIES myb_sgp4)
myb_add_fuzz(ocm)
myb_add_fuzz(sp3 SEEDS "${_fuzz_test_data}/sp3/ilrsa.lageos1.250607.excerpt.sp3")
myb_add_fuzz(tle LIBRARIES myb_sgp4 SEEDS "${_fuzz_test_data}/sgp4/SGP4-VER.TLE")
myb_add_fuzz(ionex LIBRARIES myb_ionosphere SEEDS
  "${_fuzz_test_data}/ionosphere/synthetic.inx" "${_fuzz_test_data}/ionosphere/synthetic3d.inx")
myb_add_fuzz(eop SEEDS
  "${_fuzz_test_data}/eop/eopc04-excerpt.txt" "${_fuzz_test_data}/eop/finals2000A-excerpt.txt")
myb_add_fuzz(space_weather SEEDS "${_fuzz_test_data}/weather/sw-excerpt.csv")
myb_add_fuzz(gpt3 LIBRARIES myb_troposphere)
myb_add_fuzz(blq SEEDS "${_fuzz_test_data}/site/loading.blq")
myb_add_fuzz(site_atmosphere SEEDS
  "${_fuzz_test_data}/site/atmosphere.atl" "${_fuzz_test_data}/site/atmosphere_sites.dat")
myb_add_fuzz(site_grid SEEDS
  "${_fuzz_test_data}/site/opole-excerpt.txt" "${_fuzz_test_data}/site/atmosphere_grid.dat")
myb_add_fuzz(jb2008 LIBRARIES myb_jb2008)
