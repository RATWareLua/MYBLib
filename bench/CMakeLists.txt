# Benchmarks of the hot paths (MYB_BENCH): the program myb_bench prints a table of median times
# (bench/bench.c lists the cases). With the tests, "bench_smoke" runs its short form (small
# sizes, no downloaded data, a few seconds) so that the cases keep building and running.

add_executable(myb_bench bench.c)
target_link_libraries(myb_bench PRIVATE myb)
target_include_directories(myb_bench PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_compile_definitions(myb_bench PRIVATE
  "MYB_DATA_DIR=\"${PROJECT_SOURCE_DIR}/data\""
  # C99-conformant printf (%zu, %lld) from mingw-w64 instead of msvcrt's
  $<$<BOOL:${WIN32}>:__USE_MINGW_ANSI_STDIO=1>)
myb_apply_flags(myb_bench)

if(MYB_BUILD_TESTS)
  add_test(NAME bench_smoke COMMAND myb_bench --smoke)
  set_tests_properties(bench_smoke PROPERTIES LABELS bench TIMEOUT 300)
endif()
