# replay.cmake — runs the files of directories (a corpus that libFuzzer grew, the crashes it
# kept) and single files through a harness, in any build: the replay main of a build without
# libFuzzer takes no directories (it reads only what it is named, through a list), a libFuzzer
# binary given files runs each once.
#
#   cmake -DPROGRAM=build/windows-release/tests/fuzz/fuzz_spk.exe -DINPUTS=dir1;file2
#         [-DARGS=-mutate=1000] [-DLIBFUZZER=ON] -P tests/fuzz/replay.cmake
#
# ARGS go to the program before the inputs (-mutate, -seed, -verbose of the replay main). The
# replay main gets the inputs as a list file (@list); a libFuzzer binary (LIBFUZZER=ON) gets
# them on its command line.
cmake_minimum_required(VERSION 3.28)

if(NOT PROGRAM OR NOT INPUTS)
  message(FATAL_ERROR "usage: cmake -DPROGRAM=<fuzz_name> -DINPUTS=<dirs and files> "
    "[-DARGS=...] -P replay.cmake")
endif()
set(_files "")
foreach(_input IN LISTS INPUTS)
  if(IS_DIRECTORY "${_input}")
    file(GLOB _inside LIST_DIRECTORIES false "${_input}/*")
    list(SORT _inside)
    list(APPEND _files ${_inside})
  elseif(EXISTS "${_input}")
    list(APPEND _files "${_input}")
  else()
    message(FATAL_ERROR "no such input: ${_input}")
  endif()
endforeach()
list(LENGTH _files _count)
if(_count EQUAL 0)
  message(FATAL_ERROR "no input files")
endif()
if(LIBFUZZER)
  execute_process(COMMAND "${PROGRAM}" ${ARGS} ${_files} RESULT_VARIABLE _rc)
else()
  # a list file keeps the command line short whatever the corpus
  string(RANDOM LENGTH 8 _tag)
  set(_list "${CMAKE_CURRENT_BINARY_DIR}/replay-${_tag}.list")
  list(JOIN _files "\n" _text)
  file(WRITE "${_list}" "${_text}\n")
  execute_process(COMMAND "${PROGRAM}" ${ARGS} "@${_list}" RESULT_VARIABLE _rc)
  file(REMOVE "${_list}")
endif()
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "the harness failed (${_rc}) on one of ${_count} inputs")
endif()
message(STATUS "${_count} inputs ran")
