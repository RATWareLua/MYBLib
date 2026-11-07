# expect_death.cmake — runs a program that must terminate abnormally (non-zero exit, e.g.
# __fastfail 0xC0000409 or SIGABRT) and print a line matching a regular expression on stderr.
# CTest itself reports such a process as an exception even when its output matches.
#
# Usage: cmake -DPROGRAM=<exe> -DEXPECT=<regex> -P expect_death.cmake
if(NOT PROGRAM OR NOT EXPECT)
  message(FATAL_ERROR "usage: cmake -DPROGRAM=<exe> -DEXPECT=<regex> -P expect_death.cmake")
endif()

execute_process(COMMAND "${PROGRAM}"
  OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
message(STATUS "exit: ${rc}")
message(STATUS "stderr: ${err}")
if(rc STREQUAL "0")
  message(FATAL_ERROR "the program exited normally; a failed check must terminate it")
endif()
if(NOT err MATCHES "${EXPECT}")
  message(FATAL_ERROR "stderr does not match '${EXPECT}'")
endif()
string(REGEX MATCHALL "\n" newlines "${err}")
list(LENGTH newlines line_count)
if(NOT line_count EQUAL 1)
  message(FATAL_ERROR "expected exactly one line on stderr, got ${line_count}")
endif()
message(STATUS "terminated with the expected one-line message")
