# check_stack.cmake — the stack budget of the library: every function's frame, as reported by
# Clang's -fstack-usage (one .su file per object: "file:line:column:function<TAB>bytes<TAB>
# static|dynamic|dynamic,bounded"), must be static (no VLA, no alloca) and at most BUDGET bytes.
# The library has no recursion, so the deepest stack is a finite sum of such frames.
#
# Usage: cmake -DDIRS=<dir;dir> -DBUDGET=<bytes> -P check_stack.cmake
if(NOT DIRS OR NOT BUDGET)
  message(FATAL_ERROR "usage: cmake -DDIRS=<object dirs> -DBUDGET=<bytes> -P check_stack.cmake")
endif()

set(files "")
foreach(dir IN LISTS DIRS)
  file(GLOB_RECURSE found "${dir}/*.su")
  list(APPEND files ${found})
endforeach()
list(LENGTH files nfiles)
if(nfiles EQUAL 0)
  message(FATAL_ERROR "no .su files under ${DIRS}: is MYB_STACK_USAGE on?")
endif()

set(over "")
set(largest 0)
set(largest_name "")
set(functions 0)
foreach(f IN LISTS files)
  file(STRINGS "${f}" lines)
  foreach(line IN LISTS lines)
    string(REPLACE "\t" ";" fields "${line}")
    list(LENGTH fields n)
    if(NOT n EQUAL 3)
      continue()
    endif()
    list(GET fields 0 where)
    list(GET fields 1 bytes)
    list(GET fields 2 kind)
    math(EXPR functions "${functions} + 1")
    string(REGEX REPLACE "^.*[/\\]" "" where "${where}")
    if(bytes GREATER largest)
      set(largest ${bytes})
      set(largest_name "${where}")
    endif()
    if(NOT kind STREQUAL "static" OR bytes GREATER BUDGET)
      list(APPEND over "${where}: ${bytes} bytes (${kind})")
    endif()
  endforeach()
endforeach()

message(STATUS "${functions} functions in ${nfiles} objects; largest frame ${largest} bytes "
               "(${largest_name}); budget ${BUDGET}")
if(over)
  string(REPLACE ";" "\n  " list "${over}")
  message(FATAL_ERROR "frames over the budget or dynamic:\n  ${list}")
endif()
