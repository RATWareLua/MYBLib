# check_imports.cmake — the static library must not reference libm functions whose results
# depend on the C runtime (sin from UCRT, glibc and mingw-w64 differ; mingw's fma is wrong).
# Exact operations (sqrt, floor, fabs, ...) are correctly rounded by IEEE-754 everywhere and
# are allowed, although with -fno-math-errno they are normally inlined.
#
# Usage: cmake -DNM=<llvm-nm> -DLIB=<libmyb.a> -P check_imports.cmake
if(NOT NM OR NOT LIB)
  message(FATAL_ERROR "usage: cmake -DNM=<llvm-nm> -DLIB=<library> -P check_imports.cmake")
endif()

execute_process(COMMAND "${NM}" --undefined-only --just-symbol-name "${LIB}"
  OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "llvm-nm failed (${rc}): ${err}")
endif()

string(REPLACE "\n" ";" lines "${out}")
set(symbols "")
foreach(line IN LISTS lines)
  string(STRIP "${line}" line)
  # Skip empty lines and archive member headers ("member.o:").
  if(line STREQUAL "" OR line MATCHES ":$")
    continue()
  endif()
  list(APPEND symbols "${line}")
endforeach()
list(REMOVE_DUPLICATES symbols)
list(SORT symbols)

set(forbidden_regex "^(sin|cos|tan|sincos|asin|acos|atan|atan2|sinh|cosh|tanh|asinh|acosh|atanh|exp|exp2|exp10|expm1|log|log2|log10|log1p|pow|cbrt|hypot|fma|erf|erfc|tgamma|lgamma|fmod|remainder|remquo|strtod|atof)[fl]?$")
set(found "")
foreach(symbol IN LISTS symbols)
  if(symbol MATCHES "${forbidden_regex}")
    list(APPEND found "${symbol}")
  endif()
endforeach()

list(JOIN symbols " " all)
message(STATUS "undefined symbols of ${LIB}: ${all}")
if(found)
  list(JOIN found " " found_text)
  message(FATAL_ERROR "forbidden libm imports: ${found_text}")
endif()
message(STATUS "no forbidden libm imports")
