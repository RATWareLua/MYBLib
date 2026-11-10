# Examples of myb: small programs that use only the public headers.
#
# In the tree:            cmake --preset <preset> -DMYB_BUILD_EXAMPLES=ON
# Against an installation: cmake -S examples -B build-examples -DCMAKE_PREFIX_PATH=<prefix>
# (any C17 compiler: only the library itself needs Clang). MYB_EXAMPLES_SHARED=ON links the
# shared library instead of the static one. With a C++ compiler the C++ example is built too
# (the header myb/myb.hpp, target myb::cpp).
cmake_minimum_required(VERSION 3.28)

if(NOT TARGET myb)
  project(myb_examples LANGUAGES C)
  find_package(myb 1.2 REQUIRED)
  set(_myb_examples_data "${myb_DATA_DIR}")
else()
  set(_myb_examples_data "${PROJECT_SOURCE_DIR}/data")
endif()

option(MYB_EXAMPLES_SHARED "Link the examples with the shared library" OFF)
if(MYB_EXAMPLES_SHARED)
  set(_myb_examples_library myb::myb_shared)
else()
  set(_myb_examples_library myb::myb)
endif()

foreach(_example IN ITEMS low_orbit time_scales)
  add_executable(myb_example_${_example} ${_example}.c)
  set_target_properties(myb_example_${_example} PROPERTIES
    C_STANDARD 17 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
  target_link_libraries(myb_example_${_example} PRIVATE ${_myb_examples_library})
  if(NOT WIN32)
    target_link_libraries(myb_example_${_example} PRIVATE m) # the examples' own sqrt, cos
  endif()
  target_compile_definitions(myb_example_${_example} PRIVATE
    "MYB_EXAMPLE_DATA_DIR=\"${_myb_examples_data}\"")
endforeach()

# The C++ face of the library: the header is found through myb::cpp; the library linked is
# the one chosen above.
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
  enable_language(CXX)
  add_executable(myb_example_low_orbit_cpp low_orbit.cpp)
  set_target_properties(myb_example_low_orbit_cpp PROPERTIES
    CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  get_target_property(_myb_cpp_include myb::cpp INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(myb_example_low_orbit_cpp PRIVATE ${_myb_cpp_include})
  target_link_libraries(myb_example_low_orbit_cpp PRIVATE ${_myb_examples_library})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    # the header handles no exceptions: the example is built with them switched off
    target_compile_options(myb_example_low_orbit_cpp PRIVATE -fno-exceptions)
  endif()
  if(MINGW)
    # the C++ runtime inside the program: it does not look for libstdc++ on PATH
    target_link_options(myb_example_low_orbit_cpp PRIVATE -static-libstdc++ -static-libgcc)
    target_link_libraries(myb_example_low_orbit_cpp PRIVATE
      -Wl,-Bstatic winpthread -Wl,-Bdynamic)
  endif()
endif()
