# Bindings: the library in other languages. Each is a thin face over the C interface — no
# arithmetic of its own, so a result has the bits of the same calls made in C.
#
#   cpp/     a C++17 header, myb/myb.hpp: target myb::cpp
#   python/  a package over the shared library (cffi)
#   rust/    the crates myb-sys (declarations) and myb (ownership, outcomes)

# ---- C++ ---------------------------------------------------------------------------------------
# Header-only: ownership (RAII), outcomes (myb::result<T>), lambdas for terms and events.
add_library(myb_cpp INTERFACE)
add_library(myb::cpp ALIAS myb_cpp)
set_target_properties(myb_cpp PROPERTIES EXPORT_NAME cpp)
target_include_directories(myb_cpp INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(myb_cpp INTERFACE myb)
if(CMAKE_CXX_COMPILER)
  target_compile_features(myb_cpp INTERFACE cxx_std_17)
endif()
