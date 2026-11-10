# Optional modules of empirical models of particular bodies (docs/PLAN.md §3): separate targets,
# never part of libmyb.

add_library(myb_nrlmsise00 STATIC nrlmsise00/nrlmsise00.c)
add_library(myb::nrlmsise00 ALIAS myb_nrlmsise00)
target_link_libraries(myb_nrlmsise00 PUBLIC myb)
target_include_directories(myb_nrlmsise00
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/nrlmsise00"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
myb_apply_flags(myb_nrlmsise00)

add_library(myb_sgp4 STATIC sgp4/sgp4.c)
add_library(myb::sgp4 ALIAS myb_sgp4)
target_link_libraries(myb_sgp4 PUBLIC myb)
target_include_directories(myb_sgp4
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/sgp4"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
myb_apply_flags(myb_sgp4)

add_library(myb_troposphere STATIC troposphere/troposphere.c)
add_library(myb::troposphere ALIAS myb_troposphere)
target_link_libraries(myb_troposphere PUBLIC myb)
target_include_directories(myb_troposphere
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/troposphere"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
myb_apply_flags(myb_troposphere)

add_library(myb_ionosphere STATIC ionosphere/ionosphere.c)
add_library(myb::ionosphere ALIAS myb_ionosphere)
target_link_libraries(myb_ionosphere PUBLIC myb)
target_include_directories(myb_ionosphere
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/ionosphere"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
myb_apply_flags(myb_ionosphere)
