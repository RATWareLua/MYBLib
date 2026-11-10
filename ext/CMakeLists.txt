# Optional modules of empirical models of particular bodies (docs/PLAN.en.md §3): separate targets,
# never part of libmyb.

# myb_add_ext(<name> <sources>...): the static library myb_<name> (alias and exported name
# myb::<name>) of ext/<name> and, with MYB_BUILD_SHARED, the shared one myb_<name>_shared
# (libmyb_<name>.so) on top of the shared libmyb.
function(myb_add_ext name)
  add_library(myb_${name} STATIC ${ARGN})
  add_library(myb::${name} ALIAS myb_${name})
  set_target_properties(myb_${name} PROPERTIES EXPORT_NAME ${name})
  target_link_libraries(myb_${name} PUBLIC myb)
  target_include_directories(myb_${name}
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${name}>"
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    PRIVATE "${PROJECT_SOURCE_DIR}/src")
  myb_apply_flags(myb_${name})
  if(MYB_BUILD_SHARED)
    add_library(myb_${name}_shared SHARED ${ARGN})
    add_library(myb::${name}_shared ALIAS myb_${name}_shared)
    set_target_properties(myb_${name}_shared PROPERTIES
      EXPORT_NAME ${name}_shared
      OUTPUT_NAME myb_${name}
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_VERSION_MAJOR}"
      C_VISIBILITY_PRESET hidden)
    target_link_libraries(myb_${name}_shared PUBLIC myb_shared)
    target_include_directories(myb_${name}_shared
      PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${name}>"
             "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
      PRIVATE "${PROJECT_SOURCE_DIR}/src")
    target_compile_definitions(myb_${name}_shared PRIVATE MYB_EXT_BUILDING)
    myb_apply_flags(myb_${name}_shared)
    myb_symbol_versions(myb_${name}_shared
      "${CMAKE_CURRENT_SOURCE_DIR}/${name}/myb_${name}.map")
  endif()
  set_property(GLOBAL APPEND PROPERTY MYB_EXT_MODULES ${name})
endfunction()

myb_add_ext(nrlmsise00 nrlmsise00/nrlmsise00.c)
myb_add_ext(sgp4 sgp4/sgp4.c)
myb_add_ext(troposphere troposphere/troposphere.c)
myb_add_ext(ionosphere ionosphere/ionosphere.c)
