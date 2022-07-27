include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)

#finite elements
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(
  mofem_finite_elements 
  impl/All.cpp
)
target_link_libraries(mofem_finite_elements PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_finite_elements PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_finite_elements install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_finite_elements Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_finite_elements DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
