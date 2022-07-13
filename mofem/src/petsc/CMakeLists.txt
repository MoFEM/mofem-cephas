

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)

# multi_indices
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(mofem_petsc
  impl/All.cpp
)
target_link_libraries(mofem_petsc PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_petsc PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_petsc install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_petsc Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_petsc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
