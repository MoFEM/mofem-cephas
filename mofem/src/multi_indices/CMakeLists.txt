

# multi_indices
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(mofem_multi_indices
  impl/All.cpp
)
target_link_libraries(mofem_multi_indices PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_multi_indices PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_multi_indices install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_multi_indices Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)
install(TARGETS mofem_multi_indices DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
