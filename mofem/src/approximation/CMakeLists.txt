include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

# Approximation lib
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(
  mofem_approx
  impl/All.cpp
  c_impl/base_functions.c
  c_impl/fem_tools.c
  c_impl/h1.c
  c_impl/l2.c
)
target_link_libraries(mofem_approx PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_approx PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_approx install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_approx Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_approx DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# FIXME: That is source of problems VTK out should be move to separate file.

# # Tetsing files
# if(PRECOMPILED_HEADRES)
#   set_source_files_properties(
#     impl/Hcurl.cpp
#     PROPERTIES
#     COMPILE_FLAGS "-DGENERATE_VTK_WITH_CURL_BASE -include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
#   )
  
# else(PRECOMPILED_HEADRES)
#   set_source_files_properties(
#     impl/Hcurl.cpp
#     PROPERTIES
#     COMPILE_FLAGS "-DGENERATE_VTK_WITH_CURL_BASE"
#   )
# endif(PRECOMPILED_HEADRES)

# add_executable(vtk_curl_base_on_tet impl/Hcurl.cpp)
# if(PRECOMPILED_HEADRES)
#   add_dependencies(vtk_curl_base_on_tet Includes.hpp.pch_copy)
# endif(PRECOMPILED_HEADRES)
# target_link_libraries(vtk_curl_base_on_tet
#   mofem_finite_elements
#   mofem_interfaces
#   mofem_multi_indices
#   mofem_petsc
#   mofem_approx
#   mofem_third_party
#   ${PROJECT_LIBS}
# )
