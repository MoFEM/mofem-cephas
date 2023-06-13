include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/cblas)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/post_proc)
include_directories(${PROJECT_SOURCE_DIR}/src/boundary_conditions)

if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    ScalingMethod.cpp
    PROPERTIES
    COMPILE_FLAGS 
		"-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(
  mofem_boundary_conditions
  impl/ScalingMethod.cpp
  impl/EssentialDisplacementCubitBcData.cpp
  impl/EssentialTemperatureCubitBcData.cpp
)

target_link_libraries(
	mofem_boundary_conditions PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(
	mofem_boundary_conditions PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(
	mofem_boundary_conditions install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(
		mofem_boundary_conditions Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(
	TARGETS mofem_boundary_conditions DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)