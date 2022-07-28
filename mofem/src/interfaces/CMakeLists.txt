

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

#multi_indices
set(MOFEM_INTERFACES_SRC
  impl/MatPartitioningApply_Parmetis_MoFEM.cpp
  impl/Core.cpp
  impl/DeprecatedCoreInterface.cpp
  impl/FieldCore.cpp
  impl/FECore.cpp
  impl/ProblemsCore.cpp
  impl/DeleteCore.cpp
  impl/LoopMethods.cpp
  impl/LogManager.cpp
  impl/AllCoreInterfaces.cpp
  impl/AllInterfaces.cpp
)

if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
  set_source_files_properties(
    AllInterfaces.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
  set_source_files_properties(
    LogManager.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(
  mofem_interfaces
  ${MOFEM_INTERFACES_SRC}
)
set_target_properties(mofem_interfaces PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_interfaces install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_interfaces Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
