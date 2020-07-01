# MoFEM is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# MoFEM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

#multi_indices
set(MOFEM_INTERFACES_SRC
  impl/MatPartitioningApply_Parmetis_MoFEM.cpp
  impl/All.cpp
  impl/AllInterfaces.cpp
  impl/LogManager.cpp
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
add_dependencies(mofem_interfaces install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_interfaces Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
