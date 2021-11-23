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
