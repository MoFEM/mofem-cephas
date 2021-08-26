
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

# Tetsing files
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    impl/Hcurl.cpp
    PROPERTIES
    COMPILE_FLAGS "-DGENERATE_VTK_WITH_CURL_BASE -include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
  
else(PRECOMPILED_HEADRES)
  set_source_files_properties(
    impl/Hcurl.cpp
    PROPERTIES
    COMPILE_FLAGS "-DGENERATE_VTK_WITH_CURL_BASE"
  )
endif(PRECOMPILED_HEADRES)
add_executable(vtk_curl_base_on_tet impl/Hcurl.cpp)
if(PRECOMPILED_HEADRES)
  add_dependencies(vtk_curl_base_on_tet Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)
target_link_libraries(vtk_curl_base_on_tet
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${PROJECT_LIBS}
)
