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

file(
  COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/README
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
file(
  COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/l-shape.msh
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
file(
  COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/l-shape.geo
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
file(
  COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/bc.cfg
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
file(
  COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/unsaturated.cfg
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)


include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/atom_tests)

# For convinience source files from MoFEM project are used. If user is
# implementing independent application it will be using instead project
# directories.

bfe_add_executable(
  h_adaptive_transport ${CMAKE_CURRENT_SOURCE_DIR}/h_adaptive_transport.cpp
)
target_link_libraries(h_adaptive_transport
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS}
)

if(ADOL-C_LIBRARY)

  bfe_add_executable(
    unsaturated_transport ${CMAKE_CURRENT_SOURCE_DIR}/unsaturated_transport.cpp
  )
  target_link_libraries(unsaturated_transport
    users_modules
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    mofem_cblas
    ${MoFEM_PROJECT_LIBS}
  )

endif(ADOL-C_LIBRARY)
