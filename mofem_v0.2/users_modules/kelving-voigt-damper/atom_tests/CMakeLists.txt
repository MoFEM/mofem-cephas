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

include_directories(${UM_SOURCE_DIR}/kelving-voigt-damper/src)

add_executable(damper_jacobian_test
  ${UM_SOURCE_DIR}/kelving-voigt-damper/atom_tests/damper_jacobian_test.cpp
)

target_link_libraries(damper_jacobian_test
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

cm_export_file(
  "../../basic_finite_elements/atom_tests/fluid_presssure_cub_test.h5m" export_atom_damper_tests
)

add_test(damper_jacobian_atom_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/damper_jacobian_test
  -my_file fluid_presssure_cub_test.h5m
)
