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

include_directories(${MoFEM_PROJECT_SOURCE_DIR}/users_modules/gels/src)

add_executable(gel_constitutive_equation_test
  ${UM_SOURCE_DIR}/gels/atom_tests/gel_constitutive_equation_test.cpp
)

target_link_libraries(gel_constitutive_equation_test
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

add_executable(gel_jacobian_test
  ${UM_SOURCE_DIR}/gels/atom_tests/gel_jacobian_test.cpp
)

target_link_libraries(gel_jacobian_test
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

add_test(gel_constitutive_equation_atom_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/gel_constitutive_equation_test
)

cm_export_file("../../basic_finite_elements/atom_tests/fluid_presssure_cub_test.h5m" export_atom_gel_tests)

add_test(gel_jacobian_atom_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/gel_jacobian_test
  -my_file fluid_presssure_cub_test.h5m
)
