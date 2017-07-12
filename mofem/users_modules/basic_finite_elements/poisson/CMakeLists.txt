#
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

include_directories(
  ${UM_SOURCE_DIR}/basic_finite_elements/poisson
)
include_directories(
  ${UM_SOURCE_DIR}/basic_finite_elements/poisson/src
)

bfe_add_executable(analytical_poisson
  ${CMAKE_CURRENT_SOURCE_DIR}/analytical_poisson.cpp
)
target_link_libraries(analytical_poisson
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

cm_export_file("cube_1part.h5m" export_files_poisson)
cm_export_file("cube_2part.h5m" export_files_poisson)
cm_export_file(".petscrc" export_files_poisson)

add_test(analytical_poisson
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/analytical_poisson
  -file_name cube_2part.h5m -order 3
  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps
  -ksp_atol 1e-12 -ksp_rtol 0  -ksp_monitor
  -test 1
)
