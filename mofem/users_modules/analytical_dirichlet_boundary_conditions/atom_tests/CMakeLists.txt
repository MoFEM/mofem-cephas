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

cm_export_file("sphere.jou" export_analytical_bc_atom_tests)
cm_export_file("sphere.cub" export_analytical_bc_atom_tests)

add_executable(thermal_with_analytical_bc
  ${UM_SOURCE_DIR}/analytical_dirichlet_boundary_conditions/atom_tests/thermal_with_analytical_bc.cpp)

target_link_libraries(thermal_with_analytical_bc
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

add_test(thermal_with_analytical_bc_atom ${MoFEM_MPI_RUN} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/thermal_with_analytical_bc
  -my_file sphere.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps -ksp_monitor)
