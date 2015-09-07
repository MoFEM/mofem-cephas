# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
# It can be freely used for educational and research purposes
# by other institutions. If you use this softwre pleas cite my work.
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


#Elasticity atom tests

cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/homogenisation/meshes/Cube_atom_test.cub" export_homo_atom_tests)

add_executable(homogenisation_disp_atom_test ${UM_SOURCE_DIR}/homogenisation/atom_tests/homogenisation_disp_atom_test.cpp)
target_link_libraries(homogenisation_disp_atom_test
  user_modules_homogenisation
  users_modules_obsolete
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

#displacement boundary conditions
add_test(homogenisation_disp_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/homogenisation_disp_atom_test -my_file Cube_atom_test.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps -my_order 1 -myapplied_strain 1.0,0.0,0.0,0.0,0.0,0.0)


add_executable(homogenisation_trac_atom_test ${UM_SOURCE_DIR}/homogenisation/atom_tests/homogenisation_trac_atom_test.cpp)
target_link_libraries(homogenisation_trac_atom_test
  user_modules_homogenisation
  users_modules_obsolete
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

#displacement boundary conditions
add_test(homogenisation_trac_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/homogenisation_trac_atom_test -my_file Cube_atom_test.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps -my_order 1 -myapplied_strain 1.0,0.0,0.0,0.0,0.0,0.0)

add_executable(homogenisation_periodic_atom_test ${UM_SOURCE_DIR}/homogenisation/atom_tests/homogenisation_periodic_atom_test.cpp)
target_link_libraries(homogenisation_periodic_atom_test
  user_modules_homogenisation
  users_modules_obsolete
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

#displacement boundary conditions
add_test(homogenisation_periodic_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/homogenisation_periodic_atom_test -my_file Cube_atom_test.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package mumps -my_order 1 -myapplied_strain 1.0,0.0,0.0,0.0,0.0,0.0)
