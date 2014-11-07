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
cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/meshes/elasticity_atom_test_01X.cub" export_atom_tests)
cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/meshes/elasticity_atom_test_01Y.cub" export_atom_tests)
cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/meshes/elasticity_atom_test_01Z.cub" export_atom_tests)

add_executable(elasticity_atom_test ${UM_SOURCE_DIR}/atom_tests/elasticity_atom_test.cpp)
target_link_libraries(elasticity_atom_test   
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  users_modules_obsolete
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

#elasticity
add_test(elasticity_atom_test_01X ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01X.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -my_order 5)
add_test(elasticity_atom_test_01X_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01X.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01X.cub.txt)
add_test(elasticity_atom_test_01Y ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Y.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -my_order 5)
add_test(elasticity_atom_test_01Y_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Y.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Y.cub.txt)
add_test(elasticity_atom_test_01Z ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Z.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -my_order 5)
add_test(elasticity_atom_test_01Z_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Z.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Z.cub.txt)

add_test(elasticity_atom_test_01X_pressure 
  ${MoFEM_MPI_RUN} -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01X_pressure.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -my_order 5)
add_test(elasticity_atom_test_01X_pressure_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01X_pressure.cub.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01X_pressure.cub.txt)

#nonlinear elasticty
cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/meshes/elasticity_atom_test_01X_pressure.cub" export_atom_tests)
cm_export_absolut_path_file("${MoFEM_PROJECT_SOURCE_DIR}/meshes/elasticity_atom_test_sheer_force.cub" export_atom_tests)

add_executable(nonlinear_elasticity_atom ${UM_SOURCE_DIR}/atom_tests/nonlinear_elasticity_atom.cpp)
target_link_libraries(nonlinear_elasticity_atom   
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  users_modules_obsolete
  complex_for_lazy_obsolete
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

add_test(nonlinear_elasticity_atom_force_01_X
  ${MoFEM_MPI_RUN} -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_atom 
  -my_file elasticity_atom_test_01X.cub  
  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 
  -snes_monitor -snes_type newtonls 
  -snes_linesearch_type l2 
  -snes_max_it 10 -snes_atol 1e-8 
  -snes_rtol 1e-8 -my_order 4)
add_test(nonlinear_elasticity_atom_force_01_X_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_elasticity_atom_test_01X.cub.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_elasticity_atom_test_01X.cub.txt)

add_test(nonlinear_elasticity_atom_pressure 
  ${MoFEM_MPI_RUN} -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_atom 
  -my_file elasticity_atom_test_01X_pressure.cub  
  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 
  -snes_monitor -snes_type newtonls 
  -snes_linesearch_type l2 
  -snes_max_it 10 -snes_atol 1e-8 
  -snes_rtol 1e-8 -my_order 4)
add_test(nonlinear_elasticity_atom_pressure_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_elasticity_atom_test_01X_pressure.cub.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_elasticity_atom_test_01X_pressure.cub.txt)

add_test(nonlinear_elasticity_atom_sheer_force 
  ${MoFEM_MPI_RUN} -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_atom 
  -my_file elasticity_atom_test_sheer_force.cub 
  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 
  -snes_monitor -snes_type newtonls -snes_linesearch_type l2 
  -snes_max_it 10 -snes_atol 1e-8 -snes_rtol 1e-8 -my_order 4)
add_test(nonlinear_elasticity_atom_sheer_force_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_elasticity_atom_test_sheer_force.cub.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_elasticity_atom_test_sheer_force.cub.txt)

