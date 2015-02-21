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

include_directories("${UM_SOURCE_DIR}/elasticity/src")

add_executable(elasticity ${UM_SOURCE_DIR}/elasticity/elasticity.cpp)
target_link_libraries(elasticity 
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  boost_program_options
  ${MoFEM_PROJECT_LIBS})
  
  
add_executable(dynamics ${UM_SOURCE_DIR}/elasticity/dynamics.cpp)
target_link_libraries(dynamics
  users_modules_obsolete
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

cm_export_file("LShape.cub" export_files_elasticity)
cm_export_file("data.txt" export_files_elasticity)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#testing
add_test(ElasticityTest 
  ${MoFEM_MPI_RUN} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/elasticity 
  -my_file LShape.cub 
  -ksp_type fgmres 
  -pc_type lu 
  -pc_factor_mat_solver_package mumps
  -ksp_monitor -log_summary)

add_test(DynamicsTest
  ${MoFEM_MPI_RUN} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/dynamics
  -my_file LShape.cub
  -my_time_data_file data.txt
  -ksp_type fgmres 
  -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-10 
  -snes_monitor -snes_type newtonls -snes_linesearch_type basic 
  -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 
  -ts_monitor -ts_type alpha -ts_dt 0.1 -ts_final_time 1
  -log_summary)


