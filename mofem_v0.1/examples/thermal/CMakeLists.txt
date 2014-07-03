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


cm_export_file("../../meshes/simple_thermal_problem.cub" export_thremal_example)
cm_export_file("../../meshes/heat_exchange.cub" export_thremal_example)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable(thermal_steady thermal_steady.cpp)
target_link_libraries(thermal_steady mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(thermal_steady ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_thermal_elem_atom 
#  -my_file simple_thermal_problem.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_monitor)

# add test test for temerature unsteady element
add_executable(thermal_unseady thermal_unsteady.cpp)
target_link_libraries(thermal_unseady mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_thermal_elem_unsteady_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_thermal_elem_unsteady_atom 
#  -my_file simple_thermal_problem.cub 
#  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist 
#  -snes_type newtonls -snes_linesearch_type basic 
#  -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 
#  -ts_monitor -ts_type beuler -ts_dt 1 -ts_final_time 2  -ksp_monitor -snes_max_it 3 -snes_monitor)

