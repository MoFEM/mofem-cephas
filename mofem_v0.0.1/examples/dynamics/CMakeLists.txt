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

add_executable(linear_dynamics_elasticity linear_dynamics_elasticity.cpp)
target_link_libraries(linear_dynamics_elasticity mofem_Field ${PROJECT_LIBS})

cm_export_file("../../meshes/LShape.cub" "/" export_files_linear_dynamics_elasticity)
cm_export_file("../../meshes/LShapeCoarse.cub" "/" export_linear_dynamics_elasticity)
cm_export_file("../../meshes/LShapeDense.cub" "/" export_linear_dynamics_elasticity)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#testing
add_test(LinearDynamicsElasticity ${MPI_BIN_PATH}/mpirun -np 2 ${CMAKE_CURRENT_BINARY_DIR}/linear_dynamics_elasticity
  -my_file LShapeCoarse.cub
  -ksp_type fgmres 
  -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-10 
  -snes_monitor -snes_type newtonls -snes_linesearch_type basic 
  -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 
  -ts_monitor -ts_type alpha -ts_dt 0.1 -ts_final_time 10
  -log_summary)

 

