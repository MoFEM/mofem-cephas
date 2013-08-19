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

add_executable(hook_material_forces hook_material_forces.cpp)
target_link_libraries(hook_material_forces mofem_Field ${PROJECT_LIBS})

cm_export_file("../../meshes/LShape.cub" "/" export_files_hook_material_forces)
cm_export_file("../../meshes/LShapeCoarse.cub" "/" export_files_hook_material_forces)
cm_export_file("../../meshes/LShapeDense.cub" "/" export_files_hook_material_forces)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


#testing
add_test(SimpleHookMaterialForces ${MPI_BIN_PATH}/mpirun -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/hook_material_forces
  -my_file LShapeDense.cub
  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-6 
  -snes_monitor -snes_type ls -snes_linesearch_type l2 -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8
  -log_summary)

