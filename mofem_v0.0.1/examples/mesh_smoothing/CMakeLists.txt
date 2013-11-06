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

add_executable(mesh_smoothing mesh_smoothing.cpp)
target_link_libraries(mesh_smoothing mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

cm_export_file("../../meshes/chees.cub" "/" export_files_mesh_smoothing)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


#testing
add_test(SimpleMeshSmoothingTest ${MPI_BIN_PATH}/mpirun -np 2 ${CMAKE_CURRENT_BINARY_DIR}/mesh_smoothing 
  -my_file chees.cub  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-12
  -snes_monitor -snes_type newtonls -snes_linesearch_type bt -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 -snes_converged_reason
  -my_penalty 1e12 -my_gamma 0 -my_alpha2 1 -my_meshrefine 1
  -log_summary)


   

