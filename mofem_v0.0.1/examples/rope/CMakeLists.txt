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

add_executable(wire_interface interface_fibre.cpp)
target_link_libraries(wire_interface mofem_student_body mofem_Field ${PROJECT_LIBS})
add_executable(potential_flow_rope potential_flow_rope.cpp)
target_link_libraries(potential_flow_rope mofem_student_body mofem_Field ${PROJECT_LIBS})

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/../../meshes/new_clamped_wire_coarse.cub
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#testing
add_test(ComputeFibreDirection ${MPI_RUN} -np 2 ./potential_flow_rope -my_file new_clamped_wire_coarse.cub -ksp_type fgmres -ksp_gmres_restart 1000 -pc_type asm -sub_pc_type lu -ksp_monitor -my_order 1 -log_summary)
add_test(wireTest ${MPI_RUN} -np 2 ./wire_interface -my_file solution.h5m -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_monitor -log_summary)

