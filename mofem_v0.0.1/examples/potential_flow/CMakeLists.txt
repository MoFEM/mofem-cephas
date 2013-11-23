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

add_executable(potential_flow potential_flow.cpp)
target_link_libraries(potential_flow mofem_student_body mofem_Field ${PROJECT_LIBS})

cm_export_file("../../meshes/pipe.cub" export_potential_flow)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


#testing
add_test(SimplePotentialFlowTest ${MPI_RUN} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/potential_flow
  -my_file pipe.cub -ksp_type fgmres -ksp_gmres_restart 1000 -pc_type asm -sub_pc_type lu -ksp_monitor -my_order 5 -log_summary)

install(TARGETS potential_flow DESTINATION bin)

