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

add_executable(elasticity_trans_isotropic elasticity.cpp)
target_link_libraries(elasticity_trans_isotropic mofem_Field ${PROJECT_LIBS})

add_executable(elasticity_tube elasticityTube.cpp)
target_link_libraries(elasticity_tube mofem_Field ${PROJECT_LIBS})

add_executable(elasticity_candy1 elasticityCandy1.cpp)
target_link_libraries(elasticity_candy1 mofem_Field ${PROJECT_LIBS})

#cm_export_file("../../../meshes/medium3helix.cub" "/" export_files_elasticity)
#cm_export_file("../../meshes/LShapeCoarse.cub" "/" export_files_elasticity)
#cm_export_file("../../meshes/LShapeDense.cub" "/" export_files_elasticity)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../meshes/wiremesh29_08.cub
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#testing
add_test(SimpleElasticityTest ${MPI_BIN_PATH}/mpirun -np 4 ${CMAKE_CURRENT_BINARY_DIR}/elasticity 
  -my_file LShapeCoarse.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary)




