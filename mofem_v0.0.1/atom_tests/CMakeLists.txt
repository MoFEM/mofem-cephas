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

add_executable(cubit_bc_test cubit_bc_test.cpp)
target_link_libraries(cubit_bc_test mofem_Field ${PROJECT_LIBS})

add_executable(elasticity_atom_test elasticity_atom_test.cpp)
target_link_libraries(elasticity_atom_test mofem_Field ${PROJECT_LIBS})

add_executable(l2_approximation l2_approximation.cpp)
target_link_libraries(l2_approximation mofem_Field ${PROJECT_LIBS})

cm_export_file("../meshes/disp01.cub" "/" export_atom_tests)
cm_export_file("../meshes/force01.cub" "/" export_atom_tests)
cm_export_file("../meshes/velocity01.cub" "/" export_atom_tests)
cm_export_file("../meshes/accel01.cub" "/" export_atom_tests)
cm_export_file("../meshes/temper01.cub" "/" export_atom_tests)
cm_export_file("../meshes/pressure01.cub" "/" export_atom_tests)
cm_export_file("../meshes/heatflux01.cub" "/" export_atom_tests)
cm_export_file("../meshes/comb01.cub" "/" export_atom_tests)
cm_export_file("../meshes/bcoverlap01.cub" "/" export_atom_tests)
cm_export_file("../meshes/interf01.cub" "/" export_atom_tests)
cm_export_file("../meshes/mat_elastic.cub" "/" export_atom_tests)
cm_export_file("../meshes/mat_transiso.cub" "/" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01X.cub" "/" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01Y.cub" "/" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01Z.cub" "/" export_atom_tests)

add_test(cubit_bc_test_disp01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file disp01.cub)
add_test(cubit_bc_test_disp01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/disp01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/disp01.cub.txt)

add_test(cubit_bc_test_force01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file force01.cub)
add_test(cubit_bc_test_force01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/force01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/force01.cub.txt)

add_test(cubit_bc_test_velocity01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file velocity01.cub)
add_test(cubit_bc_test_velocity01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/velocity01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/velocity01.cub.txt)

add_test(cubit_bc_test_accel01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file accel01.cub)
add_test(cubit_bc_test_accel01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/accel01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/accel01.cub.txt)

add_test(cubit_bc_test_temper01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file temper01.cub)
add_test(cubit_bc_test_temper01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/temper01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/temper01.cub.txt)

add_test(cubit_bc_test_pressure01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file pressure01.cub)
add_test(cubit_bc_test_pressure01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/pressure01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/pressure01.cub.txt)

add_test(cubit_bc_test_heatflux01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file heatflux01.cub)
add_test(cubit_bc_test_heatflux01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/heatflux01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/heatflux01.cub.txt)

add_test(cubit_bc_test_comb01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file comb01.cub)
add_test(cubit_bc_test_comb01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/comb01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/comb01.cub.txt)

add_test(cubit_bc_test_bcoverlap01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file bcoverlap01.cub)
add_test(cubit_bc_test_bcoverlap01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/bcoverlap01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/bcoverlap01.cub.txt)

add_test(cubit_bc_test_interf01 ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file interf01.cub)
add_test(cubit_bc_test_interf01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/interf01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/interf01.cub.txt)

add_test(cubit_bc_test_mat_elastic ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_elastic.cub)
add_test(cubit_bc_test_mat_elastic_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_elastic.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_elastic.cub.txt)

add_test(cubit_bc_test_mat_transiso ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_transiso.cub)
add_test(cubit_bc_test_mat_transiso_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_transiso.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_transiso.cub.txt)

#Elasticity atom tests
add_test(elasticity_atom_test_01X ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01X.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01X_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01X.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01X.cub.txt)

add_test(elasticity_atom_test_01Y ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Y.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01Y_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Y.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Y.cub.txt)

add_test(elasticity_atom_test_01Z ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Z.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01Z_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Z.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Z.cub.txt)

add_test(l2_approximation ${MPI_BIN_PATH}/mpirun -np 1 ${CMAKE_CURRENT_BINARY_DIR}/l2_approximation -my_file elasticity_atom_test_01X.cub)
add_test(l2_approximation_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/l2_approximation.txt ${CMAKE_BINARY_DIR}/atom_tests/l2_approximation.txt)







