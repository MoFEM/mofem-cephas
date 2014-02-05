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
target_link_libraries(elasticity_atom_test mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(l2_approximation l2_approximation.cpp)
target_link_libraries(l2_approximation mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(hdiv_approximation hdiv_approximation.cpp)
target_link_libraries(hdiv_approximation mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(trans_isotropic_test trans_isotropic_test.cpp)
target_link_libraries(trans_isotropic_test mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(mesh_smoothing_atom mesh_smoothing_atom.cpp)
target_link_libraries(mesh_smoothing_atom mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(mesh_insert_interface_atom mesh_insert_interface_atom.cpp)
target_link_libraries(mesh_insert_interface_atom mofem_Field ${PROJECT_LIBS})

add_executable(nonlinear_elasticity_atom nonlinear_elasticity_atom.cpp)
target_link_libraries(nonlinear_elasticity_atom mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(arc_length_nonlinear_elasticity_atom arc_length_nonlinear_elasticity_atom.cpp)
target_link_libraries(arc_length_nonlinear_elasticity_atom mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(field_axpy_atom_test field_axpy_atom_test.cpp)
target_link_libraries(field_axpy_atom_test mofem_Field ${PROJECT_LIBS})

add_executable(nonlinear_elasticity_thremal_expansion_atom nonlinear_elasticity_thremal_expansion_atom.cpp)
target_link_libraries(nonlinear_elasticity_thremal_expansion_atom mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(arc_length_nonlinear_elasticity_thermal_atom arc_length_nonlinear_elasticity_thermal_atom.cpp)
target_link_libraries(arc_length_nonlinear_elasticity_thermal_atom mofem_complex_for_lazy mofem_student_body mofem_Field ${PROJECT_LIBS})

cm_export_file("../meshes/disp01.cub" export_atom_tests)
cm_export_file("../meshes/force01.cub" export_atom_tests)
cm_export_file("../meshes/velocity01.cub" export_atom_tests)
cm_export_file("../meshes/accel01.cub" export_atom_tests)
cm_export_file("../meshes/temper01.cub" export_atom_tests)
cm_export_file("../meshes/pressure01.cub" export_atom_tests)
cm_export_file("../meshes/heatflux01.cub" export_atom_tests)
cm_export_file("../meshes/comb01.cub" export_atom_tests)
cm_export_file("../meshes/bcoverlap01.cub" export_atom_tests)
cm_export_file("../meshes/interf01.cub" export_atom_tests)
cm_export_file("../meshes/mat_elastic.cub" export_atom_tests)
cm_export_file("../meshes/mat_transiso.cub" export_atom_tests)
cm_export_file("../meshes/mat_interf.cub" export_atom_tests)
cm_export_file("../meshes/inlet_outlet.cub" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01X.cub" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01Y.cub" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01Z.cub" export_atom_tests)
cm_export_file("../meshes/trans_isotropic_cube.cub" export_atom_tests)
cm_export_file("../meshes/chees_atom.cub" export_atom_tests)
cm_export_file("../meshes/cub_with_interface.cub" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_01X_pressure.cub" export_atom_tests)
cm_export_file("../meshes/elasticity_atom_test_sheer_force.cub" export_atom_tests)
cm_export_file("../meshes/templar.cub" export_atom_tests)
cm_export_file("../meshes/cube_for_thermal_expansion.cub" export_atom_tests)


add_test(cubit_bc_atom_test_disp01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file disp01.cub)
add_test(cubit_bc_atom_test_disp01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/disp01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/disp01.cub.txt)

add_test(cubit_bc_atom_test_force01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file force01.cub)
add_test(cubit_bc_atom_test_force01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/force01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/force01.cub.txt)

add_test(cubit_bc_atom_test_velocity01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file velocity01.cub)
add_test(cubit_bc_atom_test_velocity01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/velocity01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/velocity01.cub.txt)

add_test(cubit_bc_atom_test_accel01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file accel01.cub)
add_test(cubit_bc_atom_test_accel01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/accel01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/accel01.cub.txt)

add_test(cubit_bc_atom_test_temper01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file temper01.cub)
add_test(cubit_bc_atom_test_temper01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/temper01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/temper01.cub.txt)

add_test(cubit_bc_atom_test_pressure01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file pressure01.cub)
add_test(cubit_bc_atom_test_pressure01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/pressure01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/pressure01.cub.txt)

add_test(cubit_bc_atom_test_heatflux01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file heatflux01.cub)
add_test(cubit_bc_atom_test_heatflux01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/heatflux01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/heatflux01.cub.txt)

add_test(cubit_bc_atom_test_comb01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file comb01.cub)
add_test(cubit_bc_atom_test_comb01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/comb01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/comb01.cub.txt)

add_test(cubit_bc_atom_test_bcoverlap01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file bcoverlap01.cub)
add_test(cubit_bc_atom_test_bcoverlap01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/bcoverlap01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/bcoverlap01.cub.txt)

add_test(cubit_bc_atom_test_interf01 ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file interf01.cub)
add_test(cubit_bc_atom_test_interf01_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/interf01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/interf01.cub.txt)

add_test(cubit_bc_atom_test_mat_elastic ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_elastic.cub)
add_test(cubit_bc_atom_test_mat_elastic_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_elastic.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_elastic.cub.txt)

add_test(cubit_bc_atom_test_mat_transiso ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_transiso.cub)
add_test(cubit_bc_atom_test_mat_transiso_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_transiso.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_transiso.cub.txt)

add_test(cubit_bc_atom_test_mat_interf ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_interf.cub)
add_test(cubit_bc_atom_test_mat_interf_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_interf.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_interf.cub.txt)

add_test(cubit_bc_atom_test_inlet_outlet ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file inlet_outlet.cub)
add_test(cubit_bc_atom_test_inlet_outlet_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/inlet_outlet.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/inlet_outlet.cub.txt)

#Elasticity atom tests
add_test(elasticity_atom_test_01X ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01X.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01X_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01X.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01X.cub.txt)

#Transversely Isotropic atom test
add_test(traniso_mat_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/trans_isotropic_test -my_file trans_isotropic_cube.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_monitor)
add_test(traniso_mat_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/transIso_Mat_Blessed_File.txt ${CMAKE_BINARY_DIR}/atom_tests/transIso_Mat_Blessed_File.txt)

add_test(elasticity_atom_test_01Y ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Y.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01Y_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Y.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Y.cub.txt)

add_test(elasticity_atom_test_01Z ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/elasticity_atom_test -my_file elasticity_atom_test_01Z.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -log_summary -my_order 5)
add_test(elasticity_atom_test_01Z_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/elasticity_atom_test_01Z.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/elasticity_atom_test_01Z.cub.txt)

add_test(l2_atom_approximation ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/l2_approximation -my_file elasticity_atom_test_01X.cub -log_summary)
add_test(l2_atom_approximation_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/l2_approximation.txt ${CMAKE_BINARY_DIR}/atom_tests/l2_approximation.txt)

add_test(mesh_smoothing_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_smoothing_atom -my_file chees_atom.cub  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-12 -ksp_rtol 1e-12 -snes_monitor -snes_type newtonls -snes_linesearch_type bt -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 -snes_converged_reason -my_penalty 1e10 -my_gamma 0.0 -my_alpha2 1 -log_summary)
add_test(mesh_smoothing_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_smoothing.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_smoothing.txt)

add_test(mesh_insert_interface_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file templar.cub -log_summary)
add_test(mesh_insert_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_interface.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt)

#nonlinear elasticty
add_test(nonlinear_elasticity_atom_pressure ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_atom -my_file elasticity_atom_test_01X_pressure.cub  -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 -snes_monitor -snes_type newtonls -snes_linesearch_type l2 -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 -my_order 4 -log_summary)
add_test(nonlinear_elasticity_atom_pressure_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_elasticity_atom_test_01X_pressure.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_elasticity_atom_test_01X_pressure.cub.txt)

add_test(nonlinear_elasticity_atom_sheer_force ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_atom -my_file elasticity_atom_test_sheer_force.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 -snes_monitor -snes_type newtonls -snes_linesearch_type l2 -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 -my_order 4 -log_summary)
add_test(nonlinear_elasticity_atom_sheer_force_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_elasticity_atom_test_sheer_force.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_elasticity_atom_test_sheer_force.cub.txt)

#arc length elasticty
add_test(arc_length_nonlinear_elasticity_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/arc_length_nonlinear_elasticity_atom -my_file elasticity_atom_test_01X_pressure.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_monitor ksp_monitor -ksp_atol 1e-10 -ksp_rtol 1e-10 -snes_monitor -snes_type newtonls -snes_linesearch_type l2 -snes_linesearch_monitor  -snes_max_it 16 -snes_atol 1e-8 -snes_rtol 1e-8 -my_sr 1. -my_ms 3 -my_order 5 -snes_converged_reason -log_summary)
add_test(arc_length_nonlinear_elasticity_atom_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/arc_length_nonlinear_elasticity_atom.txt ${CMAKE_BINARY_DIR}/atom_tests/arc_length_nonlinear_elasticity_atom.txt)

#hdiv approximation atom test
add_test(hdiv_atom_approximation ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/hdiv_approximation -my_file elasticity_atom_test_01X.cub -my_max_pot_proc_ref_level 2 -ksp_type fgmres -ksp_atol 1e-16 -ksp_rtol 0 -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_monitor -log_summary)
add_test(hdiv_atom_approximation_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/hdiv_approximation.txt ${CMAKE_BINARY_DIR}/atom_tests/hdiv_approximation.txt)

#field axpy
add_test(field_axpy ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/field_axpy_atom_test -my_file elasticity_atom_test_01X_pressure.cub -log_summary)
add_test(field_axpy_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/field_axpy_test.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/field_axpy_test.txt)

#thermal test
add_test(nonlinear_elasticity_atom_thermal ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elasticity_thremal_expansion_atom -my_file cube_for_thermal_expansion.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_atol 1e-10 -ksp_rtol 1e-8 -snes_monitor -snes_type newtonls -snes_linesearch_type l2 -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 -my_order 4 -log_summary)
add_test(nonlinear_elasticity_atom_thermal_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/nonlinear_thermal_expansion_cube_for_thermal_expansion.cub.txt
  ${CMAKE_BINARY_DIR}/atom_tests/nonlinear_thermal_expansion_cube_for_thermal_expansion.cub.txt)

#arc length elasticty thermal expansion
add_test(arc_length_nonlinear_elasticity_thermal_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/arc_length_nonlinear_elasticity_thermal_atom -my_file cube_for_thermal_expansion.cub -ksp_type fgmres -pc_type asm -sub_pc_type lu -ksp_monitor ksp_monitor -ksp_atol 1e-10 -ksp_rtol 1e-10 -snes_monitor -snes_type newtonls -snes_linesearch_type basic -snes_linesearch_monitor -snes_max_it 10 -snes_atol 1e-8 -snes_rtol 1e-8 -my_sr 2e-2 -my_ms 3 -ksp_monitor -snes_converged_reason -log_summary)
add_test(arc_length_nonlinear_elasticity_thermal_atom_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/arc_length_nonlinear_elasticity_thermal_atom.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/arc_length_nonlinear_elasticity_thermal_atom.txt)

