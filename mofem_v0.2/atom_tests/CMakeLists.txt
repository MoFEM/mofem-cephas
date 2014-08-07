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

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/third_party")
include_directories("${PROJECT_SOURCE_DIR}/third_party/cblas")

include_directories("${PROJECT_SOURCE_DIR}/src/approximation")
include_directories("${PROJECT_SOURCE_DIR}/src/approximation/c")
link_directories("${PROJECT_BINARY_DIR}/src/approximation")

include_directories("${PROJECT_SOURCE_DIR}/src/multi_indices")
link_directories("${PROJECT_BINARY_DIR}/src/multi_indices")

include_directories("${PROJECT_SOURCE_DIR}/src/interfaces")
link_directories("${PROJECT_BINARY_DIR}/src/interfaces")

include_directories("${PROJECT_SOURCE_DIR}/src/petsc")
link_directories("${PROJECT_BINARY_DIR}/src/petsc")

include_directories("${PROJECT_SOURCE_DIR}/src/finite_elements")
link_directories("${PROJECT_BINARY_DIR}/src/finite_elements")

include_directories("${PROJECT_SOURCE_DIR}/finite_element_library")
link_directories("${PROJECT_BINARY_DIR}/finite_element_library")

# cubit_bc_test 
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
cm_export_file("../meshes/mat_elastic_transiso.cub" export_atom_tests)
cm_export_file("../meshes/mat_interf.cub" export_atom_tests)
cm_export_file("../meshes/inlet_outlet.cub" export_atom_tests)

add_executable(cubit_bc_test cubit_bc_test.cpp)
target_link_libraries(cubit_bc_test
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

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

add_test(cubit_bc_atom_test_mat_elastic_transiso ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_elastic_transiso.cub)
add_test(cubit_bc_atom_test_mat_elastic_transiso_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_elastic_transiso.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_elastic_transiso.cub.txt)

add_test(cubit_bc_atom_test_mat_interf ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_interf.cub)
add_test(cubit_bc_atom_test_mat_interf_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_interf.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_interf.cub.txt)

add_test(cubit_bc_atom_test_inlet_outlet ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file inlet_outlet.cub)
add_test(cubit_bc_atom_test_inlet_outlet_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/inlet_outlet.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/inlet_outlet.cub.txt)

#Loop over Cubit Meshset Test
cm_export_file("../meshes/loop_atom.cub" export_atom_tests)

add_executable(loops_atom_test loops_atom_test.cpp)
target_link_libraries(loops_atom_test 
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(cubit_meshset_loop_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/loops_atom_test -my_file loop_atom.cub -log_summary)
add_test(cubit_meshset_loop_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/loop_atom.cub.txt
  ${CMAKE_BINARY_DIR}/atom_tests/loop_atom.cub.txt)

#field axpy
cm_export_file("../meshes/elasticity_atom_test_01X_pressure.cub" export_atom_tests)

add_executable(field_axpy_atom_test field_axpy_atom_test.cpp)
target_link_libraries(field_axpy_atom_test 
  mofem_cblas
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  ${PROJECT_LIBS})

add_test(field_axpy ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/field_axpy_atom_test -my_file elasticity_atom_test_01X_pressure.cub -log_summary)
add_test(field_axpy_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/field_axpy_test.txt 
  ${CMAKE_BINARY_DIR}/atom_tests/field_axpy_test.txt)

#10 node TET projection
cm_export_file("../meshes/sphere.cub" export_atom_tests)

add_executable(projection_from_10node_tet_atom projection_from_10node_tet.cpp)
target_link_libraries(projection_from_10node_tet_atom 
  mofem_cblas
  mofem_third_party
  mofem_finite_element_library
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  ${PROJECT_LIBS})

add_test(projection_from_10node_tet_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/projection_from_10node_tet_atom -my_file sphere.cub -log_summary)
add_test(projection_from_10node_tet_atom_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/10node_sphere.txt
  ${CMAKE_BINARY_DIR}/atom_tests/10node_sphere.txt)

#mesh refine
cm_export_file("../meshes/templar.cub" export_atom_tests)

add_executable(mesh_refine_atom mesh_refine_atom.cpp)
target_link_libraries(mesh_refine_atom
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(mesh_refine_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_refine_atom -my_file templar.cub -log_summary)
add_test(mesh_refine_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_refine.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_refine.txt)

#mesh_insert_interface_atom
cm_export_file("../meshes/Tinterface.cub" export_atom_tests)
cm_export_file("../meshes/4seasons.cub" export_atom_tests)

add_executable(mesh_insert_interface_atom mesh_insert_interface_atom.cpp)
target_link_libraries(mesh_insert_interface_atom 
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(mesh_insert_interface_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file templar.cub -log_summary)
add_test(mesh_insert_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_interface.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt)
add_test(mesh_insert_T_interface_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file Tinterface.cub -log_summary)
add_test(mesh_insert_T_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_T_interface.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt)
add_test(mesh_insert_T_4seasons_interface_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file 4seasons.cub -log_summary)
add_test(mesh_insert_T_4seasons_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_T_4seasons_interface.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt)

#foces and sources functions Test
add_executable(forces_and_sources_getting_orders_indices_atom_test forces_and_sources_getting_orders_indices_atom_test.cpp)
target_link_libraries(forces_and_sources_getting_orders_indices_atom_test
  mofem_third_party
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(forces_and_sources_getting_orders_indices_test 
  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_orders_indices_atom_test -my_file elasticity_atom_test_01X_pressure.cub)
add_test(forces_and_sources_getting_orders_indices_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_orders_indices_atom_test.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_orders_indices_atom_test.txt)

add_executable(forces_and_sources_getting_mult_H1_H1_atom_test forces_and_sources_getting_mult_H1_H1_atom_test.cpp)
target_link_libraries(forces_and_sources_getting_mult_H1_H1_atom_test 
  mofem_third_party
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(forces_and_sources_getting_mult_H1_H1_test 
  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_mult_H1_H1_atom_test -my_file elasticity_atom_test_01X_pressure.cub)
add_test(forces_and_sources_getting_mult_H1_H1_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_mult_H1_H1_atom_test.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_mult_H1_H1_atom_test.txt)

add_executable(forces_and_sources_calculate_jacobian forces_and_sources_calculate_jacobian.cpp)
target_link_libraries(forces_and_sources_calculate_jacobian
  mofem_third_party
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(forces_and_sources_calculate_jacobian_test
  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_calculate_jacobian -my_file elasticity_atom_test_01X_pressure.cub)
add_test(forces_and_sources_calculate_jacobian_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_calculate_jacobian.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_calculate_jacobian.txt)

add_executable(forces_and_sources_testing_volume_element forces_and_sources_testing_volume_element.cpp)
target_link_libraries(forces_and_sources_testing_volume_element
  mofem_third_party
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(forces_and_sources_testing_volume_element_test
  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_volume_element -my_file elasticity_atom_test_01X_pressure.cub)
add_test(forces_and_sources_testing_volume_element_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_volume_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_volume_element.txt)

add_executable(forces_and_sources_testing_field_approximation forces_and_sources_testing_field_approximation.cpp) 
target_link_libraries(forces_and_sources_testing_field_approximation 
  mofem_third_party
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  ${PROJECT_LIBS})

add_test(forces_and_sources_testing_field_approximation_test
  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_field_approximation -my_file elasticity_atom_test_01X_pressure.cub 
  -ksp_monitor -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist )
add_test(forces_and_sources_testing_field_approximation_test_compare ${CMAKE_COMMAND} -E compare_files 
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_field_approximation.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_field_approximation.txt)

#add_executable(forces_and_sources_body_force_atom forces_and_sources_body_force_atom_test.cpp)
#target_link_libraries(forces_and_sources_body_force_atom mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
## body forces test
#cm_export_file("../meshes/elasticity_atom_test_body_force.cub" export_atom_tests)
#add_test(forces_and_sources_body_force_atom_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_body_force_atom -my_file elasticity_atom_test_body_force.cub)
#add_test(forces_and_sources_body_force_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_body_force_atom_test.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_body_force_atom_test.txt)
#
#add_executable(forces_and_sources_getting_higher_order_skin_normals_atom forces_and_sources_getting_higher_order_skin_normals_atom_tets.cpp)
#target_link_libraries(forces_and_sources_getting_higher_order_skin_normals_atom mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_getting_higher_order_skin_normals_atom_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_higher_order_skin_normals_atom -my_file sphere.cub)
#add_test(forces_and_sources_getting_higher_order_skin_normals_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_higher_order_skin_normals_atom.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_higher_order_skin_normals_atom.txt)
#
#add_executable(forces_and_sources_testing_triangle_element forces_and_sources_testing_triangle_element.cpp)
#target_link_libraries(forces_and_sources_testing_triangle_element mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_testing_triangle_element_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_triangle_element -my_file elasticity_atom_test_01X_pressure.cub)
#add_test(forces_and_sources_testing_triangle_element_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_triangle_element.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_triangle_element.txt)
#
#add_executable(forces_and_sources_testing_edge_element forces_and_sources_testing_edge_element.cpp)
#target_link_libraries(forces_and_sources_testing_edge_element mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_testing_edge_element_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_edge_element -my_file elasticity_atom_test_01X_pressure.cub)
#add_test(forces_and_sources_testing_edge_element_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_edge_element.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_edge_element.txt)
#
#add_executable(forces_and_sources_testing_vertex_element forces_and_sources_testing_vertex_element.cpp)
#target_link_libraries(forces_and_sources_testing_vertex_element mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_testing_vertex_element_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_vertex_element -my_file elasticity_atom_test_01X_pressure.cub)
#add_test(forces_and_sources_testing_vertex_element_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_vertex_element.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_vertex_element.txt)
#
#add_executable(forces_and_sources_neumann_forces_atom forces_and_sources_neumann_forces_atom_test.cpp)
#target_link_libraries(forces_and_sources_neumann_forces_atom mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_neumann_elasticity_atom_test_01X_pressure
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_neumann_forces_atom -my_file elasticity_atom_test_01X_pressure.cub)
#add_test(forces_and_sources_neumann_elasticity_atom_test_01X_pressure_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_elasticity_atom_test_01X_pressure.cub.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_elasticity_atom_test_01X_pressure.cub.txt)
#add_test(forces_and_sources_neumann_elasticity_atom_test_01X
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_neumann_forces_atom -my_file elasticity_atom_test_01X.cub)
#add_test(forces_and_sources_neumann_elasticity_atom_test_01X_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_elasticity_atom_test_01X.cub.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_elasticity_atom_test_01X.cub.txt)
#add_test(forces_and_sources_neumann_sphre_pressure_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_neumann_forces_atom -my_file sphere.cub)
#add_test(forces_and_sources_neumann_sphre_pressure_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_sphere.cub.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_sphere.cub.txt)
#
#add_executable(forces_and_sources_neumann_forces_complex_for_lazy_atom forces_and_sources_neumann_forces_complex_for_lazy_atom_test.cpp)
#target_link_libraries(forces_and_sources_neumann_forces_complex_for_lazy_atom mofem_complex_for_lazy mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#
#add_test(forces_and_sources_neumann_forces_complex_for_lazy_atom_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_neumann_forces_complex_for_lazy_atom -my_file sphere.cub)
#add_test(forces_and_sources_neumann_forces_complex_for_lazy_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_neumann_forces_complex_for_lazy_atom.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_neumann_forces_complex_for_lazy_atom.txt)
#
## add example test for forces_and_sources_fluid_pressure_element_atom
#cm_export_file("../meshes/fluid_presssure_cub_test.cub" export_atom_tests)
#add_executable(forces_and_sources_fluid_pressure_element_atom forces_and_sources_fluid_pressure_element.cpp)
#target_link_libraries(forces_and_sources_fluid_pressure_element_atom mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_fluid_pressure_element_atom_test
#  ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_fluid_pressure_element_atom -my_file fluid_presssure_cub_test.cub)
#add_test(forces_and_sources_fluid_pressure_element_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_fluid_pressure_element.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_fluid_pressure_element.txt)
#
## add example test for temperature element
#cm_export_file("../meshes/simple_thermal_problem.cub" export_atom_tests)
#add_executable(forces_and_sources_thermal_elem_atom forces_and_sources_thermal_elem.cpp)
#target_link_libraries(forces_and_sources_thermal_elem_atom mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_thermal_elem_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_thermal_elem_atom 
#  -my_file simple_thermal_problem.cub -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_monitor)
#add_test(forces_and_sources_thermal_elem_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_thermal_elem.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_thermal_elem.txt)
#
## add test test for temerature unsteady element
#add_executable(forces_and_sources_thermal_elem_unsteady_atom forces_and_sources_thermal_elem_unsteady.cpp)
#target_link_libraries(forces_and_sources_thermal_elem_unsteady_atom mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_thermal_elem_unsteady_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_thermal_elem_unsteady_atom 
#  -my_file simple_thermal_problem.cub 
#  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist 
#  -snes_type newtonls -snes_linesearch_type basic 
#  -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8 
#  -ts_monitor -ts_type beuler -ts_dt 1 -ts_final_time 2  -ksp_monitor -snes_max_it 3 -snes_monitor)
#add_test(forces_and_sources_thermal_elem_unsteady_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_thermal_elem_unsteady.txt
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_thermal_elem_unsteady.txt)
#
## add theraml stress eleme
#cm_export_file("../meshes/elasticity_atom_test_thermal_expansion.cub" export_atom_tests)
#add_executable(forces_and_sources_thermal_stress_elem forces_and_sources_thermal_stress_elem.cpp)
#target_link_libraries(forces_and_sources_thermal_stress_elem mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_thermal_stress_elem_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_thermal_stress_elem
#  -my_file elasticity_atom_test_thermal_expansion.cub)
#add_test(forces_and_sources_thermal_stress_elem_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_thermal_stress_elem.txt
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_thermal_stress_elem.txt)
#
## add test for recording series
#add_executable(record_series_atom record_series_atom.cpp)
#target_link_libraries(record_series_atom mofem_Field ${PROJECT_LIBS})
#add_test(record_series_atom_test ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/record_series_atom
#  -my_file elasticity_atom_test_01X_pressure.cub)
#add_test(record_series_atom_test_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/record_series_atom.txt
#  ${CMAKE_BINARY_DIR}/atom_tests/record_series_atom.txt)
#
## add test for hdiv space
#add_executable(forces_and_sources_hdiv_approximation_functions forces_and_sources_hdiv_approximation_functions.cpp)
#target_link_libraries(forces_and_sources_hdiv_approximation_functions mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_hdiv_approximation_functions_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_approximation_functions)
#add_test(forces_and_sources_hdiv_approximation_functions_atom_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_hdiv_approximation_functions.txt
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_hdiv_approximation_functions.txt)
#
## add test for hdiv space - check approximation continuity
#add_executable(forces_and_sources_hdiv_continuity_check forces_and_sources_hdiv_continuity_check.cpp)
#target_link_libraries(forces_and_sources_hdiv_continuity_check mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_hdiv_continuity_check_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_continuity_check -my_file sphere.cub)
#add_test(forces_and_sources_hdiv_continuity_check_atom_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_hdiv_continuity_check.txt
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_hdiv_continuity_check.txt)
#
## add test for hdvi divergence operator
#add_executable(forces_and_sources_hdiv_divergence_operator forces_and_sources_hdiv_divergence_operator.cpp)
#target_link_libraries(forces_and_sources_hdiv_divergence_operator mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_hdiv_divergence_operator_atom ${MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_divergence_operator)
#
## add test for ultra weak finite element
#add_executable(forces_and_sources_ultra_weak_transport forces_and_sources_ultra_weak_transport.cpp)
#target_link_libraries(forces_and_sources_ultra_weak_transport mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})
#add_test(forces_and_sources_ultra_weak_transport_atom ${MPI_RUN} -np 2
#  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_ultra_weak_transport -my_file simple_thermal_problem.cub 
#  -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_package superlu_dist -ksp_monitor 
#  -my_max_post_proc_ref_level 1 -my_order 0
#  -log_summary)
#add_test(forces_and_sources_ultra_weak_transport_atom_compare ${CMAKE_COMMAND} -E compare_files 
#  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_ultra_weak_transport.txt
#  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_ultra_weak_transport.txt)
#
