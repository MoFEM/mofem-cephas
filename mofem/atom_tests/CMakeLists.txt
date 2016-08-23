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

include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/cblas)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)

include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
link_directories(${PROJECT_BINARY_DIR}/src/approximation)

include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
link_directories(${PROJECT_BINARY_DIR}/src/multi_indices)

include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
link_directories(${PROJECT_BINARY_DIR}/src/interfaces)

include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
link_directories(${PROJECT_BINARY_DIR}/src/petsc)

include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
link_directories(${PROJECT_BINARY_DIR}/src/finite_elements)

# cubit_bc_test
cm_export_file("../meshes/disp01.h5m" export_atom_tests)
cm_export_file("../meshes/force01.h5m" export_atom_tests)
cm_export_file("../meshes/velocity01.h5m" export_atom_tests)
cm_export_file("../meshes/accel01.h5m" export_atom_tests)
cm_export_file("../meshes/temper01.h5m" export_atom_tests)
cm_export_file("../meshes/pressure01.h5m" export_atom_tests)
cm_export_file("../meshes/heatflux01.h5m" export_atom_tests)
cm_export_file("../meshes/comb01.h5m" export_atom_tests)
cm_export_file("../meshes/bcoverlap01.h5m" export_atom_tests)
cm_export_file("../meshes/interf01.h5m" export_atom_tests)
cm_export_file("../meshes/mat_elastic.h5m" export_atom_tests)
cm_export_file("../meshes/mat_elastic_transiso.h5m" export_atom_tests)
cm_export_file("../meshes/mat_interf.h5m" export_atom_tests)
cm_export_file("../meshes/inlet_outlet.h5m" export_atom_tests)

set(PERCOMPILED_MOFEM_HEADER "${PROJECT_BINARY_DIR}/include/precompiled/MoFEM.hpp")

# add_executable(cubit_bc_test cubit_bc_test.cpp)
my_add_executable(cubit_bc_test cubit_bc_test.cpp)
target_link_libraries(cubit_bc_test
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(cubit_bc_atom_test_disp01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file disp01.h5m)
add_test(cubit_bc_atom_test_disp01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/disp01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/disp01.h5m.txt
)

add_test(cubit_bc_atom_test_force01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file force01.h5m)
add_test(cubit_bc_atom_test_force01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/force01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/force01.h5m.txt
)

add_test(cubit_bc_atom_test_velocity01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file velocity01.h5m)
add_test(cubit_bc_atom_test_velocity01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/velocity01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/velocity01.h5m.txt
)

add_test(cubit_bc_atom_test_accel01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file accel01.h5m)
add_test(cubit_bc_atom_test_accel01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/accel01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/accel01.h5m.txt
)

add_test(cubit_bc_atom_test_temper01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file temper01.h5m)
add_test(cubit_bc_atom_test_temper01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/temper01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/temper01.h5m.txt
)

add_test(cubit_bc_atom_test_pressure01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file pressure01.h5m)
add_test(cubit_bc_atom_test_pressure01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/pressure01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/pressure01.h5m.txt
)

add_test(cubit_bc_atom_test_heatflux01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file heatflux01.h5m)
add_test(cubit_bc_atom_test_heatflux01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/heatflux01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/heatflux01.h5m.txt
)

add_test(cubit_bc_atom_test_comb01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file comb01.h5m)
add_test(cubit_bc_atom_test_comb01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/comb01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/comb01.h5m.txt
)

add_test(cubit_bc_atom_test_bcoverlap01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file bcoverlap01.h5m)
add_test(cubit_bc_atom_test_bcoverlap01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/bcoverlap01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/bcoverlap01.h5m.txt
)

add_test(cubit_bc_atom_test_interf01 ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file interf01.h5m)
add_test(cubit_bc_atom_test_interf01_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/interf01.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/interf01.h5m.txt
)

add_test(cubit_bc_atom_test_mat_elastic ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_elastic.h5m)
add_test(cubit_bc_atom_test_mat_elastic_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_elastic.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_elastic.h5m.txt
)

add_test(cubit_bc_atom_test_mat_elastic_transiso ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_elastic_transiso.h5m)
add_test(cubit_bc_atom_test_mat_elastic_transiso_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_elastic_transiso.cub.txt
  ${CMAKE_BINARY_DIR}/atom_tests/mat_elastic_transiso.h5m.txt
)

add_test(cubit_bc_atom_test_mat_interf ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file mat_interf.h5m)
add_test(cubit_bc_atom_test_mat_interf_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mat_interf.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/mat_interf.h5m.txt
)

add_test(cubit_bc_atom_test_inlet_outlet ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test -my_file inlet_outlet.h5m)
add_test(cubit_bc_atom_test_inlet_outlet_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/inlet_outlet.cub.txt ${CMAKE_BINARY_DIR}/atom_tests/inlet_outlet.h5m.txt
)

#Loop over Cubit Meshset Test
cm_export_file("../meshes/loop_atom.h5m" export_atom_tests)

my_add_executable(loops_atom_test loops_test.cpp)
target_link_libraries(loops_atom_test
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(cubit_meshset_loop_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/loops_atom_test -my_file loop_atom.h5m)
add_test(cubit_meshset_loop_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/loop_atom.cub.txt
  ${CMAKE_BINARY_DIR}/atom_tests/loop_atom.h5m.txt
)

#field axpy
cm_export_file("../meshes/elasticity_atom_test_01X_pressure.h5m" export_atom_tests)

my_add_executable(field_axpy_atom_test field_axpy_atom_test.cpp)
target_link_libraries(field_axpy_atom_test
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(field_axpy ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/field_axpy_atom_test -my_file elasticity_atom_test_01X_pressure.h5m)
add_test(field_axpy_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/field_axpy_test.txt
  ${CMAKE_BINARY_DIR}/atom_tests/field_axpy_test.txt)

#test coord sys
my_add_executable(test_setting_coordinate_sys test_setting_coordinate_sys.cpp)
target_link_libraries(test_setting_coordinate_sys
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(test_setting_coordinate_sys ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/test_setting_coordinate_sys -my_file elasticity_atom_test_01X_pressure.h5m)

#mesh refine
cm_export_file("../meshes/templar.h5m" export_atom_tests)

my_add_executable(mesh_refine_atom mesh_refine_atom.cpp)
target_link_libraries(mesh_refine_atom
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(mesh_refine_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/mesh_refine_atom -my_file templar.h5m)
add_test(mesh_refine_atom_test_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_refine.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_refine.txt)

#mesh_insert_interface_atom
cm_export_file("../meshes/Tinterface.h5m" export_atom_tests)
cm_export_file("../meshes/4seasons.h5m" export_atom_tests)

my_add_executable(mesh_insert_interface_atom mesh_insert_interface_atom.cpp)
target_link_libraries(mesh_insert_interface_atom
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(mesh_insert_interface_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file templar.h5m
)
add_test(mesh_insert_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_interface.txt ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt
)
add_test(mesh_insert_T_interface_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file Tinterface.h5m
)
add_test(mesh_insert_T_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_T_interface.txt
  ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt
)
add_test(mesh_insert_T_4seasons_interface_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom -my_file 4seasons.h5m
)
add_test(mesh_insert_T_4seasons_interface_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/mesh_insert_T_4seasons_interface.txt
  ${CMAKE_BINARY_DIR}/atom_tests/mesh_insert_interface.txt
)

#foces and sources functions Test
my_add_executable(forces_and_sources_getting_orders_indices_atom_test forces_and_sources_getting_orders_indices_atom_test.cpp)
target_link_libraries(forces_and_sources_getting_orders_indices_atom_test
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_getting_orders_indices_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_orders_indices_atom_test -my_file elasticity_atom_test_01X_pressure.h5m)
add_test(forces_and_sources_getting_orders_indices_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_orders_indices_atom_test.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_orders_indices_atom_test.txt
)

my_add_executable(forces_and_sources_getting_mult_H1_H1_atom_test forces_and_sources_getting_mult_H1_H1_atom_test.cpp)
target_link_libraries(forces_and_sources_getting_mult_H1_H1_atom_test
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_getting_mult_H1_H1_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_mult_H1_H1_atom_test -my_file elasticity_atom_test_01X_pressure.h5m)
add_test(forces_and_sources_getting_mult_H1_H1_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_mult_H1_H1_atom_test.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_mult_H1_H1_atom_test.txt
)

my_add_executable(forces_and_sources_calculate_jacobian forces_and_sources_calculate_jacobian.cpp)
target_link_libraries(forces_and_sources_calculate_jacobian
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_calculate_jacobian_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_calculate_jacobian -my_file elasticity_atom_test_01X_pressure.h5m)
add_test(forces_and_sources_calculate_jacobian_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_calculate_jacobian.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_calculate_jacobian.txt
)

my_add_executable(forces_and_sources_testing_volume_element forces_and_sources_testing_volume_element.cpp)
target_link_libraries(forces_and_sources_testing_volume_element
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(
  forces_and_sources_testing_volume_element_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_volume_element
  -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(forces_and_sources_testing_volume_element_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_volume_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_volume_element.txt
)

my_add_executable(forces_and_sources_testing_ftensor forces_and_sources_testing_ftensor.cpp)
target_link_libraries(forces_and_sources_testing_ftensor
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(
  forces_and_sources_testing_ftensor
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_ftensor
  -my_file elasticity_atom_test_01X_pressure.h5m
)

my_add_executable(forces_and_sources_getting_higher_order_skin_normals_atom forces_and_sources_getting_higher_order_skin_normals_atom_tets.cpp)
target_link_libraries(forces_and_sources_getting_higher_order_skin_normals_atom
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

cm_export_file("../meshes/sphere.h5m" export_atom_tests)

add_test(forces_and_sources_getting_higher_order_skin_normals_atom_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_higher_order_skin_normals_atom -my_file sphere.h5m
)
add_test(forces_and_sources_getting_higher_order_skin_normals_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_getting_higher_order_skin_normals_atom.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_getting_higher_order_skin_normals_atom.txt
)

#triangle element
my_add_executable(forces_and_sources_testing_triangle_element forces_and_sources_testing_triangle_element.cpp)
target_link_libraries(forces_and_sources_testing_triangle_element
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_testing_triangle_element_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_triangle_element
  -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(forces_and_sources_testing_triangle_element_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_triangle_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_triangle_element.txt
)

#edge element
my_add_executable(forces_and_sources_testing_edge_element forces_and_sources_testing_edge_element.cpp)
target_link_libraries(forces_and_sources_testing_edge_element
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_testing_edge_element_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_edge_element
  -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(forces_and_sources_testing_edge_element_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_edge_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_edge_element.txt
)

#vertex element
my_add_executable(forces_and_sources_testing_vertex_element forces_and_sources_testing_vertex_element.cpp)
target_link_libraries(forces_and_sources_testing_vertex_element
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_testing_vertex_element_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_vertex_element
  -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(forces_and_sources_testing_vertex_element_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_vertex_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_vertex_element.txt
)

#flat prism element
my_add_executable(forces_and_sources_testing_flat_prism_element forces_and_sources_testing_flat_prism_element.cpp)
target_link_libraries(forces_and_sources_testing_flat_prism_element
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS})

add_test(forces_and_sources_testing_flat_prism_element_test
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_flat_prism_element -my_file templar.h5m
)
add_test(forces_and_sources_testing_flat_prism_element_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_testing_flat_prism_element.txt
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_testing_flat_prism_element.txt
)

# add test for recording series
my_add_executable(record_series_atom record_series_atom.cpp)
target_link_libraries(record_series_atom
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(record_series_atom_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/record_series_atom
  -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(record_series_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/record_series_atom.txt
  ${CMAKE_BINARY_DIR}/atom_tests/record_series_atom.txt
)

# add test for hcurl space
my_add_executable(forces_and_sources_hcurl_approximation_functions forces_and_sources_hcurl_approximation_functions.cpp)
target_link_libraries(forces_and_sources_hcurl_approximation_functions
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(
  forces_and_sources_hcurl_approximation_functions_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hcurl_approximation_functions
)
add_test(forces_and_sources_hcurl_approximation_functions_atom_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_hcurl_approximation_functions.txt
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_hcurl_approximation_functions.txt
)

my_add_executable(forces_and_sources_hcurl_continuity_check forces_and_sources_hcurl_continuity_check.cpp)
target_link_libraries(forces_and_sources_hcurl_continuity_check
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(forces_and_sources_hcurl_continuity_check_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_continuity_check -my_file sphere.h5m
)

# add test for hdiv space
my_add_executable(forces_and_sources_hdiv_approximation_functions forces_and_sources_hdiv_approximation_functions.cpp)
target_link_libraries(forces_and_sources_hdiv_approximation_functions
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(
  forces_and_sources_hdiv_approximation_functions_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_approximation_functions
)
add_test(forces_and_sources_hdiv_approximation_functions_atom_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_BINARY_DIR}/atom_tests/forces_and_sources_hdiv_approximation_functions.txt
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/forces_and_sources_hdiv_approximation_functions.txt
)

# add test for hdiv space - check approximation continuity
my_add_executable(forces_and_sources_hdiv_continuity_check forces_and_sources_hdiv_continuity_check.cpp)
target_link_libraries(forces_and_sources_hdiv_continuity_check
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(forces_and_sources_hdiv_continuity_check_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_continuity_check -my_file sphere.h5m
)

# add test for hdvi divergence operator
my_add_executable(forces_and_sources_hdiv_divergence_operator forces_and_sources_hdiv_divergence_operator.cpp)
target_link_libraries(forces_and_sources_hdiv_divergence_operator
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(forces_and_sources_hdiv_divergence_operator_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_divergence_operator
)

#dm
my_add_executable(dm_mofem dm_mofem.cpp)
target_link_libraries(dm_mofem
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(dm_mofem_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/dm_mofem -my_file disp01.h5m)
add_test(dm_mofem_atom_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_BINARY_DIR}/atom_tests/dm_mofem.txt
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/dm_mofem.txt
)

cm_export_file("disp01_2parts.h5m" export_atom_tests)
my_add_executable(dm_build_partitioned_mesh dm_build_partitioned_mesh.cpp)
target_link_libraries(dm_build_partitioned_mesh
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(dm_build_partitioned_mesh_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/dm_build_partitioned_mesh -my_file disp01_2parts.h5m -dm_is_partitioned
)
add_test(dm_build_partitioned_mesh_atom_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_BINARY_DIR}/atom_tests/dm_build_partitioned_mesh.txt
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/dm_build_partitioned_mesh.txt
)

#single process
my_add_executable(serial_matrix serial_matrix.cpp)
target_link_libraries(serial_matrix
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(serial_matrix ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/serial_matrix -my_file disp01.h5m)
add_test(serial_matrix_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_BINARY_DIR}/atom_tests/serial_matrix.txt
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/serial_matrix.txt
)

#test for build problems
my_add_executable(build_problems build_problems.cpp)
target_link_libraries(build_problems
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS})
add_test(build_problems_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/build_problems -my_file disp01.h5m)

#test inactivating and clearing inactive dofs
my_add_executable(activate_deactivate_dofs activate_deactivate_dofs.cpp)
target_link_libraries(activate_deactivate_dofs
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS})
add_test(activate_deactivate_dofs ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/activate_deactivate_dofs -my_file disp01.h5m)


# Projection form and to 10 node tets
my_add_executable(projection_from_10node_tet_atom projection_from_10node_tet.cpp)
target_link_libraries(projection_from_10node_tet_atom
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

add_test(
  projection_from_10node_tet_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/projection_from_10node_tet_atom -my_file sphere.h5m
)
add_test(projection_from_10node_tet_atom_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_SOURCE_DIR}/atom_tests/blessed_files/10node_sphere.txt
  ${CMAKE_BINARY_DIR}/atom_tests/10node_sphere.txt
)

if(TETGEN_LIBRARY)

  cm_export_file("../meshes/Tinterface.h5m" export_atom_tests)

  my_add_executable(tetgen_interface tetgen_interface.cpp)
  target_link_libraries(tetgen_interface
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_cblas
    mofem_third_party
    ${PROJECT_LIBS}
  )
  if(WITH_TETGEN)
    add_dependencies(tetgen_interface tetgen)
  endif(WITH_TETGEN)

  add_test(tetgen_interface_atom ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/tetgen_interface -my_file Tinterface.h5m
  )
endif(TETGEN_LIBRARY)

if(NETGEN_LIBRARY)

  my_add_executable(netgen_interface netgen_interface.cpp)
  target_link_libraries(netgen_interface
    mofem_interfaces
    mofem_multi_indices
    mofem_approx
    mofem_cblas
    mofem_third_party
    ${PROJECT_LIBS}
  )

endif(NETGEN_LIBRARY)

my_add_executable(node_merge node_merge.cpp)
target_link_libraries(node_merge
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)

cm_export_file("../meshes/node_merge_test.h5m" export_atom_tests)
add_test(node_merge ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/node_merge -my_file node_merge_test.h5m
)

my_add_executable(prisms_elements_from_surface prisms_elements_from_surface.cpp)
target_link_libraries(prisms_elements_from_surface
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
cm_export_file("../meshes/rectangle_surface_mesh.cub" export_atom_tests)
add_test(prisms_elements_from_surface ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/prisms_elements_from_surface -my_file rectangle_surface_mesh.cub
)

my_add_executable(partition_mesh partition_mesh.cpp)
target_link_libraries(partition_mesh
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
add_test(partition_mesh ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/partition_mesh -my_file disp01.h5m
)

my_add_executable(testing_base_functions testing_base_functions.cpp)
target_link_libraries(testing_base_functions
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
list(INSERT testing_base_list
  0
  legendre
  lobatto
  h1tet
  hcurltet
  hdivtet
  l2tet
  h1tri
  hdivtri
  hcurltri
  #l2tri
  h1edge
  hcurledge
)
foreach(test_base ${testing_base_list})
  add_test(testing_base_functions_${test_base}
    ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/testing_base_functions -base ${test_base}
  )
endforeach(test_base IN testing_base_list)

if(ADOL-C_LIBRARY)
  my_add_executable(adol-c_atom adol-c_atom.cpp)
  target_link_libraries(adol-c_atom ${PROJECT_LIBS})
  if(WITH_ADOL-C)
    add_dependencies(adol-c_atom adolc)
  endif(WITH_ADOL-C)
  add_test(
    adol-c_atom_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/adol-c_atom
  )
endif(ADOL-C_LIBRARY)
