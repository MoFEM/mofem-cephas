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
cm_export_file("../meshes/rectangle_tri.h5m" export_atom_tests)
cm_export_file("../meshes/rectangle_quad.h5m" export_atom_tests)
cm_export_file("../meshes/quad_mesh.h5m" export_atom_tests)
cm_export_file("../meshes/2_tets_contact.cub" export_atom_tests)
cm_export_file("../meshes/reference_rectangle.h5m" export_atom_tests)
cm_export_file("../meshes/half_of_sphere.h5m" export_atom_tests)
cm_export_file("../meshes/unit_volume_sphere.cub" export_atom_tests)
cm_export_file("../meshes/hex_cub.h5m" export_atom_tests)

set(PERCOMPILED_MOFEM_HEADER 
  "${PROJECT_BINARY_DIR}/include/precompiled/MoFEM.hpp")

function(build_test target source) 
  my_add_executable(${target} ${source})
  target_link_libraries(
    ${target}
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    mofem_matrix_function
    ${PROJECT_LIBS}
  )
endfunction(build_test)

# add_executable(cubit_bc_test cubit_bc_test.cpp)
build_test(cubit_bc_test cubit_bc_test.cpp)

list(INSERT cubit_bc_test_list
  0
  disp01.h5m
  force01.h5m
  velocity01.h5m
  accel01.h5m
  temper01.h5m
  pressure01.h5m
  heatflux01.h5m
  comb01.h5m
  bcoverlap01.h5m
  interf01.h5m
  mat_elastic.h5m
  mat_elastic_transiso.h5m
  mat_interf.h5m
  inlet_outlet.h5m
)

foreach(test_file ${cubit_bc_test_list})
  add_test(
    cubit_bc_atom_test_${test_file}
    ${CMAKE_COMMAND}
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/cubit_bc_test
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DFILE=${test_file}
    -DLOG1=${test_file}.txt
    -DLOG2=${test_file}.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
endforeach(test_file IN cubit_bc_test_list)

# testing adding cubit meshset
build_test(add_cubit_meshsets add_cubit_meshsets.cpp)

cm_export_file("add_cubit_meshsets.in" export_atom_tests)

add_test(
  add_cubit_meshsets
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/add_cubit_meshsets
  -meshsets_config add_cubit_meshsets.in -log_quiet
)

#Loop over Cubit Meshset Test
cm_export_file("../meshes/loop_atom.h5m" export_atom_tests)

build_test(loops_atom_test loops_test.cpp)

add_test(
    cubit_meshset_loop_test 
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/loops_atom_test
    -DFILE=loop_atom.h5m
    -DLOG1=loop_atom.h5m.txt
    -DLOG2=loop_atom.h5m.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#field axpy
cm_export_file("../meshes/elasticity_atom_test_01X_pressure.h5m" export_atom_tests)

build_test(field_axpy_atom_test field_axpy_atom_test.cpp)

add_test(field_axpy 
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/field_axpy_atom_test 
  -my_file elasticity_atom_test_01X_pressure.h5m
  -log_no_color
)

build_test(field_blas_set_vertex_dofs field_blas_set_vertex_dofs.cpp)
  
add_test(field_blas_set_vertex_dofs
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/field_blas_set_vertex_dofs 
  -file_name elasticity_atom_test_01X_pressure.h5m
  -log_no_color
)

#mesh refine
cm_export_file("../meshes/templar.h5m" export_atom_tests)

build_test(mesh_refine_atom mesh_refine_atom.cpp)

add_test(
    mesh_refine 
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/mesh_refine_atom
    -DFILE=templar.h5m
    -DLOG1=mesh_refine.txt
    -DLOG2=mesh_refine.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#mesh_insert_interface_atom
cm_export_file("../meshes/Tinterface.h5m" export_atom_tests)
cm_export_file("../meshes/4seasons.h5m" export_atom_tests)

build_test(mesh_insert_interface_atom mesh_insert_interface_atom.cpp)

list(INSERT testing_insert_interface_list
  0
  templar.h5m
  Tinterface.h5m
  4seasons.h5m
)
foreach(test_file ${testing_insert_interface_list})
  add_test(
    mesh_insert_interface_${test_file}
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/mesh_insert_interface_atom
    -DFILE=${test_file}
    -DLOG1=${test_file}_log
    -DLOG2=${test_file}_log
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
endforeach(test_file IN testing_insert_interface_list)

# foces and sources functions Test
build_test(forces_and_sources_getting_orders_indices_atom_test forces_and_sources_getting_orders_indices_atom_test.cpp)

add_test(
    forces_and_sources_getting_orders_indices ${test_file}
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_orders_indices_atom_test
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_getting_orders_indices_atom_test.txt
    -DLOG2=forces_and_sources_getting_orders_indices_atom_test.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(forces_and_sources_getting_mult_H1_H1_atom_test 
  forces_and_sources_getting_mult_H1_H1_atom_test.cpp)

add_test(
    forces_and_sources_getting_mult_H1_H1
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_mult_H1_H1_atom_test
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_getting_mult_H1_H1_atom_test.txt
    -DLOG2=forces_and_sources_getting_mult_H1_H1_atom_test.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(forces_and_sources_calculate_jacobian forces_and_sources_calculate_jacobian.cpp)

add_test(
    forces_and_sources_calculate_jacobian
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_calculate_jacobian
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_calculate_jacobian.txt
    -DLOG2=forces_and_sources_calculate_jacobian.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(forces_and_sources_testing_volume_element forces_and_sources_testing_volume_element.cpp)

add_test(
    forces_and_sources_testing_volume_element
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_volume_element
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_testing_volume_element.txt
    -DLOG2=forces_and_sources_testing_volume_element.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(forces_and_sources_testing_ftensor forces_and_sources_testing_ftensor.cpp)
add_test(
  forces_and_sources_testing_ftensor
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_ftensor
  -my_file elasticity_atom_test_01X_pressure.h5m
  -log_no_color
)

build_test(forces_and_sources_testing_users_base forces_and_sources_testing_users_base.cpp)
add_test(
    forces_and_sources_testing_users_base
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_users_base
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_testing_users_base.txt
    -DLOG2=forces_and_sources_testing_users_base.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(forces_and_sources_getting_higher_order_skin_normals_atom forces_and_sources_getting_higher_order_skin_normals_atom_tets.cpp)

cm_export_file("../meshes/sphere.h5m" export_atom_tests)

add_test(
    forces_and_sources_getting_higher_order_skin_normals
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_getting_higher_order_skin_normals_atom
    -DFILE=sphere.h5m
    -DLOG1=forces_and_sources_getting_higher_order_skin_normals_atom.txt
    -DLOG2=forces_and_sources_getting_higher_order_skin_normals_atom.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#triangle element
build_test(forces_and_sources_testing_triangle_element forces_and_sources_testing_triangle_element.cpp)

add_test(
    forces_and_sources_testing_triangle_element
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_triangle_element
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_testing_triangle_element.txt
    -DLOG2=forces_and_sources_testing_triangle_element.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#edge element
build_test(
  forces_and_sources_testing_edge_element 
  forces_and_sources_testing_edge_element.cpp)

add_test(
    forces_and_sources_testing_edge_element_ainsworth
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_edge_element
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DBASE=ainsworth
    -DLOG1=forces_and_sources_testing_edge_element_AINSWORTH_LEGENDRE_BASE.txt
    -DLOG2=forces_and_sources_testing_edge_element_AINSWORTH_LEGENDRE_BASE.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
add_test(
    forces_and_sources_testing_edge_element_bernstein_bezier
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_edge_element
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DBASE=bernstein_bezier
    -DLOG1=forces_and_sources_testing_edge_element_AINSWORTH_BERNSTEIN_BEZIER_BASE.txt
    -DLOG2=forces_and_sources_testing_edge_element_AINSWORTH_BERNSTEIN_BEZIER_BASE.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#vertex element
build_test(forces_and_sources_testing_vertex_element forces_and_sources_testing_vertex_element.cpp)

add_test(
    forces_and_sources_testing_vertex_element
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_vertex_element
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=forces_and_sources_testing_vertex_element.txt
    -DLOG2=forces_and_sources_testing_vertex_element.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#flat prism element
build_test(forces_and_sources_testing_flat_prism_element forces_and_sources_testing_flat_prism_element.cpp)

add_test(
    forces_and_sources_testing_flat_prism_element
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_flat_prism_element
    -DFILE=templar.h5m
    -DLOG1=forces_and_sources_testing_flat_prism_element.txt
    -DLOG2=forces_and_sources_testing_flat_prism_element.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

#contact prism element
build_test(forces_and_sources_testing_contact_prism_element forces_and_sources_testing_contact_prism_element.cpp)

add_test(
    forces_and_sources_testing_contact_prism_element_h1
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_contact_prism_element
    -DFILE=templar.h5m
    -DSPACE=h1
    -DLOG1=forces_and_sources_testing_contact_prism_element.txt
    -DLOG2=forces_and_sources_testing_contact_prism_element.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

add_test(
    forces_and_sources_testing_contact_prism_element_hdiv
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_testing_contact_prism_element
    -DFILE=templar.h5m
    -DSPACE=hdiv
    -DLOG1=forces_and_sources_testing_contact_prism_element_HDIV.txt
    -DLOG2=forces_and_sources_testing_contact_prism_element_HDIV.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

# Checking direvatives
build_test(
  check_base_functions_direvatives_on_tet
  check_base_functions_direvatives_on_tet.cpp
)

add_test(check_h1_ainsworth_base_functions_direvatives_on_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_on_tet
  -space H1TET -base ainsworth
  -log_no_color
)
add_test(check_hdiv_ainsworth_base_functions_direvatives_on_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_on_tet
  -space HDIVTET -base ainsworth
  -log_no_color
)
add_test(check_hdiv_demkowicz_base_functions_direvatives_on_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_on_tet
  -space HDIVTET -base demkowicz
  -log_no_color
)
add_test(check_hcurl_ainsworth_base_functions_direvatives_on_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_on_tet
  -space HCURLTET -base ainsworth
  -log_no_color
)
add_test(check_hcurl_demkowicz_base_functions_direvatives_on_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_on_tet
  -space HCURLTET -base demkowicz
  -log_no_color
)

build_test(
  check_base_functions_direvatives_tri
  check_base_functions_direvatives_tri.cpp
)

add_test(check_h1_ainsworth_base_functions_direvatives_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_tri 
  -space h1tri -base ainsworth
  -log_no_color
)
add_test(check_hcurl_ainsworth_base_functions_direvatives_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_tri -space hcurltri 
  -base ainsworth
  -log_no_color
)
add_test(check_hcurl_demkowicz_base_functions_direvatives_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/check_base_functions_direvatives_tri -space hcurltri 
  -base demkowicz
  -log_no_color
)

# Continuity check on tet for h1 space
build_test(forces_and_sources_h1_continuity_check 
  forces_and_sources_h1_continuity_check.cpp)

add_test(forces_and_sources_h1_continuity_check_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_h1_continuity_check 
  -my_file sphere.h5m
  -log_no_color
)

# Check approxination on H1 space in 2d and 3d
build_test(
  scalar_check_approximation_2d scalar_check_approximation.cpp)
set_target_properties(
  scalar_check_approximation_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")
build_test(
  scalar_check_approximation_3d scalar_check_approximation.cpp)
set_target_properties(
  scalar_check_approximation_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")    

add_test(h1_check_approx_2d_ainsworth_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_tri.h5m
  -base ainsworth -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_2d_ainsworth_lobatto_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_tri.h5m
  -base ainsworth_labatto -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_2d_bernstein_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_tri.h5m
  -base bernstein -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_2d_ainsworth_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_tri.h5m
  -base ainsworth -space l2 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_2d_bernstein_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_tri.h5m
  -base bernstein -space l2 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_2d_ainsworth_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_quad.h5m
  -base ainsworth -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_2d_demkowicz_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_quad.h5m
  -base demkowicz -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_2d_demkowicz_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_2d
  -file_name rectangle_quad.h5m
  -base demkowicz -space l2
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)

add_test(h1_check_approx_3d_ainsworth_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name sphere.h5m
  -base ainsworth -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_3d_ainsworth_lobatto_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name sphere.h5m
  -base ainsworth_labatto -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_3d_bernstein_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name sphere.h5m
  -base bernstein -space h1 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_3d_ainsworth_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name sphere.h5m
  -base ainsworth -space l2 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_3d_bernstein_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name sphere.h5m
  -base bernstein -space l2 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(l2_check_approx_3d_demkowicz_hex
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name hex_cub.h5m
  -base demkowicz -space l2 
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)
add_test(h1_check_approx_3d_demkowicz_hex
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/scalar_check_approximation_3d
  -file_name hex_cub.h5m
  -base demkowicz -space h1
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-16 -ksp_rtol 0
  -log_no_color
)


# Check approximation on Hcurl space in 2d
build_test(hcurl_check_approx_in_2d 
  hcurl_check_approx_in_2d.cpp)

add_test(hcurl_check_approx_in_2d_ainsworth
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_check_approx_in_2d
  -base ainsworth -ksp_type fgmres -pc_type bjacobi -ksp_atol 1e-12 -ksp_rtol 0
  -log_no_color
)

build_test(hcurl_divergence_operator_2d hcurl_divergence_operator_2d.cpp)

add_test(hcurl_divergence_operator_2d_ainsworth_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_divergence_operator_2d
  -my_file rectangle_tri.h5m
  -base ainsworth -order 5
  -log_no_color)
add_test(hcurl_divergence_operator_2d_demkowicz_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_divergence_operator_2d
  -my_file rectangle_tri.h5m
  -base demkowicz -order 5
  -log_no_color)
add_test(hcurl_divergence_operator_2d_demkowicz_on_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_divergence_operator_2d
  -my_file rectangle_quad.h5m
  -base demkowicz -order 5
  -log_no_color)

build_test(hcurl_divergence_operator_2d_embaded_in_3d 
    hcurl_divergence_operator_2d_embaded_in_3d.cpp)

add_test(hcurl_divergence_operator_2d_ainsworth_on_tri_embaded_in_3d
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_divergence_operator_2d_embaded_in_3d
  -my_file half_of_sphere.h5m
  -base ainsworth -order 5
  -log_no_color)

# add test for recording series
build_test(record_series_atom record_series_atom.cpp)

add_test(
    record_series_atom
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/record_series_atom
    -DMPI_NP=1
    -DFILE=elasticity_atom_test_01X_pressure.h5m
    -DLOG1=record_series_atom.txt
    -DLOG2=record_series_atom.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

# add test for hcurl space
build_test(forces_and_sources_hcurl_approximation_functions 
  forces_and_sources_hcurl_approximation_functions.cpp)

add_test(
    forces_and_sources_hcurl_approximation_functions
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hcurl_approximation_functions
    -DFILE=empty
    -DLOG1=forces_and_sources_hcurl_approximation_functions.txt
    -DLOG2=forces_and_sources_hcurl_approximation_functions.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(hcurl_continuity_check hcurl_continuity_check.cpp)

add_test(hcurl_ainsworth_continuity_check
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_continuity_check 
  -my_file sphere.h5m -base ainsworth
  -log_no_color
)
add_test(hcurl_demkowicz_continuity_check
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_continuity_check 
  -my_file sphere.h5m -base demkowicz
  -log_no_color
)

# add test for hdiv space
build_test(forces_and_sources_hdiv_approximation_functions 
  forces_and_sources_hdiv_approximation_functions.cpp)

add_test(
    forces_and_sources_hdiv_approximation_functions
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_hdiv_approximation_functions
    -DFILE=empty.h5m
    -DLOG1=forces_and_sources_hdiv_approximation_functions.txt
    -DLOG2=forces_and_sources_hdiv_approximation_functions.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

# add test for hdiv space - check approximation continuity
build_test(hdiv_continuity_check hdiv_continuity_check.cpp)

add_test(hdiv_continuity_check_ainsworth_base
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_continuity_check
  -my_file sphere.h5m -base hdiv_ainsworth
  -log_no_color
)
add_test(hdiv_continuity_check_demkowicz_base
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_continuity_check
  -my_file sphere.h5m -base hdiv_demkowicz
  -log_no_color
)

# add test for hdvi divergence operator
build_test(hdiv_divergence_operator hdiv_divergence_operator.cpp)

add_test(hdiv_divergence_operator_ainsworth_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_divergence_operator
  -file_name sphere.h5m -base ainsworth
  -log_no_color)
add_test(hdiv_divergence_operator_demkowicz_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_divergence_operator
  -file_name sphere.h5m -base demkowicz
  -log_no_color)
add_test(hdiv_divergence_operator_ainsworth_ho_geometry_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_divergence_operator
  -file_name sphere.h5m -base ainsworth -ho_geometry
  -log_no_color)
add_test(hdiv_divergence_operator_demkowicz_ho_geometry_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/hdiv_divergence_operator
  -file_name sphere.h5m -base demkowicz -ho_geometry
  -log_no_color)

# add test for hcurl curl operator
build_test(hcurl_curl_operator hcurl_curl_operator.cpp)

add_test(hcurl_curl_operator_ainsworth_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_curl_operator 
  -file_name sphere.h5m -base ainsworth
  -log_no_color)
add_test(hcurl_curl_operator_demkowicz_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_curl_operator 
  -file_name sphere.h5m -base demkowicz
  -log_no_color)
add_test(hcurl_curl_operator_ainsworth_ho_geometry_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_curl_operator 
  -file_name sphere.h5m -base ainsworth -ho_geometry
  -log_no_color)
add_test(hcurl_curl_operator_demkowicz_ho_geometry_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/hcurl_curl_operator 
  -file_name sphere.h5m -base demkowicz -ho_geometry
  -log_no_color)

#test integration on skeleton
build_test(
  continuity_check_on_skeleton_3d
  continuity_check_on_skeleton_3d.cpp
)

add_test(
  continuity_check_on_skeleton_3d_ainsworth_hdiv
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_3d
  -my_file sphere.h5m -base ainsworth -space hdiv
  -log_no_color
)
add_test(
  continuity_check_on_skeleton_3d_demkowicz_hdiv
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_3d
  -my_file sphere.h5m -base demkowicz -space hdiv
  -log_no_color
)
add_test(
  continuity_check_on_skeleton_3d_ainsworth_hcurl
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_3d
  -my_file sphere.h5m -base ainsworth -space hcurl
  -log_no_color
)
add_test(
  continuity_check_on_skeleton_3d_demkowicz_hcurl
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_3d
  -my_file sphere.h5m -base demkowicz -space hcurl
  -log_no_color
)

#dm
build_test(dm_mofem dm_mofem.cpp)

add_test(
    dm_mofem
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=2
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/dm_mofem
    -DFILE=disp01.h5m
    -DLOG1=dm_mofem.txt
    -DLOG2=dm_mofem.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

cm_export_file("disp01_2parts.h5m" export_atom_tests)
build_test(dm_build_partitioned_mesh dm_build_partitioned_mesh.cpp)

if(PETSc_VERSION VERSION_GREATER_EQUAL 3.14.0)
 add_test(
    dm_build_partitioned_mesh
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=2
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/dm_build_partitioned_mesh
    -DFILE=disp01_2parts.h5m 
    -DLOG1=dm_build_partitioned_mesh.txt
    -DLOG2=dm_build_partitioned_mesh.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
else()
 add_test(
    dm_build_partitioned_mesh
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=2
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/dm_build_partitioned_mesh
    -DFILE=disp01_2parts.h5m 
    -DLOG1=dm_build_partitioned_mesh_petc3.11.txt
    -DLOG2=dm_build_partitioned_mesh.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
endif()

cm_export_file("disp01_3parts.h5m" export_atom_tests)
build_test(dm_partitioned_no_field dm_partitioned_no_field.cpp)

add_test(dm_partitioned_no_field_2parts
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2 
  ${CMAKE_CURRENT_BINARY_DIR}/dm_partitioned_no_field -my_file disp01_2parts.h5m 
  -dm_is_partitioned
  -log_no_color)
add_test(dm_partitioned_no_field_3parts
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 3
  ${CMAKE_CURRENT_BINARY_DIR}/dm_partitioned_no_field -my_file disp01_3parts.h5m 
  -dm_is_partitioned
  -log_no_color)

build_test(dm_create_subdm dm_create_subdm.cpp)
add_test(
  dm_create_subdm
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/dm_create_subdm -my_file disp01_2parts.h5m 
  -dm_is_partitioned
  -log_no_color
)

#single process
build_test(serial_matrix serial_matrix.cpp)

if(PETSc_VERSION VERSION_GREATER_EQUAL 3.14.0)
  add_test(
    serial_matrix
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/serial_matrix
    -DFILE=disp01.h5m
    -DLOG1=serial_matrix.txt
    -DLOG2=serial_matrix.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
else()
  add_test(
    serial_matrix
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/serial_matrix
    -DFILE=disp01.h5m
    -DLOG1=serial_matrix_petc3.11.txt
    -DLOG2=serial_matrix.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)
endif()


#test for build problems
build_test(build_problems build_problems.cpp)
add_test(
  build_problems_tet
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/build_problems -my_file disp01.h5m
  -log_no_color
)
add_test(
  build_problems_hex
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/build_problems -my_file hex_cub.h5m
  -log_no_color
)

#test for build problems
build_test(build_large_problem build_large_problem.cpp)

add_test(
  build_large_problem
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/build_large_problem -my_file disp01.h5m
  -log_no_color 
)

# tetsing removing enetities from the problem
cm_export_file("disp01_2parts.h5m" export_atom_tests)
build_test(remove_entities_from_problem 
  remove_entities_from_problem.cpp
)

add_test(remove_entities_from_problem
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2 
  ${CMAKE_CURRENT_BINARY_DIR}/remove_entities_from_problem 
  -my_file disp01_2parts.h5m 
  -log_no_color 
)

build_test(remove_entities_from_problem_not_partitioned 
  remove_entities_from_problem_not_partitioned.cpp
)

add_test(remove_entities_from_problem_not_partitioned
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2 
  ${CMAKE_CURRENT_BINARY_DIR}/remove_entities_from_problem_not_partitioned 
  -my_file disp01_2parts.h5m 
  -log_no_color 
)

# teting build mesh on partinoned (not distributed) mesh and composite problems
#test for build problems
build_test(build_composite_problem build_composite_problem.cpp)

add_test(
  build_composite_problem_atom
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/build_composite_problem -my_file disp01.h5m
  -log_no_color
)

#test inactivating and clearing inactive dofs
build_test(activate_deactivate_dofs activate_deactivate_dofs.cpp)
add_test(
  activate_deactivate_dofs
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/activate_deactivate_dofs -my_file disp01.h5m
  -log_no_color
)

# Tetsing entity loops
build_test(loop_entities loop_entities.cpp)

add_test(loop_entities
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2 
  ${CMAKE_CURRENT_BINARY_DIR}/loop_entities -file_name disp01_2parts.h5m 
  -dm_is_partitioned
  -log_no_color
)

# Projection form and to 10 node tets
build_test(projection_from_10node_tet_atom projection_from_10node_tet.cpp)
add_test(
    projection_from_10node_tet_atom
    ${CMAKE_COMMAND}
    -DMPI_RUN=${MPI_RUN}
    -DMPI_FLAGS=${MPI_RUN_FLAGS} --oversubscribe
    -DMPI_NP=1
    -DPROG=${CMAKE_CURRENT_BINARY_DIR}/projection_from_10node_tet_atom
    -DFILE=sphere.h5m
    -DLOG1=10node_sphere.txt
    -DLOG2=10node_sphere.txt
    -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}/atom_tests/blessed_files
    -P ${PROJECT_SOURCE_DIR}/cmake/RunAndComp.cmake)

build_test(projection_from_10node_tet_vol_atom projection_from_10node_op_volume.cpp)
add_test(
    projection_from_10node_tet_vol_atom
    ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
    ${CMAKE_CURRENT_BINARY_DIR}/projection_from_10node_tet_vol_atom -my_file unit_volume_sphere.cub)

if(TETGEN_LIBRARY)

  cm_export_file("../meshes/Tinterface.h5m" export_atom_tests)

  build_test(tetgen_interface tetgen_interface.cpp)

  if(WITH_TETGEN)
    add_dependencies(tetgen_interface tetgen)
  endif(WITH_TETGEN)

  add_test(tetgen_interface_atom ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/tetgen_interface 
    -my_file Tinterface.h5m
    -log_no_color
  )
endif(TETGEN_LIBRARY)

if(NETGEN_LIBRARY)

  build_test(netgen_interface netgen_interface.cpp)


endif(NETGEN_LIBRARY)

build_test(node_merge node_merge.cpp)

cm_export_file("../meshes/node_merge_test.h5m" export_atom_tests)
add_test(node_merge ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/node_merge 
  -my_file node_merge_test.h5m
  -log_no_color
)

build_test(edge_and_bubble_shape_functions_on_quad edge_and_bubble_shape_functions_on_quad.cpp)

add_test(edge_and_bubble_shape_functions_on_quad ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/edge_and_bubble_shape_functions_on_quad 
  -my_file quad_mesh.h5m
  -log_no_color
)

build_test(gauss_points_on_outer_product gauss_points_on_outer_product.cpp)

add_test(gauss_points_on_outer_product ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/gauss_points_on_outer_product
)

build_test(partition_mesh partition_mesh.cpp)
add_test(partition_mesh ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/partition_mesh -my_file disp01.h5m
  -log_no_color
)

build_test(testing_base_functions testing_base_functions.cpp)

list(INSERT testing_base_list
  0
  legendre
  lobatto
  jacobi
  integrated_jacobi
  h1tet_ainsworth
  h1tet_bernstein_bezier 
  hcurltet_ainsworth
  hcurltet_demkowicz
  hdivtet_ainsworth
  hdivtet_demkowicz
  l2tet
  h1tri_ainsworth
  h1tri_bernstein_bezier
  h1quad
  hdivtri_ainsworth
  hdivtri_demkowicz
  hcurltri_ainsworth
  hcurltri_demkowicz
  l2tri
  h1edge_ainsworth
  h1edge_bernstein_bezier
  hcurledge_ainsworth
  hcurledge_demkowicz
)
foreach(test_base ${testing_base_list})
  add_test(testing_base_functions_${test_base}
    ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/testing_base_functions 
    -base ${test_base}
    -log_no_color
  )
endforeach(test_base IN testing_base_list)

if(ADOL-C_LIBRARY)
  build_test(adol-c_atom adol-c_atom.cpp)
  if(WITH_ADOL-C)
    add_dependencies(adol-c_atom adolc)
  endif(WITH_ADOL-C)
  add_test(
    adol-c_atom_test ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/adol-c_atom
    -log_no_color
  )
endif(ADOL-C_LIBRARY)

if(MED_LIBRARY)
  cm_export_file("../meshes/Mesh_1.med" export_atom_tests)
  build_test(reading_med_file reading_med_file.cpp)
  add_test(
    reading_med_file_test ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/reading_med_file -med_file Mesh_1.med
    -log_no_color
  )
endif(MED_LIBRARY)

# Simple interface
build_test(simple_interface simple_interface.cpp)
add_test(
  simple_interface ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/simple_interface -file_name disp01_2parts.h5m 
  -log_summary
  -log_no_color
)

build_test(test_cache_on_entities test_cache_on_entities.cpp)

add_test(
  test_cache_on_entities ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/test_cache_on_entities 
  -file_name disp01_2parts.h5m 
  -log_no_color
)

# Tetsing mesh cuttiong

cm_export_file("../meshes/test_mesh_cut.cub" export_atom_tests)
cm_export_file("../meshes/circle_cut.h5m" export_atom_tests)
cm_export_file("../meshes/cylinder_cut.h5m" export_atom_tests)

build_test(mesh_cut_test mesh_cut_test.cpp)
add_test(
  mesh_cut_circle ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/mesh_cut_test 
  -my_file circle_cut.h5m 
  -side_set 200 
  -shift -0.2,0.2,0.
  -tol 1e-4,3.3e-1,3.3e-1 -test
  -log_no_color
)

build_test(segments_distance segments_distance.cpp)
add_test(
  segments_distance ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/segments_distance 
)

build_test(find_local_coordinates_on_tet find_local_coordinates_on_tet.cpp)
add_test(
  find_local_coordinates_on_tet ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/find_local_coordinates_on_tet 
)

cm_export_file("disp01_1parts.h5m" export_atom_tests)
build_test(field_evaluator field_evaluator.cpp)
add_test(field_evaluator
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/field_evaluator -file_name disp01_1parts.h5m 
  -dm_is_partitioned
  -log_no_color
)

build_test(petsc_smart_ptr_objects petsc_smart_ptr_objects.cpp)

add_test(petsc_smart_ptr_objects
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/petsc_smart_ptr_objects
  -log_no_color
)

build_test(boundary_marker boundary_marker.cpp)

add_test(boundary_marker
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/boundary_marker 
  -file_name rectangle_tri.h5m
  -log_no_color
)

build_test(add_blockset add_blockset.cpp)

add_test(add_blockset
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 ${CMAKE_CURRENT_BINARY_DIR}/add_blockset 
  -log_no_color
)

#test integration on contact prism side element
build_test(
  continuity_check_on_contact_prism_side_ele
  continuity_check_on_contact_prism_side_ele.cpp
)

add_test(
  continuity_check_on_contact_prism_side_ele_h1
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_contact_prism_side_ele
  -my_file 2_tets_contact.cub -base ainsworth -space h1)

# continuity check for hcurl 2d

build_test(
  continuity_check_on_skeleton_with_simple_2d_for_hcurl
  continuity_check_on_skeleton_with_simple_2d_for_hcurl.cpp)
  
add_test(
  continuity_check_on_skeleton_with_simple_2d_for_hcurl_ainsworth_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hcurl 
  -file_name rectangle_tri.h5m -base ainsworth
  -log_no_color)

add_test(continuity_check_on_skeleton_with_simple_2d_for_hcurl_demkowicz_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hcurl 
  -file_name rectangle_tri.h5m -base demkowicz
  -log_no_color)

add_test(
  continuity_check_on_skeleton_with_simple_2d_for_hcurl_demkowicz_on_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hcurl 
  -file_name rectangle_quad.h5m -base demkowicz
  -log_no_color)

# continuity check for hdiv 2d

build_test(
  continuity_check_on_skeleton_with_simple_2d_for_hdiv
  continuity_check_on_skeleton_with_simple_2d_for_hdiv.cpp)
  
add_test(continuity_check_on_skeleton_with_simple_2d_for_hdiv_ainsworth_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hdiv 
  -file_name rectangle_tri.h5m -base ainsworth
  -log_no_color)

add_test(continuity_check_on_skeleton_with_simple_2d_for_hdiv_demkowicz_on_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hdiv 
  -file_name rectangle_tri.h5m -base demkowicz
  -log_no_color)

add_test(continuity_check_on_skeleton_with_simple_2d_for_hdiv_demkowicz_on_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_hdiv 
  -file_name rectangle_quad.h5m -base demkowicz
  -log_no_color)

# continuity check for h1 2d

build_test(
  continuity_check_on_skeleton_with_simple_2d_for_h1
  continuity_check_on_skeleton_with_simple_2d_for_h1.cpp)
  
add_test(continuity_check_on_skeleton_with_simple_2d_for_h1_ainsworth_tri
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_h1 
  -file_name rectangle_tri.h5m -base ainsworth
  -log_no_color)

add_test(continuity_check_on_skeleton_with_simple_2d_for_h1_debkowicz_quad
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/continuity_check_on_skeleton_with_simple_2d_for_h1 
  -file_name rectangle_quad.h5m -base demkowicz
  -log_no_color)


cm_export_file("../meshes/reference_triangle_mesh.h5m" export_atom_tests)
build_test(prism_elements_from_surface prism_elements_from_surface.cpp)

add_test(prism_elements_from_surface ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/prism_elements_from_surface 
  -my_file reference_triangle_mesh.h5m
  -log_no_color
)

build_test(
  prism_polynomial_approximation prism_polynomial_approximation.cpp)

add_test(prism_polynomial_approximation
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/prism_polynomial_approximation -pc_type lu)

build_test(
  quad_polynomial_approximation quad_polynomial_approximation.cpp)

add_test(
  quad_polynomial_approximation_h1_ainsworth
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/quad_polynomial_approximation 
  -pc_type lu
  -base ainsworth
  -space h1
  -log_no_color)

add_test(
  quad_polynomial_approximation_h1_ainsworth_lobatto
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/quad_polynomial_approximation 
  -pc_type lu
  -base ainsworth_labatto
  -space h1
  -log_no_color)

add_test(
  quad_polynomial_approximation_h1_demkowicz
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/quad_polynomial_approximation 
  -pc_type lu
  -base demkowicz
  -space h1
  -log_no_color)

add_test(
  quad_polynomial_approximation_l2_demkowicz
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/quad_polynomial_approximation 
  -pc_type lu
  -base demkowicz
  -space l2
  -log_no_color)

build_test(
  bernstein_bezier_generate_base
  bernstein_bezier_generate_base.cpp
)
add_test(bernstein_bezier_generate_base
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/bernstein_bezier_generate_base 
  -log_no_color
)

build_test(
  log
  log.cpp
)

add_test(log
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/log
  -log_severity_level noisy
  -log_scope
  -log_time
  -log_no_color
)

# Matrix function

build_test(
  matrix_function
  matrix_function.cpp
)

add_test(matrix_function
  ${MPI_RUN} ${MPI_RUN_FLAGS} --oversubscribe -np 1 
  ${CMAKE_CURRENT_BINARY_DIR}/matrix_function
  -log_severity_level inform
  -log_no_color
)
