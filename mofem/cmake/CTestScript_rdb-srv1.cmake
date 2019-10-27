set(CTEST_BUILD_OPTIONS "-DPETSC_DIR=/opt/petsc -DPETSC_ARCH=arch-linux2-c-debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMOAB_DIR=/opt/local_new_moab -DADOL-C_DIR=/opt/local_adol-c-2.5.2 -DTETGEN_DIR=/opt/tetgen1.5.0 -DMED_DIR=/opt/med -DCMAKE_INSTALL_PREFIX=/home/lukasz/tmp/cephas_users_modules -DBOOST_DIR=/opt/local_boost_1_65_1 -DSTAND_ALLONE_USERS_MODULES=1")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")
set(CTEST_BRANCH "CDashTesting")

if(NOT DASHBOARDTEST)
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

set(CTEST_CUSTOM_MEMCHECK_IGNORE
  ${CTEST_CUSTOM_MEMCHECK_IGNORE}
  #compare
  cubit_bc_atom_test_disp01_compare
  cubit_bc_atom_test_force01_compare
  cubit_bc_atom_test_velocity01_compare
  cubit_bc_atom_test_accel01_compare
  cubit_bc_atom_test_temper01_compare
  cubit_bc_atom_test_pressure01_compare
  cubit_bc_atom_test_heatflux01_compare
  cubit_bc_atom_test_comb01_compare
  cubit_bc_atom_test_bcoverlap01_compare
  cubit_bc_atom_test_interf01_compare
  cubit_bc_atom_test_mat_elastic_compare
  cubit_bc_atom_test_mat_elastic_transiso_compare
  cubit_bc_atom_test_mat_interf_compare
  cubit_bc_atom_test_inlet_outlet_compare
  cubit_meshset_loop_test_compare
  mesh_refine_atom_test_compare
  mesh_insert_interface_atom_test_compare
  mesh_insert_T_interface_atom_test_compare
  mesh_insert_T_4seasons_interface_atom_test_compare
  forces_and_sources_getting_orders_indices_test_compare
  forces_and_sources_getting_mult_H1_H1_test_compare
  forces_and_sources_calculate_jacobian_test_compare
  forces_and_sources_testing_volume_element_test_compare
  forces_and_sources_testing_users_base_compare
  forces_and_sources_getting_higher_order_skin_normals_atom_test_compare
  forces_and_sources_testing_triangle_element_test_compare
  forces_and_sources_testing_edge_element_ainsworth_test_compare
  forces_and_sources_testing_edge_element_bernstein_bezier_test_compare
  forces_and_sources_testing_vertex_element_test_compare
  forces_and_sources_testing_flat_prism_element_test_compare
  forces_and_sources_testing_contact_prism_element_test_compare
  record_series_atom_test_compare
  forces_and_sources_hcurl_approximation_functions_atom_compare
  forces_and_sources_hdiv_approximation_functions_atom_compare
  dm_mofem_atom_compare
  dm_build_partitioned_mesh_atom_compare
  serial_matrix_compare
  projection_from_10node_tet_atom_compare
)

#valgrind set up
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS
  "--trace-children=yes --quiet --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50 --verbose --demangle=yes --gen-suppressions=all")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{HOME}/tmp/cephas/source/mofem/cmake/rdb-srv1-valgrind.supp")

set(GID_SOURCE_REPO "$ENV{HOME}/tmp/cephas/source")
set(CTEST_SOURCE_DIRECTORY "${GID_SOURCE_REPO}/mofem")
set(CTEST_BINARY_DIRECTORY "$ENV{HOME}/tmp/cephas/build")

include(CTestScript.cmake)
