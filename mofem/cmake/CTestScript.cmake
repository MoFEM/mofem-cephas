set(CTEST_PROJECT_NAME "MoFEM")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "Debug")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_GIT_COMMAND NAMES git)

# MoFEM lib
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(INIT_REPOSITORY "YES")
  set(
    CTEST_CHECKOUT_COMMAND
    "${CTEST_GIT_COMMAND} clone --branch ${CTEST_BRANCH} --recurse-submodules https://bitbucket.org/likask/mofem-cephas.git ${GID_SOURCE_REPO}"
  )
else(EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} submodule update")
endif()
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "\"${CMAKE_COMMAND}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_BUILD_OPTIONS} -DWITHCOVERAGE=ON")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

#Ctest time outr
set(CTEST_TEST_TIMEOUT 1200)

# Perform the CDashTesting
ctest_start(${DASHBOARDTEST})

ctest_update(SOURCE "${GID_SOURCE_REPO}" RETURN_VALUE DOTEST)

if(INIT_REPOSITORY)
  set(DOTEST 1)
  message("Force Init Build")
else(NOT INIT_REPOSITORY)
  message ( "Found ${DOTEST} updated files." )
endif()
if(FORCETESTING)
  set(DOTEST 1)
  message ("Force build")
endif(FORCETESTING)

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
  field_axpy_compare
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
  forces_and_sources_testing_edge_element_test_compare
  forces_and_sources_testing_vertex_element_test_compare
  forces_and_sources_testing_flat_prism_element_test_compare
  record_series_atom_test_compare
  forces_and_sources_hcurl_approximation_functions_atom_compare
  forces_and_sources_hdiv_approximation_functions_atom_compare
  dm_mofem_atom_compare
  dm_build_partitioned_mesh_atom_compare
  serial_matrix_compare
  projection_from_10node_tet_atom_compare
)

if(${DOTEST} GREATER 0)
  file(WRITE ${GID_SOURCE_REPO}/has_bin_build "1")
endif(${DOTEST} GREATER 0)

if(${DOTEST} GREATER 0)
  ctest_configure()
  ctest_build()
  if(CTEST_MEMORYCHECK_COMMAND)
    ctest_memcheck()
  endif(CTEST_MEMORYCHECK_COMMAND)
  ctest_test()
  if(CTEST_COVERAGE_COMMAND)
    ctest_coverage(QUIET)
  endif(CTEST_COVERAGE_COMMAND)
  ctest_submit()
endif(${DOTEST} GREATER 0)
