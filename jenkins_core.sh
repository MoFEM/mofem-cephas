cat > script.cmake << EOF

set(GID_SOURCE_REPO "\$ENV{WORKSPACE}")
set(CTEST_SOURCE_DIRECTORY "\${GID_SOURCE_REPO}/mofem")
set(CTEST_BINARY_DIRECTORY "\$ENV{WORKSPACE}/build")

set(CTEST_SITE "Cactus")
set(CTEST_BRANCH "\$ENV{GIT_BRANCH}")
set(CTEST_BUILD_NAME "Linux-Jenkins=BranchO(\${CTEST_BRANCH})")

find_program(CTEST_CONFIGURE_COMMAND NAMES spconfig.py HINTS \${CTEST_BINARY_DIRECTORY} NO_DEFAULT_PATH)
find_program(CTEST_BUILD_COMMAND NAMES make)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

#Ctest time out
set(CTEST_TEST_TIMEOUT 12000)

ctest_start(Continuous)

ctest_configure(
  OPTIONS "-DSOURCE_DIR=\${CTEST_SOURCE_DIRECTORY} -DWITHCOVERAGE=ON -DMOFEM_BUILD_TESTS=ON")

find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS
  "--trace-children=yes --quiet --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50 --demangle=yes --gen-suppressions=all")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "\${GID_SOURCE_REPO}/mofem/cmake/jenkins-valgrind.supp")  

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
  
ctest_build()
if(CTEST_MEMORYCHECK_COMMAND)
  ctest_memcheck()
endif(CTEST_MEMORYCHECK_COMMAND)
ctest_test()
if(CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif(CTEST_COVERAGE_COMMAND)
ctest_submit()

EOF

#if [ ! $GIT_COMMIT == $GIT_PREVIOUS_COMMIT ]; then

  SPACK_ROOT=/home/jenkins/workspace/SpackInstall
  . $SPACK_ROOT/share/spack/setup-env.sh

  if [ ! -f "$SPACK_ROOT/lock_spack" ]; then

    mkdir -p build
    cd build
  
    if [ ! -f "install_hash" ]; then
    
    	spack setup mofem-cephas@develop+slepc copy_user_modules=False build_type=Debug install_id=1
        
    fi
    
  fi

#fi

$WORKSPACE/build/spconfig.py \
	-DSOURCE_DIR=$WORKSPACE/mofem \
	-DMOFEM_BUILD_TESTS=ON -DWITHCOVERAGE=ON ../mofem
        
grep CMAKE_INSTALL_PREFIX CMakeCache.txt |\
sed 's/.*mofem-cephas-develop-//' |\
tee install_hash

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus start job for branch $GIT_BRANCH\"}" \
https://hooks.slack.com/services/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

cd $WORKSPACE/build
spack load cmake
ctest -V -S $WORKSPACE/script.cmake
make install
make clean

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus end job for branch $GIT_BRANCH (See http://cdash.eng.gla.ac.uk/cdash/)\"}" \
https://hooks.slack.com/services/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