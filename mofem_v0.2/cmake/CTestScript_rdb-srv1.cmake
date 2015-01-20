set(CTEST_BUILD_OPTIONS "-DPETSC_DIR=/opt/petsc-3.5.1 -DPETSC_ARCH=arch-linux2-c-debug -DMOAB_DIR=/opt/local_new_moab -DADOL-C_DIR=/opt/local_adol-c-2.5.2 -DTETGEN_DIR=/opt/tetgen1.5.0 -DSLEPC_DIR=/opt/slepc-3.5.3 -DCMAKE_INSTALL_PREFIX=/home/lukasz/tmp/cephas_users_modules /home/lukasz/mofem-cephas/mofem_v0.2 -DCMAKE_CXX_FLAGS=\"-I/opt/local_boost_1_54_0/include -L/opt/local_boost_1_54_0/lib\"")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")
set(CTEST_BRANCH "CDashTesting")

if(NOT DASHBOARDTEST) 
  set(FORCETESTING "YES")
  set(DASHBOARDTEST "Nightly")
endif(NOT DASHBOARDTEST)

#valgrind set up
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS 
  "--trace-children=yes --quiet --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50 --verbose --demangle=yes --gen-suppressions=all")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{HOME}/tmp/cephas/source/mofem_v0.2/cmake/rdb-srv1-valgrind.supp")

include(CTestScript.cmake)


