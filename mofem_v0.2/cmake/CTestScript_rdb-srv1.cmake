set(CTEST_BUILD_OPTIONS "-DPETSC_DIR=/opt/petsc-3.5.1 -DPETSC_ARCH=arch-linux2-c-debug -DMOAB_DIR=/opt/local_new_moab -DCMAKE_INSTALL_PREFIX=/tmp/cephas/debug_examples /home/lukasz/mofem-cephas/mofem_v0.2 -DCMAKE_CXX_FLAGS=-I/opt/local_boost_1_54_0/include")

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
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{HOME}/tmp/mofem/source/mofem_v0.2/cmake/rdb-srv1-valgrind.supp")

include(CTestScript.cmake)


