set(CTEST_BUILD_OPTIONS "-DPETSC_DIR=/opt/petsc -DPETSC_ARCH=arch-linux2-c-debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMOAB_DIR=/opt/local_new_moab -DADOL-C_DIR=/opt/local_adol-c-2.5.2 -DTETGEN_DIR=/opt/tetgen1.5.0 -DMED_DIR=/opt/med -DCMAKE_INSTALL_PREFIX=/home/lukasz/tmp/cephas_users_modules -DBOOST_DIR=/opt/local_boost_1_54_0")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")
set(CTEST_BRANCH "CDashTesting")

if(NOT DASHBOARDTEST)
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

#valgrind set up
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS
  "--trace-children=yes --quiet --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50 --verbose --demangle=yes --gen-suppressions=all")
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{HOME}/tmp/cephas/source/mofem/cmake/rdb-srv1-valgrind.supp")

set(GID_SOURCE_REPO "$ENV{HOME}/tmp/cephas/source")
set(CTEST_SOURCE_DIRECTORY "${GID_SOURCE_REPO}/mofem")
set(CTEST_BINARY_DIRECTORY "$ENV{HOME}/tmp/cephas/build")

include(CTestScript.cmake)
