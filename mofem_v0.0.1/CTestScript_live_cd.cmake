set(CTEST_BUILD_OPTIONS "-DCMAKE_CXX_FLAGS=-I/opt/local/include\\ -L/opt/local/lib\\ -Wno-unused-local-typedefs -DPETSC_DIR=/opt/petsc-3.4.3 -DPETSC_ARCH=arch-linux2-c-opt -DMOAB_DIR=/opt/local-moab-4.6.0")

set(CTEST_SITE "live-cd")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST)
  # set(FORCETESTING "YES")
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

# valgrind set up
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
SET(CTEST_MEMORYCHECK_COMMAND_OPTIONS 
  "--trace-children=yes --quiet --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=50 --verbose --demangle=yes --gen-suppressions=all")
SET(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{HOME}/tmp/mofem/source/mofem_v0.0.1/live-cd-valgrind.supp")

include(CTestScript.cmake)


