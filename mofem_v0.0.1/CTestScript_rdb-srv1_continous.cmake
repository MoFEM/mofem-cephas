set(CTEST_BUILD_OPTIONS "-DPETSC_DIR=/opt/petsc-3.3-p5_openmpi-1.6.3 -DPETSC_ARCH=arch-linux2-c-opt -DMOAB_DIR=/opt/local-moab-4.6.0/ -DCMAKE_CXX_FLAGS=-I/opt/local_boost_1_54_0/include")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST) 
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


