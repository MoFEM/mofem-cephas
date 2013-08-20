set(CTEST_BUILD_OPTIONS "-DCMAKE_CXX_FLAGS=-lstdc++ -DPETSC_DIR=/opt/build_for_gcc-mp-4.4/petsc-3.3-p3 -DPETSC_ARCH=arch-darwin-c-opt -DMOAB_DIR=/opt/build_for_gcc-mp-4.4/local-moab-4.6.0/ -DBUILD_SHARED_LIBS=FALSE -DCMAKE_CXX_COMPILER=/opt/build_for_gcc-mp-4.4/local/bin/mpicxx")

set(CTEST_SITE "lukaszs-laptop.lan")
set(CTEST_BUILD_NAME "Darwin-mpicxx")
if(NOT DASHBOARDTEST) 
  set(DASHBOARDTEST "Nightly")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


