set(CTEST_BUILD_OPTIONS "-DCMAKE_CXX_FLAGS=-lstdc++ -DCMAKE_BUILD_TYPE=Debug -DPETSC_DIR=/opt/build_for_gcc-mp-4.4/petsc-3.4.3 -DPETSC_ARCH=darwin10.2.0-c-debug -DMOAB_DIR=/opt/build_for_gcc-mp-4.4/local-moab-4.6.0/ -DBUILD_SHARED_LIBS=FALSE -DCMAKE_CXX_COMPILER=/opt/build_for_gcc-mp-4.4/local/bin/mpicxx")

set(CTEST_SITE "michael_laptop")
set(CTEST_BUILD_NAME "macos-x")

if(NOT DASHBOARDTEST) 
  set(FORCETESTING "YES")
  set(DASHBOARDTEST "Experimental")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


