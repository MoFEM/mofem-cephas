set(CTEST_BUILD_OPTIONS "-DCMAKE_CXX_FLAGS=-I/opt/local_boost_1_54_0/include /home/lukasz/tmp/cephas_users_modules")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST) 
  set(FORCETESTING "YES")
  set(DASHBOARDTEST "Nightly")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


