set(CTEST_BUILD_OPTIONS "/tmp/cephas/debug_examples/users_modules")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST) 
  set(FORCETESTING "YES")
  set(DASHBOARDTEST "Nightly")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


