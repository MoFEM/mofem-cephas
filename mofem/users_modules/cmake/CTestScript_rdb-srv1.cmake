set(CTEST_BUILD_OPTIONS "-DSTAND_ALLONE_USERS_MODULES=ON -DCMAKE_CXX_FLAGS=-I/opt/local_boost_1_54_0/include -DCMAKE_EXE_LINKER_FLAGS=-L/opt/local_boost_1_54_0/lib /home/lukasz/tmp/cephas_users_modules/users_modules")

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST) 
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

include(CTestScript.cmake)


