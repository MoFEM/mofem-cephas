set(
  CTEST_BUILD_OPTIONS
  "-DSTAND_ALLONE_USERS_MODULES=ON"
  "-DCMAKE_CXX_FLAGS=-I/opt/local_boost_1_54_0/include"
  "-DCMAKE_EXE_LINKER_FLAGS=-L/opt/local_boost_1_54_0/lib"
  "/home/lukasz/tmp/cephas_users_modules/users_modules"
)

set(CTEST_SITE "rdb-srv1")
set(CTEST_BUILD_NAME "Linux-mpicxx")

if(NOT DASHBOARDTEST)
  set(DASHBOARDTEST "Continuous")
endif(NOT DASHBOARDTEST)

# modules - moisture_transport
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/users_modules/moisture_transport")
  exec_program(
    ${CTEST_GIT_COMMAND}
    "${CTEST_SOURCE_DIRECTORY}/users_modules"
    ARGS clone https://likask@bitbucket.org/likask/mofem_um_moisture_transport.git
    "${CTEST_SOURCE_DIRECTORY}/users_modules/moisture_transport"
  )
else(EXISTS "${CTEST_SOURCE_DIRECTORY}/users_modules/moisture_transport")
  exec_program(
    ${CTEST_GIT_COMMAND}
    "${CTEST_SOURCE_DIRECTORY}/users_modules/moisture_transport"
    ARGS pull
  )
endif()

# modules - ground_surface_temperature
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/users_modules/ground_surface_temperature")
  exec_program(
    ${CTEST_GIT_COMMAND}
    "${CTEST_SOURCE_DIRECTORY}/users_modules"
    ARGS clone https://likask@bitbucket.org/likask/mofem_um_ground_surface_temperature.git
    "${CTEST_SOURCE_DIRECTORY}/users_modules/ground_surface_temperature"
  )
else(EXISTS "${CTEST_SOURCE_DIRECTORY}/users_modules/ground_surface_temperature")
  exec_program(
    ${CTEST_GIT_COMMAND}
    "${CTEST_SOURCE_DIRECTORY}/users_modules/ground_surface_temperature"
    ARGS pull
  )
endif()

set(CTEST_SOURCE_DIRECTORY "/home/lukasz/tmp/cephas_users_modules/users_modules")
set(CTEST_BINARY_DIRECTORY "/home/lukasz/tmp/cephas_users_modules/build")

include(CTestScript.cmake)
