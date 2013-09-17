set(GID_SOURCE_REPO "$ENV{HOME}/tmp/mofem/source")
set(CTEST_SOURCE_DIRECTORY "${GID_SOURCE_REPO}/mofem_v0.0.1")
set(CTEST_BINARY_DIRECTORY "$ENV{HOME}/tmp/mofem/build")

set(CTEST_PROJECT_NAME "MoFEM")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "Debug")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_GIT_COMMAND NAMES git)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(INIT_REPOSITORY "YES")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone --branch CDashTesting https://bitbucket.org/likask/mofem-joseph.git ${GID_SOURCE_REPO}")
else(EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} submodule update")
endif()
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "\"${CMAKE_COMMAND}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_BUILD_OPTIONS} -DWITHCOVERAGE=ON")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

# Perform the Nightly build
ctest_start(${DASHBOARDTEST})

ctest_update(SOURCE "${GID_SOURCE_REPO}" RETURN_VALUE DOTEST)
if(INIT_REPOSITORY) 
  set(DOTEST 1)
  message("Force Init Build")
else(NOT INIT_REPOSITORY)
  message ( "Found ${DOTEST} updated files." )
endif()
if(FORCETESTING) 
  set(DOTEST 1)
  message ("Fore build")
endif(FORCETESTING)

if(${DOTEST} GREATER 0)
  ctest_configure()
  ctest_build()
  ctest_test()
  if(CTEST_COVERAGE_COMMAND)
    ctest_coverage()
  endif(CTEST_COVERAGE_COMMAND)
  ctest_submit()
endif(${DOTEST} GREATER 0)

