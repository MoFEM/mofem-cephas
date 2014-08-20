set(CTEST_SOURCE_DIRECTORY "/home/lukasz/tmp/cephas_users_modules/users_modules")
set(CTEST_BINARY_DIRECTORY "/home/lukasz/tmp/cephas_users_modules/build")

set(CTEST_PROJECT_NAME "MoFEM-UsersModules")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "Debug")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

set(CTEST_CONFIGURE_COMMAND "\"${CMAKE_COMMAND}\" -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_BUILD_OPTIONS} -DWITHCOVERAGE=ON")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

#Ctest time outr
set(CTEST_TEST_TIMEOUT 1200)

# Perform the CDashTesting
ctest_start(${DASHBOARDTEST})

if(FORCETESTING) 
  set(DOTEST 1)
  message ("Force build")
endif(FORCETESTING)

set(CTEST_CUSTOM_MEMCHECK_IGNORE
  ${CTEST_CUSTOM_MEMCHECK_IGNORE}
  #compare
)

if(${DOTEST} GREATER 0)
  ctest_configure()
  ctest_build()
  ctest_test()
  if(CTEST_COVERAGE_COMMAND)
    ctest_coverage()
  endif(CTEST_COVERAGE_COMMAND)
  ctest_submit()
endif(${DOTEST} GREATER 0)

