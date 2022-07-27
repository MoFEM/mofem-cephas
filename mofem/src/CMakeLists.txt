#include subdirectory

include_directories("${PROJECT_SOURCE_DIR}/third_party")
include_directories("${PROJECT_SOURCE_DIR}/third_party/cblas")

add_subdirectory(${PROJECT_SOURCE_DIR}/src/ftensor)

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories("${PROJECT_SOURCE_DIR}/src/approximation")
include_directories("${PROJECT_SOURCE_DIR}/src/approximation/c")
include_directories("${PROJECT_SOURCE_DIR}/src/multi_indices")
include_directories("${PROJECT_SOURCE_DIR}/src/interfaces")
include_directories("${PROJECT_SOURCE_DIR}/src/petsc")
include_directories("${PROJECT_SOURCE_DIR}/src/finite_elements")

link_directories("${PROJECT_BINARY_DIR}/src/approximation")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/approximation")

link_directories("${PROJECT_BINARY_DIR}/src/multi_indices")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/multi_indices")

link_directories("${PROJECT_BINARY_DIR}/src/interfaces")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/interfaces")

link_directories("${PROJECT_BINARY_DIR}/src/petsc")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/petsc")

link_directories("${PROJECT_BINARY_DIR}/src/finite_elements")
add_subdirectory("${PROJECT_SOURCE_DIR}/src/finite_elements")

if(STAND_ALLONE_USERS_MODULES)
  install(
    DIRECTORY 
    ${CMAKE_SOURCE_DIR}/src/
    DESTINATION 
    ${CMAKE_INSTALL_PREFIX}/include
    FILES_MATCHING 
    PATTERN "*.hpp"
    PATTERN "*.h"
    PATTERN "impl" EXCLUDE
    PATTERN "c_impl" EXCLUDE
    PATTERN "ftensor" EXCLUDE)
  else(STAND_ALLONE_USERS_MODULES) 
    file(
     GLOB DIRS
     RELATIVE ${CMAKE_SOURCE_DIR}/src
     ${CMAKE_SOURCE_DIR}/src/*
     ${CMAKE_SOURCE_DIR}/src/*/c)
    foreach(DIR ${DIRS})
      if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/${DIR})
        install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${DIR})
        file(
          GLOB HEADER_FILES
          RELATIVE ${CMAKE_SOURCE_DIR}/src/${DIR}
          ${CMAKE_SOURCE_DIR}/src/${DIR}/*.h*)
        foreach(HF ${HEADER_FILES})
          install(CODE
            "
            EXECUTE_PROCESS(
            COMMAND 
            ln -sf 
            ${CMAKE_SOURCE_DIR}/src/${DIR}/${HF} 
            ${CMAKE_INSTALL_PREFIX}/include/${DIR}/${HF}
            )
            MESSAGE(\"-- Linking ${CMAKE_INSTALL_PREFIX}/include/${DIR}/${HF}\")
            "
          )
        endforeach(HF)
      endif(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/${DIR})
    endforeach(DIR)
endif(STAND_ALLONE_USERS_MODULES)
