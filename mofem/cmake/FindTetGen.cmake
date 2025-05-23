# - Try to find TETGEN

if(NOT TETGEN_DIR)
  set(TETGEN_DIR $ENV{TETGEN_DIR})
endif(NOT TETGEN_DIR)

if(TETGEN_DIR)
  find_library(
    TETGEN_LIBRARY 
    NAMES tet 
    PATHS ${TETGEN_DIR}/lib
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    CMAKE_FIND_ROOT_PATH_BOTH
  )
  find_path(TETGEN_HEADER
    NAMES tetgen.h
    PATHS ${TETGEN_DIR}/include
  )
  message(STATUS ${TETGEN_LIBRARY})
  message(STATUS ${TETGEN_HEADER})
  if(TETGEN_LIBRARY AND TETGEN_HEADER)
    include_directories(${TETGEN_HEADER})
    add_definitions(-DWITH_TETGEN)
  endif(TETGEN_LIBRARY AND TETGEN_HEADER)
endif(TETGEN_DIR)

if(WITH_TETGEN)
  ExternalProject_Add(
    tetgen
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL https://bitbucket.org/likask/mofem-joseph/downloads/tetgen-1.5.0.tgz
    CONFIGURE_COMMAND cmake ${PROJECT_BINARY_DIR}/external/src/tetgen
    BUILD_COMMAND make
    INSTALL_COMMAND ""
  )
  execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/external/lib/libtet.a)
  add_custom_target(
    copy_tetgen_fields
    ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen/include/tetgen.h ${PROJECT_BINARY_DIR}/external/include
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen-build/libtet* ${PROJECT_BINARY_DIR}/external/lib
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen-build/tetgen ${PROJECT_BINARY_DIR}/external/bin
    DEPENDS tetgen
  )
  add_dependencies(install_prerequisites tetgen copy_tetgen_fields)
endif(WITH_TETGEN)

if(WITH_TETGEN)
  if(NOT TETGEN_LIBRARY)
    set(TETGEN_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to tetgen dir"  FORCE)
    set(TETGEN_HEADER ${TETGEN_DIR}/include)
    include_directories(${TETGEN_HEADER})
    set(TETGEN_LIBRARY ${TETGEN_DIR}/lib/libtet.a)
    add_definitions(-DWITH_TETGEN)
    message(STATUS ${TETGEN_LIBRARY})
    message(STATUS ${TETGEN_HEADER})
  endif(NOT TETGEN_LIBRARY)
endif(WITH_TETGEN)
