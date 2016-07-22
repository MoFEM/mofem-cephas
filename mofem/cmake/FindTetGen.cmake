# - Try to find TETGEN

if(NOT TETGEN_DIR)
  set(TETGEN_DIR $ENV{TETGEN_DIR})
endif(NOT TETGEN_DIR)

if(TETGEN_DIR)
  find_library(TETGEN_LIBRARY NAMES tet PATHS ${TETGEN_DIR}/lib)
  message(STATUS ${TETGEN_LIBRARY})
  if(TETGEN_LIBRARY)
    include_directories(${TETGEN_DIR}/include)
    add_definitions(-DWITH_TETGEN)
  endif(TETGEN_LIBRARY)
endif(TETGEN_DIR)

if(WITH_TETGEN)
  ExternalProject_Add(
    tetgen
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
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
    include_directories(${PROJECT_BINARY_DIR}/external/include)
    add_definitions(-DWITH_TETGEN)
    set(TETGEN_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to tetgen dir"  FORCE)
    set(TETGEN_LIBRARY IMPORTED_LOCATION ${TETGEN_DIR}/lib/libtet.a)
    message(STATUS ${TETGEN_LIBRARY})
  endif(NOT TETGEN_LIBRARY)
endif(WITH_TETGEN)
