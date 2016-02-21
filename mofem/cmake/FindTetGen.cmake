# - Try to find TETGEN

if(WITH_TETGEN)
  ExternalProject_Add(
    tetgen
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL https://bitbucket.org/likask/mofem-joseph/downloads/tetgen1.5.0.tgz
    INSTALL_COMMAND ""
  )
  add_custom_target(
    install_tetegen
    ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen/include/tetgen.h ${PROJECT_BINARY_DIR}/external/include/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen-build/libtet* ${PROJECT_BINARY_DIR}/external/lib/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_BINARY_DIR}/external/src/tetgen-build/tetgen ${PROJECT_BINARY_DIR}/external/bin/
    DEPENDS tetgen
  )
  include_directories(${PROJECT_BINARY_DIR}/external/include)
  link_directories(${PROJECT_BINARY_DIR}/external/lib)
  add_definitions(-DWITH_TETGEN)
  set(TETGEN_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to tetgen dir"  FORCE)
  set(TETGEN_LIBRARY ${TETGEN_DIR}/lib/libtet.a CACHE FILEPATH "tetgen lib"  FORCE)
  add_library(TETGEN_LIBRARY STATIC IMPORTED)
  message(STATUS ${TETGEN_LIBRARY})
elseif(TETGEN_DIR)
  find_library(TETGEN_LIBRARY NAMES tet PATHS ${TETGEN_DIR}/lib)
  message(STATUS ${TETGEN_LIBRARY})
  if(TETGEN_LIBRARY)
    include_directories(${TETGEN_DIR}/include)
    add_definitions(-DWITH_TETGEN)
  endif(TETGEN_LIBRARY)
endif(WITH_TETGEN)
