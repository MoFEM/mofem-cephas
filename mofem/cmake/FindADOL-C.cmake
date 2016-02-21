# - Try to find ADOL-C
# Once done this will define
#
#  ADOL-C_DIR - directory in which ADOL-C resides

if(WITH_ADOL-C)
  ExternalProject_Add(
    adolc
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL http://www.coin-or.org/download/source/ADOL-C/ADOL-C-2.5.2.tgz
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/external/src/adolc/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --disable-shared --prefix=${PROJECT_BINARY_DIR}/external
  )
  include_directories(${PROJECT_BINARY_DIR}/external/include)
  link_directories(${PROJECT_BINARY_DIR}/external/lib64)
  set(ADOL-C_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to adol-c" FORCE)
  set(ADOL-C_LIBRARY ${ADOL-C_DIR}/lib64/libadolc.a CACHE FILEPATH  "adol-c lib" FORCE)
  # add_library(ADOL-C_LIBRARY STATIC IMPORTED)
  add_definitions(-DWITH_ADOL_C)
  find_library(COLPACK_LIBLARY NAMES Colpack PATHS /usr/local/lib)
  list(APPEND ADOL-C_LIBRARY ${COLPACK_LIBLARY})
  message(STATUS ${ADOL-C_LIBRARY})
elseif(ADOL-C_DIR)
  find_library(ADOL-C_LIBRARY NAMES adolc PATHS ${ADOL-C_DIR}/lib)
  message(STATUS ${ADOL-C_LIBRARY})
  if(ADOL-C_LIBRARY)
    include_directories("${ADOL-C_DIR}/include")
    add_definitions(-DWITH_ADOL_C)
  endif(ADOL-C_LIBRARY)
endif(WITH_ADOL-C)
