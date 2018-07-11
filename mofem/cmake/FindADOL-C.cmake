# - Try to find ADOL-C
# Once done this will define
#
#  ADOL-C_DIR - directory in which ADOL-C resides

if(NOT ADOL-C_DIR)
  set(ADOL-C_DIR $ENV{ADOL-C_DIR})
endif(NOT ADOL-C_DIR)

if(ADOL-C_DIR)
  find_library(ADOL-C_LIBRARY NAMES adolc PATHS ${ADOL-C_DIR}/lib ${ADOL-C_DIR}/lib64)
  find_path(ADOL-C_HEADER NAMES adolc/adolc.h PATHS ${ADOL-C_DIR}/include)
  if(ADOL-C_LIBRARY AND ADOL-C_HEADER)
    find_library(
      COLPACK_LIBLARY 
      NAMES Colpack 
      PATHS ${ADOL-C_DIR}/lib ${ADOL-C_DIR}/lib64 /usr/local/lib)
    if(COLPACK_LIBLARY)
      set(ADOL-C_LIBRARY ${ADOL-C_LIBRARY} ${COLPACK_LIBLARY})
    endif(COLPACK_LIBLARY)
    include_directories(${ADOL-C_HEADER})
    add_definitions(-DWITH_ADOL_C)
  endif(ADOL-C_LIBRARY AND ADOL-C_HEADER)
  message(STATUS ${ADOL-C_LIBRARY})
  message(STATUS ${ADOL-C_HEADER})
endif(ADOL-C_DIR)

if(WITH_ADOL-C)
  ExternalProject_Add(
    adolc
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL http://bitbucket.org/likask/mofem-joseph/downloads/ADOL-C-2.5.2.tgz
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/external/src/adolc/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --disable-shared --prefix=${PROJECT_BINARY_DIR}/external
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  find_library(
    COLPACK_LIBLARY 
    NAMES Colpack 
    PATHS ${ADOL-C_DIR}/lib ${ADOL-C_DIR}/lib64 /usr/local/lib)
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/lib64)
  execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${PROJECT_BINARY_DIR}/external/lib64/libadolc.a)
  add_dependencies(install_prerequisites adolc)
endif(WITH_ADOL-C)

if(WITH_ADOL-C)
  if(NOT ADOL-C_LIBRARY)
    set(ADOL-C_DIR 
      ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to adol-c" FORCE)
    set(ADOL-C_HEADER ${ADOL-C_DIR}/include)
    set(ADOL-C_LIBRARY 
      ${ADOL-C_DIR}/lib64/libadolc.a CACHE FILEPATH "adol-c lib" FORCE)
    include_directories(${PROJECT_BINARY_DIR}/external/include)
    add_definitions(-DWITH_ADOL_C)
    if(COLPACK_LIBLARY)
      set(ADOL-C_LIBRARY ${ADOL-C_LIBRARY} ${COLPACK_LIBLARY})
    endif(COLPACK_LIBLARY)
    message(STATUS ${ADOL-C_LIBRARY})
  endif(NOT ADOL-C_LIBRARY)
  message(STATUS ${ADOL-C_LIBRARY})
  message(STATUS ${ADOL-C_HEADER})
endif(WITH_ADOL-C)
