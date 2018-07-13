# - Try to find MED

if(NOT MED_DIR)
  set(MED_DIR $ENV{MED_DIR})
endif(NOT MED_DIR)

if(MED_DIR)
  find_library(MED_LIBRARY NAMES med medC PATHS ${MED_DIR}/lib)
  find_path(MED_HEADER NAMES med.h PATHS ${MED_DIR}/include)
  message(STATUS ${MED_LIBRARY})
  message(STATUS ${MED_HEADER})
  if(MED_LIBRARY AND MED_HEADER)
    include_directories(${MED_HEADER})
    add_definitions(-DWITH_MED)
  endif(MED_LIBRARY AND MED_HEADER)
endif(MED_DIR)

if(WITH_MED)
  ExternalProject_Add(
    med
    PREFIX ${PROJECT_BINARY_DIR}/external/
    URL http://files.salome-platform.org/Salome/other/med-3.2.0.tar.gz
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/external/src/med/configure CPPFLAGS=${MOAB_CPPFLAGS} LDFLAGS=${MOAB_LDFLAGS} CC=${CMAKE_C_COMPILER} CFLAGS=${MOAB_CFLAGS} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${MOAB_CXXFLAGS} --disable-shared --disable-python --disable-fortran --prefix=${PROJECT_BINARY_DIR}/external
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  add_dependencies(install_prerequisites med)
endif(WITH_MED)

if(WITH_MED)
  if(NOT MED_LIBRARY)
    set(MED_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to med" FORCE)
    set(MED_HEADER ${MED_DIR}/include)
    set(MED_LIBRARY ${MED_DIR}/lib/libmed.a CACHE FILEPATH "med lib" FORCE)
    include_directories(${MED_HEADER})
    add_definitions(-DWITH_MED)
    message(STATUS ${MED_LIBRARY})
    message(STATUS ${MED_HEADER})
  endif(NOT MED_LIBRARY)
endif(WITH_MED)
