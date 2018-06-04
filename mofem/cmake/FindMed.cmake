# - Try to find MED

if(NOT MED_DIR)
  set(MED_DIR $ENV{MED_DIR})
endif(NOT MED_DIR)

if(MED_DIR)
  find_library(MED_LIBRARY NAMES med PATHS ${MED_DIR}/lib)
  message(STATUS ${MED_LIBRARY})
  if(MED_LIBRARY)
    include_directories(${MED_DIR}/include)
    add_definitions(-DWITH_MED)
  endif(MED_LIBRARY)
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
    include_directories(${PROJECT_BINARY_DIR}/external/include)
    set(MED_DIR ${PROJECT_BINARY_DIR}/external CACHE FILEPATH "path to med" FORCE)
    set(MED_LIBRARY ${MED_DIR}/lib/libmed.a CACHE FILEPATH "med lib" FORCE)
    add_definitions(-DWITH_MED)
    message(STATUS ${MED_LIBRARY})
  endif(NOT MED_LIBRARY)
endif(WITH_MED)
