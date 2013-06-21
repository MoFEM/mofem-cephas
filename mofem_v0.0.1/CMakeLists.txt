# Lukasz Kaczmarczyk
# likask@civil.gla.ac.uk

cmake_minimum_required (VERSION 2.6)

include(cmake/FindPETSC.cmake)
include(cmake/FindMOAB.cmake)

set(CMAKE_C_COMPILER ${PETSCVAR_CC} ) 
set(CMAKE_CXX_COMPILER ${PETSCVAR_CXX})

include(cmake/ResolveCompilerPaths.cmake)

project(123_configurational_forces C CXX)

#find_package(CGAL QUIET COMPONENTS  Core)
#if(NOT CGAL_FOUND)
#    message(FATAL_ERROR "CGAL NOT FOUND")
#endif()
#message( STATUS "CGAL_USE_FILE " ${CGAL_USE_FILE} )
#include( ${CGAL_USE_FILE} )

set(PETSC_FLAGS "${PETSCVAR_CC_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${PETSC_FLAGS} ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_MINSIZEREL "${PETSC_FLAGS} ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_MINSIZEREL}")
set(CMAKE_C_FLAGS_RELEASE "${PETSC_FLAGS} ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
#message(STATUS ${CMAKE_C_FLAGS_DEBUG})
#message(STATUS ${CMAKE_C_FLAGS_MINSIZEREL})
#message(STATUS ${CMAKE_C_FLAGS_RELEASE})

set(MOAB_FLAGS "${MOAB_CXXFLAGS} -Wno-deprecated-declarations")
string (REGEX REPLACE "-pedantic" "" "MOAB_FLAGS"  ${MOAB_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG "${MOAB_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${MOAB_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}")
set(CMAKE_CXX_FLAGS_RELEASE "${MOAB_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
#message(STATUS ${CMAKE_CXX_FLAGS_DEBUG})
#message(STATUS ${CMAKE_CXX_FLAGS_MINSIZEREL})
#message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})

find_library(MPI_F90_LIB mpi_f90 PATHS /usr/lib /usr/lib/openmpi/lib)
message(STATUS ${MPI_F90_LIB})
find_library(MPI_F77_LIB mpi_f77 PATHS /usr/lib /usr/lib/openmpi/lib)
message(STATUS ${MPI_F77_LIB})

resolve_includes(PROJECT_INCLUDE_DIRECTORIES 
  "${PETSCVAR_PACKAGES_INCLUDES} ${MOAB_INCLUDES0} ${MOAB_INCLUDES}")
#message(STATUS ${123_configurational_forces_SOURCE_DIR})
#message(STATUS ${PROJECT_INCLUDE_DIRECTORIES})
include_directories(
  ${123_configurational_forces_SOURCE_DIR}
  ${PROJECT_INCLUDE_DIRECTORIES}
  ${PETSC_DIR}/${PETSC_ARCH}/include
  ${PETSCVAR_DIR}/include
)

#SPECIFIC FROM NOW

add_library(mofem_Field 
  dataStructures.cpp
  moabField.cpp
  moabField_Core.cpp
  moabFEMethod_Core.cpp
  moabFEMethod_Student.cpp
  meshRefine.cpp
  fem.c
  H1.c
  L2.c
)
#add_library(cgal_interface 
#  cgal_interface.cpp
#)

resolve_libraries(PROJECT_LIBS "
  ${MOAB_LIBS_LINK}
  -L${PETSC_DIR}/${PETSC_ARCH}/lib 
  ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a
  ${PETSCVAR_PACKAGES_LIBS}
  ${PETSCVAR_PETSC_EXTERNAL_LIB_BASIC}
  ${MPI_F90_LIB}
  ${MPI_F77_LIB}
  ")
#message(STATUS ${PETSCVAR_BLASLAPACK_LIB})
message(STATUS ${PROJECT_LIBS})

add_executable(test_moabField test_moabField.cpp)
target_link_libraries(test_moabField mofem_Field ${PROJECT_LIBS})

add_executable(test_mesh_refinment test_mesh_refinment.cpp)
target_link_libraries(test_mesh_refinment mofem_Field ${PROJECT_LIBS})

add_executable(test_approx_function test_approx_function.cpp)
target_link_libraries(test_approx_function mofem_Field ${PROJECT_LIBS})

add_executable(test_elasticity test_elasticity.cpp)
target_link_libraries(test_elasticity mofem_Field ${PROJECT_LIBS})



#example how to add new executable
#substitute "Program_Name" by name of new program
#add_executable(Program_Name Program_Name.cpp)
#target_link_libraries(Program_Name
#  fem ${PROJECT_LIBS} ${PROJECT_LIBS})

#add_executable(test_all test_all.cpp)
#target_link_libraries(test_all cgal_interface fem  ${PROJECT_LIBS} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
#target_link_libraries(test_all fem  ${PROJECT_LIBS})




