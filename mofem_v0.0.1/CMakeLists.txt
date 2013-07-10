# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
# It can be freely used for educational and research purposes 
# by other institutions. If you use this softwre pleas cite my work. 
#
# MoFEM is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# MoFEM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MoFEM. If not, see <http://www.gnu.org/licenses/>


cmake_minimum_required (VERSION 2.6)

include(cmake/FindPETSC.cmake)
include(cmake/FindMOAB.cmake)
include(cmake/ExportFile.cmake)

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

include_directories(
  ${123_configurational_forces_SOURCE_DIR}
  ${PROJECT_INCLUDE_DIRECTORIES}
  ${PETSC_DIR}/${PETSC_ARCH}/include
  ${PETSCVAR_DIR}/include
)

#SPECIFIC FROM NOW
include_directories("${PROJECT_SOURCE_DIR}/do_not_blink")
include_directories("${PROJECT_SOURCE_DIR}/student_body")
include_directories("${PROJECT_SOURCE_DIR}/complex_for_lazy")
include_directories("${PROJECT_SOURCE_DIR}/third_party")

add_library(mofem_Field 
  third_party/tetcircumcenter.c
  do_not_blink/Core_dataStructures.cpp
  do_not_blink/Core_moabMeshMBTETRefine.cpp
  do_not_blink/moabField_Core.cpp
  complex_for_lazy/complex_for_lazy.c
  complex_for_lazy/moabFEMethod_ComplexForLazy.cpp
  student_body/moabFEMethod_LowLevelStudent.cpp
  student_body/moabFEMethod_UpLevelStudent.cpp
  moabField.cpp
  FEM.c
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

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/INSTALL
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#add executable library form subdirectory
add_subdirectory("${PROJECT_SOURCE_DIR}/third_party")
add_subdirectory("${PROJECT_SOURCE_DIR}/dirty_testing")
add_subdirectory("${PROJECT_SOURCE_DIR}/examples")


#example how to add new executable
#substitute "Program_Name" by name of new program
#add_executable(Program_Name Program_Name.cpp)
#target_link_libraries(Program_Name
#  mofem_Field ${PROJECT_LIBS} ${PROJECT_LIBS})



