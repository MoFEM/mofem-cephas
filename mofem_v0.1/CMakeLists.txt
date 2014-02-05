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

cmake_minimum_required (VERSION 2.8)

include(cmake/FindPETSC.cmake)
include(cmake/FindMOAB.cmake)
include(cmake/ExportFile.cmake)
include(cmake/FindTAO.cmake)

message(${PETSCVAR_CXX})

set(CMAKE_C_COMPILER ${PETSCVAR_CC} ) 
string(REGEX MATCH  "^.*mpicxx" PETSCVAR_CXX_ "${PETSCVAR_CXX}")
message(${PETSCVAR_CXX_})
set(CMAKE_CXX_COMPILER ${PETSCVAR_CXX_})
#get mpipath bin directory
string(REGEX REPLACE "mpicc" "" "MPI_BIN_PATH" ${CMAKE_C_COMPILER})
find_program(MPI_RUN mpirun HINTS ${MPI_BIN_PATH} PATH ${MPI_BIN_PATH})
include(cmake/ResolveCompilerPaths.cmake)

project(MoFEM C Fortran CXX)

set(MoFEM_VERSION_MAJOR 0)
set(MoFEM_VERSION_MINOR 1)
set(MoFEM_VERSION_BUILD 1)

#git revision
include(cmake/GetGitRevisionDescription.cmake)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
#add_definitions("-DGIT_SHA1=${GIT_SHA1}")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
configure_file(
  "${PROJECT_SOURCE_DIR}/version.h.in"
  "${PROJECT_BINARY_DIR}/version.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

#add coverage
if(WITHCOVERAGE) 
  set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
  add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")
endif(WITHCOVERAGE) 

#add check for verbose testing
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

#find_package(CGAL QUIET COMPONENTS  Core)
#if(NOT CGAL_FOUND)
#    message(FATAL_ERROR "CGAL NOT FOUND")
#endif()
#message( STATUS "CGAL_USE_FILE " ${CGAL_USE_FILE} )
#include( ${CGAL_USE_FILE} )

resolve_includes(PROJECT_INCLUDE_DIRECTORIES 
  "${PETSCVAR_PACKAGES_INCLUDES} ${MOAB_INCLUDES0} ${MOAB_INCLUDES}"
  )

include_directories(
  ${123_configurational_forces_SOURCE_DIR}
  ${PROJECT_INCLUDE_DIRECTORIES}
  ${PETSC_DIR}/${PETSC_ARCH}/include
  ${PETSCVAR_DIR}/include
  )

if(TAO_LIBRARY) 
  include_directories(${TAO_DIR}/include)
  add_definitions(-DTAO)
endif(TAO_LIBRARY)

include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/do_not_blink")
include_directories("${PROJECT_SOURCE_DIR}/student_body")
include_directories("${PROJECT_SOURCE_DIR}/projections")
include_directories("${PROJECT_SOURCE_DIR}/complex_for_lazy")
include_directories("${PROJECT_SOURCE_DIR}/third_party")
include_directories("${PROJECT_SOURCE_DIR}/third_party/cblas")
link_directories("${PROJECT_BINARY_DIR}/third_party/cblas")

add_library(mofem_Field 
  FieldInterface.cpp
  SnesCtx.cpp
  third_party/tricircumcenter.c
  third_party/tetcircumcenter.c
  do_not_blink/CoreDataStructures.cpp
  do_not_blink/TETRefine.cpp
  do_not_blink/FieldCore.cpp
  FEM.c
  H1.c
  Hdiv.c
  L2.c
)

add_library(mofem_student_body
  student_body/FEMethod_LowLevelStudent.cpp
  student_body/FEMethod_UpLevelStudent.cpp
  student_body/DirihletBC.cpp
  student_body/petscShellMATs_ConstrainsByMarkAinsworth.cpp
  student_body/ArcLengthTools.cpp
)

add_library(mofem_complex_for_lazy
  complex_for_lazy/complex_for_lazy.c
  complex_for_lazy/FEMethod_ComplexForLazy.cpp
  complex_for_lazy/FEMethod_SurfaceConstrains.cpp
  complex_for_lazy/ConfigurationalFractureMechanics.cpp
)

add_subdirectory("${PROJECT_SOURCE_DIR}/scripts")
add_subdirectory("${PROJECT_SOURCE_DIR}/third_party")
add_subdirectory("${PROJECT_SOURCE_DIR}/third_party/cblas")

resolve_libraries(PROJECT_LIBS "
  ${MOAB_LIBS_LINK}
  -L${PETSC_DIR}/${PETSC_ARCH}/lib 
  ${PETSCVAR_PETSC_WITH_EXTERNAL_LIB}
  ${TAO_LIBRARY}
  ${MPI_F90_LIB}
  ${MPI_F77_LIB}")
#message(STATUS ${PETSCVAR_BLASLAPACK_LIB})
set(PROJECT_LIBS "mofem_cblas ${PROJECT_LIBS}")
message(STATUS ${PROJECT_LIBS})

file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/INSTALL
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

#add executable library form subdirectory
enable_testing()
include(CTest)
add_subdirectory("${PROJECT_SOURCE_DIR}/atom_tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/examples")
add_subdirectory("${PROJECT_SOURCE_DIR}/benchmarks")

install(DIRECTORY "${PROJECT_SOURCE_DIR}/benchmarks" DESTINATION mofem)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  )
set(DOXYGEN_LINE_COMMAND 
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
add_custom_target(doc
  ${DOXYGEN_LINE_COMMAND}
  COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)


