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

cmake_minimum_required (VERSION 2.8.7)

#petsc and moab configuration are used to confiugre mofem
include(cmake/FindPETSC.cmake)
include(cmake/FindMOAB.cmake)

message(${PETSCVAR_CXX})

set(CMAKE_C_COMPILER ${PETSCVAR_CC} )
string(REGEX MATCH  "^.*mpicxx" PETSCVAR_CXX_ "${PETSCVAR_CXX}")
set(CMAKE_CXX_COMPILER ${PETSCVAR_CXX_})

#get mpipath bin directory
string(REGEX REPLACE "mpicc" "" "MPI_BIN_PATH" ${CMAKE_C_COMPILER})
find_program(MPI_RUN mpirun HINTS ${MPI_BIN_PATH} /usr/bin PATH ${MPI_BIN_PATH})
include(cmake/ResolveCompilerPaths.cmake)

project(MoFEM C Fortran CXX)

set(MoFEM_VERSION_MAJOR 0)
set(MoFEM_VERSION_MINOR 3)
set(MoFEM_VERSION_BUILD 32)

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES(
  "int main() { double * __restrict__ x;}\n"
  SYMM_RESTRICT
)
if(SYMM_RESTRICT)
  set(DEFINE_RESTRICT "-Drestrict=__restrict__")
else(SYMM_RESTRICT)
  set(DEFINE_RESTRICT "-Drestrict=")
  add_definitions("-Drestrict=")
endif(SYMM_RESTRICT)
add_definitions(${DEFINE_RESTRICT})

# Installing extermal packages
include(ExternalProject)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/include
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/lib
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/lib64
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/bin
)
add_custom_target(
  install_prerequisites
  COMMAND ${CMAKE_COMMAND} -E touch_nocreate ${PROJECT_BINARY_DIR}/CMakeCache.txt
  COMMENT "Install prerequisites compiled with MoFEM"
)

#SIGMA tools
include(cmake/FindiMesh.cmake)
include(cmake/FindCGM.cmake)
include(cmake/FindLASSO.cmake)
include(cmake/FindMESHKIT.cmake)

#PETSC tools
include(cmake/FindTAO.cmake)
include(cmake/FindSLEPC.cmake)

#Other tools
include(cmake/FindADOL-C.cmake)
include(cmake/FindTetGen.cmake)
include(cmake/FindNetGen.cmake)
include(cmake/FindTriangle.cmake)

#CGAL tools
#find_package(CGAL HINTS ${CGAL_DIR}/lib/CGAL)
#if(CGAL_CONFIG_LOADED)
#  include_directories(${CGAL_INCLUDE_DIRS})
#  link_directories(${CGAL_LIBRARIES_DIR})
#  add_definitions( -DWITH_CGAL )
#  message(STATUS ${CGAL_Core_LIBRARY})
#endif(CGAL_CONFIG_LOADED)

#other
include(cmake/ExportFile.cmake)

#git revision
include(cmake/GetGitRevisionDescription.cmake)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
#add_definitions("-DGIT_SHA1=${GIT_SHA1}")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/include/config.h.in"
  "${PROJECT_SOURCE_DIR}/include/config.h"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/include/version.h.in"
  "${PROJECT_SOURCE_DIR}/include/version.h"
)
include_directories("${PROJECT_SOURCE_DIR}/include")

#add coverage
if(WITHCOVERAGE)
  set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
  add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")
endif(WITHCOVERAGE)

#add check for verbose testing
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

resolve_includes(PROJECT_INCLUDE_DIRECTORIES
  "${PETSCVAR_PACKAGES_INCLUDES} ${MOAB_INCLUDES0} ${MOAB_INCLUDES1}"
)

include_directories(
  ${PROJECT_INCLUDE_DIRECTORIES}
  ${PETSC_DIR}/${PETSC_ARCH}/include
  ${PETSCVAR_DIR}/include
)

if(TAO_LIBRARY)
  include_directories(${TAO_DIR}/include)
  add_definitions(-DTAO)
endif(TAO_LIBRARY)

resolve_libraries(PROJECT_LIBS "
  ${SLEPC_LIBRARY}
  ${ADOL-C_LIBRARY}
  ${TETGEN_LIBRARY}
  ${NETGEN_LIBRARY}
  ${MOAB_LIBS_LINK}
  -L${PETSC_DIR}/${PETSC_ARCH}/lib
  ${PETSCVAR_PETSC_WITH_EXTERNAL_LIB}
  ${TAO_LIBRARY}
  ${MPI_F90_LIB}
  ${MPI_F77_LIB}"
)

message(STATUS ${PROJECT_LIBS})
file(COPY
  ${PROJECT_SOURCE_DIR}/INSTALL
  DESTINATION ${PROJECT_BINARY_DIR}
)

#add executable library form subdirectory
enable_testing()
include(CTest)

#add generation of examplaes
configure_file(
  ${PROJECT_SOURCE_DIR}/MoFEMConfig.cmake.in
  ${PROJECT_BINARY_DIR}/MoFEMConfig.cmake
)

configure_file(
  ${PROJECT_SOURCE_DIR}/MoFEMConfig-version.cmake.in
  ${PROJECT_BINARY_DIR}/MoFEMConfig-version.cmake
)
install(FILES ${PROJECT_BINARY_DIR}/MoFEMConfig.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${PROJECT_BINARY_DIR}/MoFEMConfig-version.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})

add_subdirectory(${PROJECT_SOURCE_DIR}/scripts)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party)
add_subdirectory(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/atom_tests)
add_subdirectory(${PROJECT_SOURCE_DIR}/users_modules.in)
add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
