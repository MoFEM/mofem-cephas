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

cmake_minimum_required (VERSION 3.1.3)

project(MoFEM-UserModules C CXX)

#serch inlude directories before system directories
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

include(cmake/GetGitRevisionDescription.cmake)
git_describe(GIT_VERSION "--tags")
if(GIT_VERSION)
  string(
    REGEX 
    REPLACE "^v([0-9]+)\\..*" "\\1" GIT_MoFEM_VERSION_MAJOR "${GIT_VERSION}"
  )
  string(
    REGEX 
    REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" GIT_MoFEM_VERSION_MINOR "${GIT_VERSION}"
  )
  string(
    REGEX 
    REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" GIT_MoFEM_VERSION_BUILD "${GIT_VERSION}"
  )
  find_package(
    MoFEM 
    ${GIT_MoFEM_VERSION_MAJOR}.${GIT_MoFEM_VERSION_MINOR}.${GIT_MoFEM_VERSION_MINOR}
    REQUIRED HINTS ${MOFEM_DIR} ../ ../../ $ENV{MOFEM_DIR}
  )  
else(GIT_VERSION)
  find_package(
    MoFEM REQUIRED HINTS ${MOFEM_DIR} ../ ../../ $ENV{MOFEM_DIR}
  )  
endif(GIT_VERSION)

message(STATUS "MoFEM version ${MoFEM_VERSION}")
add_definitions(${MoFEM_DEFINE_RESTRICT})

macro(CHECK_MOFEM_VERSION VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

  if(${MoFEM_VERSION} VERSION_LESS ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
    message(
      FATAL_ERROR
      "(Update MoFEM) Wrong MoFEM version: "
      "${MoFEM_VERSION} < ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}"
    )
  endif(${MoFEM_VERSION} VERSION_LESS ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

endmacro(CHECK_MOFEM_VERSION)

find_package(Git)
# find_package(PythonInterp)
# find_package(PythonLibs)
# if(PYTHON_INCLUDE_DIRS)
#   include_directories ( ${PYTHON_INCLUDE_DIRS} )
# endif(PYTHON_INCLUDE_DIRS)

find_program(
  MBCONVERT NAMES mbconvert
  PATH ${MoFEM_PETSC_DIR}/${MoFEM_PETSC_ARCH}/bin
)

include(ExternalProject)
include(${MoFEM_PROJECT_SOURCE_DIR}/cmake/ExportFile.cmake)

set(CMAKE_C_COMPILER ${MoFEM_CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MoFEM_CMAKE_CXX_COMPILER})
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
  ${MoFEM_PROJECT_INCLUDE_DIRECTORIES}
  ${MoFEM_PETSC_DIR}/include
  ${MoFEM_PETSC_DIR}/${PETSC_ARCH}/include
  ${MoFEM_PETSC_DIR}/${MoFEM_PETSC_ARCH}/include
  ${MoFEM_PETSCVAR_DIR}/include
)

add_definitions(${MOAB_DEFINITIONS})

include_directories(${MoFEM_PROJECT_SOURCE_DIR}/include)
include_directories(${MoFEM_INSTALL_DIR}/include)

include_directories(${MoFEM_PROJECT_SOURCE_DIR}/third_party)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/third_party/cblas)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/petsc)
include_directories(${MoFEM_PROJECT_SOURCE_DIR}/src/finite_elements)
# boost dir
if(MoFEM_BOOST_INCLUDE_DIR)
  include_directories(${MoFEM_BOOST_INCLUDE_DIR})
endif(MoFEM_BOOST_INCLUDE_DIR)
# cblas 
if(MoFEM_CBLAS_INCLUDE_DIR) 
  include_directories(${MoFEM_CBLAS_INCLUDE_DIR})
endif(MoFEM_CBLAS_INCLUDE_DIR)

link_directories(${MoFEM_INSTALL_DIR}/lib)
link_directories(${PROJECT_BINARY_DIR})

#add executable library form subdirectory
enable_testing()
include(CTest)

#add coverage
if(WITHCOVERAGE)
  set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
  add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")
endif(WITHCOVERAGE)

# source path
set(UM_SOURCE_DIR ${MoFEM_PROJECT_SOURCE_DIR}/users_modules)
if(STAND_ALLONE_USERS_MODULES)
  set(UM_SOURCE_DIR ${PROJECT_SOURCE_DIR})
endif(STAND_ALLONE_USERS_MODULES)

# third party libraries
if(TETGEN_LIBRARY)
  include_directories("${TETGEN_DIR}/include")
  add_definitions( -DWITH_TETGEN )
  message(STATUS "WITH TETEGEN")
endif(TETGEN_LIBRARY)

if(NETGEN_LIBRARY)
  include_directories("${NETGEN_DIR}/include")
  include_directories("${NETGEN_SRC_DIR}/libsrc/include")
  add_definitions( -DWITH_NETGEN )
  message(STATUS "WITH NETEGEN")
endif(NETGEN_LIBRARY)

if(CGAL_DIR)
  find_package(CGAL HINTS ${CGAL_DIR}/lib/CGAL)
  if(CGAL_CONFIG_LOADED)
    include_directories(${CGAL_INCLUDE_DIRS})
    link_directories(${CGAL_LIBRARIES_DIR})
    add_definitions( -DWITH_CGAL )
    message(STATUS "WITH CGAL")
  endif(CGAL_CONFIG_LOADED)
endif(CGAL_DIR)

if(ADOL-C_LIBRARY)
  include_directories("${ADOL-C_DIR}/include")
  add_definitions( -DWITH_ADOL_C )
  message(STATUS "WITH ADLOC-C")
endif(ADOL-C_LIBRARY)

if(MED_LIBRARY)
  include_directories("${MED_DIR}/include")
  add_definitions(-DWITH_MED )
  message(STATUS "WITH MED")
endif(MED_LIBRARY)

if(SLEPC_LIBRARY)
  include_directories("${SLEPC_DIR}/${PETSC_ARCH}/include")
  include_directories("${SLEPC_DIR}/include")
  add_definitions( -DWITH_SLEPC_DIR )
  message(STATUS "WITH SLEPC")
endif(SLEPC_LIBRARY)

# Add modules
include(ModulesList.cmake)

# Complle tools
add_subdirectory(${PROJECT_SOURCE_DIR}/tools)


