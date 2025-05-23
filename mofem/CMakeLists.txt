cmake_minimum_required (VERSION 3.1.3)
#cmake_policy(SET CMP0144 NEW)

# Serch includes directories before system directories
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

# Petsc (parse config files)
include(cmake/FindPETSC.cmake)
# MoAB (parse config files)
include(cmake/FindMOAB.cmake)

# Test compilers like it is set by PETSc
set(CMAKE_C_COMPILER ${PETSCVAR_CC})
string(REGEX MATCH  "^.*mpic.?.?" PETSCVAR_CXX_ "${PETSCVAR_CXX}")
set(CMAKE_CXX_COMPILER ${PETSCVAR_CXX_})
if(PETSCVAR_FC)
  set(CMAKE_Fortran_COMPILER ${PETSCVAR_FC})
endif(PETSCVAR_FC)

project(MoFEM 
  VERSION 0.14.1
  DESCRIPTION "Finite element framework" 
  LANGUAGES C CXX)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
endif()


if(PETSCVAR_CUDA_INCLUDE)
  message(STATUS "PETSc is compiled with CUDA")
  set(PETSC_HAVE_CUDA TRUE)
  add_definitions(-DPETSC_HAVE_CUDA=1)
endif(PETSCVAR_CUDA_INCLUDE)

# Get mpirun path bin directory
string(REGEX REPLACE "mpicc" "" "MPI_BIN_PATH" ${CMAKE_C_COMPILER})
find_program(MPI_RUN mpirun HINTS ${MPI_BIN_PATH} /usr/bin PATH ${MPI_BIN_PATH})
include(cmake/ResolveCompilerPaths.cmake)


# Set fallback version. The fallback version is a version of MoFEM which is used
# as the default version if the version from the git tags can not be extracted. For example, 
# that is key when someone installs MoFEM from tarball.
set(MoFEM_FALLBACK_VERSION "${CMAKE_PROJECT_VERSION}" CACHE STRING "Fallback version" FORCE)
message(STATUS "Fallback version v${MoFEM_FALLBACK_VERSION}")

# Options
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(STAND_ALLONE_USERS_MODULES 
  "If if ON copy files, otherwise link" OFF)
option(WITHCOVERAGE "Add gcc coverage compile tags" OFF)
option(MOFEM_BUILD_TESTS "If is ON enable testing" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# If Fortran is set by PETSc or the user
if(CMAKE_Fortran_COMPILER)
  # Optional language
  enable_language(Fortran)
endif(CMAKE_Fortran_COMPILER)

# Add options unique to compilation with Spack
include(cmake/WithSpack.cmake)

# http://cellperformance.beyond3d.com/articles/2006/05/demystifying-the-restrict-keyword.html
# The restrict keyword can be considered an extension of the strict aliasing rule.
# It allows the programmer to declare that pointers that share the same type (or
# were otherwise validly created) do not alias each other. By using restrict the
# programmer can declare that any loads and stores through the qualified pointer
# (or through another pointer copied either directly or indirectly from the
# restricted pointer) are the only loads and stores to the same address during the
# lifetime of the pointer. In other words, the pointer is not aliased by any
# pointers other than its own copies.
# Restrict is a "no data hazards will be generated" contract between the
# programmer and the compiler. The compiler relies on this information to make
# optimizations. If the data is aliased, the results are undefined and a
# programmer should not expect the compiler to output a warning. The compiler
# assumes the programmer is not lying.

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
add_definitions("-D_LIBCPP_DISABLE_DEPRECATION_WARNINGS")

# Installing external packages It is used to install prerequisites, ADOL-C,
# TetGen, MED can be installed like that
include(ExternalProject)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/include
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/lib
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/lib64
  COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/external/bin)
add_custom_target(
  install_prerequisites
  # COMMAND ${CMAKE_COMMAND} -E touch_nocreate ${PROJECT_BINARY_DIR}/CMakeCache.txt
  COMMENT "Install prerequisites compiled with MoFEM")

include(FindPkgConfig)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)
include(FindPkgConfig)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)
set(
  CMAKE_PREFIX_PATH
  ${BLAS_DIR}
  ${CBLAS_DIR}
  ${PETSC_DIR}
  ${PETSC_DIR}/${PETSC_ARCH} 
  ${SLEPC_DIR}
  ${SLEPC_DIR}/${PETSC_ARCH} 
  ${CMAKE_PREFIX_PATH})  

#Find PETSc
pkg_check_modules(PETSc REQUIRED PETSc)
message(STATUS "PETSc version " ${PETSc_VERSION})
message(STATUS "PETSc lib dirs " ${PETSc_LIBRARY_DIRS})
message(STATUS "PETSc includ " ${PETSc_INCLUDE_DIRS})

# PETSC tools
pkg_check_modules(SLEPC slepc)
message(STATUS "SLEPc version " ${SLEPC_VERSION})
message(STATUS "SLEPc lib dirs " ${SLEPC_LIBRARY_DIRS})
message(STATUS "SLEPc include " ${SLEPC_INCLUDE_DIRS})

if(SLEPC_FOUND)

  # For some unknown reason pkg-config returning empty SLEPC_LINK_LIBRARIES,
  # following config force to find SLEPc lib 
  if(NOT SLEPC_LINK_LIBRARIES)
    find_library(
      SLEPC_LINK_LIBRARIES 
      NAMES slepc 
      PATHS 
      ${SLEPC_LIBRARY_DIRS}
      NO_DEFAULT_PATH)
  endif(NOT SLEPC_LINK_LIBRARIES)

  if(SLEPC_LINK_LIBRARIES)
    add_definitions(-DWITH_SLEPC)
    include_directories(${SLEPC_INCLUDE_DIRS})
  endif(SLEPC_LINK_LIBRARIES)
  
endif(SLEPC_FOUND)

# ccache (https://ccache.dev) 
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  message(STATUS "CCache found:" ${CCACHE_PROGRAM} )
endif()

# SIGMA tools
include(cmake/FindCGM.cmake)
include(cmake/FindMESHKIT.cmake)

# Find boost
include(cmake/FindBoost.cmake)

# Python
message(STATUS "PETSCVAR_PYTHON ${PETSCVAR_PYTHON}")
execute_process( COMMAND ${PETSCVAR_PYTHON} -c 
  "import sys; print(sys.exec_prefix)" 
  OUTPUT_VARIABLE Python_ROOT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process( COMMAND ${Python_ROOT_DIR}/bin/python3-config --embed --includes 
  OUTPUT_VARIABLE Python_INCLUDES_CONFIG)
execute_process( COMMAND ${Python_ROOT_DIR}/bin/python3-config --embed --ldflags 
  OUTPUT_VARIABLE Python_LIBRARY_CONFIG)

# Other tools
include(cmake/FindADOL-C.cmake)
include(cmake/FindTetGen.cmake)
include(cmake/FindTriangle.cmake)
include(cmake/FindMed.cmake)

# Other
include(cmake/ExportFile.cmake)

# Git revision and mofem version
find_package(Git)
include(cmake/GetGitRevisionSimple.cmake)
get_git_hash(${PROJECT_SOURCE_DIR} GIT_SHA1)
message(STATUS "MoFEM GIT_SHA1 ${GIT_SHA1}")
get_git_tag(${PROJECT_SOURCE_DIR} "v${MoFEM_FALLBACK_VERSION}" VERSION)
message(STATUS "MoFEM GIT_TAG ${VERSION}")
get_git_version(
  ${VERSION} MoFEM_VERSION_MAJOR MoFEM_VERSION_MINOR MoFEM_VERSION_BUILD)
set(CMAKE_PROJECT_VERSION 
  ${MoFEM_VERSION_MAJOR}.${MoFEM_VERSION_MINOR}.${MoFEM_VERSION_BUILD})

message(STATUS "MoFEM Version v${CMAKE_PROJECT_VERSION}")
add_definitions(
  -DMoFEM_VERSION_MAJOR=${MoFEM_VERSION_MAJOR} 
  -DMoFEM_VERSION_MINOR=${MoFEM_VERSION_MINOR} 
  -DMoFEM_VERSION_BUILD=${MoFEM_VERSION_BUILD})
add_definitions(-DGIT_SHA1_NAME="${GIT_SHA1}")

# Add coverage
if(WITHCOVERAGE)
  set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -fprofile-arcs -ftest-coverage")
  add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
  set(CMAKE_EXE_LINKER_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")
endif(WITHCOVERAGE)

resolve_includes(PROJECT_INCLUDE_DIRECTORIES 
  "${PETSCVAR_PACKAGES_INCLUDES}
  ${MOAB_INCLUDES0} ${MOAB_INCLUDES1} ${MOAB_CPPFLAGS}
  ${Python_INCLUDES_CONFIG}")

# Check CBlas
find_path(
  CBLAS_INCLUDE_DIR 
  cblas.h 
  HINTS 
  ${BLAS_DIR}/include
  ${CBLAS_DIR}/include 
  ${PETSC_DIR}/${PETSC_ARCH}/include
  /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers)
message(STATUS "CBLAS include found: " ${CBLAS_INCLUDE_DIR})
if(CBLAS_INCLUDE_DIR)
  include_directories(${CBLAS_INCLUDE_DIR})
endif(CBLAS_INCLUDE_DIR)

include_directories(
  ${PROJECT_INCLUDE_DIRECTORIES}
  ${MOAB_INCLUDE_DIR}
  ${PETSC_DIR}/include
  ${PETSC_DIR}/${PETSC_ARCH}/include
  ${PETSCVAR_DIR}/include
  ${BOOST_INCLUDE_DIR})

# Check Metis install with PETSc
if(PETSCVAR_METIS_INCLUDE)
  resolve_includes(METIS_INCLUDE_DIRECTORIES
    "${PETSCVAR_METIS_INCLUDE}")
  include_directories(${METIS_INCLUDE_DIRECTORIES})
  add_definitions(-DMETIS)
endif(PETSCVAR_METIS_INCLUDE)

# Check ParMeris install with PETSc
if(PETSCVAR_PARMETIS_INCLUDE)
  resolve_includes(PARMETIS_INCLUDE_DIRECTORIES
    "${PETSCVAR_PARMETIS_INCLUDE}")
  include_directories(${PARMETIS_INCLUDE_DIRECTORIES})
  add_definitions(-DPARMETIS)
endif(PETSCVAR_PARMETIS_INCLUDE)

# Find jupytext
find_file(JUPYTEXT jupytext)
message(STATUS ${JUPYTEXT})

# System project libraries
resolve_libraries(SYSTEM_PROJECT_LIBS "
# MoAB
  ${MESHKIT_LIBS_LINK}
  ${MOAB_LIBS_LINK}
# PETSc
  ${SLEPC_LINK_LIBRARIES}
  -L${PETSC_DIR}/${PETSC_ARCH}/lib
  ${PETSCVAR_PETSC_WITH_EXTERNAL_LIB}
# Python
  ${Python_LIBRARY_CONFIG}
# Boost
  ${Boost_LIBRARIES}
  ${Boost_FILESYSTEM_LIBRARY_RELEASE}
  ${Boost_LOG_LIBRARY_RELEASE}
  ${Boost_LOG_SETUP_LIBRARY_RELEASE}
  ${Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE}
  ${Boost_SYSTEM_LIBRARY_RELEASE}
  ${Boost_THREAD_LIBRARY_RELEASE}
# Other
  ${MPI_F90_LIB}
  ${MPI_F77_LIB}")

# Set project libs compiled with the MoFEM as external libs
set(OPTIONAL_PROJECT_LIBS
  ${TETGEN_LIBRARY}
  ${MED_LIBRARY}
  ${ADOL-C_LIBRARY})
set(PROJECT_LIBS 
  ${OPTIONAL_PROJECT_LIBS} ${SYSTEM_PROJECT_LIBS} ${PETSCVAR_BLASLAPACK_LIB})
message(STATUS ${PROJECT_LIBS})

# Add executable library from a subdirectory
if(MOFEM_BUILD_TESTS)
  enable_testing()
  include(CTest)
endif(MOFEM_BUILD_TESTS)

# Precompiled headers
include(cmake/PrecompiledHeaders.cmake)

add_subdirectory(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(${PROJECT_SOURCE_DIR}/cmake)
add_subdirectory(${PROJECT_SOURCE_DIR}/scripts)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/tools)
add_subdirectory(${PROJECT_SOURCE_DIR}/tutorials)
add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
if(MOFEM_BUILD_TESTS)
  add_subdirectory(${PROJECT_SOURCE_DIR}/atom_tests)
endif(MOFEM_BUILD_TESTS)

# Generate config files
configure_file(
  ${PROJECT_SOURCE_DIR}/MoFEMConfig-version.cmake.in
  ${PROJECT_BINARY_DIR}/MoFEMConfig-version.cmake)
configure_file(
  ${PROJECT_SOURCE_DIR}/MoFEMConfig.cmake.in
  ${PROJECT_BINARY_DIR}/MoFEMConfig.cmake)

# Install project files
install(
  FILES ${PROJECT_BINARY_DIR}/MoFEMConfig.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX})
install(
  FILES ${PROJECT_BINARY_DIR}/MoFEMConfig-version.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX})

# Install cmake files
install(
  DIRECTORY
  ${PROJECT_SOURCE_DIR}/cmake/
  DESTINATION
  ${CMAKE_INSTALL_PREFIX}/lib/cmake/mofem
  FILES_MATCHING
  PATTERN "*.cmake"
  PATTERN "*.cmake.in")

# Install external directory
install(
  DIRECTORY
  ${PROJECT_BINARY_DIR}/external
  DESTINATION
  ${CMAKE_INSTALL_PREFIX})

# copy scripts to html directory created by doxygen
if(DOXYGEN_FOUND)
  add_custom_target(doxygen_copy_scripts
    ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/scripts ${PROJECT_BINARY_DIR}/html/scripts
  )
  add_dependencies(doc doxygen_copy_scripts)
endif(DOXYGEN_FOUND)

add_definitions(-Wfatal-errors)