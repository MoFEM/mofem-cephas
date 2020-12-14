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

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

if(MOFEM_BUILD_TESTS)
  file(GLOB_RECURSE FTENSOR_CONFORMANCE_TEST_FILES
    ${PROJECT_SOURCE_DIR}/src/ftensor/tests/conformance/test_*.cpp
  )
  file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/All.cpp_tmp)
  foreach(LINE ${FTENSOR_CONFORMANCE_TEST_FILES})
    set(TO_INCLUDE "#include \"${LINE}\"\n")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/All.cpp_tmp ${TO_INCLUDE})
  endforeach(LINE ${TENSOR_CONFORMANCE_TEST_FILES})
  add_custom_target(
    all_cpp_generate_for_ftensor
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/All.cpp_tmp
    ${CMAKE_CURRENT_BINARY_DIR}/All.cpp
    COMMENT
    "Copy gnerated ftensor All.cpp file"
  )
  set_source_files_properties(All.cpp PROPERTIES GENERATED TRUE)
  add_executable(
    ftensor_conformance_test
    ${PROJECT_SOURCE_DIR}/src/ftensor/tests/conformance/main.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/All.cpp
  )
  add_dependencies(ftensor_conformance_test all_cpp_generate_for_ftensor)
  add_dependencies(ftensor_conformance_test install_prerequisites)
  add_test(
    ftensor_conformance_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/ftensor_conformance_test
  )
  add_subdirectory(${PROJECT_SOURCE_DIR}/src/ftensor/tests/speed)
endif(MOFEM_BUILD_TESTS)

if(ADOL-C_LIBRARY)
  if(MOFEM_BUILD_TESTS)
    find_path(TAPLESS_ADOL-C adtl.h HINTS ${ADOL-C_DIR}/include/adolc)
    # message(STATUS ${TAPLESS_ADOL-C})
    # Older versions of ADOL-C have don't have tapeless functionality
    if(TAPLESS_ADOL-C)
      add_executable(
        ftensor_and_adolc_tapeless_test
        ${PROJECT_SOURCE_DIR}/src/ftensor/tests/ftensor_and_adolc_tapeless.cpp
      )
      add_dependencies(ftensor_and_adolc_tapeless_test install_prerequisites)
      target_link_libraries(ftensor_and_adolc_tapeless_test ${ADOL-C_LIBRARY})
      add_test(
        ftensor_adloc_tapeless_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
        ${CMAKE_CURRENT_BINARY_DIR}/ftensor_and_adolc_tapeless_test
      )
    endif(TAPLESS_ADOL-C)
    add_executable(
      ftensor_and_adolc_test
      ${PROJECT_SOURCE_DIR}/src/ftensor/tests/ftensor_and_adolc.cpp
    )
    add_dependencies(ftensor_and_adolc_test install_prerequisites)
    target_link_libraries(ftensor_and_adolc_test ${ADOL-C_LIBRARY})
    add_test(
      ftensor_adloc_test ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
      ${CMAKE_CURRENT_BINARY_DIR}/ftensor_and_adolc_test
    )
  endif(MOFEM_BUILD_TESTS)
endif(ADOL-C_LIBRARY)

function(link_header_for_ftensor source destination)
  file(
    GLOB HEADER_FILES
    RELATIVE ${PROJECT_SOURCE_DIR}/src/${source}
    ${PROJECT_SOURCE_DIR}/src/${source}/*.h*)
  foreach(HF ${HEADER_FILES})
    install(CODE
      "
      EXECUTE_PROCESS(
      COMMAND 
      ln -sf 
      ${PROJECT_SOURCE_DIR}/src/${source}/${HF} 
      ${CMAKE_INSTALL_PREFIX}/include/${destination}/${HF}
      )
      MESSAGE(\"-- Linking ${CMAKE_INSTALL_PREFIX}/include/${destination}/${HF}\")
      "
    )
  endforeach(HF)
endfunction(link_header_for_ftensor)

# Install header files
if(STAND_ALLONE_USERS_MODULES)
  install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}/include/ftensor)
  install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/src/ftensor/src/
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}/include/ftensor
    FILES_MATCHING PATTERN "*.hpp")  
else(STAND_ALLONE_USERS_MODULES)
  install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ftensor)
  install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ftensor/FTensor)
  link_header_for_ftensor("ftensor/src" "ftensor")
  link_header_for_ftensor("ftensor/src/FTensor" "ftensor/FTensor")
  # Install FTensor directories
  file(
    GLOB DIRS
    RELATIVE ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor
    ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor/*)
  foreach(DIR ${DIRS})
    if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor/${DIR})
      install(DIRECTORY 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ftensor/FTensor/${DIR})
      file(
        GLOB HEADER_FILES
        RELATIVE ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor/${DIR}
        ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor/${DIR}/*.h*)
      link_header_for_ftensor(
        "ftensor/src/FTensor/${DIR}" "ftensor/FTensor/${DIR}")  
    endif(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/src/ftensor/src/FTensor/${DIR})
  endforeach(DIR)
endif(STAND_ALLONE_USERS_MODULES)

add_library(
  mofem_matrix_function
  src/impl/MatrixFunction.cpp)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_matrix_function Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_matrix_function DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)