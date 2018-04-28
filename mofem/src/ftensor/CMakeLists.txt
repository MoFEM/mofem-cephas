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

include_directories("${PROJECT_SOURCE_DIR}/src/ftensor/src")

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

if(ADOL-C_LIBRARY)

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

endif(ADOL-C_LIBRARY)
