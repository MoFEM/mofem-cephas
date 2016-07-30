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
add_executable(
  ftensor_conformance_test
  ${PROJECT_SOURCE_DIR}/src/ftensor/tests/conformance/main.cpp
  ${FTENSOR_CONFORMANCE_TEST_FILES}
)
add_test(
  ftensor_conformance_test ${CMAKE_CURRENT_BINARY_DIR}/ftensor_conformance_test
)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/ftensor/tests/speed)

if(ADOL-C_LIBRARY)

  add_executable(
    ftensor_and_adolc_tapeless_test
    ${PROJECT_SOURCE_DIR}/src/ftensor/tests/ftensor_and_adolc_tapeless.cpp
  )
  target_link_libraries(ftensor_and_adolc_tapeless_test ${ADOL-C_LIBRARY})
  add_test(
    ftensor_adloc_tapeless_test ${CMAKE_CURRENT_BINARY_DIR}/ftensor_and_adolc_tapeless_test
  )

  add_executable(
    ftensor_and_adolc_test
    ${PROJECT_SOURCE_DIR}/src/ftensor/tests/ftensor_and_adolc.cpp
  )
  target_link_libraries(ftensor_and_adolc_test ${ADOL-C_LIBRARY})
  add_test(
    ftensor_adloc_test ${CMAKE_CURRENT_BINARY_DIR}/ftensor_and_adolc_test
  )


endif(ADOL-C_LIBRARY)
