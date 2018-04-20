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

file(GLOB_RECURSE FTENSOR_ONE_OVER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/one_over_*.cpp
)
foreach(ONE_OVER_SRC_FILE ${FTENSOR_ONE_OVER_FILES})
  string(REGEX REPLACE ".cpp" "" ONE_OVER ${ONE_OVER_SRC_FILE})
  string(REGEX REPLACE ".*/" "" ONE_OVER ${ONE_OVER})
  # message(STATUS ${ONE_OVER})
  add_executable(${ONE_OVER} ${ONE_OVER_SRC_FILE})
  add_dependencies(${ONE_OVER} install_prerequisites)
endforeach(ONE_OVER_SRC_FILE ${FTENSOR_ONE_OVER_FILES})

# Copy script doing tests
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/one_over_script
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
