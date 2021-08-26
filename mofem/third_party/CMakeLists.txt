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

#include subdirectory
include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)

add_library(mofem_third_party
  IntegrationRules.cpp
  gm_rule.c
  tetcircumcenter.c
  tricircumcenter.c
  phg-quadrule/quad.c)
target_link_libraries(mofem_third_party PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_third_party PROPERTIES VERSION ${PROJECT_VERSION})
install(TARGETS mofem_third_party DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# FIXME: Install directory, and other stuff, like kiss lib,
# that should not be there.  It is for back compatibility.
install(
  DIRECTORY
  ${CMAKE_SOURCE_DIR}/third_party
  DESTINATION
  ${CMAKE_INSTALL_PREFIX})
