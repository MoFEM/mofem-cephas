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

# List of sources for user_modules libaries.
set(UM_LIB_SOURCES "")

include(${UM_SOURCE_DIR}/basic_finite_elements/AddModule.cmake)
include(${UM_SOURCE_DIR}/analytical_dirichlet_boundary_conditions/AddModule.cmake)
include(${UM_SOURCE_DIR}/convective_mass_element/AddModule.cmake)
include(${UM_SOURCE_DIR}/obsolete/AddModule.cmake)

# Users modules library, common for all programs
add_library(users_modules ${UM_LIB_SOURCES})

# Atom test user modules
add_subdirectory(
  ${UM_SOURCE_DIR}/basic_finite_elements/atom_tests
  ${PROJECT_BINARY_DIR}/basic_finite_elements/atom_tests
)
add_subdirectory(
  ${UM_SOURCE_DIR}/analytical_dirichlet_boundary_conditions/atom_tests
  ${PROJECT_BINARY_DIR}/analytical_dirichlet_boundary_conditions/atom_tests
)
add_subdirectory(
  ${UM_SOURCE_DIR}/convective_mass_element/atom_tests
  ${PROJECT_BINARY_DIR}/convective_mass_element/atom_tests
)
add_subdirectory(
  ${UM_SOURCE_DIR}/gels/atom_tests
  ${PROJECT_BINARY_DIR}/gels/atom_tests
)
add_subdirectory(
  ${UM_SOURCE_DIR}/obsolete/atom_tests
  ${PROJECT_BINARY_DIR}/obsolete/atom_tests
)

file(
  GLOB_RECURSE INSTLLED_MODULES
  FOLLOW_SYMLINKS
  ?*/InstalledAddModule.cmake
)

foreach(LOOP_MODULE ${INSTLLED_MODULES})
  message(STATUS "Add module ... ${LOOP_MODULE}")
  include(${LOOP_MODULE})
endforeach(LOOP_MODULE)

# Users Programs
add_subdirectory(
  ${UM_SOURCE_DIR}/thermal
  ${PROJECT_BINARY_DIR}/thermal
)
add_subdirectory(
  ${UM_SOURCE_DIR}/ultraweak
  ${PROJECT_BINARY_DIR}/ultraweak
)
add_subdirectory(
  ${UM_SOURCE_DIR}/elasticity
  ${PROJECT_BINARY_DIR}/elasticity
)
add_subdirectory(
  ${UM_SOURCE_DIR}/cohesive_interface
  ${PROJECT_BINARY_DIR}/cohesive_interface
)
add_subdirectory(
  ${UM_SOURCE_DIR}/nonlinear_elasticity
  ${PROJECT_BINARY_DIR}/nonlinear_elasticity
)
add_subdirectory(
  ${UM_SOURCE_DIR}/helmholtz
  ${PROJECT_BINARY_DIR}/helmholtz
)
add_subdirectory(
  ${UM_SOURCE_DIR}/gels
  ${PROJECT_BINARY_DIR}/gels
)
