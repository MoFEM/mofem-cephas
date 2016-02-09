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

file(
  GLOB_RECURSE INSTLLED_MODULES
  FOLLOW_SYMLINKS
  ?*/InstalledAddModule.cmake
)

# git pull for all users modules
add_custom_target(
  update_users_modules
  COMMENT "Update all modules ..." VERBATIM
)
foreach(LOOP_MODULE ${INSTLLED_MODULES})
  # message(STATUS "${LOOP_MODULE}")
  string(REGEX REPLACE
    "/+InstalledAddModule.cmake" ""
    MODULE_DIRECTORY ${LOOP_MODULE}
  )
  # message(STATUS "${MODULE_DIRECTORY}")
  string(REGEX REPLACE
    ".*/+" ""
    UPDATE_MODULE_NAME ${MODULE_DIRECTORY}
  )
  message(STATUS "Add module ... ${MODULE_DIRECTORY}")
  include(${LOOP_MODULE})
  message(STATUS "Add custom target ... update_${UPDATE_MODULE_NAME}")
  add_custom_target(
    update_${UPDATE_MODULE_NAME}
    COMMAND git pull
    WORKING_DIRECTORY ${MODULE_DIRECTORY}
    COMMENT "Update module ... ${UPDATE_MODULE_NAME}" VERBATIM
  )
  add_dependencies(update_users_modules update_${UPDATE_MODULE_NAME})
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
  ${UM_SOURCE_DIR}/reliability
  ${PROJECT_BINARY_DIR}/reliability
)
