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

if(WITH_MODULE_OBSOLETE)
  if(NOT EXISTS ${UM_SOURCE_DIR}/obsolete)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_obsolete.git obsolete
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/obsolete)
endif(WITH_MODULE_OBSOLETE)

if(WITH_MODULE_HOMOGENISATION)
  if(NOT EXISTS ${UM_SOURCE_DIR}/homogenisation)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_homogenisation homogenisation
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/homogenisation)
endif(WITH_MODULE_HOMOGENISATION)


if(WITH_MODULE_FRACTURE_MECHANICS)
  if(NOT EXISTS ${UM_SOURCE_DIR}/fracture_mechanics)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_fracture_mechanics fracture_mechanics
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/fracture_mechanics)
endif(WITH_MODULE_FRACTURE_MECHANICS)

if(WITH_MODULE_GELS)
  if(NOT EXISTS ${UM_SOURCE_DIR}/gels)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_gels gels
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/gels)
endif(WITH_MODULE_GELS)

if(WITH_MODULE_STRAIN_PLASTICITY)
  if(NOT EXISTS ${UM_SOURCE_DIR}/strain_plasticity)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_small_strain_plasticity small_strain_plasticity
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/strain_plasticity)
endif(WITH_MODULE_STRAIN_PLASTICITY)

if(WITH_MODULE_SOLID_SHELL_PRISM_ELEMENT)
  if(NOT EXISTS ${UM_SOURCE_DIR}/solid_shell_prism_element)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_solid_shell_prism_element solid_shell_prism_element
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/solid_shell_prism_element)
endif(WITH_MODULE_SOLID_SHELL_PRISM_ELEMENT)

if(WITH_MODULE_MINIMAL_SURFACE_EQUATION)
  if(NOT EXISTS ${UM_SOURCE_DIR}/minimal_surface_equation)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_minimal_surface_equation minimal_surface_equation
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/minimal_surface_equation)
endif(WITH_MODULE_MINIMAL_SURFACE_EQUATION)

if(WITH_MODULE_MINIMAL_HELMHOLTZ)
  if(NOT EXISTS ${UM_SOURCE_DIR}/helmholtz)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} clone https://bitbucket.org/likask/mofem_um_helmholtz helmholtz
      WORKING_DIRECTORY ${UM_SOURCE_DIR}
    )
  endif(NOT EXISTS ${UM_SOURCE_DIR}/helmholtz)
endif(WITH_MODULE_MINIMAL_HELMHOLTZ)

file(
  GLOB_RECURSE INSTLLED_MODULES
  FOLLOW_SYMLINKS
  ?*/InstalledAddModule.cmake
)

# inattal modules && git pull for all users modules
add_custom_target(
  update_users_modules
  COMMENT "Update all modules ..." VERBATIM
)
add_custom_target(
  checkout_CDashTesting
  COMMENT "Checkout CDashTesting branch ..." VERBATIM
)
add_custom_target(
  checkout_master
  COMMENT "Checkout master branch ..." VERBATIM
)
add_custom_target(
  merge_CDashTesting
  COMMENT "Make merge CDashTesting branch ..." VERBATIM
)
foreach(LOOP_MODULE ${INSTLLED_MODULES})
  string(REGEX REPLACE
    "/+InstalledAddModule.cmake" ""
    MODULE_DIRECTORY ${LOOP_MODULE}
  )
  string(REGEX REPLACE
    ".*/+" ""
    MODULE_NAME ${MODULE_DIRECTORY}
  )
  string(TOUPPER ${MODULE_NAME} MODULE_NAME)
  message(STATUS "Add definitions to the compiler command -DWITH_MODULE_${MODULE_NAME}")
  add_definitions(-DWITH_MODULE_${MODULE_NAME})
endforeach(LOOP_MODULE)
foreach(LOOP_MODULE ${INSTLLED_MODULES})
  # message(STATUS "${LOOP_MODULE}")
  string(REGEX REPLACE
    "/+InstalledAddModule.cmake" ""
    MODULE_DIRECTORY ${LOOP_MODULE}
  )
  message(STATUS "Add module ... ${MODULE_DIRECTORY}")
  string(REGEX REPLACE
    ".*/+" ""
    MODULE_NAME ${MODULE_DIRECTORY}
  )
  message(STATUS "Add custom targets for ${MODULE_NAME}")
  add_custom_target(
    update_${MODULE_NAME}
    COMMAND ${GIT_EXECUTABLE} pull
    WORKING_DIRECTORY ${MODULE_DIRECTORY}
    COMMENT "Update module ... ${MODULE_NAME}" VERBATIM
  )
  add_dependencies(update_users_modules update_${MODULE_NAME})
  add_custom_target(
    checkout_CDashTesting_${MODULE_NAME}
    COMMAND ${GIT_EXECUTABLE} checkout CDashTesting
    WORKING_DIRECTORY ${MODULE_DIRECTORY}
    COMMENT "Checkout CDashTesting baranch for module ${MODULE_NAME}" VERBATIM
  )
  add_dependencies(checkout_CDashTesting checkout_CDashTesting_${MODULE_NAME})
  add_custom_target(
    checkout_master_${MODULE_NAME}
    COMMAND ${GIT_EXECUTABLE} checkout master
    WORKING_DIRECTORY ${MODULE_DIRECTORY}
    COMMENT "Checkout master baranch for module ${MODULE_NAME}" VERBATIM
  )
  add_dependencies(checkout_master checkout_master_${MODULE_NAME})
  add_custom_target(
    merge_CDashTesting_${MODULE_NAME}
    COMMAND ${GIT_EXECUTABLE} merge --ff CDashTesting
    WORKING_DIRECTORY ${MODULE_DIRECTORY}
    COMMENT "Make merge CDashTesting branch for module ${MODULE_NAME}" VERBATIM
  )
  add_dependencies(merge_CDashTesting merge_CDashTesting_${MODULE_NAME})
  # include module
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
