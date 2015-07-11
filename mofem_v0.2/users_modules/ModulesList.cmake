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

include(basic_finite_elements/AddModule.cmake)
include(analytical_dirichlet_boundary_conditions/AddModule.cmake)
include(convective_mass_element/AddModule.cmake)
include(obsolete/AddModule.cmake)

# Users modules library, common for all programs
add_library(users_modules ${UM_LIB_SOURCES})

# Atom test user modules
add_subdirectory(basic_finite_elements/atom_tests)
add_subdirectory(analytical_dirichlet_boundary_conditions/atom_tests)
add_subdirectory(convective_mass_element/atom_tests)
add_subdirectory(obsolete/atom_tests)

# Users Programs
add_subdirectory(thermal)
add_subdirectory(ultraweak)
add_subdirectory(ground_surface_temperature)
add_subdirectory(elasticity)
add_subdirectory(cohesive_interface)
add_subdirectory(nonlinear_elasticity)
add_subdirectory(fracture_mechanics)
add_subdirectory(helmholtz)
add_subdirectory(homogenisation)
add_subdirectory(moisture_transport)
add_subdirectory(degradation_model)
add_subdirectory(composite_laminates)