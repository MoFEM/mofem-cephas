# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
# It can be freely used for educational and research purposes 
# by other institutions. If you use this softwre pleas cite my work. 
#
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

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# For convinience source files from MoFEM project are used. If user is
# implementing independent application it will be using instead project
# directories.

add_executable(thermal_steady ${MoFEM_PROJECT_SOURCE_DIR}/examples/thermal/thermal_steady.cpp)
target_link_libraries(thermal_steady 
  mofem_cblas
  mofem_third_party
  mofem_finite_element_library
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  ${MoFEM_PROJECT_LIBS})

add_executable(thermal_unseady ${MoFEM_PROJECT_SOURCE_DIR}/examples/thermal/thermal_unsteady.cpp)
target_link_libraries(thermal_unseady 
  mofem_cblas
  mofem_third_party
  mofem_finite_element_library
  mofem_approx
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  ${MoFEM_PROJECT_LIBS})


