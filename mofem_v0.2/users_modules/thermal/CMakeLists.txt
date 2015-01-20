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
    ${UM_SOURCE_DIR}/thermal/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

include_directories("${UM_SOURCE_DIR}/ground_surface_temperature/src")
include_directories("${UM_SOURCE_DIR}/ground_surface_temperature/third_party")

# For convinience source files from MoFEM project are used. If user is
# implementing independent application it will be using instead project
# directories.

add_executable(thermal_steady ${UM_SOURCE_DIR}/thermal/thermal_steady.cpp)
target_link_libraries(thermal_steady 
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

add_executable(thermal_unseady ${UM_SOURCE_DIR}/thermal/thermal_unsteady.cpp)
target_link_libraries(thermal_unseady 
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

if(ADOL-C_LIBRARY)

  add_executable(dm_thermal_unseady 
    ${UM_SOURCE_DIR}/thermal/dm_thermal_unsteady.cpp
    ${UM_SOURCE_DIR}/ground_surface_temperature/third_party/spa.c)

else(ADOL-C_LIBRARY) 

  add_executable(dm_thermal_unseady ${UM_SOURCE_DIR}/thermal/dm_thermal_unsteady.cpp)

endif(ADOL-C_LIBRARY)

target_link_libraries(dm_thermal_unseady 
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  boost_program_options
  ${MoFEM_PROJECT_LIBS})


