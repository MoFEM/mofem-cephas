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

if(ADOL-C_LIBRARY)

  add_executable(ground_surface_temperature
    ${UM_SOURCE_DIR}/ground_surface_temperature/atom_tests/ground_surface_temperature.cpp 
    ${UM_SOURCE_DIR}/ground_surface_temperature/third_party/spa.c)

  target_link_libraries(ground_surface_temperature
    mofem_finite_element_library
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    users_modules_obsolete
    complex_for_lazy_obsolete
    mofem_approx
    mofem_third_party
    mofem_cblas
    ${MoFEM_PROJECT_LIBS})

  cm_export_file("Valey2.cub" export_files_ground_surface_temperature)
  add_test(ground_surface_temperature ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/ground_surface_temperature -my_file Valey2.cub)
  add_test(ground_surface_temperature_compare ${CMAKE_COMMAND} -E compare_files 
    ${CMAKE_SOURCE_DIR}/ground_surface_temperature/atom_tests/blessed_files/ground_surface_temperature.txt 
    ${CMAKE_BINARY_DIR}/ground_surface_temperature/atom_tests/ground_surface_temperature.txt)

endif(ADOL-C_LIBRARY)


