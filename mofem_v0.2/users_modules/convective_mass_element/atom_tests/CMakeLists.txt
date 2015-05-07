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

  find_program(GSED NAMES gsed sed)
  if(NOT GSED)
    message(FATAL_ERROR gsed/sed not found)
  endif (NOT GSED)

  add_executable(convective_matrix
    ${UM_SOURCE_DIR}/convective_mass_element/atom_tests/convective_matrix.cpp)

  target_link_libraries(convective_matrix
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

  configure_file(
    "${UM_SOURCE_DIR}/convective_mass_element/atom_tests/convective_matrix.sh.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/convective_matrix.sh"
    @ONLY
  )
  file(COPY convective_matrix.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)

  cm_export_file("LShape.cub" export_files_convective_matrix)
  add_test(convective_matrix convective_matrix.sh)
  add_test(convective_matrix_compare ${CMAKE_COMMAND} -E compare_files
    ${CMAKE_SOURCE_DIR}/convective_mass_element/atom_tests/blessed_files/convective_matrix.txt
    ${CMAKE_BINARY_DIR}/convective_mass_element/atom_tests/convective_matrix.txt)

endif(ADOL-C_LIBRARY)
