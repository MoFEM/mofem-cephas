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
 

include_directories(${UM_SOURCE_DIR}/fracture_mechanics)
include_directories(${UM_SOURCE_DIR}/elasticity/src)
include_directories(${UM_SOURCE_DIR}/nonlinear_elasticity/src)
include_directories(${UM_SOURCE_DIR}/convective_mass_element/src)

add_library(fracture_mechanics
  ${UM_SOURCE_DIR}/fracture_mechanics/FaceSplittingTool.cpp 
  ${UM_SOURCE_DIR}/fracture_mechanics/ConfigurationalFractureMechanics.cpp
)

 
add_executable(spatial_problem 
  ${UM_SOURCE_DIR}/fracture_mechanics/spatial_problem.cpp 
)
target_link_libraries(spatial_problem
  fracture_mechanics
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

add_executable(material_force 
  ${UM_SOURCE_DIR}/fracture_mechanics/material_force.cpp
)
target_link_libraries(material_force 
  fracture_mechanics
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

add_executable(arc_length_material_coupled 
  ${UM_SOURCE_DIR}/fracture_mechanics/arc_length_material_coupled.cpp
)
target_link_libraries(arc_length_material_coupled 
  fracture_mechanics
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

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

cm_export_file("crack1.cub" export_files_material_mechanics)

configure_file(
  "${UM_SOURCE_DIR}/fracture_mechanics/convergence_study.sh.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/convergence_study.sh"
  @ONLY
  )

configure_file(
  "${UM_SOURCE_DIR}/fracture_mechanics/arc_length.sh.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/arc_length.sh"
  @ONLY
  )

file(COPY convergence_study.sh 
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)
file(COPY arc_length.sh
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)

#testing
add_test(ConvergenceStudy_sh ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.sh
  -f crack1.cub -l 1e-4 -o 2 -r 1 -d TEST_SIF -p 2)
add_test(ArcLength_sh ${CMAKE_CURRENT_BINARY_DIR}/arc_length.sh
  -f TEST_SIF/out_spatial_2_1.h5m -a 2e-3 -n 1 -g 1e-6 -d TEST_ARC -p 2
  --snes_linesearch_type=basic --my_nb_sub_steps=1)

if(ADOL-C_LIBRARY)

  add_executable(dynamic_material_coupled 
    ${UM_SOURCE_DIR}/fracture_mechanics/dynamic_material_coupled.cpp
  )
  target_link_libraries(dynamic_material_coupled
    fracture_mechanics
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

endif(ADOL-C_LIBRARY)
