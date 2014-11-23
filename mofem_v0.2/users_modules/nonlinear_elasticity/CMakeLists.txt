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

#include subdirectory
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/nonlinear_elasticity")

add_executable(arc_length_nonlinear_elasticity 
  ${UM_SOURCE_DIR}/nonlinear_elasticity/arc_length_nonlinear_elasticity.cpp)
target_link_libraries(arc_length_nonlinear_elasticity 
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

cm_export_file("dome_thin_symm.cub" export_files_arc_length_nonlinearelasticity)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/arc_length_nonlinear_elasticity.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/arc_length_nonlinear_elasticity.sh"
  @ONLY
  )

#testing
add_test(ArcLengthNonLinearElasticityTest 
  ${MoFEM_MPI_RUN} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/arc_length_nonlinear_elasticity
  -my_file dome_thin_symm.cub -ksp_type fgmres 
  -pc_type lu -pc_factor_mat_solver_package mumps
  -ksp_atol 1e-10 -ksp_rtol 1e-10 
  -snes_monitor 
  -snes_type newtonls -snes_linesearch_type basic -snes_linesearch_monitor 
  -snes_max_it 16 -snes_atol 1e-8 -snes_rtol 1e-8 
  -my_sr 1e-6 -my_ms 3 -snes_converged_reason
  -log_summary)

if(ADOL-C_LIBRARY)

  add_executable(nonlinear_dynamics
    ${UM_SOURCE_DIR}/nonlinear_elasticity/nonlinear_dynamics.cpp)
  target_link_libraries(nonlinear_dynamics
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


