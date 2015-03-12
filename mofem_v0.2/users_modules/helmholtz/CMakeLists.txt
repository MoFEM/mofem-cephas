# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)                                                                                                                   
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.•
# It can be freely used for educational and research purposes•
# by other institutions. If you use this softwre pleas cite my work.•
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
file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/meshes/cube_No_reflecting_BC_Problem.cub
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
 file(COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/meshes/impinging_sphere.cub
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
# For convinience source files from MoFEM project are used. If user is
# implementing independent application it will be using instead project
# directories.

include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/helmholtz/src")
#add_subdirectory(atom_tests)
#add_subdirectory(meshes)

add_executable(acoustic_steady ${UM_SOURCE_DIR}/helmholtz/acoustic_steady.cpp)
target_link_libraries(acoustic_steady
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})
#add_executable(incident_wave_acoustic ${UM_SOURCE_DIR}/helmholtz/incident_wave_acoustic.cpp)
#target_link_libraries(incident_wave_acoustic
#  mofem_finite_element_library
#  mofem_finite_elements
#  mofem_interfaces
#  mofem_multi_indices
#  mofem_petsc
#  mofem_approx
#  mofem_third_party
#  mofem_cblas
#  ${MoFEM_PROJECT_LIBS})

add_executable(error_norm ${UM_SOURCE_DIR}/helmholtz/error_norm.cpp)
target_link_libraries(error_norm
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})


add_executable(best_approximation ${UM_SOURCE_DIR}/helmholtz/best_approximation.cpp)
target_link_libraries(best_approximation
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

#add_executable(exact_solution ${UM_SOURCE_DIR}/helmholtz/exact_solution.cpp)
#target_link_libraries(exact_solution
#  mofem_finite_element_library
#  mofem_finite_elements
#  mofem_interfaces
#  mofem_multi_indices
#  mofem_petsc
#  mofem_approx
#  mofem_third_party
#  mofem_cblas
#  ${MoFEM_PROJECT_LIBS})


#add_executable(thermal_unseady ${UM_SOURCE_DIR}/thermal/thermal_unsteady.cpp)
#target_link_libraries(thermal_unseady•
#  mofem_finite_element_library
#  mofem_finite_elements
#  mofem_interfaces
#  mofem_multi_indices
#  mofem_petsc
#  mofem_approx
#  mofem_third_party
#  mofem_cblas
#  ${MoFEM_PROJECT_LIBS})

#cm_export_file("cube_No_reflecting_BC_Problem.cub" export_files_acoustic_steady)

#testing
#add_test(Helmholtz_Operator_test•
#  ${MoFEM_MPI_RUN} -np 4 ${CMAKE_CURRENT_BINARY_DIR}/acoustic_steady
#  -my_file cube_No_reflecting_BC_Problem.cub -ksp_type fgmres•
#  -pc_type lu -pc_factor_mat_solver_package superlu_dist
#  -ksp_atol 1e-10 -ksp_rtol 1e-10•
#  -ksp_monitor•
#  -my_order 1 
#  -my_max_post_proc_ref_level 0
#  -log_summary)

