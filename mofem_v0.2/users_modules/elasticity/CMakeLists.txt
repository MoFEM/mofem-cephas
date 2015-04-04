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

include_directories("${UM_SOURCE_DIR}/elasticity/src")
include_directories("${UM_SOURCE_DIR}/nonlinear_elastic_materials/src")

add_executable(elasticity 
  ${UM_SOURCE_DIR}/elasticity/elasticity.cpp
  ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
)
target_link_libraries(elasticity 
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

cm_export_file("LShape.cub" export_files_elasticity)
cm_export_file("data.txt" export_files_elasticity)

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

#testing
add_test(ElasticityTest 
  ${MoFEM_MPI_RUN} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/elasticity 
  -my_file LShape.cub 
  -ksp_type fgmres 
  -pc_type lu 
  -pc_factor_mat_solver_package mumps
  -ksp_monitor -log_summary)

add_test(MultiGridViaOrderSolver
  ${MoFEM_MPI_RUN} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/elasticity 
  -my_file LShape.cub 
  -ksp_type fgmres -ksp_monitor 
  -my_order 4
  -pc_type mg 
  -mofem_mg_verbose 1 -mofem_mg_coarse_order 2 -mofem_mg_levels 3
  -mg_coarse_ksp_type preonly -mg_coarse_pc_type lu -mg_coarse_pc_factor_mat_solver_package mumps
  -mg_levels_1_ksp_type gmres -mg_levels_1_ksp_max_it 30
  -mg_levels_2_ksp_type gmres -mg_levels_2_ksp_max_it 30
  -log_summary)

