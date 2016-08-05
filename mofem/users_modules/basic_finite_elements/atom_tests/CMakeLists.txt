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

set(UM_CURRENT_SOURCE_DIR ${UM_SOURCE_DIR}/basic_finite_elements/atom_tests)

# body forces test
bfe_add_executable(body_force_atom ${UM_CURRENT_SOURCE_DIR}/body_force_atom_test.cpp)
target_link_libraries(body_force_atom
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)

cm_export_file("elasticity_atom_test_body_force.h5m" export_atom_tests_basic_fe)
add_test(body_force_atom_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/body_force_atom -my_file elasticity_atom_test_body_force.h5m
)
add_test(forces_and_sources_body_force_atom_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_CURRENT_SOURCE_DIR}/blessed_files/body_force_atom_test.txt
  ${CMAKE_CURRENT_BINARY_DIR}/body_force_atom_test.txt
)

# neumann forces
cm_export_file("elasticity_atom_test_01X_pressure.h5m" export_atom_tests_basic_fe)

bfe_add_executable(neumann_forces_atom ${UM_CURRENT_SOURCE_DIR}/neumann_forces.cpp)
target_link_libraries(neumann_forces_atom
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS})

cm_export_file("elasticity_atom_test_01X.h5m" export_atom_tests_basic_fe)
cm_export_file("sphere.h5m" export_atom_tests_basic_fe)

add_test(neumann_elasticity_atom_test_01X_pressure
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/neumann_forces_atom -my_file elasticity_atom_test_01X_pressure.h5m
)
add_test(neumann_elasticity_atom_test_01X_pressure_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_CURRENT_SOURCE_DIR}/blessed_files/forces_and_sources_elasticity_atom_test_01X_pressure.cub.txt
  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_elasticity_atom_test_01X_pressure.h5m.txt
)
add_test(neumann_elasticity_atom_test_01X
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/neumann_forces_atom -my_file elasticity_atom_test_01X.h5m
)
add_test(neumann_elasticity_atom_test_01X_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_CURRENT_SOURCE_DIR}/blessed_files/forces_and_sources_elasticity_atom_test_01X.cub.txt
  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_elasticity_atom_test_01X.h5m.txt
)
add_test(neumann_sphre_pressure_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/neumann_forces_atom -my_file sphere.h5m
)
add_test(neumann_sphre_pressure_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_CURRENT_SOURCE_DIR}/blessed_files/forces_and_sources_sphere.cub.txt
  ${CMAKE_CURRENT_BINARY_DIR}/forces_and_sources_sphere.h5m.txt
)

# add example test for forces_and_sources_fluid_pressure_element_atom
cm_export_file("fluid_presssure_cub_test.h5m" export_atom_tests_basic_fe)
bfe_add_executable(fluid_pressure_element_atom ${UM_CURRENT_SOURCE_DIR}/fluid_pressure_element.cpp)
target_link_libraries(fluid_pressure_element_atom
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)

add_test(fluid_pressure_element_atom_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/fluid_pressure_element_atom -my_file fluid_presssure_cub_test.h5m
)

# add example test for temperature element
cm_export_file("simple_thermal_problem.h5m" export_atom_tests_basic_fe)
bfe_add_executable(thermal_elem_atom ${UM_CURRENT_SOURCE_DIR}/thermal_elem.cpp)
target_link_libraries(thermal_elem_atom
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)

add_test(thermal_elem_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/thermal_elem_atom
-my_file simple_thermal_problem.h5m -ksp_type fgmres -pc_type bjacobi -ksp_atol 0 -ksp_rtol 1e-12 -ksp_monitor)

# add test test for temerature unsteady element
bfe_add_executable(thermal_elem_unsteady_atom ${UM_CURRENT_SOURCE_DIR}/thermal_elem_unsteady.cpp)
target_link_libraries(thermal_elem_unsteady_atom
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS})

add_test(thermal_elem_unsteady_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/thermal_elem_unsteady_atom
  -my_file simple_thermal_problem.h5m
  -ksp_type fgmres -pc_type bjacobi -ksp_atol 0 -ksp_rtol 1e-12
  -snes_type newtonls -snes_linesearch_type basic
  -snes_max_it 100 -snes_atol 1e-8 -snes_rtol 1e-8
  -ts_monitor -ts_type beuler -ts_dt 1 -ts_final_time 2  -ksp_monitor -snes_max_it 3 -snes_monitor
)

# add theraml stress eleme
cm_export_file("elasticity_atom_test_thermal_expansion.h5m" export_atom_tests_basic_fe)
bfe_add_executable(thermal_stress_elem ${UM_CURRENT_SOURCE_DIR}/thermal_stress_elem.cpp)
target_link_libraries(thermal_stress_elem
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)

add_test(thermal_stress_elem_atom_test ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/thermal_stress_elem
  -my_file elasticity_atom_test_thermal_expansion.h5m
)

# add test for ultra weak finite element
bfe_add_executable(ultra_weak_transport ${UM_CURRENT_SOURCE_DIR}/ultra_weak_transport.cpp)
target_link_libraries(ultra_weak_transport
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS})

add_test(ultra_weak_transport_atom ${MoFEM_MPI_RUN} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/ultra_weak_transport -my_file simple_thermal_problem.h5m
  -ksp_type fgmres -pc_type bjacobi -ksp_monitor -ksp_atol 0 -ksp_rtol 1e-12
  -my_max_post_proc_ref_level 0 -my_order 0
)

# Testing field approximation
bfe_add_executable(
  field_approximation ${UM_CURRENT_SOURCE_DIR}/field_approximation.cpp
)
target_link_libraries(field_approximation
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)

add_test(field_approximation_test
  ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/field_approximation
  -my_file elasticity_atom_test_01X_pressure.h5m
  -ksp_monitor -ksp_type fgmres -pc_type bjacobi -ksp_atol 0 -ksp_rtol 1e-12
)
add_test(field_approximation_test_compare ${CMAKE_COMMAND} -E compare_files
  ${CMAKE_CURRENT_SOURCE_DIR}/blessed_files/field_approximation.txt
  ${CMAKE_CURRENT_BINARY_DIR}/field_approximation.txt
)

if(ADOL-C_LIBRARY)

  cm_export_file("LShape.h5m" export_atom_tests_basic_fe)

  bfe_add_executable(nonlinear_elastic_basic_test ${UM_CURRENT_SOURCE_DIR}/nonlinear_elastic.cpp)
  target_link_libraries(nonlinear_elastic_basic_test
    users_modules
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_cblas
    mofem_third_party
    ${MoFEM_PROJECT_LIBS}
  )

  add_test(nonlinear_elastic_basic_test ${MoFEM_MPI_RUN} -np 1
    ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_basic_test -my_file LShape.h5m
  )

endif(ADOL-C_LIBRARY)

# kelving voigt damper
if(ADOL-C_LIBRARY)

  bfe_add_executable(damper_jacobian_test
    ${UM_CURRENT_SOURCE_DIR}/damper_jacobian_test.cpp
  )

  target_link_libraries(damper_jacobian_test
    users_modules
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    mofem_cblas
    ${MoFEM_PROJECT_LIBS}
  )

  add_test(damper_jacobian_atom_test
    ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/damper_jacobian_test
    -my_file fluid_presssure_cub_test.h5m
  )

endif(ADOL-C_LIBRARY)


# Convective element

if(ADOL-C_LIBRARY)

  bfe_add_executable(convective_matrix
    ${UM_CURRENT_SOURCE_DIR}/convective_matrix.cpp
  )

  target_link_libraries(convective_matrix
    users_modules
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    mofem_cblas
    ${MoFEM_PROJECT_LIBS}
  )

  cm_export_file("LShape.cub" export_atom_tests_basic_fe)
  add_test(convective_matrix
    ${MoFEM_MPI_RUN} -np 1 ${CMAKE_CURRENT_BINARY_DIR}/convective_matrix -my_file LShape.cub
  )

endif(ADOL-C_LIBRARY)
