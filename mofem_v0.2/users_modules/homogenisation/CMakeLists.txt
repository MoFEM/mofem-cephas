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


include_directories(${UM_SOURCE_DIR}/nonlinear_elasticity)
include_directories(${UM_SOURCE_DIR}/nonlinear_elasticity/src)
include_directories(${UM_SOURCE_DIR}/nonlinear_elastic_materials/src)
include_directories(${UM_SOURCE_DIR}/elasticity/src)

include_directories(${UM_SOURCE_DIR}/homogenisation/src)
include_directories(${UM_SOURCE_DIR}/moisture_transport/src)

add_library(user_modules_homogenisation
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Disp.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Homogenized_Stress_Disp.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Traction.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_RigidBodyTranslation.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_RigidBodyRotation.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Homogenized_Stress_Traction.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Periodic.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Homogenized_Stress_Periodic.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Disp_Multi_Rhs.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Traction_Multi_Rhs.cpp
  ${UM_SOURCE_DIR}/homogenisation/src/impl/ElasticFE_RVELagrange_Periodic_Multi_Rhs.cpp
)

add_subdirectory(atom_tests)

# ===================== Simple cube ============================================

add_executable(Cube_Disp ${UM_SOURCE_DIR}/homogenisation/Cube_Disp.cpp)
target_link_libraries(Cube_Disp
  user_modules_homogenisation
  users_modules_obsolete
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas ${MoFEM_PROJECT_LIBS}
)

add_executable(Cube_Trac ${UM_SOURCE_DIR}/homogenisation/Cube_Trac.cpp)
target_link_libraries(Cube_Trac
  user_modules_homogenisation
  users_modules_obsolete
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas ${MoFEM_PROJECT_LIBS}
)

add_executable(Cube_Periodic ${UM_SOURCE_DIR}/homogenisation/Cube_Periodic.cpp)
target_link_libraries(Cube_Periodic
  user_modules_homogenisation
  users_modules_obsolete
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

# ===================== Simple cube with high-order geometry============================================

add_executable(Cube_Disp_HOG ${UM_SOURCE_DIR}/homogenisation/Cube_Disp_HOG.cpp)
target_link_libraries(Cube_Disp_HOG
user_modules_homogenisation
users_modules_obsolete
users_modules
mofem_finite_elements
mofem_interfaces
mofem_multi_indices
mofem_petsc
mofem_approx
mofem_third_party
mofem_cblas ${MoFEM_PROJECT_LIBS}
)

add_executable(Cube_Trac_HOG ${UM_SOURCE_DIR}/homogenisation/Cube_Trac_HOG.cpp)
target_link_libraries(Cube_Trac_HOG
user_modules_homogenisation
users_modules_obsolete
users_modules
mofem_finite_elements
mofem_interfaces
mofem_multi_indices
mofem_petsc
mofem_approx
mofem_third_party
mofem_cblas ${MoFEM_PROJECT_LIBS}
)


# =================== Simple cube with reinforcement ================================

add_executable(Cube_With_Inclusion_Disp ${UM_SOURCE_DIR}/homogenisation/Cube_With_Inclusion_Disp.cpp)
target_link_libraries(Cube_With_Inclusion_Disp
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(Cube_With_Inclusion_Traction ${UM_SOURCE_DIR}/homogenisation/Cube_With_Inclusion_Traction.cpp)
target_link_libraries(Cube_With_Inclusion_Traction
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(Cube_With_Inclusion_Periodic ${UM_SOURCE_DIR}/homogenisation/Cube_With_Inclusion_Periodic.cpp)
target_link_libraries(Cube_With_Inclusion_Periodic
  user_modules_homogenisation
  users_modules_obsolete
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

# =================== Adding textile RVE homogenisation (without interface elements)================================
add_executable(TransIso_Potential_Flow_NoInt ${UM_SOURCE_DIR}/homogenisation/TransIso_Potential_Flow_NoInt.cpp)
target_link_libraries(TransIso_Potential_Flow_NoInt
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(TransIso_Disp_NoInt ${UM_SOURCE_DIR}/homogenisation/TransIso_Disp_NoInt.cpp)
target_link_libraries(TransIso_Disp_NoInt
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(TransIso_Trac_NoInt ${UM_SOURCE_DIR}/homogenisation/TransIso_Trac_NoInt.cpp)
target_link_libraries(TransIso_Trac_NoInt
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(TransIso_Periodic_NoInt ${UM_SOURCE_DIR}/homogenisation/TransIso_Periodic_NoInt.cpp)
target_link_libraries(TransIso_Periodic_NoInt
  user_modules_homogenisation
  users_modules_obsolete
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


# ===================== Simple cube or cube with reinformenct for moisture transport ============================================
add_executable(MoistureTrans_Disp ${UM_SOURCE_DIR}/homogenisation/MoistureTrans_Disp.cpp)
target_link_libraries(MoistureTrans_Disp
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(MoistureTrans_Trac ${UM_SOURCE_DIR}/homogenisation/MoistureTrans_Trac.cpp)
target_link_libraries(MoistureTrans_Trac
  user_modules_homogenisation
  users_modules_obsolete
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

add_executable(MoistureTrans_Periodic ${UM_SOURCE_DIR}/homogenisation/MoistureTrans_Periodic.cpp)
target_link_libraries(MoistureTrans_Periodic
  user_modules_homogenisation
  users_modules_obsolete
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

# ===================== FE2 elasticity  ============================================
add_executable(FE2_Elasticity_Disp ${UM_SOURCE_DIR}/homogenisation/FE2_Elasticity_Disp.cpp)
target_link_libraries(FE2_Elasticity_Disp
  user_modules_homogenisation
  users_modules_obsolete
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

add_executable(FE2_Elasticity_Trac ${UM_SOURCE_DIR}/homogenisation/FE2_Elasticity_Trac.cpp)
target_link_libraries(FE2_Elasticity_Trac
user_modules_homogenisation
users_modules_obsolete
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

add_executable(FE2_Elasticity_Periodic ${UM_SOURCE_DIR}/homogenisation/FE2_Elasticity_Periodic.cpp)
target_link_libraries(FE2_Elasticity_Periodic
user_modules_homogenisation
users_modules_obsolete
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
