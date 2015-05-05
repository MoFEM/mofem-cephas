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



include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/homogenisation/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/stochastic/src")

# ===================== Transversely isotropic RVE ============================================

# 1. Potential flow problem

#  add_executable(trans_iso_potential_flow ${UM_SOURCE_DIR}/stochastic/trans_iso_potential_flow.cpp)
#target_link_libraries(trans_iso_potential_flow
#  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
#  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})

# 2. Linear displacement boundary condition

  add_executable(trans_iso_RVE_disp_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_disp_7vars.cpp)
target_link_libraries(trans_iso_RVE_disp_7vars
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})

# 3. Periodic boundary condition

  add_executable(trans_iso_RVE_Periodic_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Periodic_7vars.cpp)
target_link_libraries(trans_iso_RVE_Periodic_7vars
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})  

# 4. Constant traction boundary condition
  
  add_executable(trans_iso_RVE_Traction_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Traction_7vars.cpp)
target_link_libraries(trans_iso_RVE_Traction_7vars
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})  
  
# 5. Waviness - Linear displacement boundary condition

  add_executable(trans_iso_RVE_Geom_disp ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Geom_disp.cpp)
target_link_libraries(trans_iso_RVE_Geom_disp
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})   
