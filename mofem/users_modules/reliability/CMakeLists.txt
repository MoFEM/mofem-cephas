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

include_directories(${UM_SOURCE_DIR}/obsolete/src)
include_directories(${UM_SOURCE_DIR}/obsolete/src/c)
include_directories(${UM_SOURCE_DIR}/obsolete/src/c_impl)

include_directories(${PROJECT_SOURCE_DIR}/nonlinear_elasticity)
include_directories(${PROJECT_SOURCE_DIR}/nonlinear_elasticity/src)
include_directories(${PROJECT_SOURCE_DIR}/nonlinear_elastic_materials/src)
include_directories(${PROJECT_SOURCE_DIR}/elasticity/src)

include_directories(${PROJECT_SOURCE_DIR}/homogenisation/src)
include_directories(${PROJECT_SOURCE_DIR}/moisture_transport/src)
include_directories("${UM_SOURCE_DIR}/cohesive_interface/src")
include_directories(${UM_SOURCE_DIR}/small_strain_plasticity/src)

include_directories("${UM_SOURCE_DIR}/reliability/src")
include_directories("${UM_SOURCE_DIR}/reliability/third_party")
#include_directories("${UM_SOURCE_DIR}/composite_laminates/src")
#include_directories("${UM_SOURCE_DIR}/degradation_model/src")

include_directories("${UM_SOURCE_DIR}/basic_finite_elements/src")


add_library(user_modules_reliability
  ${UM_SOURCE_DIR}/reliability/third_party/ltqnorm.c
)	   

add_executable(psfe_rve_mechanical ${PROJECT_SOURCE_DIR}/reliability/psfe_rve_mechanical.cpp)
target_link_libraries(psfe_rve_mechanical
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

# 2. Linear displacement boundary condition

  add_executable(trans_iso_RVE_disp_7vars ${PROJECT_SOURCE_DIR}/reliability/trans_iso_RVE_disp_7vars.cpp)
target_link_libraries(
  trans_iso_RVE_disp_7vars
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
  
# 3. Periodic boundary condition

  add_executable(trans_iso_RVE_Periodic_7vars ${UM_SOURCE_DIR}/reliability/trans_iso_RVE_Periodic_7vars.cpp)
target_link_libraries(
  trans_iso_RVE_Periodic_7vars
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

# 4. Constant traction boundary condition
  
  add_executable(trans_iso_RVE_Traction_7vars ${UM_SOURCE_DIR}/reliability/trans_iso_RVE_Traction_7vars.cpp)
target_link_libraries(
  trans_iso_RVE_Traction_7vars
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
  
# 5. Waviness - Linear displacement boundary condition

  add_executable(trans_iso_RVE_Geom_disp ${UM_SOURCE_DIR}/reliability/trans_iso_RVE_Geom_disp.cpp)
target_link_libraries(
  trans_iso_RVE_Geom_disp
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

# =======================================================================================#
#                                                                                        #
#                            MsFE_Reliability_xxx                                        #
#  Multi-scale finite element (MsFE) based reliability analysis for single layer         #
#  (0) CompStruct_Rel: reliability analysis with analytical limit-state function         #
#  (1) MsFE_Reliability_Disp: MsFE based reliability analysis for RVE with linear        #
#      displacement BC                                                                   #
#                                                                                        #
# ========================================================================================

#------------------------------ 
# Case 0: Reliability analysis 
#------------------------------  
add_executable(CompStruct_Rel 
    ${UM_SOURCE_DIR}/reliability/CompStruct_Rel.cpp
    )
target_link_libraries(
  CompStruct_Rel
  user_modules_reliability
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


#---------------------------------------- 
# Case 2: Reliability analysis - Laminate
#---------------------------------------- 
# ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/PCMGSetUpViaApproxOrders.cpp
add_executable(MsFE_Reliability_Laminate_FORM_Disp 
               ${UM_SOURCE_DIR}/reliability/MsFE_Reliability_Laminate_FORM_Disp.cpp
               )
target_link_libraries(
  MsFE_Reliability_Laminate_FORM_Disp
  user_modules_reliability
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
  
  
  #---------------------------------------- 
# Case 3: Reliability analysis - Laminate
#---------------------------------------- 
# ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/PCMGSetUpViaApproxOrders.cpp
add_executable(FE_Reliability_Laminate_FORM_Disp 
               ${UM_SOURCE_DIR}/reliability/FE_Reliability_Laminate_FORM_Disp.cpp
               )
target_link_libraries(
  FE_Reliability_Laminate_FORM_Disp
  user_modules_reliability
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
  