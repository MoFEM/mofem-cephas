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


include_directories(${UM_SOURCE_DIR}/elasticity/src)
include_directories(${UM_SOURCE_DIR}/nonlinear_elastic_materials/src)

include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/homogenisation/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/stochastic/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/stochastic/third_party")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/composite_laminates/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/moisture_transport/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/degradation_model/src")


add_library(user_modules_reliability
  ${UM_SOURCE_DIR}/stochastic/third_party/ltqnorm.c
)

# ===================== Transversely isotropic RVE ============================================

# 1. Potential flow problem

#  add_executable(trans_iso_potential_flow ${UM_SOURCE_DIR}/stochastic/trans_iso_potential_flow.cpp)
#target_link_libraries(trans_iso_potential_flow
#  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
#  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})

# 2. Linear displacement boundary condition

  add_executable(trans_iso_RVE_disp_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_disp_7vars.cpp)
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

  add_executable(trans_iso_RVE_Periodic_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Periodic_7vars.cpp)
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
  
  add_executable(trans_iso_RVE_Traction_7vars ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Traction_7vars.cpp)
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

  add_executable(trans_iso_RVE_Geom_disp ${UM_SOURCE_DIR}/stochastic/trans_iso_RVE_Geom_disp.cpp)
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
  
# 6. MCS - Linear displacement boundary condition

  add_executable(MCS_RVE_Disp ${UM_SOURCE_DIR}/stochastic/MCS_RVE_Disp.cpp)
target_link_libraries(
  MCS_RVE_Disp
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
#                                     FE2 & PSFE2                                        #
#  Perturbation-based Stochastic Finite Element 2                                        #
#                                                                                        #
# ========================================================================================

add_executable(FE2_TransIso_Disp_REL 
               ${UM_SOURCE_DIR}/stochastic/FE2_TransIso_Disp_REL.cpp
               )
target_link_libraries(
  FE2_TransIso_Disp_REL
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

#add_executable(PSFE2_TransIso_Disp 
#               ${UM_SOURCE_DIR}/stochastic/PSFE2_TransIso_Disp.cpp
#               )
#target_link_libraries(
#  PSFE2_TransIso_Disp
#  user_modules_homogenisation
#  users_modules_obsolete
#  users_modules
#  mofem_finite_elements 
#  mofem_interfaces 
#  mofem_multi_indices 
#  mofem_petsc 
#  mofem_approx 
#  mofem_third_party 
#  mofem_cblas 
#  ${MoFEM_PROJECT_LIBS}
#  )
  
add_executable(PSFE2_Composite_Laminate
${UM_SOURCE_DIR}/stochastic/PSFE2_Composite_Laminate.cpp
${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
)
target_link_libraries(
	PSFE2_Composite_Laminate
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
	boost_program_options
	${MoFEM_PROJECT_LIBS}
	)
	
	
add_executable(FE2_Degradation_TransIso_Disp_REL 
               ${UM_SOURCE_DIR}/stochastic/FE2_Degradation_TransIso_Disp_REL.cpp
               )
target_link_libraries(
  FE2_Degradation_TransIso_Disp_REL
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
  
  
add_executable(PSFE2_Degradation_TransIso_Disp 
               ${UM_SOURCE_DIR}/stochastic/PSFE2_Degradation_TransIso_Disp.cpp
               )
target_link_libraries(
  PSFE2_Degradation_TransIso_Disp
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
    ${UM_SOURCE_DIR}/stochastic/CompStruct_Rel.cpp
    )
target_link_libraries(
  CompStruct_Rel
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

#------------------------------------------ 
# Case 1: Reliability analysis - Single Ply
#------------------------------------------ 
#add_executable(MsFE_Reliability_Disp 
#               ${UM_SOURCE_DIR}/stochastic/MsFE_Reliability_Disp.cpp
#               )
#target_link_libraries(
#  MsFE_Reliability_Disp
#  user_modules_reliability
#  user_modules_homogenisation
#  users_modules_obsolete 
#  users_modules
#  mofem_finite_elements 
#  mofem_interfaces 
#  mofem_multi_indices 
#  mofem_petsc  
#  mofem_approx 
#  mofem_third_party 
#  mofem_cblas 
#  ${MoFEM_PROJECT_LIBS}
#  )

#---------------------------------------- 
# Case 2: Reliability analysis - Laminate
#---------------------------------------- 
add_executable(MsFE_Reliability_Laminate_Disp 
               ${UM_SOURCE_DIR}/stochastic/MsFE_Reliability_Laminate_Disp.cpp
               ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
               )
target_link_libraries(
  MsFE_Reliability_Laminate_Disp
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

add_executable(MsFE_Reliability_Laminate_FORM_Disp 
               ${UM_SOURCE_DIR}/stochastic/MsFE_Reliability_Laminate_FORM_Disp.cpp
               ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
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
  
add_executable(MsFE_Reliability_Laminate_Periodic 
               ${UM_SOURCE_DIR}/stochastic/MsFE_Reliability_Laminate_Periodic.cpp
               ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
               )
target_link_libraries(
  MsFE_Reliability_Laminate_Periodic
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
    
add_executable(MsFE_Reliability_Laminate_FORM_Periodic 
               ${UM_SOURCE_DIR}/stochastic/MsFE_Reliability_Laminate_FORM_Periodic.cpp
               ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
               )
target_link_libraries(
  MsFE_Reliability_Laminate_FORM_Periodic
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
#         Degradation
#---------------------------------------- 
add_executable(MsFE_Degradation_Reliability_Laminate_Disp 
               ${UM_SOURCE_DIR}/stochastic/MsFE_Degradation_Reliability_Laminate_Disp.cpp
               ${UM_SOURCE_DIR}/elasticity/src/impl/PCMGSetUpViaApproxOrders.cpp
               )
target_link_libraries(
  MsFE_Degradation_Reliability_Laminate_Disp
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
  
#------------------------------ 
# Case 1: Reliability analysis - MCS 
#------------------------------ 
add_executable(FE2_TransIso_Reliability_MCS_Disp 
               ${UM_SOURCE_DIR}/stochastic/FE2_TransIso_Reliability_MCS_Disp.cpp
               )
target_link_libraries(
  FE2_TransIso_Reliability_MCS_Disp
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
  

add_executable(Reliability_MCS_Laminate_Disp 
               ${UM_SOURCE_DIR}/stochastic/Reliability_MCS_Laminate_Disp.cpp
               )
target_link_libraries(
  Reliability_MCS_Laminate_Disp
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
  
  
# =======================================================================================#
#                                                                                        #
#                            MsFE_Reliability_xxx                                        #
#  Multi-scale finite element (MsFE) based reliability analysis for multi-layer laminate #
#  (0) CompStruct_Rel: reliability analysis with analytical limit-state function         #
#  (1) MsFE_Reliability_Disp: MsFE based reliability analysis for RVE with linear        #
#      displacement BC                                                                   #
#                                                                                        #
# ========================================================================================