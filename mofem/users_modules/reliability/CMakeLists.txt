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


#include_directories(${UM_SOURCE_DIR}/elasticity/src)
#include_directories(${UM_SOURCE_DIR}/nonlinear_elastic_materials/src)

#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/homogenisation/src")
#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/reliability/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/reliability/third_party")
#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/composite_laminates/src")
#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/moisture_transport/src")
#include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/degradation_model/src")


add_library(user_modules_reliability
  ${UM_SOURCE_DIR}/reliability/third_party/ltqnorm.c
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
  boost_program_options
  ${MoFEM_PROJECT_LIBS}
)