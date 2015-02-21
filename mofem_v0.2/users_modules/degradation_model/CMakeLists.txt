# Copyright (C) 2013, Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
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
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/moisture_transport/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/degradation_model/src")
include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/obsolete/src")


# ===================== Thermal homogenisaiton ============================================
add_executable(Thermal_Disp_MRhs ${UM_SOURCE_DIR}/degradation_model/Thermal_Disp_MRhs.cpp)
target_link_libraries(Thermal_Disp_MRhs
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})

# ===================== moisture transport homogenisaiton ============================================
add_executable(MoistureTrans_Disp_MRhs ${UM_SOURCE_DIR}/degradation_model/MoistureTrans_Disp_MRhs.cpp)
target_link_libraries(MoistureTrans_Disp_MRhs
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})

# ===================== moisture transport homogenisaiton ============================================
add_executable(MoistureTrans_Disp_MRhs_with_incl ${UM_SOURCE_DIR}/degradation_model/MoistureTrans_Disp_MRhs_with_incl.cpp)
target_link_libraries(MoistureTrans_Disp_MRhs_with_incl
  mofem_finite_element_library mofem_finite_elements mofem_interfaces mofem_multi_indices mofem_petsc users_modules_obsolete user_modules_homogenisation
  mofem_approx mofem_third_party mofem_cblas ${MoFEM_PROJECT_LIBS})



# ===================== Diffusion unsteady ============================================
add_executable(thermal_unsteady_inputDmat ${UM_SOURCE_DIR}/degradation_model/thermal_unsteady_inputDmat.cpp)
target_link_libraries(thermal_unsteady_inputDmat
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

# ===================== Diffusion unsteady ============================================
add_executable(diffusion_unsteady_inputDmat ${UM_SOURCE_DIR}/degradation_model/diffusion_unsteady_inputDmat.cpp)
target_link_libraries(diffusion_unsteady_inputDmat
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})

# ===================== Degredation ============================================
add_executable(wt_unsteady_macro ${UM_SOURCE_DIR}/degradation_model/wt_unsteady_macro.cpp)
target_link_libraries(wt_unsteady_macro
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})


# ===================== Mechanical with isomatric materials ============================================
add_executable(mech_macro ${UM_SOURCE_DIR}/degradation_model/mech_macro_iso.cpp)
target_link_libraries(mech_macro
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  users_modules_obsolete
  user_modules_homogenisation
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})


  configure_file(
    "${UM_SOURCE_DIR}/degradation_model/do_vtk1.sh.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/do_vtk1.sh"
    @ONLY
  )

  file(COPY do_vtk1.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)


# ===================== Mechanical with trans-iso fibres ============================================
add_executable(mech_macro_trans_iso ${UM_SOURCE_DIR}/degradation_model/mech_macro_trans_iso.cpp)
target_link_libraries(mech_macro_trans_iso
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  users_modules_obsolete
  user_modules_homogenisation
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})
  
  
  
add_executable(Strain_calculation ${UM_SOURCE_DIR}/degradation_model/Strain_calculation.cpp)
target_link_libraries(Strain_calculation
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  users_modules_obsolete
  user_modules_homogenisation
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})




