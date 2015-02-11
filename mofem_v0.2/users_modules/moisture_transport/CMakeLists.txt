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

include_directories("${MoFEM_PROJECT_SOURCE_DIR}/users_modules/moisture_transport/src")


# ===================== Diffusion steady ============================================
add_executable(diffusion_steady ${UM_SOURCE_DIR}/moisture_transport/diffusion_steady.cpp)
target_link_libraries(diffusion_steady
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
add_executable(diffusion_unsteady ${UM_SOURCE_DIR}/moisture_transport/diffusion_unsteady.cpp)
target_link_libraries(diffusion_unsteady
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})


# ===================== Darceys Flow steady ============================================
add_executable(darceys_flow_steady ${UM_SOURCE_DIR}/moisture_transport/darceys_flow_steady.cpp)
target_link_libraries(darceys_flow_steady
  mofem_finite_element_library
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  ${MoFEM_PROJECT_LIBS})
