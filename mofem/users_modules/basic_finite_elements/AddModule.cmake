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

include_directories(${UM_SOURCE_DIR}/basic_finite_elements/src)

set(UM_LIB_SOURCES
  ${UM_LIB_SOURCES}
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/ConstrainMatrixCtx.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/ArcLengthTools.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/DirichletBC.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/ThermalElement.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/NodeForce.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/EdgeForce.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/SurfacePressure.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/PostProcOnRefMesh.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/PCMGSetUpViaApproxOrders.cpp
  ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/ConvectiveMassElement.cpp
)

if(ADOL-C_LIBRARY)
  set(UM_LIB_SOURCES
    ${UM_LIB_SOURCES}
    ${UM_SOURCE_DIR}/basic_finite_elements/src/impl/NonLinearElasticElement.cpp
  )
endif(ADOL-C_LIBRARY)

add_subdirectory(${PROJECT_SOURCE_DIR}/basic_finite_elements/elasticity)
add_subdirectory(${PROJECT_SOURCE_DIR}/basic_finite_elements/nonlinear_elasticity)
add_subdirectory(${PROJECT_SOURCE_DIR}/basic_finite_elements/ultraweak)
add_subdirectory(${PROJECT_SOURCE_DIR}/basic_finite_elements/thermal)
add_subdirectory(${PROJECT_SOURCE_DIR}/basic_finite_elements/cohesive_interface)
