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

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/third_party")
include_directories("${PROJECT_SOURCE_DIR}/third_party/cblas")
include_directories("${PROJECT_SOURCE_DIR}/src/approximation")
include_directories("${PROJECT_SOURCE_DIR}/src/approximation/c")
include_directories("${PROJECT_SOURCE_DIR}/src/multi_indices")
include_directories("${PROJECT_SOURCE_DIR}/src/interfaces")
include_directories("${PROJECT_SOURCE_DIR}/src/petsc")
include_directories("${PROJECT_SOURCE_DIR}/src/finite_elements")
include_directories("${PROJECT_SOURCE_DIR}/finite_element_library")

include_directories("${PROJECT_BINARY_DIR}/include")

link_directories("${PROJECT_BINARY_DIR}/src/approximation")
link_directories("${PROJECT_BINARY_DIR}/src/multi_indices")
link_directories("${PROJECT_BINARY_DIR}/src/interfaces")
link_directories("${PROJECT_BINARY_DIR}/src/petsc")
link_directories("${PROJECT_BINARY_DIR}/src/finite_elements")
link_directories("${PROJECT_BINARY_DIR}/finite_element_library")

add_library(mofem_finite_element_library
  ${PROJECT_SOURCE_DIR}/finite_element_library/impl/DirichletBC.cpp
  ${PROJECT_SOURCE_DIR}/finite_element_library/impl/ThermalElement.cpp
  ${PROJECT_SOURCE_DIR}/finite_element_library/impl/NonLienarElasticElement.cpp
  ${PROJECT_SOURCE_DIR}/finite_element_library/impl/ConstrainMatrixCtx.cpp
)

