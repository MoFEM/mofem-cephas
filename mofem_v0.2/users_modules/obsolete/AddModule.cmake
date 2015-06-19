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

# Obsolete, i.e. implementation is obsolete and need to be changed or replaced
# by alternative classes or functions. Use of oboslete implementation should be
# avoided. This code will be removed in futuer versions of MoFEM

include_directories(${UM_SOURCE_DIR}/obsolete)
include_directories(${UM_SOURCE_DIR}/obsolete/c)
include_directories(${UM_SOURCE_DIR}/obsolete/c_impl)

add_library(users_modules_obsolete
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_SurfaceConstrains.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_LowLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_UpLevelStudent.cpp
  ${UM_SOURCE_DIR}/obsolete/c_impl/complex_for_lazy.c
  ${UM_SOURCE_DIR}/obsolete/impl/FEMethod_ComplexForLazy.cpp
  ${UM_SOURCE_DIR}/obsolete/impl/SurfacePressureComplexForLazy.cpp
)
