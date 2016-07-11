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

#multi_indices
set(MOFEM_INTERFACES_SRC
  impl/LoopMethods.cpp
  #Core
  impl/Core.cpp
  impl/CreateMat.cpp
  impl/Vectors.cpp
  impl/FieldBlas.cpp
  impl/GetAdjacancies.cpp
  impl/FieldCore.cpp
  impl/FECore.cpp
  impl/ProblemsCore.cpp
  impl/CommCore.cpp
  impl/DeleteCore.cpp
  impl/MeshRefinment.cpp
  impl/PrismInterfaceCore.cpp
  impl/SeriesRecorderCore.cpp
  impl/NodeMerger.cpp
  impl/BitLevelCoupler.cpp
  impl/PrismsFromSurfaceInterface.cpp
)

if(TETGEN_LIBRARY)
  set(MOFEM_INTERFACES_SRC ${MOFEM_INTERFACES_SRC} impl/TetGenInterface.cpp)
endif(TETGEN_LIBRARY)

if(NETGEN_LIBRARY)
  set(MOFEM_INTERFACES_SRC ${MOFEM_INTERFACES_SRC} impl/NetGenInterface.cpp)
endif(NETGEN_LIBRARY)

add_library(mofem_interfaces
  ${MOFEM_INTERFACES_SRC}
)

install(TARGETS mofem_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
