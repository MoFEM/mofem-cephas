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

#multi_indices
add_library(mofem_interfaces
  impl/LoopMethods.cpp
  #Core
  impl/Core.cpp
  impl/FieldInterfaceCore.cpp
  impl/MeshRefinmentCore.cpp
  impl/PrismInterfaceCore.cpp
  impl/SeriesRecorderCore.cpp
)


#install(FILES   
#  Core.hpp
#  FieldInterface.hpp
#  FieldUnknownInterface.hpp
#  LoopMethods.hpp
#  MeshRefinment.hpp
#  PrismInterface.hpp
#  SeriesRecorder.hpp
#  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS mofem_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

