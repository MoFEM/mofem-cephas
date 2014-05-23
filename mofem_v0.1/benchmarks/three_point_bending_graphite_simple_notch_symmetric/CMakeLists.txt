# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
# It can be freely used for educational and research purposes 
# by other institutions. If you use this software pleas cite my work. 
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

cm_export_file("../../meshes/T3_horiz.cub" export_files_3pb_graphite)
cm_export_file("../../meshes/T3_horiz.jou" export_files_3pb_graphite)

install(FILES ${CMAKE_BINARY_DIR}/examples/material_forces/arc_length.sh 
  DESTINATION mofem/benchmarks/three_point_bending_graphite_simple_notch_symmetric
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)
install(FILES ${CMAKE_BINARY_DIR}/examples/material_forces/convergence_study.sh 
  DESTINATION mofem/benchmarks/three_point_bending_graphite_simple_notch_symmetric
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/T3_horiz.cub
  DESTINATION mofem/benchmarks/three_point_bending_graphite_simple_notch_symmetric
  PERMISSIONS OWNER_READ OWNER_WRITE )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/T3_horiz.jou
  DESTINATION mofem/benchmarks/three_point_bending_graphite_simple_notch_symmetric
  PERMISSIONS OWNER_READ OWNER_WRITE )





