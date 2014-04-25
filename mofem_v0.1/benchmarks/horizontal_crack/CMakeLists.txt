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

cm_export_file("../../meshes/radwanska_benchmark.cub" export_files_plate_stability)
cm_export_file("../../meshes/randwanska_benchmark.jou" export_files_plate_stability)
cm_export_file("../../meshes/radwanska_benchmark2.cub" export_files_plate_stability)
cm_export_file("../../meshes/radwanska_benchmark2.jou" export_files_plate_stability)

install(FILES ${CMAKE_BINARY_DIR}/examples/arc_length_nonlinear_elasticity/arc_length_nonlinear_elasticity.sh
  DESTINATION mofem/benchmarks/plate_stability
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/radwanska_benchmark.cub
  DESTINATION mofem/benchmarks/plate_stability
  PERMISSIONS OWNER_READ OWNER_WRITE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/radwanska_benchmark2.cub
  DESTINATION mofem/benchmarks/plate_stability
  PERMISSIONS OWNER_READ OWNER_WRITE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/randwanska_benchmark.jou
  DESTINATION mofem/benchmarks/plate_stability
  PERMISSIONS OWNER_READ OWNER_WRITE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/radwanska_benchmark2.jou
  DESTINATION mofem/benchmarks/plate_stability
  PERMISSIONS OWNER_READ OWNER_WRITE)




