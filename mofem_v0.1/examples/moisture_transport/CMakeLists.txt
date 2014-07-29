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


add_executable(moisture_trans moisture_trans.cpp)
target_link_libraries(moisture_trans mofem_student_body mofem_forces_and_sources mofem_Field ${PROJECT_LIBS})


add_executable(RVE_moisture_trans_disp RVE_moisture_trans_disp.cpp)
target_link_libraries(RVE_moisture_trans_disp mofem_student_body mofem_forces_and_sources mofem_homogenisation mofem_Field ${PROJECT_LIBS})

add_executable(RVE_moisture_trans_trac RVE_moisture_trans_trac.cpp)
target_link_libraries(RVE_moisture_trans_trac mofem_student_body mofem_forces_and_sources mofem_homogenisation mofem_Field ${PROJECT_LIBS})
