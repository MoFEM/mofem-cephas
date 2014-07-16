# Copyright (C) 2013, Zahur Ullah (Zahur.Ullah AT Glasgow.ac.uk)
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

add_executable(RVE_Cube_Disp RVE_Cube_Disp.cpp)
target_link_libraries(RVE_Cube_Disp mofem_student_body mofem_homogenisation mofem_Field ${PROJECT_LIBS})

add_executable(RVE_Cube_Traction RVE_Cube_Traction.cpp)
target_link_libraries(RVE_Cube_Traction mofem_student_body mofem_homogenisation mofem_Field ${PROJECT_LIBS})

add_executable(RVE_Cube_Periodic RVE_Cube_Periodic.cpp)
target_link_libraries(RVE_Cube_Periodic mofem_student_body mofem_homogenisation mofem_Field ${PROJECT_LIBS})

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


#add_executable(RVE_Cube_Disp_transoport RVE_Cube_Disp_transoport.cpp)
#target_link_libraries(RVE_Cube_Disp_transoport mofem_student_body mofem_homogenisation mofem_Field ${PROJECT_LIBS})


