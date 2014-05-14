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

# ================= Potential flow for fibre directions ==================================
add_executable(potential_flow_all potential_flow_all.cpp)
target_link_libraries(potential_flow_all mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Linear displacement BCs ========================================
add_executable(RVE_disp RVE_disp.cpp)
target_link_libraries(RVE_disp mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Traction BCs ========================================
add_executable(RVE_trac RVE_trac.cpp)
target_link_libraries(RVE_trac mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Periodic BCs ========================================
add_executable(RVE_periodic RVE_periodic.cpp)
target_link_libraries(RVE_periodic mofem_student_body mofem_Field ${PROJECT_LIBS})




# ================= Potential flow for fibre directions with no inteface elements ==================================
add_executable(potential_flow_all_no_interface potential_flow_all_no_interface.cpp)
target_link_libraries(potential_flow_all_no_interface mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Linear displacement BCs with no interface element ========================================
add_executable(RVE_disp_no_interface RVE_disp_no_interface.cpp)
target_link_libraries(RVE_disp_no_interface mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Traction BCs with no interface element ========================================
add_executable(RVE_trac_no_interface RVE_trac_no_interface.cpp)
target_link_libraries(RVE_trac_no_interface mofem_student_body mofem_Field ${PROJECT_LIBS})

# ======================= Periodic BCs with no interface element ========================================
add_executable(RVE_periodic_no_interface RVE_periodic_no_interface.cpp)
target_link_libraries(RVE_periodic_no_interface mofem_student_body mofem_Field ${PROJECT_LIBS})
