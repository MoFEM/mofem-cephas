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

add_executable(spatial_problem spatial_problem.cpp)
target_link_libraries(spatial_problem mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(material_force material_force.cpp)
target_link_libraries(material_force mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(material_problem material_problem.cpp)
target_link_libraries(material_problem mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(arc_length_material_coupled arc_length_material_coupled.cpp)
target_link_libraries(arc_length_material_coupled mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(spatial_problem_with_thermal_field spatial_problem_with_thermal_field.cpp)
target_link_libraries(spatial_problem_with_thermal_field mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

add_executable(split_faces split_faces.cpp)
target_link_libraries(split_faces mofem_complex_for_lazy mofem_face_splitting mofem_student_body mofem_Field ${PROJECT_LIBS})

file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/README
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

install(TARGETS spatial_problem DESTINATION bin)
install(TARGETS material_force DESTINATION bin)
install(TARGETS arc_length_material_coupled DESTINATION bin)

cm_export_file("../../meshes/crack1.cub" export_files_material_mechanics)
cm_export_file("../../meshes/T3_horiz.cub" export_files_material_mechanics)
cm_export_file("../../meshes/brick_slice_thermal.cub" export_files_material_mechanics)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/convergence_study.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/convergence_study.sh"
  @ONLY
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/arc_length.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/arc_length.sh"
  @ONLY
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/arc_length.sh DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.sh DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE )

#testing
add_test(ConvergenceStudy_sh ${CMAKE_CURRENT_BINARY_DIR}/convergence_study.sh -f crack1.cub -l 1e-4 -o 2 -r 1 -d TEST_SIF)
add_test(ArcLength_sh ${CMAKE_CURRENT_BINARY_DIR}/arc_length.sh -f TEST_SIF/out_spatial_2_1.h5m 
  -a 2e-3 -n 1 -g 1e-6 -d TEST_ARC --snes_linesearch_type=basic --my_nb_sub_steps=1)

