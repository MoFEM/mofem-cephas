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

add_executable(test_moabField test_moabField.cpp)
target_link_libraries(test_moabField mofem_Field ${PROJECT_LIBS})

add_executable(test_mesh_refinment test_mesh_refinment.cpp)
target_link_libraries(test_mesh_refinment mofem_Field ${PROJECT_LIBS})

add_executable(test_approx_function test_approx_function.cpp)
target_link_libraries(test_approx_function mofem_Field ${PROJECT_LIBS})

add_executable(test_elasticity test_elasticity.cpp)
target_link_libraries(test_elasticity mofem_Field ${PROJECT_LIBS})

