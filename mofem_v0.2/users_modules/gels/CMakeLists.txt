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

include_directories(${UM_SOURCE_DIR}/gels/src)
include_directories(${UM_SOURCE_DIR}/nonlinear_elasticity/src)
include_directories(${UM_SOURCE_DIR}/nonlinear_elastic_materials/src)

cm_export_file("meshes/brick_of_gel.cub" export_files_gels)
cm_export_file("meshes/brick_of_gel.jou" export_files_gels)
cm_export_file("gel_config.in" export_files_gels)

add_executable(gel_analysis
  ${UM_SOURCE_DIR}/gels/gel_analysis.cpp
)

target_link_libraries(gel_analysis
  users_modules
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  mofem_cblas
  boost_program_options
  ${MoFEM_PROJECT_LIBS}
)
