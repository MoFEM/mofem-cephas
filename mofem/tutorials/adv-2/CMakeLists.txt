include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../adv-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-2/src)

set(permissions_default 
  OWNER_WRITE 
  OWNER_READ
  GROUP_READ)

set(permissions_execute 
  OWNER_WRITE 
  OWNER_READ
  GROUP_READ
  OWNER_EXECUTE 
  GROUP_EXECUTE)

add_custom_target(
  adv_2_copy_mesh_files_directory
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files ${CMAKE_CURRENT_BINARY_DIR}/mesh_files
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Copying mesh_files directory" VERBATIM)

tutorials_copy_and_install("*.msh" "${permissions_default}")
tutorials_copy_and_install("*.cfg" "${permissions_default}")
# tutorials_copy_and_install("*.jou" "${permissions_default}")
tutorials_copy_and_install("*.cub" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.geo" "${permissions_execute}")
tutorials_copy_and_install("*.png" "${permissions_execute}")
tutorials_copy_and_install("*.md" "${permissions_execute}")
tutorials_copy_and_install("*.cfg" "${permissions_execute}")
tutorials_copy_and_install("*.jou" "${permissions_execute}")
tutorials_copy_and_install("meshes/*" "${permissions_execute}")

tutorials_build_and_install(
  thermo_elastic_2d ${CMAKE_CURRENT_SOURCE_DIR}/thermo_elastic.cpp)
set_target_properties(
  thermo_elastic_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")

tutorials_build_and_install(
  thermo_elastic_3d ${CMAKE_CURRENT_SOURCE_DIR}/thermo_elastic.cpp)
set_target_properties(
  thermo_elastic_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")

#################
##### TESTS #####
#################

### atom_test = 1 : 3D isotropic thermal expansion != 2D plane strain
### atom_test = 2 : 3D no thermal expansion in Z = 2D plane strain
### atom_test = 2 : 2D plane stress

add_test(
  convection_and_radiation_3d_exp_iso
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  thermo_elastic_3d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/conv_rad_3d_exp_iso_4p.h5m
  -order 2
  -pc_type lu
  -pc_factor_mat_solver_type mumps  
  -ts_type beuler
  -ts_max_time 10 
  -ts_dt 1
  -init_temp 293 
  -ref_temp 293
  -field_eval_coords 0.75,0.75,0.25
  -atom_test 1
)

add_test(
  convection_and_radiation_3d_exp_no_z
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  thermo_elastic_3d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/conv_rad_3d_exp_no_z_4p.h5m
  -order 2
  -pc_type lu
  -pc_factor_mat_solver_type mumps  
  -ts_type beuler
  -ts_max_time 10 
  -ts_dt 1
  -init_temp 293 
  -ref_temp 293
  -field_eval_coords 0.75,0.75,0.25
  -atom_test 2
)

add_test(
  convection_and_radiation_2d_plane_strain
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  thermo_elastic_2d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/conv_rad_2d_4p.h5m
  -order 2
  -pc_type lu
  -pc_factor_mat_solver_type mumps  
  -ts_type beuler
  -ts_max_time 10 
  -ts_dt 1
  -init_temp 293 
  -ref_temp 293
  -field_eval_coords 0.75,0.75
  -plane_strain
  -atom_test 2
)

add_test(
  convection_and_radiation_2d_plane_stress
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  thermo_elastic_2d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/conv_rad_2d_4p.h5m
  -order 2
  -pc_type lu
  -pc_factor_mat_solver_type mumps  
  -ts_type beuler
  -ts_max_time 10 
  -ts_dt 1
  -init_temp 293 
  -ref_temp 293
  -field_eval_coords 0.75,0.75
  -atom_test 3
)