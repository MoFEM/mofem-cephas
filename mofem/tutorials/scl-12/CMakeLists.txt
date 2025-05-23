include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

include_directories(
  ${UM_SOURCE_DIR}/tutorials/cor-2to5/src
)

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

tutorials_copy_and_install("*.jou" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.md" "${permissions_default}")
tutorials_copy_and_install("*.cfg" "${permissions_default}")
tutorials_copy_and_install("*.cub" "${permissions_default}")


tutorials_build_and_install(
    electrostatics_3d ${CMAKE_CURRENT_SOURCE_DIR}/electrostatics.cpp)
  set_target_properties(
    electrostatics_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")
  
tutorials_build_and_install(
    electrostatics_2d ${CMAKE_CURRENT_SOURCE_DIR}/electrostatics.cpp)
  set_target_properties(
    electrostatics_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")

  add_custom_target(
      tut_scl_12_copy_mesh_files_directory
      ALL
      COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files ${CMAKE_CURRENT_BINARY_DIR}/mesh_files
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Copying mesh_files directory" VERBATIM)

#################
##### TESTS #####
#################

# sigma=+-10.0 on the electrode; ep_air:1, ep1: 2, ep2: 3;
  add_test(
    electrostatic_2d
    ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
    electrostatics_2d -order 2 -geom_order 1
    -file_name ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files/elec_2D_4p.h5m
    -atom_test 1)

  add_test(
    electrostatic_3d
    ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
    electrostatics_3d -order 2 -geom_order 1
    -file_name ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files/elec_3D_4p.h5m
    -atom_test 1)

  # sigma=+10.0 on the electrode 1; ep 1; 
  # roughness on electrode 1, Amplitude= 0.01; lamda=1.
  add_test(
    electrostatic_wavy_3d
    ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4 
    electrostatics_3d -order 2 -geom_order 2
    -file_name ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files/elec_wavy_4p.h5m 
    -atom_test 2)