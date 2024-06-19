# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
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

tutorials_copy_and_install("*.msh" "${permissions_default}")
tutorials_copy_and_install("*.cfg" "${permissions_default}")
# tutorials_copy_and_install("*.jou" "${permissions_default}")
tutorials_copy_and_install("meshes/*.cub" "${permissions_default}")
tutorials_copy_and_install("*.cub" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.geo" "${permissions_execute}")

tutorials_build_and_install(
  nonlinear_elastic ${CMAKE_CURRENT_SOURCE_DIR}/nonlinear_elastic.cpp)
set_target_properties(
  nonlinear_elastic PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")

tutorials_build_and_install(
  nonlinear_elastic_2d ${CMAKE_CURRENT_SOURCE_DIR}/nonlinear_elastic.cpp)
set_target_properties(
  nonlinear_elastic_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")

add_test(
nonlinear_elastic
${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic
-file_name ${CMAKE_CURRENT_SOURCE_DIR}/beam_3D.cub -order 2
-ts_max_time 0.3 -ts_dt 0.1 -log_no_color -test 1)

add_test(
  test_rotation_0_2d
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_D0_2D.cub 
  -order 2 -ts_max_time 1.0 -ts_dt 0.05
  -log_no_color -test 2)

add_test(
  test_rotation_1_2d
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_D1_2D.cub -order 2 -ts_max_time 1.0 -ts_dt 0.05
  -log_no_color -test 3)

add_test(
  test_nelastic_mpc_links
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d -ts_max_steps 25
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_MPC_link.cub
  -order 2 -ts_max_time 1.0 -ts_dt 0.05 
  -log_no_color -test 4 -pc_type lu)

  # add_test(
  #   test_nelastic_mpc_3_part
  #   ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 3
  #   ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d -ts_max_steps 25
  #   -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_MPC_link_3_part.h5m -order 2 -ts_max_time 1.0 -ts_dt 0.05 
  #   -log_no_color -test 4 -pc_type lu)
  
  add_test(
  test_nelastic_mpc_master
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d -ts_max_steps 25
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_MPC_link_master.cub
  -order 2 -ts_max_time 1.0 -ts_dt 0.05 
  -log_no_color -test 5 -pc_type lu)


add_test(
  test_nelastic_mpc_master_swap
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d -ts_max_steps 25
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_MPC_link_master_swap.cub -order 2 -ts_max_time 1.0 -ts_dt 0.05
  -log_no_color -test 6 -pc_type lu)

add_test(
  test_nelastic_mpc_link_Y
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic_2d -ts_max_steps 25
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/LShape_Rot_MPC_link_Y.cub -order 2 -ts_max_time 1.0 -ts_dt 0.05
  -log_no_color -test 7 -pc_type lu)

add_test(
  test_nelastic_mpc_link_3D_repr
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/nonlinear_elastic -ts_max_steps 25
  -file_name ${CMAKE_CURRENT_SOURCE_DIR}/meshes/beam3D_MPCs.cub -order 1 -ts_dt 0.1 -ts_max_time 1.0 
  -pc_type lu -log_no_color -test 8 -pc_type lu)

