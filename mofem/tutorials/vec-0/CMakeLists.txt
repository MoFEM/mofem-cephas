include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

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
tutorials_copy_and_install("*.cub" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.geo" "${permissions_execute}")

tutorials_copy_and_install("../vec-2/beam_2D.cub" "${permissions_default}")
tutorials_copy_and_install("../vec-2/beam_3D.cub" "${permissions_default}")
tutorials_copy_and_install("../vec-2/meshes/LShape_Rot_MPC_link.cub" "${permissions_default}")
tutorials_copy_and_install("../vec-2/meshes/LShape_Rot_MPC_link_master.cub" "${permissions_default}")

tutorials_build_and_install(
  elastic_2d ${CMAKE_CURRENT_SOURCE_DIR}/elastic.cpp)
set_target_properties(
  elastic_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1")

tutorials_build_and_install(
  elastic_3d ${CMAKE_CURRENT_SOURCE_DIR}/elastic.cpp)
set_target_properties(
  elastic_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1")

  tutorials_build_and_install(
    elastic_2d_petsc ${CMAKE_CURRENT_SOURCE_DIR}/elastic.cpp)
  set_target_properties(
    elastic_2d_petsc PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=0")
  
  tutorials_build_and_install(
    elastic_3d_petsc ${CMAKE_CURRENT_SOURCE_DIR}/elastic.cpp)
  set_target_properties(
    elastic_3d_petsc PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=0")  

add_test(
  elastic
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_2d
  -file_name part.h5m -order 4
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -ksp_type fgmres -test -log_no_color)

add_test(
  elastic_mpc_XYZ
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_2d_petsc
  -file_name LShape_Rot_MPC_link.cub -order 2 
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -test -log_no_color -pc_type lu)

add_test(
  elastic_mpc_master_only
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_2d_petsc
  -file_name LShape_Rot_MPC_link_master.cub -order 2 
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -test -log_no_color -pc_type lu)

add_test(
  elastic_3d_rotation0
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_3d
  -file_name LShape_Rot_D0.cub -order 4
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -test -log_no_color)

add_test(
  elastic_3d_rotation_blockset
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_3d
  -file_name LShape_Rot_D1.cub -order 4
  -fieldsplit_1_pc_type cg -fieldsplit_1_pc_type lu 
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -test -log_no_color)

add_test(
  elastic_3d_dam_fluid_pressure
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/elastic_3d
  -file_name dam.cub -order 4
  -fieldsplit_1_pc_type cg -fieldsplit_1_pc_type lu 
  -ksp_rtol 1e-12 -ksp_atol 1e-12 -test -log_no_color)
