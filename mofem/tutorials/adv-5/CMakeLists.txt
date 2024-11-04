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

#tutorials_build_and_install(
  #seepage_2d ${CMAKE_CURRENT_SOURCE_DIR}/seepage.cpp)
#set_target_properties(
  #seepage_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")

tutorials_build_and_install(
  seepage_2d_1 ${CMAKE_CURRENT_SOURCE_DIR}/seepage.cpp)
set_target_properties(
  seepage_2d_1 PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")

tutorials_build_and_install(
  seepage_3d ${CMAKE_CURRENT_SOURCE_DIR}/seepage.cpp)
set_target_properties(
  seepage_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")

add_test(
  seepage_2d
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  seepage_2d_1
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/test_seepage.cub
  -ts_type beuler
  -ts_max_time 1 
  -ts_max_steps 2
  -ts_dt 0.2
  -ksp_monitor 
  -young_modulus 0.0009
  -poisson_ratio 0.2
  -conductivity 0.00012
  -biot_constant 1
  -storage 1
  -plane strain
  -order 2 
  -field_eval_coords 0 0 
  -log_sl inform 
  -log_no_color)





