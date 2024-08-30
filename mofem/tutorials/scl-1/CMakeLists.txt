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

tutorials_copy_and_install("*.jou" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.md" "${permissions_default}")
tutorials_copy_and_install("*.cub" "${permissions_default}")

tutorials_build_and_install(
  poisson_2d_homogeneous ${CMAKE_CURRENT_SOURCE_DIR}/poisson_2d_homogeneous.cpp)
set_target_properties(
  poisson_2d_homogeneous PROPERTIES COMPILE_FLAGS 
  "-DEXECUTABLE_DIMENSION=2")

tutorials_build_and_install(
  poisson_3d_homogeneous ${CMAKE_CURRENT_SOURCE_DIR}/poisson_2d_homogeneous.cpp)
set_target_properties(
  poisson_3d_homogeneous PROPERTIES COMPILE_FLAGS 
  "-DEXECUTABLE_DIMENSION=3")

  add_test(
    poisson_2d_homo ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/poisson_2d_homogeneous
    -file_name ${CMAKE_CURRENT_SOURCE_DIR}/square_2parts.h5m 
    -order 2
    -atom_test 1
  )