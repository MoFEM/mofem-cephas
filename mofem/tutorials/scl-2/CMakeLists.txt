#

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
  poisson_2d_nonhomogeneous
  ${CMAKE_CURRENT_SOURCE_DIR}/poisson_2d_nonhomogeneous.cpp)

  add_test(
    poisson_2d_nonhomo ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/poisson_2d_nonhomogeneous
      -file_name ${CMAKE_CURRENT_SOURCE_DIR}/mesh2d_2parts.h5m
      -order 2
      -atom_test 1
    )