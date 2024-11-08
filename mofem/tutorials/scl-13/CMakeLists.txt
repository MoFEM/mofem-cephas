include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../scl-1/src)

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
tutorials_copy_and_install("*.cub" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.md" "${permissions_default}")
tutorials_copy_and_install("*.ipynb" "${permissions_default}")
tutorials_copy_and_install("meshes/*.cub" "${permissions_default}")
tutorials_copy_and_install("meshes/*.stl" "${permissions_default}")

tutorials_build_and_install(
  poisson_cutfem ${CMAKE_CURRENT_SOURCE_DIR}/poisson_3d_cutfem.cpp)

add_test(
  poisson_cutfem_test
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 3
  poisson_cutfem
  -file_name mesh_dense_p3.h5m
  -order 2
  -penalty 1
  -phi -1 
  -nitsche 1 
  -log_sl inform
  -is_test)
