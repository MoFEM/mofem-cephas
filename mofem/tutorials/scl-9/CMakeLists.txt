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
tutorials_copy_and_install("*.md" "${permissions_execute}")

tutorials_build_and_install(
  heat_method ${CMAKE_CURRENT_SOURCE_DIR}/heat_method.cpp)

# FIXME Add test
# add_test(
#   heat_method
#   ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
#   ${CMAKE_CURRENT_BINARY_DIR}/heat_method
#   -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh2d.h5m -log_no_color)