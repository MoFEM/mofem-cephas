
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-2/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../adv-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../adv-1/src)

include_directories(${MoFEM_INSTALL_DIR}/third_party/phg-quadrule)

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
tutorials_copy_and_install("*.py" "${permissions_execute}")

tutorials_build_and_install(
  incompressible_2d_schur ${CMAKE_CURRENT_SOURCE_DIR}/incompressible_elasticity.cpp)
tutorials_build_and_install(
  incompressible_3d_schur ${CMAKE_CURRENT_SOURCE_DIR}/incompressible_elasticity.cpp)

set_target_properties(
  incompressible_2d_schur PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1")
set_target_properties(
  incompressible_3d_schur PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1")

tutorials_build_and_install(
  incompressible_2d_petsc ${CMAKE_CURRENT_SOURCE_DIR}/incompressible_elasticity.cpp)
tutorials_build_and_install(
  incompressible_3d_petsc ${CMAKE_CURRENT_SOURCE_DIR}/incompressible_elasticity.cpp)

set_target_properties(
  incompressible_2d_petsc PROPERTIES COMPILE_FLAGS
  "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=0")
set_target_properties(
  incompressible_3d_petsc PROPERTIES COMPILE_FLAGS
  "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=0") 

# add_test(
#   incompressible_2d_schur
#   ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
#   incompressible_2d_schur
#   -file_name ${CMAKE_CURRENT_BINARY_DIR}/sphere_2d.h5m
#   -ts_type beuler
#   -ts_max_time 1 -ts_dt 0.2 
#   -spring_stiffness 0.1
#   -log_no_color)

add_test(
  cookes_membranme_2d
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  incompressible_2d_petsc
  -order 2
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_2d.h5m
  -is_discontinuous_pressure 1
  -ksp_monitor
  -snes_monitor
  -young_modulus 240.5656
  -poisson_ratio 0.4999
  -post_proc_vol 1 
  -post_proc_skin 0
  -eps_stab 0.
  -ts_max_time 1 -ts_dt 0.2
  -is_a_t_petsc_fieldsplit 1
  -pc_type lu
  -field_eval_coords 48,60
  -plane_strain 1 
  -atom_test 1
  )
  