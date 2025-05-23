include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-2/src)
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
tutorials_copy_and_install("*.cub" "${permissions_default}")
tutorials_copy_and_install("*.h5m" "${permissions_default}")
tutorials_copy_and_install("*.petsc" "${permissions_default}")
tutorials_copy_and_install("*.sh" "${permissions_execute}")
tutorials_copy_and_install("*.geo" "${permissions_execute}")
tutorials_copy_and_install("*.png" "${permissions_execute}")
tutorials_copy_and_install("*.md" "${permissions_execute}")
tutorials_copy_and_install("*.cfg" "${permissions_execute}")
tutorials_copy_and_install("*.jou" "${permissions_execute}")
tutorials_copy_and_install("*.py" "${permissions_execute}")
tutorials_copy_and_install("meshes/*" "${permissions_execute}")

if(Boost_NUMPY_LIBRARY_RELEASE) 
  
  tutorials_build_and_install(
    plastic_2d ${CMAKE_CURRENT_SOURCE_DIR}/plastic.cpp)
  tutorials_build_and_install(
    plastic_3d ${CMAKE_CURRENT_SOURCE_DIR}/plastic.cpp)

  set_target_properties(
    plastic_3d PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1 -DPYTHON_SDF -DADD_CONTACT")
  set_target_properties(
    plastic_2d PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1 -DPYTHON_SDF -DADD_CONTACT")

  target_link_libraries(plastic_2d ${Boost_LIBRARIES})
  target_link_libraries(plastic_3d ${Boost_LIBRARIES})

else (Boost_NUMPY_LIBRARY_RELEASE)

  tutorials_build_and_install(
    plastic_2d ${CMAKE_CURRENT_SOURCE_DIR}/plastic.cpp)
  tutorials_build_and_install(
    plastic_3d ${CMAKE_CURRENT_SOURCE_DIR}/plastic.cpp)

  set_target_properties(
    plastic_3d PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1")
  set_target_properties(
    plastic_2d PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1")

endif(Boost_NUMPY_LIBRARY_RELEASE)

add_test(
  plastic_3d_necking_test_operators
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  plastic_3d 
  -file_name necking_3D_p2.h5m
  -ts_dt 0.01
  -ts_max_time 1.2
  -ts_max_steps 1
  -ts_type theta
  -snes_atol 1e-9
  -order 1
  -large_strains 1
  -scale 1
  -hardening_viscous 0
  -cn0 1 
  -zeta 0.1
  -ksp_monitor 
  -test_operators
  -log_no_color)

add_test(
  plastic_3d_necking
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  plastic_3d 
  -file_name necking_3D_p2.h5m
  -ts_dt 0.01
  -ts_max_time 1.2
  -ts_max_steps 1
  -ts_type theta
  -snes_atol 1e-9
  -order 2
  -large_strains 1
  -scale 1
  -hardening_viscous 0
  -cn0 1 -ksp_monitor 
  -log_view
  -log_no_color)

add_test(
  plastic_3d_necking_petsc
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  plastic_3d 
  -file_name necking_3D_p2.h5m
  -ts_dt 0.01
  -ts_max_time 1.2
  -ts_max_steps 1
  -ts_type theta
  -snes_atol 1e-9
  -order 2
  -large_strains 1
  -scale 1
  -hardening_viscous 0
  -cn0 1 -ksp_monitor 
  -pc_type lu
  -log_view
  -log_no_color)

  add_test(
  plastic_2d_uniaxial_tension
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 2
  plastic_2d 
  -file_name coarse_plastic_square_tension_2p.h5m
  -ts_dt 0.02
  -ts_max_time 0.5
  -ts_max_steps 25
  -ts_type theta
  -snes_atol 1e-9
  -order 2
  -large_strains 1
  -scale 1
  -hardening_viscous 0
  -cn0 1 -ksp_monitor 
  -field_eval_coords 0,0,0
  -atom_test 1
  -log_view
  -log_no_color)