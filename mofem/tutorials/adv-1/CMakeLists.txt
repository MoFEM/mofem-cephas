
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-0/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../vec-2/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../adv-0/src)

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

set(BOOST_ROOT ${MoFEM_BOOST_DIR})
message(STATUS ${BOOST_ROOT})
find_package(
  Boost
  OPTIONAL_COMPONENTS
  numpy
)

if(Boost_NUMPY_LIBRARY_RELEASE) 

  message(STATUS "Boost numpy found - SDF in pythoon enabled")

  tutorials_build_and_install(
    contact_2d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  tutorials_build_and_install(
    contact_2d_schur ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d_schur ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  set_target_properties(
    contact_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DPYTHON_SDF") 
  set_target_properties(
    contact_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DPYTHON_SDF") 

  set_target_properties(
    contact_2d_schur PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1 -DPYTHON_SDF") 
  set_target_properties(
    contact_3d_schur PROPERTIES COMPILE_FLAGS 
    "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1 -DPYTHON_SDF") 

  target_link_libraries(contact_2d ${Boost_LIBRARIES})
  target_link_libraries(contact_3d ${Boost_LIBRARIES})

  target_link_libraries(contact_2d_schur ${Boost_LIBRARIES})
  target_link_libraries(contact_3d_schur ${Boost_LIBRARIES})

else (Boost_NUMPY_LIBRARY_RELEASE)

  message(STATUS "Boost numpy NOT found - SDF is hardcodes in C++")

  tutorials_build_and_install(
    contact_2d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  tutorials_build_and_install(
    contact_2d_schur ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d_schur ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  set_target_properties(
    contact_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")
  set_target_properties(
    contact_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")  
  
  set_target_properties(
    contact_2d_schur PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DSCHUR_ASSEMBLE=1")
  set_target_properties(
    contact_3d_schur PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DSCHUR_ASSEMBLE=1")  

endif(Boost_NUMPY_LIBRARY_RELEASE)

# set(IS_MFRONT_MOD false)  

# Add mfront interface
# if(MFRONT_INTERFACE_INSTALL_PREFIX)
#   include(${MFRONT_INTERFACE_INSTALL_PREFIX}/mfront_interface/mfront_interface_all_targets.cmake)
#   include_directories(${MFRONT_INTERFACE_INSTALL_PREFIX}/mfront_interface/src)
#   set(IS_MFRONT_MOD true)  
#   else(MFRONT_INTERFACE_INSTALL_PREFIX)
# if(EXISTS ${PROJECT_SOURCE_DIR}/users_modules/mfront_interface/src)
#   include_directories(${PROJECT_SOURCE_DIR}/users_modules/mfront_interface/src)
#   set(IS_MFRONT_MOD true)  
#   endif()
# endif(MFRONT_INTERFACE_INSTALL_PREFIX)

# if(IS_MFRONT_MOD)
#   message("WITH MFRONT")
#   target_link_libraries(contact_2d
#     mfront_lib)
#   target_link_libraries(contact_3d
#     mfront_lib)
#   target_link_libraries(contact_2d_schur
#     mfront_schur_lib)
#   target_link_libraries(contact_3d_schur
#     mfront_schur_lib)
# endif(IS_MFRONT_MOD)

# Copy scripts to build a directory
add_custom_target(
  um_copy_sdf_files_directory
  ALL
  COMMAND ${CMAKE_COMMAND} -E copy_directory 
  ${CMAKE_CURRENT_SOURCE_DIR}/sdf_files ${CMAKE_CURRENT_BINARY_DIR}/sdf_files
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Copying sdf_files directory" VERBATIM)

add_custom_target(
    um_copy_mesh_files_directory
    ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory 
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh_files ${CMAKE_CURRENT_BINARY_DIR}/mesh_files
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying mesh_files directory" VERBATIM)

#################
##### TESTS #####
#################

if(Boost_NUMPY_LIBRARY_RELEASE) 

add_test(
  hertz_2d_plane_stress
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_2d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/cyl_2d_4p.h5m
  -pc_type lu
  -ts_type beuler
  -ts_max_time 1 -ts_dt 1
  -sdf_file ./sdf_files/sdf_hertz_2d_plane.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -atom_test 1
)

add_test(
  hertz_2d_plane_stress_schur
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_2d_schur
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/cyl_2d_4p.h5m
  -pc_type fieldsplit
  -mofem_mg_levels 2
  -order 2 
  -eps_stab 1e-4
  -ts_type beuler
  -ts_max_time 1 
  -ts_max_steps 10
  -ts_dt 1
  -sdf_file ./sdf_files/sdf_hertz_2d_plane.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -is_large_strain 1
  -atom_test 1
)

add_test(
  hertz_2d_plane_stress_schur_hookes ## Tests for hookeOps with is_large_strain = 0
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_2d_schur
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/cyl_2d_4p.h5m
  -pc_type fieldsplit
  -mofem_mg_levels 2
  -order 2 
  -eps_stab 1e-4
  -ts_type beuler
  -ts_max_time 1 
  -ts_max_steps 10
  -ts_dt 1
  -sdf_file ./sdf_files/sdf_hertz_2d_plane.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -is_large_strain 0
  -atom_test 1
)

add_test(
  hertz_3D
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_3d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_3d_4p.h5m
  -pc_type lu
  -ts_type beuler
  -ts_max_time 1 
  -ts_dt 1
  -sdf_file ./sdf_files/sdf_hertz_3d.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -atom_test 3
)

add_test(
  hertz_3d_schur
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_3d_schur
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_3d_4p.h5m
  -snes_max_it 40
  -ksp_atol 1e-10
  -ksp_rtol 1e-8
  -pc_type fieldsplit
  -fieldsplit_1_ksp_type gmres
  -fieldsplit_1_ksp_rtol 1e-2
  -fieldsplit_1_pc_type mg
  -mofem_mg_levels 2
  -order 2 
  -eps_stab 1e-5
  -ts_type beuler
  -ts_max_time 1 
  -ts_max_steps 10
  -ts_dt 1
  -sdf_file ./sdf_files/sdf_hertz_3d.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -atom_test 3
)

add_test(
  wavy_3D
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_3d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/wavy_3d_4p.h5m
  -is_large_strain 1
  -base demkowicz
  -pc_type lu
  -ts_type beuler
  -ts_max_time 1 -ts_dt 1
  -sdf_file ./sdf_files/sdf_plane_3d.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -atom_test 7
)

endif(Boost_NUMPY_LIBRARY_RELEASE)

add_test(
  wavy_3d_hookes
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
  contact_3d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/wavy_3d_4p.h5m
  -is_large_strain 0
  -base demkowicz
  -pc_type lu
  -ts_type beuler
  -ts_max_time 1 -ts_dt 1
  -sdf_file ./sdf_files/sdf_plane_3d.py
  -young_modulus 1e3
  -poisson_ratio 0.4
  -atom_test 7
)
# if(IS_MFRONT_MOD)
#   if(TFEL_DIR)
#     message("TFEL DIR FOUND")

#     add_test(adv_1_compile_example_behaviour_test
#       ${TFEL_DIR}/bin/mfront --obuild --interface=generic ${PROJECT_SOURCE_DIR}/mfront_interface/behaviours/LinearElasticity.mfront ${PROJECT_SOURCE_DIR}/mfront_interface/behaviours/SignoriniHyperElasticity.mfront)

#     if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#       set(DEFAULT_LIB_EXTENSION "dylib")
#     else()
#       set(DEFAULT_LIB_EXTENSION "so")
#     endif()

#     add_test(
#       wavy_2d_mfront
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_2d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/wavy_2d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_plane_2d.py
#       -use_mfront 1
#       -is_axisymmetric 0
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 6
#     )

#     add_test(
#       hertz_mfront_2d_plane_strain
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_2d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/cyl_2d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_hertz_2d_plane.py
#       -use_mfront 1
#       -is_axisymmetric 0
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 2
#     )

#     add_test(
#       hertz_mfront_2d_axisymm_flat
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_2d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_2d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1 
#       -sdf_file ./sdf_files/sdf_hertz_2d_axisymm_sphere.py
#       -use_mfront 1
#       -is_axisymmetric 1
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 5
#     )

#     add_test(
#       hertz_mfront_2d_axisymm_sphere
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_2d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/cyl_2d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1 
#       -sdf_file ./sdf_files/sdf_hertz_2d_axisymm_plane.py
#       -use_mfront 1
#       -is_axisymmetric 1
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 4
#     )

#     add_test(
#       hertz_mfront_3d
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_3d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_3d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_hertz_3d.py
#       -use_mfront 1
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 3
#     )

#     add_test(
#       wavy_mfront_half_space_2d
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_2d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/half_space_2d_4p.h5m
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_wavy_2d.py
#       -use_mfront 1
#       -is_axisymmetric 0
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity     
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 6
#     )
    
#     add_test(
#       hertz_mfront_3d_schur
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_3d_schur
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_3d_4p.h5m
#       -pc_type fieldsplit
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_hertz_3d.py
#       -use_mfront 1
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 3
#     )


#     add_test(
#       wavy_3d_mfront #hex
#       ${MPI_RUN} ${MPI_RUN_FLAGS} -np 4
#       contact_3d
#       -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/wavy_3d_4p.h5m
#       -base demkowicz
#       -pc_type lu
#       -ts_type beuler
#       -ts_max_time 1 -ts_dt 1
#       -sdf_file ./sdf_files/sdf_plane_3d.py
#       -use_mfront 1
#       -is_axisymmetric 0
#       -mi_lib_path_10 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
#       -mi_block_10 LinearElasticity
#       -mi_param_10_0 1e3
#       -mi_param_10_1 0.4
#       -atom_test 7
#     )

#   endif(TFEL_DIR)
# endif()