
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

# cmake_policy(SET CMP0144 NEW)
set(BOOST_ROOT ${MoFEM_BOOST_DIR})
message(STATUS ${BOOST_ROOT})
find_package(
  Boost
  OPTIONAL_COMPONENTS
  numpy
)

if(Python_LIBRARY AND Boost_NUMPY_LIBRARY_RELEASE) 

  include_directories(${_Python_INCLUDE_DIR})

  tutorials_build_and_install(
    contact_2d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  set_target_properties(
    contact_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2 -DPYTHON_SFD") 
  set_target_properties(
    contact_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3 -DPYTHON_SFD") 

  target_link_libraries(contact_2d ${Python_LIBRARY} ${Boost_LIBRARIES})
  target_link_libraries(contact_3d ${Python_LIBRARY} ${Boost_LIBRARIES})

else (Python_LIBRARY AND Boost_NUMPY_LIBRARY_RELEASE)

  tutorials_build_and_install(
    contact_2d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)
  tutorials_build_and_install(
    contact_3d ${CMAKE_CURRENT_SOURCE_DIR}/contact.cpp)

  set_target_properties(
    contact_2d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=2")
  set_target_properties(
    contact_3d PROPERTIES COMPILE_FLAGS "-DEXECUTABLE_DIMENSION=3")  

endif(Python_LIBRARY AND Boost_NUMPY_LIBRARY_RELEASE)

set(IS_MFRONT_MOD false)  

# Add mfront interface
if(MFRONT_INTERFACE_INSTALL_PREFIX)
  include(${MFRONT_INTERFACE_INSTALL_PREFIX}/mfront_interface/mfront_interface_all_targets.cmake)
  include_directories(${MFRONT_INTERFACE_INSTALL_PREFIX}/mfront_interface/src)
  set(IS_MFRONT_MOD true)  
  else(MFRONT_INTERFACE_INSTALL_PREFIX)
if(EXISTS ${PROJECT_SOURCE_DIR}/mfront_interface/src)
  include_directories(${PROJECT_SOURCE_DIR}/mfront_interface/src)
  set(IS_MFRONT_MOD true)  
  endif()
endif(MFRONT_INTERFACE_INSTALL_PREFIX)

if(IS_MFRONT_MOD)
  message("WITH MFRONT")
  target_link_libraries(contact_2d
    mfront_lib_all)
  target_link_libraries(contact_3d
    mfront_lib_all)
endif(IS_MFRONT_MOD)

# Copy scripts to build directory
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
  
add_test(
  contact_2d
  ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  contact_2d
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/sphere_2d.h5m
  -ts_type beuler
  -ts_max_time 1 -ts_dt 0.2 
  -spring_stiffness 0.1
  -log_no_color
  -sdf_file ${CMAKE_CURRENT_BINARY_DIR}/sdf_files/sdf_ydirection.py
)

if(IS_MFRONT_MOD)
  if(TFEL_DIR)
    message("TFEL DIR FOUND")

    add_test(adv_1_compile_example_behaviour_test
      ${TFEL_DIR}/bin/mfront --obuild --interface=generic ${PROJECT_SOURCE_DIR}/mfront_interface/behaviours/LinearElasticity.mfront)

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      target_compile_definitions(mfront_lib_all PUBLIC DEFAULT_LIB_EXTENSION="dylib")
      set(DEFAULT_LIB_EXTENSION "dylib")
    else()
      target_compile_definitions(mfront_lib_all PUBLIC DEFAULT_LIB_EXTENSION="so")
      set(DEFAULT_LIB_EXTENSION "so")
    endif()

    add_test(
      hertz_mfront_axisymm
      ${MoFEM_MPI_RUN} ${MPI_RUN_FLAGS} -np 4
      contact_2d
      -file_name ${CMAKE_CURRENT_BINARY_DIR}/mesh_files/hertz_axisymm_4p.h5m
      -ts_max_time 1 -ts_dt 0.2 
      -sdf_file ./sdf_files/sdf_hertz.py
      -use_mfront 1
      -is_axisymmetric 1
      -mi_lib_path_2 ${CMAKE_CURRENT_BINARY_DIR}/src/libBehaviour.${DEFAULT_LIB_EXTENSION}
      -mi_block_2 LinearElasticity
      -mi_param_2_0 1e3
      -mi_param_2_1 0.4
      -atom_test 1
    )

  endif(TFEL_DIR)
endif()