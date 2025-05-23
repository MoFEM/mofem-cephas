#

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)
# include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)
include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/cblas)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/post_proc)
include_directories(${PROJECT_SOURCE_DIR}/src/boundary_conditions)

include_directories(${PROJECT_SOURCE_DIR}/tools)
include_directories(${PROJECT_SOURCE_DIR}/tools/src)

set(OPTIONAL_PROJECT_LIBS
  ${TETGEN_LIBRARY}
  ${ADOL-C_LIBRARY}
  ${MED_LIBRARY}
)

set(MoFEM_PROJECT_LIBS 
  ${OPTIONAL_PROJECT_LIBS} ${MoFEM_SYSTEM_PROJECT_LIBS}  ${PROJECT_LIBS})
message(STATUS ${MoFEM_PROJECT_LIBS})


function(tools_copy_and_install extension file_permissions)
  file(GLOB files ${CMAKE_CURRENT_SOURCE_DIR}/${extension})
  foreach(file ${files})
    file(
      COPY ${file}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      FILE_PERMISSIONS ${file_permissions})     
    install(FILES ${file}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      PERMISSIONS ${file_permissions})
  endforeach(file)
endfunction()

set(permissions_default 
  OWNER_WRITE 
  OWNER_READ
  GROUP_READ)

file(
  COPY meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(
  DIRECTORY meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(MED_LIBRARY)
  add_executable(read_med ${CMAKE_CURRENT_SOURCE_DIR}/read_med.cpp)
  target_link_libraries(read_med
    mofem_interfaces
    mofem_finite_elements
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    ${MoFEM_PROJECT_LIBS}
  )
  install(TARGETS read_med DESTINATION bin)
  
  add_executable(write_med ${CMAKE_CURRENT_SOURCE_DIR}/write_med.cpp)
  target_link_libraries(write_med
    mofem_interfaces  
    mofem_finite_elements
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    ${MoFEM_PROJECT_LIBS}
  )
  install(TARGETS write_med DESTINATION bin)
endif(MED_LIBRARY)

add_executable(add_meshsets ${CMAKE_CURRENT_SOURCE_DIR}/add_meshsets.cpp)
target_link_libraries(add_meshsets
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS add_meshsets DESTINATION bin)

add_executable(mofem_part ${CMAKE_CURRENT_SOURCE_DIR}/mofem_part.cpp)
target_link_libraries(mofem_part
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS mofem_part DESTINATION bin)

add_executable(split_sideset ${CMAKE_CURRENT_SOURCE_DIR}/split_sideset.cpp)
target_link_libraries(split_sideset
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS split_sideset DESTINATION bin)

add_executable(mesh_cut ${CMAKE_CURRENT_SOURCE_DIR}/mesh_cut.cpp)
target_link_libraries(mesh_cut
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS mesh_cut DESTINATION bin)

add_executable(field_to_vertices ${CMAKE_CURRENT_SOURCE_DIR}/field_to_vertices.cpp)
target_link_libraries(field_to_vertices
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS field_to_vertices DESTINATION bin)

add_executable(remove_mofem_meshsets
  ${CMAKE_CURRENT_SOURCE_DIR}/remove_mofem_meshsets.cpp)
target_link_libraries(remove_mofem_meshsets
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS remove_mofem_meshsets DESTINATION bin)

add_executable(delete_ho_nodes ${CMAKE_CURRENT_SOURCE_DIR}/delete_ho_nodes.cpp)
target_link_libraries(delete_ho_nodes
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS delete_ho_nodes DESTINATION bin)

add_executable(meshset_to_vtk ${CMAKE_CURRENT_SOURCE_DIR}/meshset_to_vtk.cpp)
target_link_libraries(meshset_to_vtk
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS meshset_to_vtk DESTINATION bin)

add_executable(uniform_mesh_refinement 
  ${CMAKE_CURRENT_SOURCE_DIR}/uniform_mesh_refinement.cpp)
target_link_libraries(uniform_mesh_refinement
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS uniform_mesh_refinement DESTINATION bin)

add_executable(extrude_prisms ${CMAKE_CURRENT_SOURCE_DIR}/extrude_prisms.cpp)
target_link_libraries(extrude_prisms
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS extrude_prisms DESTINATION bin)

add_executable(get_skin ${CMAKE_CURRENT_SOURCE_DIR}/get_skin.cpp)
target_link_libraries(get_skin
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS get_skin DESTINATION bin)

add_executable(fix_file ${CMAKE_CURRENT_SOURCE_DIR}/fix_file.cpp)
target_link_libraries(fix_file
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS fix_file DESTINATION bin)

add_executable(read_vtk ${CMAKE_CURRENT_SOURCE_DIR}/read_vtk.cpp)
target_link_libraries(read_vtk
  mofem_interfaces
  mofem_finite_elements
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_third_party
  ${MoFEM_PROJECT_LIBS}
)
install(TARGETS read_vtk DESTINATION bin)

cm_export_file("convert.py" export_tools)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/convert.py
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
  GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

cm_export_file("make_mesh_png.py" export_tools)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/make_mesh_png.py
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE 
  GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ADOL-C_LIBRARY)

  include_directories(
    ${PROJECT_SOURCE_DIR}/basic_finite_elements/nonlinear_elastic_materials/src
  )

  # FIXME:  
  # add_executable(mesh_smoothing ${CMAKE_CURRENT_SOURCE_DIR}/mesh_smoothing.cpp)
  # target_link_libraries(mesh_smoothing
  #   # users_modules
  #   mofem_interfaces
  #   mofem_finite_elements
  #   mofem_multi_indices
  #   mofem_petsc
  #   mofem_approx
  #   mofem_third_party
  #   ${MoFEM_PROJECT_LIBS}
  # )
  # install(TARGETS mesh_smoothing DESTINATION bin)

endif(ADOL-C_LIBRARY)


add_test(read_vtk
  ${MPI_RUN} ${MPI_RUN_FLAGS} -np 1
  ${CMAKE_CURRENT_BINARY_DIR}/read_vtk
  -file_name ${CMAKE_CURRENT_BINARY_DIR}/meshes/box_layer_heter.vtk
  -log_no_color
  -output_file out_test_read_vtk.h5m
  -atom_test 1
)