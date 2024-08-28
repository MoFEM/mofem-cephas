
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/tutorials/src)

include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)

include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
link_directories(${PROJECT_BINARY_DIR}/src/approximation)

include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
link_directories(${PROJECT_BINARY_DIR}/src/multi_indices)

include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
link_directories(${PROJECT_BINARY_DIR}/src/interfaces)

include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
link_directories(${PROJECT_BINARY_DIR}/src/petsc)

include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
link_directories(${PROJECT_BINARY_DIR}/src/finite_elements)

include_directories(${PROJECT_SOURCE_DIR}/src/post_proc)
include_directories(${PROJECT_SOURCE_DIR}/src/boundary_conditions)

set(OPTIONAL_PROJECT_LIBS
  ${TETGEN_LIBRARY}
  ${ADOL-C_LIBRARY}
  ${MED_LIBRARY}
)
set(MoFEM_PROJECT_LIBS 
  ${OPTIONAL_PROJECT_LIBS} ${MoFEM_SYSTEM_PROJECT_LIBS} ${PROJECT_LIBS})
message(STATUS ${MoFEM_PROJECT_LIBS})


set(TUTORIALS_INSTALL_DIR "tutorials") 

# Install the headers
install(DIRECTORY "${PROJECT_SOURCE_DIR}/${TUTORIALS_INSTALL_DIR}/"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${TUTORIALS_INSTALL_DIR}"
        FILES_MATCHING PATTERN "*.hpp")

function(tutorials_add_executable target source)
  if(MoFEM_PRECOMPILED_HEADRES)
    set_source_files_properties(${source} 
    PROPERTIES COMPILE_FLAGS "-include ${PERCOMPILED_HEADER}")
  endif(MoFEM_PRECOMPILED_HEADRES)
  add_executable(${target} ${source})
endfunction(tutorials_add_executable)

function(tutorials_install tutorials_mod_dir)
  install(
    DIRECTORY
    ${TUTORIALS_INSTALL_DIR} 
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}/tutorials/${TUTORIALS_MOD_DIR}
    FILES_MATCHING
    PATTERN "README"
    PATTERN "*.h5m"
    PATTERN "*.cub"
    PATTERN "*.jou"
    PATTERN "*.msh"
    PATTERN "*.cfg"
    PATTERN "*.geo"
    PATTERN "src" EXCLUDE
    PATTERN "doc" EXCLUDE
    PATTERN "atom_tests" EXCLUDE)
endfunction(tutorials_install)

function(tutorials_build_and_install target source)
  tutorials_add_executable(${target} ${source})
  target_link_libraries(${target}
    # users_modules
    mofem_finite_elements
    mofem_interfaces
    mofem_multi_indices
    mofem_petsc
    mofem_approx
    mofem_third_party
    mofem_matrix_function
    mofem_post_proc
    mofem_boundary_conditions
    ${MoFEM_PROJECT_LIBS})
  install(TARGETS ${target} 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/tutorials/${TUTORIALS_MOD_DIR})

endfunction(tutorials_build_and_install)

function(tutorials_copy_and_install extension file_permissions)
  file(GLOB files ${CMAKE_CURRENT_SOURCE_DIR}/${extension})
  foreach(file ${files})
    file(
      COPY ${file}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      FILE_PERMISSIONS ${file_permissions})     
    install(FILES ${file}
      DESTINATION ${TUTORIALS_INSTALL_DIR}/${TUTORIALS_MOD_DIR}
      PERMISSIONS ${file_permissions})
  endforeach(file)
endfunction()



# FIND JUPYTEXT
execute_process(
  COMMAND which jupytext
  RESULT_VARIABLE JUPYTEXT_NOT_FOUND
  OUTPUT_VARIABLE JUPYTEXT_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check the jupytext_found
if(JUPYTEXT_NOT_FOUND EQUAL 0)
  message(STATUS "Jupytext path: ${JUPYTEXT_PATH}")
else()
  message(WARNING "Jupytext not found")
endif()

  add_custom_target(notebook
  COMMENT "Making jupyter notebooks" VERBATIM
)


# Function to copy and convert markdown files to Jupyter notebooks
function(copy_and_convert_md_to_ipynb source_dir build_dir)
  get_filename_component(TUTORIAL_DIR_NAME ${build_dir} NAME)

  # Define a custom target for the sub-directory
  add_custom_target(notebook_${TUTORIAL_DIR_NAME}
    COMMENT "Building Jupyter notebooks for ${TUTORIAL_DIR_NAME}"
    VERBATIM
  )

  # Initialize the list of IPYNB outputs
  # set(IPYNB_OUTPUTS "")

  file(GLOB MD_FILES "${source_dir}/*.md")

  foreach(MD_FILE ${MD_FILES})
  # WE mean without extension.
    get_filename_component(FILENAME_WE ${MD_FILE} NAME_WE)
    if(FILENAME_WE STREQUAL "README")
      # Skip README.md files
      continue()
    endif()

    set(SOURCE_MD_PATH ${MD_FILE})
    set(DESTINATION_MD_PATH ${build_dir}/${FILENAME_WE}.md)
    set(DESTINATION_IPYNB_PATH ${build_dir}/${FILENAME_WE}.ipynb)

    # Add a custom command to copy the markdown file
    add_custom_command(
      OUTPUT ${DESTINATION_MD_PATH}
      COMMAND ${CMAKE_COMMAND} -E copy ${SOURCE_MD_PATH} ${DESTINATION_MD_PATH}
      DEPENDS ${SOURCE_MD_PATH}
      COMMENT "Copying ${SOURCE_MD_PATH} to build directory"
    )

    # Add a custom command to convert the markdown file to a Jupyter notebook
    if(JUPYTEXT_NOT_FOUND EQUAL 0)
    add_custom_command(
      OUTPUT ${DESTINATION_IPYNB_PATH}
      COMMAND ${JUPYTEXT_PATH} --to notebook ${DESTINATION_MD_PATH}
      DEPENDS ${DESTINATION_MD_PATH}
      COMMENT "Converting ${DESTINATION_MD_PATH} to ${DESTINATION_IPYNB_PATH}"
    )
    else()
    message(WARNING "Jupytext not found. Unable to convert ${DESTINATION_MD_PATH} to a Jupyter notebook.")
  endif()

    # Add a custom target
    add_custom_target(notebook_${FILENAME_WE} ALL
      DEPENDS ${DESTINATION_IPYNB_PATH}
    )
    add_dependencies(notebook notebook_${FILENAME_WE})
    # Make the sub-directory custom target depend on all notebook targets
    add_dependencies(notebook_${TUTORIAL_DIR_NAME} notebook_${FILENAME_WE})
  endforeach()

endfunction(copy_and_convert_md_to_ipynb)

set(TUTORIALS_DIRS
  fun-0
  fun-1
  fun-2
  msh-1
  msh-2
  scl-0
  scl-1
  scl-2
  scl-3
  scl-4
  scl-5
  scl-6
  scl-7
  scl-8
  scl-9
  scl-10
  scl-11
  scl-12
  clx-0
  vec-0
  vec-1
  vec-2
  vec-3
  vec-4
  vec-5
  vec-6
  mix-0
  mix-1
  mix-2
  max-0
  max-1
  adv-0
  adv-1
  adv-2
  adv-3
  adv-4
  adv-5
  cor-0to1
  cor-2to5
  cor-6
  cor-7
  cor-8
  cor-9
  cor-10
)

option(BUILD_TUTORIALS "Compile tutorials" ON)
message(STATUS "Build tutorials ${BUILD_TUTORIALS}")

foreach(LOOP_DIR ${TUTORIALS_DIRS})

  string(REGEX REPLACE "-" "" OP_NAME ${LOOP_DIR})
  string(TOUPPER ${OP_NAME} OP_NAME)
  set(VAR_OP_NAME "BUILD_TUT_${OP_NAME}")

  option("${VAR_OP_NAME}" "Build tutorial {OP_NAME}" ${BUILD_TUTORIALS})
  set(VAR_OP_NAME "${${VAR_OP_NAME}}")

  message(STATUS "Build status ${OP_NAME} ${VAR_OP_NAME}")

  if(VAR_OP_NAME)

    set(TUTORIALS_MOD_DIR ${LOOP_DIR})
    add_subdirectory(${LOOP_DIR})
    tutorials_install(${LOOP_DIR})
    # TODO: we need a little help with this, maybe you can contribute?
    # get_property(test_names DIRECTORY ${LOOP_DIR} PROPERTY TESTS)
    # set_property(TEST ${test_names} PROPERTY LABELS "um-tutorials-ctests")

  endif(VAR_OP_NAME)

endforeach(LOOP_DIR)