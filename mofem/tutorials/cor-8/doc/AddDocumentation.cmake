# copy dox/figures to html directory created by doxygen
add_custom_target(spring_element
  ${CMAKE_COMMAND} -E copy_directory
  ${ADD_DOC_DIRECTORY}/figures ${PROJECT_BINARY_DIR}/html)
add_dependencies(doc spring_element)