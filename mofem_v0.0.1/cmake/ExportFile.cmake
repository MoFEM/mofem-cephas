#
# export file: copy it to the build tree on every build invocation and add rule for installation
#
function(cm_export_file FILE target_name)
  if(NOT TARGET ${target_name})
    add_custom_target(${target_name} ALL COMMENT "Exporting files into build tree")
  endif (NOT TARGET ${target_name})
  get_filename_component(FILENAME "${FILE}" NAME)
  add_custom_command(TARGET ${target_name} COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
endfunction (cm_export_file)

#example of usage
#cm_export_file("API/someHeader0.hpp")
#cm_export_file("API/someHeader1.hpp")
