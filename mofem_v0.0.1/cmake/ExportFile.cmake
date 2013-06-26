#
# export file: copy it to the build tree on every build invocation and add rule for installation
#
function    (cm_export_file FILE DEST)
  if    (NOT TARGET export-files)
    add_custom_target(export-files ALL COMMENT "Exporting files into build tree")
  endif (NOT TARGET export-files)
  get_filename_component(FILENAME "${FILE}" NAME)
  add_custom_command(TARGET export-files COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${DEST}/${FILENAME}")
  install(FILES "${FILE}" DESTINATION "${DEST}")
endfunction (cm_export_file)

#example of usage
#cm_export_file("API/someHeader0.hpp" "include/API/")
#cm_export_file("API/someHeader1.hpp" "include/API/")
