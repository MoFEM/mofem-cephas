find_package(Git)
find_package(Doxygen)

if(DOXYGEN_FOUND)

  # running doxygen
  configure_file(
    ${PROJECT_SOURCE_DIR}/Doxyfile.in
    ${PROJECT_BINARY_DIR}/Doxyfile
  )
  
  # set(DOXYGEN_GENERATE_HEADER
  #   ${DOXYGEN_EXECUTABLE} -w html
  #   ${PROJECT_SOURCE_DIR}/doc/header.html
  #   ${PROJECT_SOURCE_DIR}/doc/footer.html
  #   ${PROJECT_SOURCE_DIR}/doc/style.html
  #   ${PROJECT_BINARY_DIR}/Doxyfile
  # )
  # add_custom_target(doxygen_generating_header
  #   ${DOXYGEN_GENERATE_HEADER}
  #   COMMAND cat ${PROJECT_SOURCE_DIR}/doc/google_analytics.html | tee -a ${PROJECT_SOURCE_DIR}/doc/header.html
  #   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  #   COMMENT "Generating doxygen header" VERBATIM
  # )
  #
  # add_custom_target(doxygen_generating_footer
  #   ${DOXYGEN_GENERATE_HEADER}
  #   COMMAND cat ${PROJECT_SOURCE_DIR}/doc/chatillo.html | tee -a ${PROJECT_SOURCE_DIR}/doc/footer.html
  #   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  #   COMMENT "Generating doxygen header" VERBATIM
  # )

  set(DOXYGEN_LINE_COMMAND
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
  )
  add_custom_target(doxygen_generating_documenation
    ${DOXYGEN_LINE_COMMAND}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating documentation with Doxygen" VERBATIM
  )

  add_custom_target(doc
    DEPENDS doxygen_generating_documenation
    COMMENT "Making documentation" VERBATIM
  )

  # copy dox/figures to html directory created by doxygen
  add_custom_target(doxygen_copy_figures
    ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/doc/figures ${PROJECT_BINARY_DIR}/html
  )
  add_dependencies(doc doxygen_copy_figures)

  # copy slider to html directory created by doxygen
  add_custom_target(doxygen_copy_slider
    ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/doc/slider ${PROJECT_BINARY_DIR}/html
  )
  add_dependencies(doc doxygen_copy_slider)

  # copy google site verification file
  add_custom_target(doxygen_google_site_verification
    ${CMAKE_COMMAND} -E copy
    ${PROJECT_SOURCE_DIR}/doc/google94051b0650bf905b.html ${PROJECT_BINARY_DIR}/html
  )
  add_dependencies(doc doxygen_google_site_verification)

  # copy script files 
  add_custom_target(doxygen_copy_script_files
    ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/scripts ${PROJECT_BINARY_DIR}/html/scripts
  )
  add_dependencies(doc doxygen_copy_script_files)

  # copy mastadon files
  add_custom_target(doxygen_mastadoc_files
    ${CMAKE_COMMAND} -E copy
    ${PROJECT_SOURCE_DIR}/doc/mastodon* ${PROJECT_BINARY_DIR}/html
  )
  add_dependencies(doc doxygen_mastadoc_files)

  # users modules docs
  # FIXME: This only works if modules are in ${PROJECT_SOURCE_DIR}/users_modules
  file(
    GLOB_RECURSE INSTLLED_DOC_MODULES
    FOLLOW_SYMLINKS
    ${PROJECT_SOURCE_DIR}/users_modules/*AddDocumentation.cmake
  )
  if(EXTERNAL_MODULE_SOURCE_DIRS)
    foreach(LOOP_DIR ${EXTERNAL_MODULE_SOURCE_DIRS})
      message(STATUS "Search module directory: " ${LOOP_DIR})
      file(
        GLOB_RECURSE EXTERNAL_INSTLLED_DOC_MODULES
        FOLLOW_SYMLINKS
        ${LOOP_DIR}/*AddDocumentation.cmake
      )    
      message(STATUS "Found: " ${EXTERNAL_INSTLLED_DOC_MODULES})
      set(INSTLLED_DOC_MODULES ${INSTLLED_DOC_MODULES} ${EXTERNAL_INSTLLED_DOC_MODULES})
    endforeach(LOOP_DIR)
  endif(EXTERNAL_MODULE_SOURCE_DIRS)

  foreach(LOOP_DOC_MODULE ${INSTLLED_DOC_MODULES})
    message(STATUS "Add module documentation ... ${LOOP_DOC_MODULE}")
    string(REGEX REPLACE
      "/+AddDocumentation.cmake" "" ADD_DOC_DIRECTORY ${LOOP_DOC_MODULE})
    include(${LOOP_DOC_MODULE})
  endforeach(LOOP_DOC_MODULE)

  find_program(RSYNC NAMES rsync)
  if(RSYNC)

    # This work that script publish_doc_likask.sh is created first in source
    # directory, then copied to binary direcyory with set premissions to execute.
    # Script publish_doc_likask.sh, run sftp, which itself run its on script
    # from PROJECT_SOURCE_DIR/publish_doc_likask

    configure_file(
      ${PROJECT_SOURCE_DIR}/scripts/publish_doc_cdash.sh.in
      ${PROJECT_BINARY_DIR}/scripts/tmp/publish_doc_cdash.sh
      @ONLY
    )
    file(COPY ${PROJECT_BINARY_DIR}/scripts/tmp/publish_doc_cdash.sh
      DESTINATION ${PROJECT_BINARY_DIR}/scripts/
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
    )
    add_custom_target(doc_publish_cdash
      ${PROJECT_BINARY_DIR}/scripts/publish_doc_cdash.sh
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      COMMENT "Publish (cdash.eng.gla.ac) doxygen documentation" VERBATIM
    )

  endif(RSYNC)

endif(DOXYGEN_FOUND)
