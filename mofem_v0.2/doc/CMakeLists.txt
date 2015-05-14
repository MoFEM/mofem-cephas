# MoFEM is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# MoFEM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

if(DOXYGEN_FOUND)

  # copy dox/figures to html direcyory created by doxygen
  add_custom_target(doxygen_copy_figures
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/doc/figures ${PROJECT_BINARY_DIR}/html
    DEPENDS doxygen_generating_documenation
  )
  add_dependencies(doc doxygen_copy_figures)

  #copy figures form users UM documentation
  add_custom_target(doxygen_copy_figures_from_user_modules
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/users_modules/doc/figures ${PROJECT_BINARY_DIR}/html
    DEPENDS doxygen_generating_documenation
  )
  add_dependencies(doc doxygen_copy_figures_from_user_modules)

  find_program(SFTP NAMES sftp)
  if(SFTP)

    # This work that script publish_doc.sh is created first in source
    # directory, then copied to binary direcyory with set premissions to execute.
    # Script publish_doc.sh, run sftp, which itself run its on script
    # from PROJECT_SOURCE_DIR/publish_doc

    configure_file(
      ${PROJECT_SOURCE_DIR}/scripts/publish_doc.sh.in
      ${PROJECT_SOURCE_DIR}/scripts/publish_doc.sh
      @ONLY
    )
    file(COPY ${PROJECT_SOURCE_DIR}/scripts/publish_doc.sh
      DESTINATION ${PROJECT_BINARY_DIR}/scripts/
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
    )
    add_custom_target(doc_publish_likask
      ${PROJECT_BINARY_DIR}/scripts/publish_doc.sh
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      COMMENT "Publish (likask) doxygen documentation" VERBATIM
    )

  endif(SFTP)

endif(DOXYGEN_FOUND)
