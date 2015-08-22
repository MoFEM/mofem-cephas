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

find_package(Doxygen)

if(DOXYGEN_FOUND)

  configure_file(
    ${PROJECT_SOURCE_DIR}/Doxyfile.in
    ${PROJECT_BINARY_DIR}/Doxyfile
  )
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
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/doc/figures ${PROJECT_BINARY_DIR}/html_users_modules
    DEPENDS doxygen_generating_documenation
  )
  add_dependencies(doc doxygen_copy_figures)

  #copy pdfs form users UM documentation
  add_custom_target(doxygen_copy_pdfs_from_user_modules
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/doc/pdfs ${PROJECT_BINARY_DIR}/html_users_modules
    DEPENDS doxygen_generating_documenation
  )
  add_dependencies(doc doxygen_copy_pdfs_from_user_modules)

endif(DOXYGEN_FOUND)
