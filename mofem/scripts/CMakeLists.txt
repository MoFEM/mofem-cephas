# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
# It can be freely used for educational and research purposes
# by other institutions. If you use this softwre pleas cite my work.
#
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


configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_build.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_build.sh"
  @ONLY
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_fast_check.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_fast_check.sh"
  @ONLY
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_update.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_update.sh"
  @ONLY
  )

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_build.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_fast_check.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_update.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
