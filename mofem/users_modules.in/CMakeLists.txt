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

if(STAND_ALLONE_USERS_MODULES) 
  install(
    DIRECTORY 
    ${PROJECT_SOURCE_DIR}/users_modules 
    DESTINATION ${CMAKE_INSTALL_PREFIX})
else(STAND_ALLONE_USERS_MODULES) 
  if(CMAKE_HOST_UNIX)
    install(CODE
      "
      EXECUTE_PROCESS(
        COMMAND ln -sf ${PROJECT_SOURCE_DIR}/users_modules ${CMAKE_INSTALL_PREFIX}/
      )
      "
    )
  else(CMAKE_HOST_UNIX)
    install(
      DIRECTORY ${PROJECT_SOURCE_DIR}/users_modules 
      DESTINATION ${CMAKE_INSTALL_PREFIX}
     )
  endif(CMAKE_HOST_UNIX)
endif(STAND_ALLONE_USERS_MODULES)
