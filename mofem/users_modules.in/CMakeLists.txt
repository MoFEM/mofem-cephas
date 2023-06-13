# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
# It can be freely used for educational and research purposes
# by other institutions. If you use this softwre pleas cite my work.
#

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
      MESSAGE(\"-- Linking ${PROJECT_SOURCE_DIR}/users_modules\")
      "
    )
  else(CMAKE_HOST_UNIX)
    install(
      DIRECTORY ${PROJECT_SOURCE_DIR}/users_modules 
      DESTINATION ${CMAKE_INSTALL_PREFIX}
     )
  endif(CMAKE_HOST_UNIX)
endif(STAND_ALLONE_USERS_MODULES)
