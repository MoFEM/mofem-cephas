# User modules are installed from separate repositories. The codel below is obsolete.

# if(STAND_ALLONE_USERS_MODULES) 
#   install(
#     DIRECTORY 
#     ${PROJECT_SOURCE_DIR}/users_modules 
#     DESTINATION ${CMAKE_INSTALL_PREFIX})
# else(STAND_ALLONE_USERS_MODULES) 
#   if(CMAKE_HOST_UNIX)
#     install(CODE
#       "
#       EXECUTE_PROCESS(
#         COMMAND ln -sf ${PROJECT_SOURCE_DIR}/users_modules ${CMAKE_INSTALL_PREFIX}/
#       )
#       MESSAGE(\"-- Linking ${PROJECT_SOURCE_DIR}/users_modules\")
#       "
#     )
#   else(CMAKE_HOST_UNIX)
#     install(
#       DIRECTORY ${PROJECT_SOURCE_DIR}/users_modules 
#       DESTINATION ${CMAKE_INSTALL_PREFIX}
#      )
#   endif(CMAKE_HOST_UNIX)
# endif(STAND_ALLONE_USERS_MODULES)
