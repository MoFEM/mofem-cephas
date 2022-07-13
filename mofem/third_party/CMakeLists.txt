

#include subdirectory
include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)

add_library(mofem_third_party
  IntegrationRules.cpp
  gm_rule.c
  tetcircumcenter.c
  tricircumcenter.c
  phg-quadrule/quad.c)
target_link_libraries(mofem_third_party PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_third_party PROPERTIES VERSION ${PROJECT_VERSION})
install(TARGETS mofem_third_party DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# FIXME: Install directory, and other stuff, like kiss lib,
# that should not be there.  It is for back compatibility.
install(
  DIRECTORY
  ${CMAKE_SOURCE_DIR}/third_party
  DESTINATION
  ${CMAKE_INSTALL_PREFIX})
