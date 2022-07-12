# MIT License
#
# Copyright (c) 2022 
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
