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

# multi_indices
if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(mofem_multi_indices
  impl/All.cpp
)
target_link_libraries(mofem_multi_indices PUBLIC ${MoFEM_PROJECT_LIBS})
set_target_properties(mofem_multi_indices PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_multi_indices install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_multi_indices Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)
install(TARGETS mofem_multi_indices DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
