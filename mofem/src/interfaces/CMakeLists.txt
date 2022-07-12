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

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

#multi_indices
set(MOFEM_INTERFACES_SRC
  impl/MatPartitioningApply_Parmetis_MoFEM.cpp
  impl/Core.cpp
  impl/DeprecatedCoreInterface.cpp
  impl/FieldCore.cpp
  impl/FECore.cpp
  impl/ProblemsCore.cpp
  impl/DeleteCore.cpp
  impl/LoopMethods.cpp
  impl/LogManager.cpp
  impl/AllCoreInterfaces.cpp
  impl/AllInterfaces.cpp
)

if(PRECOMPILED_HEADRES)
  set_source_files_properties(
    All.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
  set_source_files_properties(
    AllInterfaces.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
  set_source_files_properties(
    LogManager.cpp
    PROPERTIES
    COMPILE_FLAGS "-include ${PROJECT_BINARY_DIR}/include/precompiled/Includes.hpp"
  )
endif(PRECOMPILED_HEADRES)
add_library(
  mofem_interfaces
  ${MOFEM_INTERFACES_SRC}
)
set_target_properties(mofem_interfaces PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(mofem_interfaces install_prerequisites)
if(PRECOMPILED_HEADRES)
  add_dependencies(mofem_interfaces Includes.hpp.pch_copy)
endif(PRECOMPILED_HEADRES)

install(TARGETS mofem_interfaces DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
