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

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/third_party/cblas)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)

if(PRECOMPILED_HEADRES)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(OUT_PCH_SUFFIX "pch")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(OUT_PCH_SUFFIX "gch")
  endif()

  # Includes.hpp
  set_source_files_properties(
    Includes.hpp
    PROPERTIES
    LANGUAGE CXX
    COMPILE_FLAGS "-x c++-header"
    OBJECT_OUTPUTS "Includes.hpp.pch"
  )
  add_library(Includes.hpp.pch OBJECT Includes.hpp)
  add_custom_target(
    Includes.hpp.pch_copy
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Includes.hpp.pch.dir/Includes.hpp.o
    ${CMAKE_CURRENT_SOURCE_DIR}/Includes.hpp.${OUT_PCH_SUFFIX}
    COMMENT
    "Copy precompiled Includes.hpp header"
  )
  add_dependencies(Includes.hpp.pch_copy Includes.hpp.pch)

  # MoFEM.hpp
  set_source_files_properties(
    MoFEM.hpp
    PROPERTIES
    LANGUAGE CXX
    COMPILE_FLAGS "-x c++-header"
    OBJECT_OUTPUTS "MoFEM.hpp.pch"
  )
  add_library(MoFEM.hpp.pch OBJECT MoFEM.hpp)
  add_dependencies(MoFEM.hpp.pch Includes.hpp.pch_copy)
  add_custom_target(
    MoFEM.hpp.pch_copy
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/MoFEM.hpp.pch.dir/MoFEM.hpp.o
    ${CMAKE_CURRENT_SOURCE_DIR}/MoFEM.hpp.${OUT_PCH_SUFFIX}
    COMMENT
    "Copy precompiled MoFEM.hpp header"
  )
  add_dependencies(MoFEM.hpp.pch_copy MoFEM.hpp.pch)

endif(PRECOMPILED_HEADRES)
