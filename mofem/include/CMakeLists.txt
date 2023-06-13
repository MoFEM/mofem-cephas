include_directories(${PROJECT_SOURCE_DIR}/third_party)
include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/third_party/phg-quadrule)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation)
include_directories(${PROJECT_SOURCE_DIR}/src/approximation/c)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/post_proc)
include_directories(${PROJECT_SOURCE_DIR}/src/boundary_conditions)

# Configure a header file to pass some of the CMake settings to the source code
configure_file(
  ${PROJECT_SOURCE_DIR}/include/config.h.in
  ${PROJECT_BINARY_DIR}/include/config.h
)
set_source_files_properties(
  ${PROJECT_BINARY_DIR}/include/config.h PROPERTIES GENERATED TRUE
)
include_directories(${PROJECT_BINARY_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

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
  )
  add_library(Includes.hpp.pch OBJECT Includes.hpp)
  add_dependencies(Includes.hpp.pch install_prerequisites)
  add_custom_target(
    Includes.hpp.pch_copy
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Includes.hpp.pch.dir/Includes.hpp.o
    ${CMAKE_CURRENT_BINARY_DIR}/precompiled/Includes.hpp.${OUT_PCH_SUFFIX}
    COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/Includes.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/precompiled/Includes.hpp
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
  )
  add_library(MoFEM.hpp.pch OBJECT MoFEM.hpp)
  add_dependencies(MoFEM.hpp.pch install_prerequisites)
  add_custom_target(
    MoFEM.hpp.pch_copy
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/MoFEM.hpp.pch.dir/MoFEM.hpp.o
    ${CMAKE_CURRENT_BINARY_DIR}/precompiled/MoFEM.hpp.${OUT_PCH_SUFFIX}
    COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_SOURCE_DIR}/MoFEM.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/precompiled/MoFEM.hpp
    COMMENT
    "Copy precompiled MoFEM.hpp header"
  )
  add_dependencies(MoFEM.hpp.pch_copy MoFEM.hpp.pch)

endif(PRECOMPILED_HEADRES)

if(STAND_ALLONE_USERS_MODULES)
  install(
    DIRECTORY 
    ${CMAKE_SOURCE_DIR}/include/
    DESTINATION 
    ${CMAKE_INSTALL_PREFIX}/include
    FILES_MATCHING 
    PATTERN "*.hpp"
    PATTERN "*.h")
else(STAND_ALLONE_USERS_MODULES)
  install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
  file(
    GLOB HEADER_FILES
    RELATIVE ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/*.h*)
  foreach(HF ${HEADER_FILES})
    install(CODE
      "
      EXECUTE_PROCESS(
      COMMAND 
      ln -sf 
      ${CMAKE_SOURCE_DIR}/include/${HF} 
      ${CMAKE_INSTALL_PREFIX}/include/${HF}
      )
      MESSAGE(\"-- Linking ${CMAKE_INSTALL_PREFIX}/include/${HF}\")
      "
    )
  endforeach(HF)
endif(STAND_ALLONE_USERS_MODULES)
