

file(GLOB_RECURSE FTENSOR_ONE_OVER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/one_over_*.cpp
)
foreach(ONE_OVER_SRC_FILE ${FTENSOR_ONE_OVER_FILES})
  string(REGEX REPLACE ".cpp" "" ONE_OVER ${ONE_OVER_SRC_FILE})
  string(REGEX REPLACE ".*/" "" ONE_OVER ${ONE_OVER})
  # message(STATUS ${ONE_OVER})
  add_executable(${ONE_OVER} ${ONE_OVER_SRC_FILE})
  add_dependencies(${ONE_OVER} install_prerequisites)
endforeach(ONE_OVER_SRC_FILE ${FTENSOR_ONE_OVER_FILES})

# Copy script doing tests
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/one_over_script
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/
)
