# - Try to find triangle

# If unset, try environment

if(NOT TRIANGLE_DIR)
  set(TRIANGLE_DIR $ENV{TRIANGLE_DIR})
endif(NOT TRIANGLE_DIR)

find_library(TRIANGLE_LIBRARY NAMES triangle PATHS "${TRIANGLE_DIR}/lib")
message(STATUS ${TRIANGLE_LIBRARY})

if(TRIANGLE_LIBRARY) 
  include_directories("${TRIANGLE_DIR}/include")
endif(TRIANGLE_LIBRARY)




