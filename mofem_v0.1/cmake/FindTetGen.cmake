# - Try to find TETGEN

# If unset, try environment

if(NOT TETGEN_DIR)
  set(TETGEN_DIR $ENV{TETGEN_DIR})
endif(NOT TETGEN_DIR)

find_library(TETGEN_LIBRARY NAMES tet PATHS "${TETGEN_DIR}/lib")
message(STATUS ${TETGEN_LIBRARY})

if(TETGEN_LIBRARY) 
  include_directories("${TETGEN_DIR}/include")
endif(TETGEN_LIBRARY)




