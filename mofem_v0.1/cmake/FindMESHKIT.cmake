# - Try to find MESHKIT

# If unset, try environment

if(NOT MESHKIT_DIR)
  set(MESHKIT_DIR $ENV{MESHKIT_DIR})
endif(NOT MESHKIT_DIR)

find_library(MESHKIT_ALGS_LIBRARY NAMES MKalgs PATHS "${MESHKIT_DIR}/lib")
find_library(MESHKIT_UTILS_LIBRARY NAMES MKutils PATHS "${MESHKIT_DIR}/lib")
message(STATUS ${MESHKIT_ALGS_LIBRARY})
message(STATUS ${MESHKIT_UTILS_LIBRARY})

if(MESHKIT_ALGS_LIBRARY) 
  include_directories("${MESHKIT_DIR}/inlude")
endif(MESHKIT_ALGS_LIBRARY)




