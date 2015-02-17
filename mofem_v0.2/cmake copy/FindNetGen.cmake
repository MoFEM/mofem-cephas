# - Try to find NETGEN

# If unset, try environment

if(NOT NETGEN_DIR)
  set(NETGEN_DIR $ENV{NETGEN_DIR})
endif(NOT NETGEN_DIR)

if(NOT NETGEN_SRC_DIR)
  set(NETGEN_SRC_DIR $ENV{NETGEN_SRC_DIR})
endif(NOT NETGEN_SRC_DIR)


find_library(NETGEN_LIBRARY NAMES nglib PATHS "${NETGEN_DIR}/lib")
message(STATUS ${NETGEN_LIBRARY})

if(NETGEN_LIBRARY) 
  include_directories("${NETGEN_DIR}/include")
  include_directories("${NETGEN_SRC_DIR}/libsrc/include")
  add_definitions( -DWITH_NETGEN )
endif(NETGEN_LIBRARY)


