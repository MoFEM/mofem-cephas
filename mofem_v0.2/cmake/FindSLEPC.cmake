# - Try to find SLEPC
# Once done this will define
#
#  SLEPC_DIR - directory in which MOAB resides

# If unset, try environment
if(SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif(SLEPC_DIR)

find_library(SLEPC_LIBRARY NAMES slepc PATHS "${SLEPC_DIR}/lib")
message(STATUS ${ADOL-C_LIBRARY})

if(SLEPC_LIBRARY) 
  include_directories("${SLEPC_DIR}/include")
  add_definitions( -DWITH_SLEPC )
endif(SLEPC_LIBRARY)






