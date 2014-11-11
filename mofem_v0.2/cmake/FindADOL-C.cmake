# - Try to find ADOL-C
# Once done this will define
#
#  ADOL-C_DIR - directory in which MOAB resides

# If unset, try environment
if(ADOL-C_DIR)
  set(ADOL-C_DIR $ENV{ADOL-C_DIR})
endif(ADOL-C_DIR)

find_library(ADOL-C_LIBRARY NAMES adolc PATHS "${ADOL-C_DIR}/lib")
message(STATUS ${ADOL-C_LIBRARY})

if(ADOL-C_LIBRARY) 
  include_directories("${ADOL-C_DIR}/include")
  add_definitions( -DWITH_ADOL_C )
endif(ADOL-C_LIBRARY)






