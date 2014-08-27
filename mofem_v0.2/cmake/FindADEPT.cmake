# - Try to find ADEPT
# Once done this will define
#
#  ADEPT_DIR - directory in which MOAB resides

# If unset, try environment
if(ADEPT_DIR)
  set(ADEPT_DIR $ENV{ADEPT_DIR})
endif(ADEPT_DIR)

find_library(ADEPT_LIBRARY NAMES adept PATHS "${ADEPT_DIR}/lib")
message(STATUS ${ADEPT_LIBRARY})

if(ADEPT_LIBRARY) 
  include_directories("${ADEPT_DIR}/include")
endif(ADEPT_LIBRARY)






