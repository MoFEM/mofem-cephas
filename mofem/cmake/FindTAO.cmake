# - Try to find TAO
# Once done this will define
#
#  TAO_DIR - directory in which MOAB resides

# If unset, try environment
# If unset, try environment
if(TAO_DIR)
  set(TAO_DIR $ENV{TAO_DIR})
endif(TAO_DIR)

if(TAO_DIR)
  find_library(TAO_LIBRARY NAMES tao PATHS "${TAO_DIR}/${PETSC_ARCH}/lib")
  message(STATUS ${TAO_LIBRARY})
endif(TAO_DIR)
