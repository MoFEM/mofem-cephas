# - Try to find TAO
# Once done this will define
#
#  TAO_DIR - directory in which MOAB resides

# If unset, try environment
if(TAO_DIR)
  set(TAO_DIR $ENV{TAO_DIR})
  message(${TAO_DIR}/${PETSC_ARCH}/lib)
  find_library(TAO_LIBRARY NAMES tao PATHS "${TAO_DIR}/${PETSC_ARCH}/lib")
endif(TAO_DIR)




