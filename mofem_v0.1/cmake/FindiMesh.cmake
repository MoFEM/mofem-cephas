# - Try to find MESHKIT

# If unset, try environment

find_library(iMESH_LIBRARY NAMES iMesh PATHS "${MOAB_DIR}/lib")
message(STATUS ${iMESH_LIBRARY})




