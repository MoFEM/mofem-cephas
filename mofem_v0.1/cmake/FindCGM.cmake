# - Try to find CGM

# If unset, try environment

if(NOT CGM_DIR)
  set(CGM_DIR $ENV{CGM_DIR})
endif(NOT CGM_DIR)


find_library(CGM_LIBRARY NAMES cgm PATHS "${CGM_DIR}/lib")
find_library(CGM_iGEOM_LIBRARY NAMES iGeom  PATHS "${CGM_DIR}/lib")
message(STATUS ${CGM_LIBRARY})
message(STATUS ${CGM_iGEOM_LIBRARY})

if(CGM_LIBRARY) 
  include_directories("${CGM_DIR}/include")
endif(CGM_LIBRARY) 




