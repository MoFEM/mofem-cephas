# - Try to find MED

if(NOT MED_DIR)
  set(MED_DIR $ENV{MED_DIR})
endif(NOT MED_DIR)

if(MED_DIR)
  find_library(MED_LIBRARY NAMES med medC PATHS ${MED_DIR}/lib)
  find_path(MED_HEADER NAMES med.h PATHS ${MED_DIR}/include)
  message(STATUS "MED LIB ${MED_LIBRARY}")
  message(STATUS "MED HEADER ${MED_HEADER}")
  if(MED_LIBRARY AND MED_HEADER)
    include_directories(${MED_HEADER})
    add_definitions(-DWITH_MED)
  endif(MED_LIBRARY AND MED_HEADER)
endif(MED_DIR)

