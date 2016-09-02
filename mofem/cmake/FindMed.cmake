# - Try to find MED

if(NOT MED_DIR)
  set(MED_DIR $ENV{MED_DIR})
endif(NOT MED_DIR)

if(MED_DIR)
  find_library(MED_LIBRARY NAMES med PATHS ${MED_DIR}/lib)
  message(STATUS ${MED_LIBRARY})
  if(MED_LIBRARY)
    include_directories(${MED_DIR}/include)
    add_definitions(-DWITH_MED)
  endif(MED_LIBRARY)
endif(MED_DIR)
