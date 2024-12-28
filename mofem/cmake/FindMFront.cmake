# - Try to find MFront

if(NOT MGIS_DIR)
  message( WARNING "MGIS_DIR is not defined! MFront Interface will not be configured." )
endif(NOT MGIS_DIR)

if(MGIS_DIR)
  find_library(MGIS_LIBRARY NAMES libMFrontGenericInterface.dylib libMFrontGenericInterface.so PATHS ${MGIS_DIR}/lib)
  message(STATUS "MGIS_LIBRARY ${MGIS_LIBRARY}")
  # find_path(MGIS_HEADER PATHS ${MGIS_DIR}/include/MGIS/)
  add_library(MFrontGenericInterface SHARED IMPORTED)
  set_target_properties(MFrontGenericInterface PROPERTIES IMPORTED_LOCATION ${MGIS_LIBRARY})
  if(MGIS_LIBRARY)
    # include_directories(${MGIS_HEADER})
    include_directories(${MGIS_DIR}/include)
    include_directories(${MGIS_DIR}/lib)
    add_definitions(-DWITH_MGIS)
  endif(MGIS_LIBRARY)
endif(MGIS_DIR)

# if(NOT TFEL_DIR)
#   message( WARNING "TFEL_DIR is not defined! Unable to test behavior compilation." )
# endif(NOT TFEL_DIR)

# if(TFEL_DIR)
#   message(STATUS "TFEL_DIR FOUND")
# endif(TFEL_DIR)
