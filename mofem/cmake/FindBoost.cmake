
# Add BOOST directory
if(NOT BOOST_DIR)
  set(BOOST_DIR $ENV{BOOST_DIR})
endif(NOT BOOST_DIR)

find_library(
  BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options
  HINTS
  ${BOOST_DIR}/lib
  ${PETSC_DIR}/${PETSC_ARCH}/lib
)
message(STATUS "Boost found: " ${BOOST_PROGRAM_OPTIONS_LIB})
if(NOT BOOST_PROGRAM_OPTIONS_LIB)
  message(FATAL_ERROR "boost program options library not found")
endif(NOT BOOST_PROGRAM_OPTIONS_LIB)

find_package(
  Boost COMPONENTS 
  program_options log log_setup thread system filesystem REQUIRED)

find_path(
  BOOST_INCLUDE_DIR
  NAMES boost/multi_index_container.hpp
  HINTS
  ${BOOST_DIR}/include
  ${PETSC_DIR}/${PETSC_ARCH}/include
)
message(STATUS "Boost include found: " ${BOOST_INCLUDE_DIR})
if(NOT BOOST_INCLUDE_DIR)
  message(FATAL_ERROR "boost program include dir not found")
endif(NOT BOOST_INCLUDE_DIR)
include_directories(${BOOST_INCLUDE_DIR})
