
# Add BOOST directory
if(NOT BOOST_DIR)
  set(BOOST_DIR $ENV{BOOST_DIR})
endif(NOT BOOST_DIR)

find_package(
  Boost COMPONENTS 
  program_options log log_setup thread system filesystem REQUIRED)

if(NOT Boost_LIBRARIES)
  message(FATAL_ERROR "boost libraries not found")
endif(NOT Boost_LIBRARIES)

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
