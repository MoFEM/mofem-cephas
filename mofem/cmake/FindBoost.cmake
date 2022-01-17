# Add BOOST directory
set(BOOST_ROOT "${BOOST_DIR}")

find_package(
  Boost 
  REQUIRED COMPONENTS
  program_options log log_setup thread system filesystem)

message(STATUS "Boost_LIBRARIES ${Boost_LIBRARIES}")
message(STATUS "Boost_LIBRARY_DIRS ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_INCLUDE_DIRS ${Boost_INCLUDE_DIRS}")

if(NOT Boost_LIBRARIES)
  message(FATAL_ERROR "boost libraries not found")
endif(NOT Boost_LIBRARIES)

find_path(
  BOOST_INCLUDE_DIR
  NAMES boost/multi_index_container.hpp
  HINTS
  ${Boost_INCLUDE_DIRS}
  ${BOOST_DIR}/include
  ${PETSC_DIR}/${PETSC_ARCH}/include
)
message(STATUS "Boost include found: " ${BOOST_INCLUDE_DIR})
if(NOT BOOST_INCLUDE_DIR)
  message(FATAL_ERROR "Boost include dir not found")
endif(NOT BOOST_INCLUDE_DIR)
include_directories(${BOOST_INCLUDE_DIR})

add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
