# - Try to find MOAB
# Once done this will define
#
#  MOAB_DIR - directory in which MOAB resides

# If unset, try environment
if(NOT MOAB_DIR)
  set(MOAB_DIR $ENV{MOAB_DIR})
endif(NOT MOAB_DIR)

find_file(MOAB_VARIABLES_FILE moab.make HINTS ${MOAB_DIR}/lib)

if(NOT MOAB_VARIABLES_FILE)
  message(FATAL_ERROR ${MOAB_VARIABLES_FILE})
endif(NOT MOAB_VARIABLES_FILE)

find_file (MBCONVERT NAMES mbconvert mbconvert.exe HINTS ${MOAB_DIR}/bin)
if(NOT MBCONVERT)
  message(FATAL_ERROR ${MBCONVERT})
endif(NOT MBCONVERT)

set(MOAB_INCLUDES_COUNTER 0)

file(STRINGS ${MOAB_VARIABLES_FILE} MOAB_VARIABLES)
foreach(LINE ${MOAB_VARIABLES})
  if(NOT LINE MATCHES "^#.*")
    string(REGEX REPLACE "=" ";" FIELDS ${LINE})
    list(GET FIELDS 0 VAR)
    string(REGEX REPLACE " " "" VARSTRIP ${VAR})

    if(${VARSTRIP} STREQUAL MOAB_INCLUDES)
      set(VARSTRIPEXT "${VARSTRIP}${MOAB_INCLUDES_COUNTER}")
      string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
      set("${VARSTRIPEXT}" ${VAL} CACHE INTERNAL "moab varible")
      set(MOAB_INCLUDES_COUNTER 1)
      # message(STATUS ${VARSTRIPEXT})
      # message(STATUS ${VAL})
    else (${VARSTRIP} STREQUAL MOAB_INCLUDES)
      string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
      # message(STATUS ${VARSTRIP})
      # message(STATUS ${VAL})
      set(${VARSTRIP} ${VAL} CACHE INTERNAL "moab varible")
    endif(${VARSTRIP} STREQUAL MOAB_INCLUDES)

  endif(NOT LINE MATCHES "^#.*")
endforeach(LINE ${MOAB_VARIABLES})

# Add moab definitions
if(MOAB_DEFINITIONS)
  resolve_definitions(MOAB_DEFINITIONS ${MOAB_CPPFLAGS})
  message(STATUS ${MOAB_DEFINITIONS})
  add_definitions(${MOAB_DEFINITIONS})
endif(MOAB_DEFINITIONS)

if(MOAB_HDF5_ENABLED)
  add_definitions("-DMOAB_HDF5_ENABLED")
  set(MOAB_DEFINITIONS "-DMOAB_HDF5_ENABLED")
endif(MOAB_HDF5_ENABLED)
