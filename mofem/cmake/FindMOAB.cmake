# - Try to find MOAB
# Once done this will define
#
#  MOAB_DIR - directory in which MOAB resides

# If unset, try environment
if(NOT MOAB_DIR)
  set(MOAB_DIR $ENV{MOAB_DIR})
endif(NOT MOAB_DIR)

find_file (MOAB_VARIBLES_FILE moab.make
  HINTS ${MOAB_DIR}/lib)
if(NOT MOAB_VARIBLES_FILE)
  message(FATAL_ERROR ${MOAB_VARIBLES_FILE})
endif(NOT MOAB_VARIBLES_FILE)

find_file (MBCONVERT mbconvert HINTS ${MOAB_DIR}/bin)
if(NOT MBCONVERT)
  message(FATAL_ERROR ${MBCONVERT})
endif(NOT MBCONVERT)

set(MOAB_INCLUDES_COUNTER 0)

file(STRINGS ${MOAB_VARIBLES_FILE} MOAB_VARIBLES)
foreach(LINE ${MOAB_VARIBLES})
  if(NOT LINE MATCHES "^#.*")
    string(REGEX REPLACE "=" ";" FIELDS ${LINE})
    list(GET FIELDS 0 VAR)
    string(REGEX REPLACE " " "" VARSTRIP ${VAR})
    list(GET FIELDS 1 VAL)
    set("${VARSTRIP}" ${VAL} CACHE INTERNAL "moab varible")
    #message(STATUS ${VARSTRIP})
    #message(STATUS ${VAL})
    if(${VARSTRIP} STREQUAL "${MOAB_INCLUDES}")
      set(VARSTRIPEXT "${VARSTRIP}${MOAB_INCLUDES_COUNTER}")
      set("${VARSTRIPEXT}" ${VAL} CACHE INTERNAL "moab varible")
      set(MOAB_INCLUDES_COUNTER 1)
      #message(STATUS ${VARSTRIPEXT})
      #message(STATUS ${VAL})
    endif(${VARSTRIP} STREQUAL "${MOAB_INCLUDES}")
  endif(NOT LINE MATCHES "^#.*")
endforeach(LINE ${MOAB_VARIBLES})
