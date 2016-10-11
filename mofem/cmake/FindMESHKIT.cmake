# - Try to find MESHKIT

# If unset, try environment

if(NOT MESHKIT_DIR)
  set(MESHKIT_DIR $ENV{MESHKIT_DIR})
endif(NOT MESHKIT_DIR)

find_file(MESHKIT_VARIBLES_FILE meshkit.make HINTS ${MESHKIT_DIR}/lib)

if(MESHKIT_VARIBLES_FILE)

  set(MESHKIT_INCLUDES_COUNTER 0)

  file(STRINGS ${MESHKIT_VARIBLES_FILE} MESHKIT_VARIBLES)
  foreach(LINE ${MESHKIT_VARIBLES})
    if(NOT LINE MATCHES "^#.*")
      string(REGEX REPLACE "=" ";" FIELDS ${LINE})
      list(GET FIELDS 0 VAR)
      string(REGEX REPLACE " " "" VARSTRIP ${VAR})

      if(${VARSTRIP} STREQUAL MESHKIT_INCLUDES)
        set(VARSTRIPEXT "${VARSTRIP}${MESHKIT_INCLUDES_COUNTER}")
        string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
        set("${VARSTRIPEXT}" ${VAL} CACHE INTERNAL "moab varible")
        set(MESHKIT_INCLUDES_COUNTER 1)
        # message(STATUS ${VARSTRIPEXT})
        # message(STATUS ${VAL})
      else (${VARSTRIP} STREQUAL MESHKIT_INCLUDES)
        string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
        # message(STATUS ${VARSTRIP})
        # message(STATUS ${VAL})
        set(${VARSTRIP} ${VAL} CACHE INTERNAL "moab varible")
      endif(${VARSTRIP} STREQUAL MESHKIT_INCLUDES)

    endif(NOT LINE MATCHES "^#.*")
  endforeach(LINE ${MESHKIT_VARIBLES})

  # Add moab definitions
  if(MESHKIT_DEFINITIONS)
    resolve_definitions(MESHKIT_DEFINITIONS ${MESHKIT_CPPFLAGS})
    message(STATUS ${MESHKIT_DEFINITIONS})
    add_definitions(${MESHKIT_DEFINITIONS})
  endif(MESHKIT_DEFINITIONS)

endif(MESHKIT_VARIBLES_FILE)
