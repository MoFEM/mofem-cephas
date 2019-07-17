# - Try to find CGM

# If unset, try environment

if(NOT CGM_DIR)
  set(CGM_DIR $ENV{CGM_DIR})
endif(NOT CGM_DIR)

find_file(CGM_VARIABLES_FILE cgm.make HINTS ${CGM_DIR}/lib)

if(CGM_VARIABLES_FILE)

  set(CGM_INCLUDES_COUNTER 0)

  file(STRINGS ${CGM_VARIABLES_FILE} CGM_VARIABLES)
  foreach(LINE ${CGM_VARIABLES})
    if(NOT LINE MATCHES "^#.*")
      string(REGEX REPLACE "=" ";" FIELDS ${LINE})
      list(GET FIELDS 0 VAR)
      string(REGEX REPLACE " " "" VARSTRIP ${VAR})

      if(${VARSTRIP} STREQUAL CGM_INCLUDES)
        set(VARSTRIPEXT "${VARSTRIP}${CGM_INCLUDES_COUNTER}")
        string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
        set("${VARSTRIPEXT}" ${VAL} CACHE INTERNAL "moab varible")
        set(CGM_INCLUDES_COUNTER 1)
        # message(STATUS ${VARSTRIPEXT})
        # message(STATUS ${VAL})
      else (${VARSTRIP} STREQUAL CGM_INCLUDES)
        string(REGEX REPLACE "${VARSTRIP} *=" "" VAL ${LINE})
        # message(STATUS ${VARSTRIP})
        # message(STATUS ${VAL})
        set(${VARSTRIP} ${VAL} CACHE INTERNAL "moab varible")
      endif(${VARSTRIP} STREQUAL CGM_INCLUDES)

    endif(NOT LINE MATCHES "^#.*")
  endforeach(LINE ${CGM_VARIABLES})

  # Add moab definitions
  if(CGM_DEFINITIONS)
    resolve_definitions(CGM_DEFINITIONS ${CGM_CPPFLAGS})
    message(STATUS ${CGM_DEFINITIONS})
    add_definitions(${CGM_DEFINITIONS})
  endif(CGM_DEFINITIONS)

endif(CGM_VARIABLES_FILE)
