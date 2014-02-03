# - Try to find GRUMMP
# Once done this will define
#
#  GRUMMP_DIR - directory in which GRUMMP resides

# If unset, try environment
if(NOT GRUMMP_DIR)
  set(GRUMMP_DIR $ENV{GRUMMP_DIR})
endif(NOT GRUMMP_DIR)

if(NOT GRUMMP_ARCH)
  set(GRUMMP_ARCH $ENV{GRUMMP_ARCH})
endif(NOT GRUMMP_ARCH)

#message(STATUS ${GRUMMP_DIR})
#message(STATUS ${GRUMMP_ARCH})

find_file(GRUMMP_LIB_PATH lib/${GRUMMP_ARCH}
  HINTS ${GRUMMP_DIR})
if(NOT GRUMMP_LIB_PATH)
  message(FATAL_ERROR ${GRUMMP_LIB_PATH})
endif(NOT GRUMMP_LIB_PATH)

find_file (GRUMP_DEFAULT_RULES_FILE src/conf/default-rules
  HINTS ${GRUMMP_DIR})
if(NOT GRUMP_DEFAULT_RULES_FILE)
  message(FATAL_ERROR ${GRUMP_DEFAULT_RULES_FILE})
endif(NOT GRUMP_DEFAULT_RULES_FILE)

file(STRINGS ${GRUMP_DEFAULT_RULES_FILE} GRUMP_DEFAULT_RULES) 
foreach(LINE ${GRUMP_DEFAULT_RULES})
  #message(STATUS ${LINE})
  if(LINE MATCHES "^EXTERNLIBS.*") 
    string(REGEX REPLACE "=" ";" FIELDS ${LINE})
    list(GET FIELDS 0 VAR)
    string(REGEX REPLACE " " "" VARSTRIP ${VAR})
    list(GET FIELDS 1 VAL)
    set("GRUMMP_${VARSTRIP}" ${VAL} CACHE INTERNAL "grummp varible")
  endif(LINE MATCHES "^EXTERNLIBS.*")   
  if(LINE MATCHES "^INCLUDES.*") 
    string(REGEX REPLACE "-I../.." "-I${GRUMMP_DIR}" LINE ${LINE})
    string(REGEX REPLACE "=" ";" FIELDS ${LINE})
    list(GET FIELDS 0 VAR)
    string(REGEX REPLACE " " "" VARSTRIP ${VAR})
    list(GET FIELDS 1 VAL)
    set("GRUMMP_${VARSTRIP}" ${VAL} CACHE INTERNAL "grummp varible")
  endif(LINE MATCHES "^INCLUDES.*")   
endforeach(LINE ${GRUMP_DEFAULT_RULES})


