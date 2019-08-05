# - Try to find MOAB
# Once done this will define
#
#  PETSC_DIR - directory in which PETSc resides
#  PETSC_ARCH - build architecture
#

# If unset, try environment
if(NOT PETSC_DIR)
  set(PETSC_DIR $ENV{PETSC_DIR})
endif(NOT PETSC_DIR)
if(NOT PETSC_ARCH)
  set(PETSC_ARCH $ENV{PETSC_ARCH})
endif(NOT PETSC_ARCH)

find_file(PETSC_VARIABLES_FILE petscvariables
  HINTS
  ${PETSC_DIR}/${PETSC_ARCH}/conf
  ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf
  ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc-conf
)
if(NOT PETSC_VARIABLES_FILE)
  message(FATAL_ERROR ${PETSC_VARIABLES_FILE})
endif(NOT PETSC_VARIABLES_FILE)

file(STRINGS ${PETSC_VARIABLES_FILE} PETSC_VARIABLES)
foreach(LINE ${PETSC_VARIABLES})
  string(REGEX REPLACE " = " ";" FIELDS ${LINE})
  list(LENGTH FIELDS LISTLEN)
  if(LISTLEN EQUAL 2)
    list(GET FIELDS 0 VAR)
    list(GET FIELDS 1 VAL)
    set("PETSCVAR_${VAR}" ${VAL} CACHE INTERNAL "petsc varible")
    #message(STATUS PETSCVAR_${VAR})
  endif(LISTLEN EQUAL 2)
endforeach(LINE ${PETSC_VARIABLES})
