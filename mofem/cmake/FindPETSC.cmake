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

find_file(PETSC_VARIBLES_FILE petscvariables
  HINTS
  ${PETSC_DIR}/${PETSC_ARCH}/conf
  ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf
  ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc-conf
)
if(NOT PETSC_VARIBLES_FILE)
  message(FATAL_ERROR ${PETSC_VARIBLES_FILE})
endif(NOT PETSC_VARIBLES_FILE)

file(STRINGS ${PETSC_VARIBLES_FILE} PETSC_VARIBLES)
foreach(LINE ${PETSC_VARIBLES})
  string(REGEX REPLACE " = " ";" FIELDS ${LINE})
  list(LENGTH FIELDS LISTLEN)
  if(LISTLEN EQUAL 2)
    list(GET FIELDS 0 VAR)
    list(GET FIELDS 1 VAL)
    set("PETSCVAR_${VAR}" ${VAL} CACHE INTERNAL "petsc varible")
    #message(STATUS PETSCVAR_${VAR})
  endif(LISTLEN EQUAL 2)
endforeach(LINE ${PETSC_VARIBLES})
