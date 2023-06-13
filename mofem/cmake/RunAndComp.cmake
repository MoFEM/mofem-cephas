if(BASE AND SPACE)
  message(
    "${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -base ${BASE} -space ${SPACE} -log_no_color")
  execute_process(
    COMMAND 
    ${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -base ${BASE} -space ${SPACE} -log_no_color
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
elseif(BASE)
  message(
    "${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -base ${BASE} -log_no_color")
  execute_process(
    COMMAND 
    ${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -base ${BASE} -log_no_color
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
elseif(SPACE)
  message(
    "${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -space ${SPACE} -log_no_color")
  execute_process(
    COMMAND 
    ${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -space ${SPACE} -log_no_color
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
else()
  message(
    "${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -log_no_color")
  execute_process(
    COMMAND 
    ${MPI_RUN} ${MPI_FLAGS} -np ${MPI_NP}
    ${PROG} -my_file ${FILE} -log_no_color
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
endif()

if(CMD_RESULT)
  message(FATAL_ERROR "Error running ${PROG} for file ${FILE}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files 
    ${SOURCE_DIR}/${LOG1}
    ${BINARY_DIR}/${LOG2}
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
if(CMD_RESULT)
    message(
      FATAL_ERROR 
      "Difftent files computed/blessed 
      ${SOURCE_DIR}/${LOG1} ${BINARY_DIR}/${LOG2}")
endif()