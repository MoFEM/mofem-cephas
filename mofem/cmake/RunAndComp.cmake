macro(EXEC_CHECK 
    mpi_run  np 
    cmd arg base space log1 log2 source_dir binary_dir)
    if(base AND space)
        execute_process(
           COMMAND 
           ${mpi_run} ${mpi_flg} -np ${np}
           ${cmd} -my_file ${arg} -base ${base} -space ${space} -log_no_color
           WORKING_DIRECTORY ${binary_dir}
           RESULT_VARIABLE CMD_RESULT)
    elseif(base) 
        execute_process(
           COMMAND 
           ${mpi_run} ${mpi_flg} -np ${np}
           ${cmd} -my_file ${arg} -base ${base} -log_no_color
           WORKING_DIRECTORY ${binary_dir}
           RESULT_VARIABLE CMD_RESULT)
    elseif(space)
        execute_process(
           COMMAND 
           ${mpi_run} ${mpi_flg} -np ${np}
           ${cmd} -my_file ${arg} -space ${space} -log_no_color
           WORKING_DIRECTORY ${binary_dir}
           RESULT_VARIABLE CMD_RESULT)
    else()
        execute_process(
           COMMAND            
           ${mpi_run} ${mpi_flg} -np ${np} 
           ${cmd} -my_file ${arg} -log_no_color
           WORKING_DIRECTORY ${binary_dir}
           RESULT_VARIABLE CMD_RESULT)
    endif()
    if(CMD_RESULT)
        message(FATAL_ERROR "Error running ${CMD} for file ${ARG}")
    endif()
endmacro()

if(BASE AND SPACE)
    exec_check(
        ${MPI_RUN} ${MPI_RUN_FLAGS} ${MPI_NP} 
        ${PROG} ${FILE} ${BASE} ${SPACE} ${LOG1} ${LOG2} ${SOURCE_DIR} ${BINARY_DIR})
elseif(BASE)
    exec_check(
        ${MPI_RUN} ${MPI_RUN_FLAGS} ${MPI_NP} 
        ${PROG} ${FILE} ${BASE} NO ${LOG1} ${LOG2} ${SOURCE_DIR} ${BINARY_DIR})
elseif(SPACE)
    exec_check(
        ${MPI_RUN} ${MPI_RUN_FLAGS} ${MPI_NP} 
        ${PROG} ${FILE} NO ${SPACE} ${LOG1} ${LOG2} ${SOURCE_DIR} ${BINARY_DIR}) 
else()
    exec_check(
        ${MPI_RUN} ${MPI_RUN_FLAGS} ${MPI_NP} 
        ${PROG} ${FILE} NO NO ${LOG1} ${LOG2} ${SOURCE_DIR} ${BINARY_DIR}) 
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files 
    ${SOURCE_DIR}/atom_tests/blessed_files/${LOG1}
    ${BINARY_DIR}/${LOG2}
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
if(CMD_RESULT)
    message(FATAL_ERROR "Error comparison for ${CMD} for file ${ARG}")
endif()