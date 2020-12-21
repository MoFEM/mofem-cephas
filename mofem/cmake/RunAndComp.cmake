# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
# It can be freely used for educational and research purposes
# by other institutions. If you use this softwre pleas cite my work.
#
# MoFEM is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# MoFEM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

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
    ${SOURCE_DIR}/${LOG1}
    ${BINARY_DIR}/${LOG2}
    WORKING_DIRECTORY ${BINARY_DIR}
    RESULT_VARIABLE CMD_RESULT)
if(CMD_RESULT)
    message(FATAL_ERROR "Error comparison for ${CMD} for file ${ARG}")
endif()