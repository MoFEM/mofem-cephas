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