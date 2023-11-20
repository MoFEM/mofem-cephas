#!/bin/bash

echo MOFEM_SRC_DIR $MOFEM_SRC_DIR
echo MOFEM_INSTALL_DIR $MOFEM_INSTALL_DIR
echo MOFEM_BUILD_DIR $MOFEM_BUILD_DIR

source ${MOFEM_SRC_DIR}/scripts/docker_build_core_lib_script.sh
source ${MOFEM_SRC_DIR}/scripts/docker_build_users_modules_script.sh
