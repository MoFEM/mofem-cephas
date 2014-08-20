#!/bin/sh

CTEST_SCRIPTS_FILE_PATH=/tmp/cephas/debug_examples/users_modules/scripts
CTSET_SCRIPT=CTestScript_rdb-srv1.cmake
CWD=`pwd`

cd $CTEST_SCRIPTS_FILE_PATH
/usr/bin/ctest -VV --http1.0 -S $CTSET_SCRIPT >> /home/lukasz/tests_users_modules.log 2>&1
cd $CWD
