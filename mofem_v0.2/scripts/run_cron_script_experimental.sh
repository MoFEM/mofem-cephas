#!/bin/sh

CTEST_SCRIPTS_FILE_PATH=/home/lukasz/bitbucket/mofem-joseph/mofem_v0.1/cmake
CTSET_SCRIPT=CTestScript_rdb-srv1_experimental.cmake
CWD=`pwd`

cd $CTEST_SCRIPTS_FILE_PATH
/usr/bin/ctest -VV --http1.0 -S $CTSET_SCRIPT >> /home/lukasz/tests.log 2>&1
cd $CWD
