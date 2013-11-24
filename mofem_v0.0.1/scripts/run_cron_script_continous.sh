#!/bin/sh

CTEST_SCRIPTS_FILE_PATH=/home/lukasz/bitbucket/mofem-joseph/mofem_v0.0.1/cmake
CTSET_SCRIPT=CTestScript_rdb-srv1_continous.cmake

if test -f /home/lukasz/tests.lock
then
  echo "lock"
else
  touch /home/lukasz/tests.lock
  /usr/bin/ctest -VV --http1.0 -S $CTEST_SCRIPTS_FILE_PATH/$CTSET_SCRIPT >> /home/lukasz/tests.log 2>&1
  rm -v /home/lukasz/tests.lock
fi
