#!/bin/sh

CTEST_SCRIPTS_FILE_PATH=/home/lukasz/mofem-cephas/mofem_v0.2/cmake
CTSET_SCRIPT=CTestScript_rdb-srv1_continous.cmake
CWD=`pwd`

if test -f /home/lukasz/tests.lock
then
  echo "lock"
else
  touch /home/lukasz/tests.lock
  cd $CTEST_SCRIPTS_FILE_PATH
  /usr/bin/ctest -VV --http1.0 -S $CTSET_SCRIPT >> /home/lukasz/tests.log 2>&1
  cd $CWD
  rm -v /home/lukasz/tests.lock
fi
