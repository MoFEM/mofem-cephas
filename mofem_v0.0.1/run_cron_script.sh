#!/bin/sh

if test -f /home/lukasz/tests.lock
then
  echo "lock"
else
  touch /home/lukasz/tests.lock
  cd /home/lukasz/bitbucket/mofem-joseph/mofem_v0.0.1/
  /usr/bin/ctest -VV --http1.0 -S CTestScript_rdb-srv1.cmake >> /home/lukasz/tests.log 2>&1
  rm -v /home/lukasz/tests.lock
fi
