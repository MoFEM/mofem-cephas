#!/bin/sh

cd /home/lukasz/bitbucket/mofem-joseph/mofem_v0.0.1/
/usr/bin/ctest -VV --http1.0 -S CTestScript_rdb-srv1_experimental.cmake  >> /home/lukasz/tests.log 2>&1

