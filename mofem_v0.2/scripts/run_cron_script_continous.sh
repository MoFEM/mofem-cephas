#!/bin/sh

if test -f /home/lukasz/tests_cephas.lock
then
  echo "lock"
else
  CWD=`pwd`
  CTEST_SCRIPTS_FILE_PATH=/home/lukasz/mofem-cephas/mofem_v0.2/cmake
  CTEST_USER_MODULES_PATH=/home/lukasz/tmp/cephas_users_modules/users_modules/
  CTSET_SCRIPT=CTestScript_rdb-srv1_continous.cmake
  BUILD_DIR=/home/lukasz/tmp/cephas/build
  touch /home/lukasz/tests_cephas.lock
  cd $CTEST_SCRIPTS_FILE_PATH
  /usr/bin/ctest -VV --http1.0 -S $CTSET_SCRIPT >> /home/lukasz/tests_cephas.log 2>&1
  if [ -f /home/lukasz/tmp/cephas/source/has_bin_build ]; then
    rm /home/lukasz/tmp/cephas/source/has_bin_build
    cd $BUILD_DIR
    /usr/bin/make install
    chmod u+x $CTEST_USER_MODULES_PATH/scripts/run_cron_script.sh
    $CTEST_USER_MODULES_PATH/scripts/run_cron_script.sh
  fi  
  cd $CWD
  rm -v /home/lukasz/tests_cephas.lock
fi
