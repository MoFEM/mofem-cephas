#!/bin/sh

export LD_LIBRARY_PATH=/opt/local_boost_1_65_1/lib

if [ -e $HOME/tests_cephas.lock ]
then
  echo "lock"
else
  CWD=`pwd`
  CTEST_SCRIPTS_FILE_PATH=$HOME/mofem-cephas/mofem/cmake
  CTEST_USER_MODULES_PATH=$HOME/tmp/cephas_users_modules/users_modules/
  CTSET_SCRIPT=CTestScript_rdb-srv1.cmake
  BUILD_DIR=$HOME/tmp/cephas/build
  touch $HOME/tests_cephas.lock
  cd $CTEST_SCRIPTS_FILE_PATH
  /opt/local/bin/ctest -VV --http1.0 -S $CTSET_SCRIPT >> $HOME/tests_cephas.log 2>&1
  if [ -e $HOME/tmp/cephas/source/has_bin_build ]; then
    rm $HOME/tmp/cephas/source/has_bin_build
    cd $BUILD_DIR
    /usr/bin/make install
    chmod u+x $CTEST_USER_MODULES_PATH/scripts/run_cron_script.sh
    $CTEST_USER_MODULES_PATH/scripts/run_cron_script.sh
  fi
  cd $CWD
  rm -v $HOME/tests_cephas.lock
fi
