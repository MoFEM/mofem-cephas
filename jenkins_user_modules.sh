cat > script.cmake << EOF

set(GID_SOURCE_REPO "\$ENV{WORKSPACE}")
set(CTEST_SOURCE_DIRECTORY "\${GID_SOURCE_REPO}")
set(CTEST_BINARY_DIRECTORY "\$ENV{WORKSPACE}/build")

set(CTEST_SITE "Cactus")
set(CTEST_BRANCH "\$ENV{GIT_BRANCH}")
set(CTEST_BUILD_NAME "Linux-Jenkins-Branch(\${CTEST_BRANCH})")

find_program(CTEST_CONFIGURE_COMMAND NAMES spconfig.py HINTS \${CTEST_BINARY_DIRECTORY}  NO_DEFAULT_PATH)
find_program(CTEST_BUILD_COMMAND NAMES make)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

#Ctest time out
set(CTEST_TEST_TIMEOUT 12000)

ctest_start(Continuous)

ctest_configure(
  OPTIONS "-DSOURCE_DIR=\${CTEST_SOURCE_DIRECTORY} -DWITHCOVERAGE=ON -DMOFEM_UM_BUILD_TESTS=ON -DMOFEM_DIR=${WORKSPACE}/build/um_view")

ctest_build()
#if(CTEST_MEMORYCHECK_COMMAND)
#  ctest_memcheck()
#endif(CTEST_MEMORYCHECK_COMMAND)
ctest_test()
if(CTEST_COVERAGE_COMMAND)
  ctest_coverage()
endif(CTEST_COVERAGE_COMMAND)
ctest_submit()

EOF

#if [ ! $GIT_COMMIT == $GIT_PREVIOUS_COMMIT ]; then

  SPACK_ROOT=/home/jenkins/workspace/SpackInstall
  . $SPACK_ROOT/share/spack/setup-env.sh
  
  CORE_WORKSPACE=/home/jenkins/workspace/CoreDevelop
  
  if [ ! -f "$SPACK_ROOT/lock_spack" ]; then
  
    mkdir -p build
    cd build
    
    if [ ! -f "install_hash" ]; then
          
  	  	CORE_INSTALL_HASH=`cat $CORE_WORKSPACE/build/install_hash`

        cd $WORKSPACE/build
  		  spack view --verbose symlink -i ${WORKSPACE}/build/um_view /$CORE_INSTALL_HASH
        spack setup mofem-users-modules@develop copy_user_modules=False build_type=Debug install_id=2 ^/$CORE_INSTALL_HASH

    fi

  fi

#fi

cd $WORKSPACE/build
$WORKSPACE/build/spconfig.py \
	-DSOURCE_DIR=$WORKSPACE \
	-DMOFEM_UM_BUILD_TESTS=ON -DWITHCOVERAGE=ON \
	-DMOFEM_DIR=$WORKSPACE/build/um_view $WORKSPACE
        
grep CMAKE_INSTALL_PREFIX CMakeCache.txt |\
sed 's/.*mofem-users-modules-develop-//' |\
tee install_hash

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus start job for users modules branch $GIT_BRANCH\"}" \
https://hooks.slack.com/services/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

cd $WORKSPACE/build
spack load cmake
ctest -V -S $WORKSPACE/script.cmake
make install
make clean

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus end job for users modules branch $GIT_BRANCH (See http://cdash.eng.gla.ac.uk/cdash/)\"}" \
https://hooks.slack.com/services/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

