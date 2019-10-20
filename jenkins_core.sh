cat > script.cmake << EOF

set(GID_SOURCE_REPO "\$ENV{WORKSPACE}")
set(CTEST_SOURCE_DIRECTORY "\${GID_SOURCE_REPO}/mofem")
set(CTEST_BINARY_DIRECTORY "\$ENV{WORKSPACE}/build")

set(CTEST_SITE "Jenkins")
set(CTEST_BRANCH "\$ENV{GIT_BRANCH}")
set(CTEST_BUILD_NAME "Linux-Jenkins=BranchO(\${CTEST_BRANCH})")

find_program(CTEST_CONFIGURE_COMMAND NAMES spconfig.py HINTS \${CTEST_BINARY_DIRECTORY} NO_DEFAULT_PATH)
find_program(CTEST_BUILD_COMMAND NAMES make)
find_program(CTEST_COVERAGE_COMMAND NAMES gcov)

#Ctest time out
set(CTEST_TEST_TIMEOUT 12000)

ctest_start(Continuous)

ctest_configure(
  OPTIONS "-DSOURCE_DIR=\${CTEST_SOURCE_DIRECTORY} -DWITHCOVERAGE=ON -DMOFEM_BUILD_TESTS=ON")
  
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

  SPACK_ROOT=/var/lib/jenkins/workspace/SpackBuild
  . $SPACK_ROOT/share/spack/setup-env.sh

  if [ ! -f "$SPACK_ROOT/lock_spack" ]; then

    mkdir -p build
    cd build
  
    if [ ! -f "install_hash" ]; then
    
    	spack setup mofem-cephas@develop+slepc copy_user_modules=False build_type=Debug
      
      $WORKSPACE/build/spconfig.py \
        	-DSOURCE_DIR=$WORKSPACE/mofem \
        	-DMOFEM_BUILD_TESTS=ON -DWITHCOVERAGE=ON ../mofem
        
    	grep CMAKE_INSTALL_PREFIX CMakeCache.txt |\
    	sed 's/.*mofem-cephas-develop-//' |\
    	tee install_hash
        
    fi
    

  fi

#fi

