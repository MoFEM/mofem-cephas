cat > script.cmake << EOF

set(CTEST_PROJECT_NAME "Cephas-usersmodules")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.eng.gla.ac.uk")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=Cephas-usersmodules")
set(CTEST_DROP_SITE_CDASH TRUE)

set(GID_SOURCE_REPO "\$ENV{WORKSPACE}")
set(CTEST_SOURCE_DIRECTORY "\${GID_SOURCE_REPO}")
set(CTEST_BINARY_DIRECTORY "\$ENV{WORKSPACE}/build")

set(CTEST_SITE "Cactus")
set(CTEST_BRANCH "\$ENV{GIT_BRANCH}")
set(CTEST_BUILD_NAME "Fracture module (\${CTEST_BRANCH})")

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
#if(CTEST_COVERAGE_COMMAND)
#  ctest_coverage(QUIET)
#endif(CTEST_COVERAGE_COMMAND)
ctest_submit()
EOF


SPACK_ROOT=/home/jenkins/workspace/SpackInstall
. $SPACK_ROOT/share/spack/setup-env.sh

spack uninstall -y --dependents mofem-users-modules@lukasz build_type=RelWithDebInfo
spack install --no-cache mofem-users-modules@lukasz build_type=RelWithDebInfo
rm -rf  $WORKSPACE/um_view 
spack \
  view symlink -i $WORKSPACE/um_view \
  mofem-users-modules@lukasz build_type=RelWithDebInfo 

mkdir -p build
cd build

spack setup mofem-fracture-module@lukasz \
  copy_user_modules=False build_type=RelWithDebInfo \
  ^mofem-users-modules@lukasz build_type=RelWithDebInfo

./spconfig.py \
  -DMOFEM_UM_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DEXTERNAL_MODULE_SOURCE_DIRS=$WORKSPACE/extra_module \
  -DMOFEM_DIR=$WORKSPACE/um_view \
  $WORKSPACE/um_view/users_modules/

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus start job for fracture module branch $GIT_BRANCH\"}" \
https://hooks.slack.com/services/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

spack load cmake
ctest -V -S $WORKSPACE/script.cmake

make install
make clean

curl -X POST -H 'Content-type: application/json' \
--data "{\"text\":\"Jenkins cactus end job for fracture module branch $GIT_BRANCH (See http://cdash.eng.gla.ac.uk/cdash/)\"}" \
https://hooks.slack.com/services/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




