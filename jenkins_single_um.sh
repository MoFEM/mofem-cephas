SPACK_ROOT=/var/lib/jenkins/workspace/SpackBuild
. $SPACK_ROOT/share/spack/setup-env.sh

spack install mofem-users-modules@develop build_type=RelWithDebInfo
spack \
  view symlink -i $WORKSPACE/um_view \
  mofem-users-modules@develop build_type=RelWithDebInfo 

mkdir -p build
cd build

spack setup mofem-fracture-module@develop \
  copy_user_modules=False build_type=RelWithDebInfo \
  ^mofem-users-modules@develop build_type=RelWithDebInfo

./spconfig.py \
  -DMOFEM_UM_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DEXTERNAL_MODULE_SOURCE_DIRS=$WORKSPACE/extra_module \
  -DMOFEM_DIR=$WORKSPACE/um_view \
  $WORKSPACE/um_view/users_modules/

spack load cmake
ctest -D Experimental

make install

