#!/bin/bash

set -e

echo MOFEM_SRC_DIR $MOFEM_SRC_DIR
echo MOFEM_INSTALL_DIR $MOFEM_INSTALL_DIR
echo MOFEM_BUILD_DIR $MOFEM_BUILD_DIR
mkdir -p $MOFEM_BUILD_DIR
cd $MOFEM_BUILD_DIR

# Determine number of cores to compile code

NBCORES=$(cat /proc/cpuinfo | grep processor | wc -l)
NB_MAX=12

if (($NBCORES > $NB_MAX)); then
  NB=$NB_MAX
else
  NB=$NBCORES
fi
echo Nb. of cores $NBCORES and nb. of cores used to compilation $NB

echo "Configure"
/opt/local/bin/cmake \
  -DCMAKE_BUILD_TYPE=MinSizeRel \
  -DCMAKE_CXX_FLAGS="-Wall" \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DMOAB_DIR=$PETSC_DIR/$PETSC_ARCH \
  -DBUILD_SHARED_LIBS=yes \
  -DWITH_ADOL-C=1 \
  -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR $MOFEM_SRC_DIR
echo "Build"
make -j $NB install
ctest --output-on-failure
make clean

echo "Configure users modules"
cd $MOFEM_INSTALL_DIR
/opt/local/bin/cmake -DBUILD_SHARED_LIBS=yes -DCMAKE_CXX_FLAGS="-Wall" users_modules
echo "Build users modules"
cd $MOFEM_INSTALL_DIR/basic_finite_elements/atom_tests
make -j $NB
ctest --output-on-failure
make clean
echo "All done"
