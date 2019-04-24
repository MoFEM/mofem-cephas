#!/bin/bash

# Determine number of cores to compile code
NBCORES=$(cat /proc/cpuinfo | grep processor | wc -l)
NB_MAX=12
if (($NBCORES > $NB_MAX)); then
  NB=$NB_MAX
else
  NB=$NBCORES
fi
echo Nb. of cores $NBCORES and nb. of cores used to compilation $NB

set -e

echo "Configure"
# Make build directory
mkdir -p $MOFEM_BUILD_DIR

cd $MOFEM_BUILD_DIR
# Configure mofem core lib
/opt/local/bin/cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-Wall -Wno-sign-compare" \
  -DMPI_RUN_FLAGS="--allow-run-as-root" \
  -DPETSC_DIR=$PETSC_DIR \
  -DPETSC_ARCH=$PETSC_ARCH \
  -DMOAB_DIR=$PETSC_DIR/$PETSC_ARCH \
  -DMED_DIR=/opt/med \
  -DADOL-C_DIR=$PETSC_DIR/$PETSC_ARCH \
  -DTETGEN_DIR=/opt/tetgen1.5.0 \
  -DBUILD_SHARED_LIBS=yes \
  -DCMAKE_INSTALL_PREFIX=$MOFEM_INSTALL_DIR $MOFEM_SRC_DIR

# Install and build core library
echo "Build"
make -k -j $NB install

# Run tests and send results to CDash
ctest --output-on-failure -D Experimental
make clean

#-DCGM_DIR=$PETSC_DIR/$PETSC_ARCH \
#-DMESHKIT_DIR=$PETSC_DIR/$PETSC_ARCH \
